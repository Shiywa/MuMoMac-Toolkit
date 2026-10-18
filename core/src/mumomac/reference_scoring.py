"""Unified reference scores; ranks shared across all signatures and versions."""
from importlib import metadata
import hashlib,json,time
import numpy as np,pandas as pd,anndata as ad
from .validation import prepare_input
from .reference import REFERENCE_ORDER,reference_signatures,load_reference
from .reliability import input_fingerprint,symbol_fingerprint

def _plans(p,species,versions,coverage_policy):
    if coverage_policy not in ('reference','complete','available'):raise ValueError('coverage_policy must be reference, complete or available')
    present=set(p.symbols);index={g:i for i,g in enumerate(p.symbols)};rows=[];plans={};signatures={}
    for version in versions:
      for pid,s in reference_signatures(version,species).items():
        full=s['genes'];found=[g for g in full if g in present];absent=[g for g in full if g not in present]
        og=[g for g in s['original_genes'] if g in present];ag=[g for g in s['added_genes'] if g in present]
        cov=len(found)/len(full) if full else 0.;ret=len(full)/len(s['mouse_genes'])
        inherited=pid.startswith('S');component=s['scoring_formula'].startswith('0.5_')
        if coverage_policy=='reference':
            eligible=(not absent and bool(found)) if inherited else cov>=.8 and ret>=.8 and len(found)>=3
        elif coverage_policy=='complete':eligible=not absent and bool(found)
        else:eligible=bool(found)
        if component:eligible=eligible and bool(og) and bool(ag)
        eligible=bool(eligible and not p.zero_cells.all());key=f'{pid}__{version}'
        used=og+ag if component else found
        counts=np.asarray((p.matrix[:,[index[g] for g in found]]>0).sum(1)).ravel() if found else np.zeros(p.adata.n_obs,dtype=int)
        rows.append(dict(feature_id=pid,display_name=s['display_name'],feature_class=s['feature_class'],signature_version=version,gene_set_version=s['gene_set_version'],gene_set_hash=s['gene_set_hash'],reference_version='1.0.0',scoring_formula=s['scoring_formula'],n_reference_genes=len(s['mouse_genes']),n_mapped_genes=len(full),n_present_genes=len(found),mapping_retention=ret,coverage=cov,coverage_policy=coverage_policy,scoring_eligible=eligible,applicability_status='DESCRIPTIVE_MONOCYTE_ONLY; HISTORICAL_MACROPHAGE_SCOPE' if pid.startswith('L') else 'CONTINUOUS_ONLY' if pid=='A01' else 'EXPLORATORY_RNA',median_cell_detected_members=float(np.median(counts)),fraction_cells_zero_members=float(np.mean(counts==0)),present_genes=';'.join(found),absent_genes=';'.join(absent),unmapped_genes=';'.join(s['unmapped_genes']),developmental_candidate=s['developmental_candidate'],construct_validity=s['construct_validity'],default_scope=s['default_scope'],biological_limits=s['biological_limits']))
        plans[key]={'eligible':eligible,'component':component,'counts':counts,'spec':s}
        if eligible:
            if component:signatures[key+'_seed']=og;signatures[key+'_added']=ag
            else:signatures[key]=found
    return pd.DataFrame(rows),plans,signatures

def get_program_coverage(adata,*,species='mouse',signature_version='both',layer='counts',expression_type=None,gene_symbol_key=None,duplicate_policy='error',coverage_policy='reference',chunk_size=500):
    versions=['original','expanded'] if signature_version=='both' else [signature_version]
    p=prepare_input(adata,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key,duplicate_policy=duplicate_policy,chunk_size=chunk_size)
    return _plans(p,species,versions,coverage_policy)[0]

def score_signature_versions(adata,*,species='mouse',signature_version='both',layer='counts',expression_type=None,gene_symbol_key=None,duplicate_policy='error',coverage_policy='reference',max_rank=1500,chunk_size=500,copy=True,overwrite=False):
    if signature_version not in ('original','expanded','both'):raise ValueError('Choose original, expanded or both')
    if metadata.version('pyucell')!='0.8.0':raise RuntimeError('Pinned pyucell 0.8.0 required')
    start=time.perf_counter();p=prepare_input(adata,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key,duplicate_policy=duplicate_policy,chunk_size=chunk_size)
    if p.adata.isbacked:raise ValueError('Full workflow loads h5ad in memory; backed streaming is not implemented')
    if 'mumomac_reference' in p.adata.uns and not overwrite:raise ValueError('Reference output exists; explicitly overwrite')
    versions=['original','expanded'] if signature_version=='both' else [signature_version]
    coverage,plans,sigs=_plans(p,species,versions,coverage_policy)
    if not isinstance(max_rank,int) or max_rank<1 or (sigs and max(map(len,sigs.values()))>max_rank):raise ValueError('Invalid max_rank for largest component')
    from pyucell import get_rankings,compute_scores_from_ranks
    values={v:np.full((p.adata.n_obs,22),np.nan,np.float32) for v in versions}
    detected={v:np.zeros((p.adata.n_obs,22),dtype=np.int16) for v in versions}
    for v in versions:
      for j,pid in enumerate(REFERENCE_ORDER):detected[v][:,j]=plans[f'{pid}__{v}']['counts']
    for lo in range(0,p.adata.n_obs,chunk_size):
      hi=min(p.adata.n_obs,lo+chunk_size)
      if not sigs:break
      b=ad.AnnData(p.matrix[lo:hi,:].copy(),var=pd.DataFrame(index=pd.Index(p.symbols)))
      ranks=get_rankings(b,max_rank=max_rank,ties_method='average',device='cpu')
      compute_scores_from_ranks(b,ranks,sigs,max_rank=max_rank,missing_genes='skip',suffix='')
      for v in versions:
       for j,pid in enumerate(REFERENCE_ORDER):
        key=f'{pid}__{v}';plan=plans[key]
        if plan['eligible']:
          if plan['component']:
            original=b.obs[key+'_seed'].to_numpy(np.float32);added=b.obs[key+'_added'].to_numpy(np.float32)
            values[v][lo:hi,j]=.5*(original+added)
          else:values[v][lo:hi,j]=b.obs[key].to_numpy(np.float32)
    for v in versions:values[v][p.zero_cells,:]=np.nan
    out=p.adata.copy() if copy else p.adata
    for v in versions:
      out.obsm[f'X_mumomac_reference_{v}']=values[v]
      out.obsm[f'X_mumomac_detected_{v}']=detected[v]
    config=dict(reference_version='1.0.0',species=species,signature_version=signature_version,versions=versions,input_layer=p.provenance['layer'],expression_type=p.expression_type,gene_symbol_key=gene_symbol_key or '',duplicate_policy=duplicate_policy,coverage_policy=coverage_policy,max_rank=max_rank,ties_method='average',scorer='pyucell0.8.0/publicCPU',gene_universe_hash=hashlib.sha256(json.dumps(p.symbols).encode()).hexdigest())
    confighash=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    out.uns['mumomac_reference']=dict(configuration=config,score_config_hash=confighash,feature_order=list(REFERENCE_ORDER),gene_coverage=coverage,input_provenance=p.provenance,input_fingerprint=input_fingerprint(p.adata,config['input_layer']),symbol_fingerprint=symbol_fingerprint(p.adata,gene_symbol_key),score_hashes={v:hashlib.sha256(values[v].tobytes()).hexdigest() for v in versions},timing_seconds=time.perf_counter()-start,positive_calls_enabled=False,interpretation='Continuous RNA features; L04 transition-associated compatibility feature, L05 auxiliary, A01 covariate; no function superiority certification')
    return out

def score_historical_raw(adata,*,layer='X',expression_type='log1p',**kwargs):
    """Explicit legacy raw mean only. Adjusted scores/calls require saved calibration.

    Input must be the retained-gene log1pCP10k expression used historically.
    This function never learns q20/GMM thresholds from a new query.
    """
    if expression_type!='log1p':raise ValueError('Historical raw requires declared log1pCP10k')
    p=prepare_input(adata,layer=layer,expression_type=expression_type,**kwargs)
    index={g:i for i,g in enumerate(p.symbols)};rows={}
    for pid,s in reference_signatures('original','mouse').items():
        if not set(s['genes'])<=set(index):rows[pid]=np.full(p.adata.n_obs,np.nan)
        else:rows[pid]=np.asarray(p.matrix[:,[index[g] for g in s['genes']]].mean(1)).ravel()
    out=pd.DataFrame(rows,index=p.adata.obs_names)
    out.attrs['historical_scope']='Raw original-member mean only; normalization asserted by caller; no q20/GMM calibration or positive calls'
    return out
