"""Small adapter to pinned public pyUCell APIs; no UCell reimplementation."""
from importlib import metadata
from concurrent.futures import ThreadPoolExecutor
import datetime,hashlib,json,time,warnings
import numpy as np,pandas as pd
import anndata as ad
from .validation import prepare_input,_program_plan
from .signatures import PROGRAM_ORDER,signature_manifest
from .reliability import input_fingerprint,symbol_fingerprint

def score_programs(adata,*,species='mouse',signature_version='original',scoring_method='ucell',
                   layer='counts',expression_type=None,gene_symbol_key=None,
                   duplicate_policy='error',coverage_policy='complete',min_genes=1,
                   max_rank=1500,ties_method='average',chunk_size=500,n_jobs=1,
                   copy=True,overwrite=False,add_obs=False,obs_prefix='mumomac_',feature_scope='s16'):
    """Compute definition-transfer RNA scores, with frozen resources and explicit loss."""
    if feature_scope=='reference' or signature_version=='both':
        from .reference_scoring import score_signature_versions
        if scoring_method!='ucell' or ties_method!='average':raise ValueError('Reference scoring freezes public UCell average ties')
        if n_jobs!=1:raise ValueError('Full reference uses shared-rank single-worker chunks; n_jobs is currently legacy S16 only')
        if add_obs:raise ValueError('Reference API exports version-specific obsm matrices; add_obs is legacy S16 only')
        return score_signature_versions(adata,species=species,signature_version=signature_version,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key,duplicate_policy=duplicate_policy,coverage_policy='reference' if coverage_policy=='complete' else coverage_policy,max_rank=max_rank,chunk_size=chunk_size,copy=copy,overwrite=overwrite)
    if feature_scope!='s16':raise ValueError('feature_scope must be s16 or reference')
    if scoring_method!='ucell':raise ValueError('v0.1 supports only the frozen portable ucell method.')
    if ties_method!='average':raise ValueError('v0.1 freezes ties_method=average; other conventions are unsupported.')
    if not isinstance(max_rank,int) or max_rank<1:raise ValueError('max_rank must be a positive integer.')
    if not isinstance(n_jobs,int) or n_jobs<1:raise ValueError('n_jobs must be a positive worker count.')
    started=time.perf_counter()
    p=prepare_input(adata,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key,
                    duplicate_policy=duplicate_policy,chunk_size=chunk_size)
    if not overwrite and ('X_mumomac_s16' in p.adata.obsm or 'mumomac' in p.adata.uns):raise ValueError('MuMoMac output already exists; explicitly set overwrite=True for reserved outputs.')
    cols=[obs_prefix+pid for pid in PROGRAM_ORDER]
    if add_obs and any(c in p.adata.obs for c in cols):raise ValueError('Requested obs output would overwrite an existing field; use another prefix.')
    report,plan=_program_plan(p,species,signature_version,coverage_policy,min_genes)
    signatures={}
    for pid in PROGRAM_ORDER:
        x=plan[pid]
        if x['eligible']:
            signatures['O_'+pid]=x['original']
            if x['has_added']:signatures['A_'+pid]=x['added']
    if signatures and max(map(len,signatures.values()))>max_rank:
        raise ValueError('max_rank must be at least the largest scored component signature to keep the truncated UCell definition valid.')
    if species=='human':warnings.warn('Human scoring is an Experimental one-to-one definition extension; mapping release is Missing and cross-species absolute calibration is unverified.',UserWarning,stacklevel=2)
    if not report.scoring_eligibility.all():warnings.warn(f'{int((~report.scoring_eligibility).sum())} programs are Not Evaluable; see gene coverage.',UserWarning,stacklevel=2)
    if metadata.version('pyucell')!='0.8.0':raise RuntimeError('This scorer adapter was verified with pyucell==0.8.0; incompatible installed version.')
    from pyucell import get_rankings,compute_scores_from_ranks
    prepared_time=time.perf_counter()-started;matrix=np.full((p.adata.n_obs,16),np.nan,dtype=np.float32)
    def compute(start,end):
        b=p.matrix[start:end,:].copy()
        t=ad.AnnData(b,var=pd.DataFrame(index=pd.Index(p.symbols)),obs=pd.DataFrame(index=[str(i) for i in range(end-start)]))
        ranks=get_rankings(t,max_rank=max_rank,ties_method=ties_method,device='cpu')
        compute_scores_from_ranks(t,ranks,signatures,max_rank=max_rank,missing_genes='skip',suffix='')
        out=np.full((end-start,16),np.nan,np.float32)
        for j,pid in enumerate(PROGRAM_ORDER):
            x=plan[pid]
            if x['eligible']:
                out[:,j]=t.obs['O_'+pid].to_numpy(dtype=np.float32)
                if x['has_added']:out[:,j]=.5*(out[:,j]+t.obs['A_'+pid].to_numpy(dtype=np.float32))
        out[p.zero_cells[start:end],:]=np.nan
        return start,end,out
    compute_start=time.perf_counter()
    if signatures:
        with ThreadPoolExecutor(max_workers=n_jobs) as pool:
            # Bound pending chunks to worker count rather than materializing all ranks.
            for base in range(0,p.adata.n_obs,chunk_size*n_jobs):
                jobs=[pool.submit(compute,s,min(s+chunk_size,p.adata.n_obs)) for s in range(base,min(base+chunk_size*n_jobs,p.adata.n_obs),chunk_size)]
                for f in jobs:s,e,value=f.result();matrix[s:e]=value
    scoring_time=time.perf_counter()-compute_start
    if np.isinf(matrix).any():raise RuntimeError('Unexpected infinite scorer output.')
    for j,pid in enumerate(PROGRAM_ORDER):
        if plan[pid]['eligible'] and not np.isfinite(matrix[~p.zero_cells,j]).all():raise RuntimeError(f'Unexpected missing scorer output for {pid}.')
    out=p.adata.copy() if copy else p.adata
    manifest=signature_manifest();universe_hash=hashlib.sha256(json.dumps(p.symbols,ensure_ascii=False).encode()).hexdigest()
    config={'package_version':metadata.version('mumomac-toolkit'),'scorer_version':'pyucell0.8.0/publicCPU','signature_version':signature_version,
            'definition_version':manifest['definition_versions'][signature_version],'species':species,'scoring_method':scoring_method,
            'max_rank':max_rank,'ties_method':ties_method,'coverage_policy':coverage_policy,'min_genes_engineering':min_genes,
            'effective_gene_universe_sha256':universe_hash,'input_layer':p.provenance['layer'],'expression_type':p.expression_type,
            'duplicate_policy':duplicate_policy,'expanded_weights':'0.5 Original +0.5 Added; no-added=>Original'}
    config_hash=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    versions={name:metadata.version(name) for name in ['pyucell','anndata','numpy','scipy','pandas']}
    out.obsm['X_mumomac_s16']=matrix
    out.uns['mumomac']={'schema_version':'1','program_order':list(PROGRAM_ORDER),'configuration':config,
                       'input_fingerprint':input_fingerprint(p.adata,config['input_layer']),
                       'symbol_metadata_fingerprint':symbol_fingerprint(p.adata,gene_symbol_key),
                       'score_config_hash':config_hash,'software_versions':versions,'gene_coverage':report,
                       'execution':{'chunk_size':chunk_size,'n_jobs':n_jobs,'device':'cpu','copy':copy},
                       'definition_resource_hashes':manifest['resource_sha256'],
                       'species_mapping':manifest['source_files']['human_mapping'] if species=='human' else {'status':'Same-species mouse definitions'},
                       'input_provenance':p.provenance,'zero_expression_cell_indices':np.flatnonzero(p.zero_cells),
                       'timing_seconds':{'validation_preparation':prepared_time,'scoring':scoring_time,'total_before_output_copy':time.perf_counter()-started},
                       'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                       'interpretation':'Exploratory RNA definition transfer; not full Atlas model, function measurement, identity classifier or cross-study numerical calibration.'}
    if add_obs:
        for j,c in enumerate(cols):out.obs[c]=matrix[:,j]
    return out
