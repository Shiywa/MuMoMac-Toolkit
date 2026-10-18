"""User-identity, biological-unit comparisons and explicit exploratory combinations."""
from importlib import resources
import itertools,json,math
import numpy as np,pandas as pd
from scipy import stats
from .io import read_adata
from .reference import REFERENCE_ORDER,load_reference,reference_signatures
from .input import validate_input_metadata
from .comparison import compare_programs,_adjust

def assess_identity(adata,*,identity_key,marker_sets=None,layer='counts',expression_type=None,gene_symbol_key=None):
    a=read_adata(adata)
    if identity_key not in a.obs:raise KeyError(identity_key)
    r=json.loads(resources.files('mumomac.resources').joinpath('identity_reference_v1.json').read_text())
    canonical=set(r['identities']);aliases={'Kupffer cells':'Kupffer cell','Alveolar Macrophages':'Alveolar macrophage','Cross-tissue shared macrophages':'Macrophage, tissue specialization unresolved'}
    rows=[];marker_sets=marker_sets or {}
    if marker_sets:
      from .validation import prepare_input
      p=prepare_input(a,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key);index={g:i for i,g in enumerate(p.symbols)}
    for label,n in a.obs[identity_key].value_counts(dropna=False).items():
        c=aliases.get(label,label);known=c in canonical
        markers=marker_sets.get(label,marker_sets.get(c,[]));found=[g for g in markers if g in index] if markers else [];fraction=np.nan
        if found:
            cells=np.flatnonzero(a.obs[identity_key].eq(label).to_numpy());fraction=float(np.mean(np.asarray((p.matrix[cells,:][:,[index[g] for g in found]]>0).sum(1)).ravel()>0))
        overlaps=sorted(set(markers)&set(g for f in load_reference()['features'] for g in f['signatures']['original']['genes']))
        rows.append(dict(input_label=str(label),reference_label=c if known else 'unknown',display_label=r['display_mapping'].get(c,c) if known else str(label),n_cells=int(n),mapping_evidence='User-provided exact label or explicit spelling alias' if known else 'User identity outside frozen vocabulary; not force-mapped',applicable_reference_class='monocyte' if known and 'monocyte' in c.lower() else 'macrophage' if known else 'unknown',marker_check_status='MEASURED_DESCRIPTIVE' if markers else 'NOT_TESTED_NO_DECLARED_MARKER_SET',marker_genes=';'.join(markers),present_marker_genes=';'.join(found),fraction_cells_any_declared_marker=fraction,marker_program_gene_overlap=';'.join(overlaps),uncertainty='Marker check is descriptive; overlapping markers cannot independently validate programs. Fine identity not independently certified; no classifier run',shared_category_interpretation='Cross-tissue shared macrophages is a broad annotation category, not common ontogeny' if c in r['display_mapping'] else ''))
    return pd.DataFrame(rows)

def aggregate_program_scores(adata,*,sample_key,identity_key=None,condition_key=None,biological_unit_key=None,study_key=None,tissue_key=None,min_cells=1):
    a=read_adata(adata)
    if 'mumomac_reference' not in a.uns:raise ValueError('Run score_signature_versions first')
    if min_cells<1:raise ValueError('min_cells must be positive')
    validate_input_metadata(a,sample_key=sample_key,biological_unit_key=biological_unit_key,condition_key=condition_key,identity_key=identity_key,study_key=study_key)
    keys=list(dict.fromkeys(k for k in [sample_key,study_key,identity_key,tissue_key,condition_key,biological_unit_key] if k))
    for k in keys:
      if k not in a.obs:raise KeyError(k)
      if a.obs[k].isna().any():raise ValueError(f'Missing grouping metadata {k}; formal biological mapping cannot be inferred')
    rows=[];meta=a.uns['mumomac_reference'];versions=list(meta['configuration']['versions'])
    for key,ids in a.obs.groupby(keys,observed=True,dropna=False,sort=False).indices.items():
      key=key if isinstance(key,tuple) else (key,);base=dict(zip(keys,key));ids=np.array(ids)
      for v in versions:
       x=np.asarray(a.obsm[f'X_mumomac_reference_{v}'])[ids]
       for j,pid in enumerate(REFERENCE_ORDER):
        value=x[:,j];finite=value[np.isfinite(value)]
        rows.append(dict(**base,feature_id=pid,signature_version=v,n_cells=len(ids),n_valid_cells=len(finite),mean_score=float(finite.mean()) if len(finite) else np.nan,median_score=float(np.median(finite)) if len(finite) else np.nan,score_sd=float(finite.std(ddof=1)) if len(finite)>1 else np.nan,available_program=bool(len(finite) and len(ids)>=min_cells),biological_unit_valid=bool(biological_unit_key),score_config_hash=meta['score_config_hash']))
    out=pd.DataFrame(rows)
    out.attrs['reference_summary']={'sample_key':sample_key,'biological_unit_key':biological_unit_key or '', 'condition_key':condition_key or '', 'identity_key':identity_key or '', 'study_key':study_key or '', 'tissue_key':tissue_key or '', 'strata_keys':list(dict.fromkeys(k for k in [study_key,identity_key,tissue_key] if k)), 'minimum_cells':min_cells,'biological_assertion':'Caller declares actual donor/animal; software verifies consistency, not biological truth'}
    return out

def biological_unit_scores(summary):
    m=summary.attrs.get('reference_summary',{});unit=m.get('biological_unit_key','')
    if not unit:raise ValueError('Actual donor/animal key missing')
    keys=list(dict.fromkeys([unit]+m['strata_keys']+([m['condition_key']] if m['condition_key'] else [])+['feature_id','signature_version']))
    z=summary[summary.available_program & summary.mean_score.notna()]
    out=z.groupby(keys,observed=True,dropna=False).agg(mean_score=('mean_score','mean'),n_libraries=(m['sample_key'],'nunique'),n_cells=('n_cells','sum')).reset_index()
    out.attrs.update(summary.attrs);return out

def compare_conditions(summary,*,reference,case,condition_key=None,biological_unit_key=None,paired=False,n_bootstrap=2000,n_permutations=9999,seed=20261014,correction='bh'):
    meta=summary.attrs.get('reference_summary',{})
    if not meta:raise ValueError('Missing traceable reference summary metadata')
    group=condition_key or meta.get('condition_key');unit=biological_unit_key or meta.get('biological_unit_key') or None
    if not group:raise ValueError('condition_key is required')
    if unit and unit!=meta.get('biological_unit_key'):raise ValueError('Use the declared actual donor/animal key, not a replacement library key')
    if not paired:
      z=summary.rename(columns={'feature_id':'program_id'}).copy()
      if unit and meta.get('study_key'):
          z[unit]=z[meta['study_key']].astype(str)+'::'+z[unit].astype(str)
      z.attrs['mumomac_summary']={'sample_key':meta['sample_key'],'strata_keys':meta['strata_keys']+['signature_version'],'group_key':group,'biological_unit_key':unit or ''}
      out=compare_programs(z,group_key=group,reference=reference,case=case,biological_unit_key=unit,n_bootstrap=n_bootstrap,n_permutations=n_permutations,seed=seed,correction=correction)
      out=out.rename(columns={'program_id':'feature_id'});out.attrs['comparison_design']='Independent units, study-specific; equal library means within each unit'
      return out
    if not unit:raise ValueError('Paired design requires actual donor/animal IDs')
    if reference==case:raise ValueError('Groups must differ')
    if n_bootstrap<1 or n_permutations<1:raise ValueError('Resampling counts must be positive')
    z=biological_unit_scores(summary);by=meta['strata_keys']+['signature_version','feature_id'];rng=np.random.default_rng(seed);rows=[]
    for keys,d in z.groupby(by,observed=True,dropna=False):
      keys=keys if isinstance(keys,tuple) else (keys,);base=dict(zip(by,keys))
      p=d[d[group].isin([reference,case])].pivot(index=unit,columns=group,values='mean_score')
      if reference not in p or case not in p:continue
      p=p[[reference,case]].dropna();delta=(p[case]-p[reference]).to_numpy();n=len(delta)
      row=dict(**base,reference=reference,case=case,n_paired_units=n,mean_difference=float(delta.mean()) if n else np.nan,standardized_effect=float(delta.mean()/delta.std(ddof=1)) if n>1 and delta.std(ddof=1)>0 else np.nan,ci_low=np.nan,ci_high=np.nan,p_value=np.nan,p_adjusted=np.nan,inference_status='Limited: fewer than 3 complete pairs',statistical_unit='Declared biological paired unit')
      if n>=3:
       boots=delta[rng.integers(0,n,(n_bootstrap,n))].mean(1);row['ci_low'],row['ci_high']=np.quantile(boots,[.025,.975])
       if n<=16:signs=np.array(list(itertools.product([-1,1],repeat=n)));vals=(signs*delta).mean(1);pv=float(np.mean(abs(vals)>=abs(delta.mean())-1e-12));method='Exact paired sign permutation'
       else:signs=rng.choice([-1,1],(n_permutations,n));vals=(signs*delta).mean(1);pv=(1+int(np.sum(abs(vals)>=abs(delta.mean())-1e-12)))/(n_permutations+1);method='Monte Carlo paired sign permutation'
       row.update(p_value=pv,p_method=method,inference_status='Paired RNA association; symmetry/exchangeability assumed; limited sample precision')
      rows.append(row)
    out=pd.DataFrame(rows)
    if len(out):
      p=out.p_value.to_numpy(float);q=_adjust(np.where(np.isfinite(p),p,1.),correction);q[~np.isfinite(p)]=np.nan;out['p_adjusted']=q
    out.attrs['comparison_design']='Paired biological-unit comparison within each study/identity/tissue/version';return out

def compare_within_identity(summary,**kwargs):
    if not summary.attrs.get('reference_summary',{}).get('identity_key'):raise ValueError('identity_key required')
    return compare_conditions(summary,**kwargs)

def summarize_composition(adata,*,sample_key,identity_key,condition_key=None,biological_unit_key=None,study_key=None):
    a=read_adata(adata);keys=list(dict.fromkeys(k for k in [study_key,sample_key,condition_key,biological_unit_key] if k))
    validate_input_metadata(a,sample_key=sample_key,biological_unit_key=biological_unit_key,condition_key=condition_key,identity_key=identity_key,study_key=study_key)
    if a.obs[keys+[identity_key]].isna().any().any():raise ValueError('Composition metadata contain missing values')
    n=a.obs.groupby(keys+[identity_key],observed=True).size().rename('n_cells').reset_index()
    n['fraction_of_input_sample']=n.n_cells/n.groupby(keys,observed=True).n_cells.transform('sum')
    n.attrs['scope']='Composition of input subset only; not whole tissue abundance; distinct from within-identity expression effects'
    return n

def _correlation(x,y):
    ok=np.isfinite(x)&np.isfinite(y)
    if ok.sum()<3 or np.ptp(x[ok])<1e-12 or np.ptp(y[ok])<1e-12:return np.nan
    return float(stats.spearmanr(x[ok],y[ok]).statistic)

def analyze_combinations(summary,*,adjust_activity=True,qc_covariates=None):
    meta=summary.attrs.get('reference_summary',{})
    if not meta:raise ValueError('Summary provenance required')
    unit=meta.get('biological_unit_key') or meta['sample_key'];bio=bool(meta.get('biological_unit_key'))
    keys=list(dict.fromkeys(meta['strata_keys']+([meta['condition_key']] if meta['condition_key'] else [])+['signature_version']))
    groupkeys=list(dict.fromkeys([unit]+keys))
    x=summary[summary.available_program].groupby(groupkeys+['feature_id'],observed=True).mean_score.mean().unstack('feature_id').reset_index()
    qc_covariates=list(qc_covariates or [])
    if qc_covariates:
      for cov in qc_covariates:
        if cov not in summary:raise KeyError(cov)
      qs=summary.groupby(groupkeys,observed=True)[qc_covariates].mean().reset_index();x=x.merge(qs,on=groupkeys,how='left',validate='one_to_one')
    registry=load_reference();rows=[]
    for k,z in x.groupby(keys,observed=True,dropna=False):
      k=k if isinstance(k,tuple) else (k,);base=dict(zip(keys,k));version=base['signature_version'];sigs=reference_signatures(version)
      for cname,(f1,f2) in registry['combinations'].items():
        v1=z[f1].to_numpy(float) if f1 in z else np.full(len(z),np.nan);v2=z[f2].to_numpy(float) if f2 in z else np.full(len(z),np.nan)
        rho=_correlation(v1,v2);adjusted=np.nan;qc_adjusted=np.nan
        if adjust_activity and 'A01' in z:
          activity=z.A01.to_numpy(float);ok=np.isfinite(v1)&np.isfinite(v2)&np.isfinite(activity)
          if ok.sum()>=5:
           a=np.c_[np.ones(ok.sum()),stats.rankdata(activity[ok])];xx=stats.rankdata(v1[ok]);yy=stats.rankdata(v2[ok]);rx=xx-a@np.linalg.lstsq(a,xx,rcond=None)[0];ry=yy-a@np.linalg.lstsq(a,yy,rcond=None)[0]
           if np.std(rx)>1e-12 and np.std(ry)>1e-12:adjusted=float(np.corrcoef(rx,ry)[0,1])
        if qc_covariates:
          qc=z[qc_covariates].to_numpy(float);valid=np.isfinite(v1)&np.isfinite(v2)&np.isfinite(qc).all(1)
          if valid.sum()>=len(qc_covariates)+4:
            design=np.c_[np.ones(valid.sum()),np.column_stack([stats.rankdata(qc[valid,k]) for k in range(qc.shape[1])])]
            if len(design)>np.linalg.matrix_rank(design)+2:
              xx=stats.rankdata(v1[valid]);yy=stats.rankdata(v2[valid]);rx=xx-design@np.linalg.lstsq(design,xx,rcond=None)[0];ry=yy-design@np.linalg.lstsq(design,yy,rcond=None)[0]
              if np.std(rx)>1e-12 and np.std(ry)>1e-12:qc_adjusted=float(np.corrcoef(rx,ry)[0,1])
        overlap=sorted(set(sigs[f1]['genes'])&set(sigs[f2]['genes']))
        rows.append(dict(**base,combination_id=cname,feature1=f1,feature2=f2,n_units=int((np.isfinite(v1)&np.isfinite(v2)).sum()),unit_type='Declared biological unit' if bio else 'Technical sample descriptive',spearman_rho=rho,partial_rank_rho_activity=adjusted,partial_rank_rho_technical_qc=qc_adjusted,qc_covariates=';'.join(qc_covariates),gene_overlap_n=len(overlap),gene_overlap=';'.join(overlap),p_value=np.nan,status='EXPLORATORY_ASSOCIATION' if np.isfinite(rho) else 'NOT_EVALUABLE',binary_calls_used=False,interpretation='Residence-associated co-variation; no proven independent interaction gain' if cname=='C02' else 'Continuous RNA co-variation; technical confounding possible'))
    return pd.DataFrame(rows)

def historical_combinations(calls,*,provenance):
    """Use existing frozen calls only. No new threshold calibration is performed."""
    if provenance.get('scoring_definition')!='historical_v4' or not provenance.get('threshold_hash') or not provenance.get('cell_alignment_verified'):
        raise ValueError('Frozen historical thresholds and verified cell alignment required')
    if provenance['threshold_hash']!=load_reference()['historical_threshold_hash']:
        raise ValueError('Historical threshold hash differs from the frozen v4 resource')
    if calls.index.duplicated().any():raise ValueError('Historical cell IDs must be unique')
    if 'A01' in calls and calls.A01.isin([1,2,3]).any():raise ValueError('A01 is continuous and cannot have positive binary calls')
    rows=[]
    for name,(f1,f2) in load_reference()['combinations'].items():
      if f1 not in calls or f2 not in calls:raise KeyError('Missing original combination member')
      values=calls[[f1,f2]].to_numpy(float)
      if not np.isfinite(values).all() or not np.isin(values,[-2,-1,0,1,2,3]).all():raise ValueError('Historical calls must use frozen integer call codes, not UCell scores')
      z=calls[[f1,f2]];eligible=z.ge(0).all(axis=1);z=z[eligible];a=z[f1].eq(2);b=z[f2].eq(2)
      pa=a.mean();pb=b.mean();joint=(a&b).mean();lift=joint/(pa*pb) if pa*pb>0 else np.nan
      rows.append(dict(combination_id=name,n_eligible_cells=len(z),n_joint=int((a&b).sum()),lift=lift,statistical_scope='Descriptive historical calls; cells not biological units',threshold_hash=provenance['threshold_hash']))
    return pd.DataFrame(rows)
