"""Transparent current-input technical diagnosis, never silent signature selection."""
import hashlib,json
import numpy as np,pandas as pd
from scipy import sparse,stats
from .validation import prepare_input
from .reference import REFERENCE_ORDER
from .reference_scoring import score_signature_versions
from .reliability import input_fingerprint,symbol_fingerprint
from .analysis import _correlation

ADVISOR_POLICY={'version':'1.0.0','retention':[.75,.5,.25],'seeds':[101,203,307],
                'min_biological_units':3,'min_rank_samples':3,'minimum_own_sd':1e-8,
                'practical_standardized_difference':.02,'bootstrap':2000,
                'unit_resampling_seed':20261014,'depth_strata':'input within-library UMI tertiles; no response labels',
                'detection_strata':'Original present-member detection fraction tertiles; no response labels',
                'no_auto_switch':True}

def _unit_mapping(obs,sample_key,biological_unit_key,condition_key,study_key):
    if sample_key not in obs:raise KeyError(sample_key)
    keys=[sample_key]+([study_key] if study_key and study_key!=sample_key else [])
    for k in [biological_unit_key,condition_key]:
      if k:
       if k not in obs:raise KeyError(k)
       if obs[k].isna().any():raise ValueError(f'Missing actual {k} values')
       if obs.groupby(keys,observed=True,dropna=False)[k].nunique(dropna=False).gt(1).any():raise ValueError(f'Technical sample has multiple {k} assignments')
    if obs[sample_key].isna().any():raise ValueError('Missing sample identifiers')
    unit=obs[biological_unit_key].astype(str) if biological_unit_key else obs[sample_key].astype(str)
    if study_key:unit=obs[study_key].astype(str)+'::'+unit
    return unit.to_numpy()

def _tertiles(values):
    v=np.asarray(values,float);q=np.quantile(v,[1/3,2/3]);labels=np.where(v<=q[0],'lower',np.where(v<=q[1],'intermediate','higher'))
    return labels,q.tolist()

def _delta_ci(delta,rng,nboot):
    a=np.asarray(delta,float);a=a[np.isfinite(a)]
    if len(a)<3:return np.nan,np.nan,np.nan
    boots=a[rng.integers(0,len(a),(nboot,len(a)))].mean(1)
    lo,hi=np.quantile(boots,[.025,.975]);return float(a.mean()),float(lo),float(hi)

def _effect(scores,obs,units,condition_key,reference,case):
    if not condition_key or not reference or not case:return np.nan,0,0
    d=pd.DataFrame({'unit':units,'condition':obs[condition_key].to_numpy(),'score':scores})
    if d.groupby('unit').condition.nunique().gt(1).any():return np.nan,0,0
    d=d.groupby(['unit','condition'],observed=True).score.mean().reset_index()
    a=d[d.condition.eq(case)].score.dropna();b=d[d.condition.eq(reference)].score.dropna()
    return float(a.mean()-b.mean()) if len(a) and len(b) else np.nan,len(a),len(b)

def advise_signature_version(scored,*,mode='coverage',sample_key=None,biological_unit_key=None,condition_key=None,study_key=None,identity_key=None,reference=None,case=None,umi_counts_verified=False,chunk_size=500):
    if mode not in ('coverage','umi'):raise ValueError('Advisor mode must be coverage or umi')
    if 'mumomac_reference' not in scored.uns:raise ValueError('Reference scores are required')
    m=scored.uns['mumomac_reference'];config=m['configuration'];versions=list(config['versions'])
    if set(versions)!={'original','expanded'}:
        rows=[dict(feature_id=pid,recommendation='NOT_EVALUABLE',measurement_eligibility='MISSING_VERSION_BASELINE',technical_stability_evidence='NOT_TESTED',phenotype_response_evidence='NOT_TESTED',biological_validation_status='NOT_ESTABLISHED_BY_THIS_INPUT',reason='Both Original and Expanded complete baseline scores are required; no missing baseline is invented',missing_evidence='Missing signature-version baseline',automatic_version_switch=False,reference_scope='Current scored versions only') for pid in REFERENCE_ORDER]
        return {'metrics':m['gene_coverage'].copy(),'robustness':pd.DataFrame(),'interaction':pd.DataFrame(),'recommendations':pd.DataFrame(rows),'policy':dict(ADVISOR_POLICY),'test_status':'NOT_TESTED','limitations':'Missing version baseline'}
    if input_fingerprint(scored,config['input_layer'])!=m['input_fingerprint'] or symbol_fingerprint(scored,config['gene_symbol_key'] or None)!=m['symbol_fingerprint']:
        raise ValueError('Scored input changed; recompute both-version baseline before direct testing')
    for v in versions:
      if hashlib.sha256(np.asarray(scored.obsm[f'X_mumomac_reference_{v}']).tobytes()).hexdigest()!=m['score_hashes'][v]:raise ValueError('Baseline score changed; cannot certify direct test')
    coverage=m['gene_coverage'].copy();base={v:np.asarray(scored.obsm[f'X_mumomac_reference_{v}'],dtype=float).copy() for v in versions}
    eligibility={(r.feature_id,r.signature_version):bool(r.scoring_eligible) for r in coverage.itertuples()}
    covered=coverage[['feature_id','signature_version','n_reference_genes','n_mapped_genes','n_present_genes','coverage','mapping_retention','median_cell_detected_members','fraction_cells_zero_members','scoring_eligible','developmental_candidate','default_scope']].copy()
    recommendations=[]
    def empty(reason):
      for pid in REFERENCE_ORDER:
        ok=all(eligibility.get((pid,v),False) for v in versions)
        recommendations.append(dict(feature_id=pid,recommendation='INCONCLUSIVE' if mode=='coverage' and ok else 'NOT_EVALUABLE',measurement_eligibility='ELIGIBLE_BOTH' if ok else 'INCOMPLETE_COVERAGE_OR_MAPPING',technical_stability_evidence='NOT_TESTED',phenotype_response_evidence='NOT_TESTED',biological_validation_status='NOT_ESTABLISHED_BY_THIS_INPUT',reason=reason,missing_evidence='Trusted UMI assay, actual unit IDs and paired perturbation' if mode=='coverage' else reason,automatic_version_switch=False,reference_scope='Full Expanded L/A are development candidates'))
      return {'metrics':covered,'robustness':pd.DataFrame(),'interaction':pd.DataFrame(),'recommendations':pd.DataFrame(recommendations),'policy':dict(ADVISOR_POLICY),'test_status':'NOT_TESTED','limitations':reason}
    if mode=='coverage':return empty('Coverage measured; no perturbation performed. No technical or biological winner.')
    if config['expression_type']!='counts' or not umi_counts_verified:return empty('Trusted raw UMI counts not explicitly declared; integer values alone do not establish UMI assay')
    if sample_key is None:return empty('Technical sample metadata required for paired sample stability')
    obs=scored.obs;units=_unit_mapping(obs,sample_key,biological_unit_key,condition_key,study_key)
    sample=obs[sample_key].astype(str).to_numpy()
    if study_key:sample=obs[study_key].astype(str).to_numpy()+'::'+sample
    p=prepare_input(scored,layer=config['input_layer'],expression_type='counts',gene_symbol_key=config['gene_symbol_key'] or None,duplicate_policy=config['duplicate_policy'])
    x=scored.X if config['input_layer']=='X' else scored.layers[config['input_layer']]
    x=sparse.csr_matrix(x)
    depth=np.asarray(p.matrix.sum(1)).ravel();depth_labels=np.empty(len(depth),object);boundaries={}
    for s in pd.unique(sample):
        mask=sample==s;labels,q=_tertiles(depth[mask]);depth_labels[mask]=labels;boundaries[s]=q
    detected={v:np.asarray(scored.obsm[f'X_mumomac_detected_{v}']) for v in versions}
    rawrows=[];unitrows=[];stratumrows=[];baselines={};rng_ci=np.random.default_rng(ADVISOR_POLICY['unit_resampling_seed'])
    for v in versions:
     for j,pid in enumerate(REFERENCE_ORDER):
      values=base[v][:,j].astype(float);finite=values[np.isfinite(values)];sd=float(np.std(finite,ddof=1)) if len(finite)>1 else np.nan
      sm=pd.DataFrame({'sample':sample,'unit':units,'score':values}).groupby(['sample','unit']).score.mean().reset_index();sample_sd=sm.score.std(ddof=1)
      effect,nc,nr=_effect(values,obs,units,condition_key,reference,case)
      baselines[(pid,v)]=(sd,sample_sd,sm,effect,nc,nr)
    # Identical molecule draws for the two versions: one thinned matrix per retention/seed.
    for retention in ADVISOR_POLICY['retention']:
     for seed in ADVISOR_POLICY['seeds']:
      rng=np.random.default_rng(seed);thin=x.copy();thin.data=rng.binomial(x.data.astype(np.int64),retention);thin.eliminate_zeros()
      query=scored.copy()
      for key in list(query.obsm):
       if key.startswith('X_mumomac_reference_') or key.startswith('X_mumomac_detected_'):del query.obsm[key]
      del query.uns['mumomac_reference']
      if config['input_layer']=='X':query.X=thin
      else:query.layers[config['input_layer']]=thin
      result=score_signature_versions(query,species=config['species'],signature_version='both',layer=config['input_layer'],expression_type='counts',gene_symbol_key=config['gene_symbol_key'] or None,duplicate_policy=config['duplicate_policy'],coverage_policy=config['coverage_policy'],max_rank=config['max_rank'],chunk_size=chunk_size,copy=False)
      for v in versions:
       new=np.asarray(result.obsm[f'X_mumomac_reference_{v}'],dtype=float)
       for j,pid in enumerate(REFERENCE_ORDER):
        before=base[v][:,j];after=new[:,j];ok=np.isfinite(before)&np.isfinite(after);sd,sample_sd,sm,effect,nc,nr=baselines[(pid,v)]
        diff=np.abs(after-before);mae=float(np.mean(diff[ok])) if ok.any() else np.nan
        normalized=mae/sd if sd>ADVISOR_POLICY['minimum_own_sd'] else np.nan
        ns=pd.DataFrame({'sample':sample,'unit':units,'score':after}).groupby(['sample','unit']).score.mean().reset_index()
        paired=sm.merge(ns,on=['sample','unit'],suffixes=('_base','_thin'));drift=np.abs(paired.score_thin-paired.score_base)
        rank=_correlation(paired.score_base.to_numpy(),paired.score_thin.to_numpy()) if len(paired)>=3 else np.nan
        thin_effect,_,_=_effect(after,obs,units,condition_key,reference,case)
        rawrows.append(dict(feature_id=pid,signature_version=v,retention=retention,technical_seed=seed,n_cells=int(ok.sum()),cell_mae=mae,own_baseline_sd=sd,cell_standardized_mae=normalized,cell_spearman=_correlation(before,after),sample_mean_drift=float(drift.mean()) if len(drift) else np.nan,sample_standardized_drift=float(drift.mean()/sample_sd) if sample_sd>1e-8 else np.nan,sample_rank_spearman=rank,n_samples=len(paired),n_biological_units=len(pd.unique(units)) if biological_unit_key else 0,baseline_condition_effect=effect,thinned_condition_effect=thin_effect,effect_direction_retained=bool(effect*thin_effect>0) if np.isfinite(effect) and abs(effect)>1e-12 else None,n_case_units=nc if biological_unit_key else 0,n_reference_units=nr if biological_unit_key else 0,non_evaluable_fraction=float(1-ok.mean()),status='TESTED' if ok.any() else 'NOT_EVALUABLE',score_scale='Relative to own fixed baseline SD; raw score sizes not comparable'))
        for unit in pd.unique(units):
          ix=units==unit;valid=ix&ok;sp=paired[paired.unit.eq(unit)];cellstd=float(np.mean(diff[valid])/sd) if valid.any() and sd>1e-8 else np.nan
          samplestd=float(np.abs(sp.score_thin-sp.score_base).mean()/sample_sd) if len(sp) and sample_sd>1e-8 else np.nan
          unitrows.append(dict(feature_id=pid,signature_version=v,retention=retention,technical_seed=seed,unit=unit,cell_standardized_mae=cellstd,sample_standardized_drift=samplestd,unit_is_biological=bool(biological_unit_key)))
        original_detection=detected['original'][:,j]/max(1,int(coverage[(coverage.feature_id==pid)&coverage.signature_version.eq('original')].n_present_genes.iloc[0]));detect_labels,detect_q=_tertiles(original_detection)
        dimensions={'UMI':depth_labels,'Original_detection':detect_labels}
        if identity_key and identity_key in obs:dimensions['identity']=obs[identity_key].astype(str).to_numpy()
        for dimension,labels in dimensions.items():
         for label in pd.unique(labels):
          for unit in pd.unique(units):
           ix=(labels==label)&(units==unit)&ok
           if not ix.any():continue
           value=float(np.mean(diff[ix])/sd) if sd>1e-8 else np.nan
           stratumrows.append(dict(feature_id=pid,signature_version=v,retention=retention,technical_seed=seed,dimension=dimension,stratum=label,unit=unit,n_cells=int(ix.sum()),cell_standardized_mae=value,mean_detected_members=float(detected[v][ix,j].mean()),unit_is_biological=bool(biological_unit_key),detection_boundaries=json.dumps(detect_q) if dimension=='Original_detection' else ''))
    robustness=pd.DataFrame(rawrows);u=pd.DataFrame(unitrows);strata=pd.DataFrame(stratumrows);interactions=[]
    # Average technical seeds within actual unit before bootstrap or recommendations.
    metrics=['cell_standardized_mae','sample_standardized_drift'];paired_results=[]
    for pid in REFERENCE_ORDER:
      for retention in ADVISOR_POLICY['retention']:
        for metric in metrics:
          z=u[(u.feature_id==pid)&u.retention.eq(retention)].groupby(['unit','signature_version'])[metric].mean().unstack('signature_version')
          delta=(z['original']-z['expanded']).dropna() if {'original','expanded'}<=set(z.columns) else pd.Series(dtype=float)
          mean,lo,hi=_delta_ci(delta,rng_ci,ADVISOR_POLICY['bootstrap']) if biological_unit_key else (float(delta.mean()) if len(delta) else np.nan,np.nan,np.nan)
          paired_results.append(dict(feature_id=pid,retention=retention,metric=metric,n_independent_units=len(delta) if biological_unit_key else 0,delta_original_minus_expanded=mean,ci_low=lo,ci_high=hi,bootstrap_unit='Declared biological unit; technical seeds averaged first' if biological_unit_key else 'No inferential CI; technical libraries only'))
    paired_metrics=pd.DataFrame(paired_results)
    for (pid,retention,dimension,label),z in strata.groupby(['feature_id','retention','dimension','stratum']):
      d=z.groupby(['unit','signature_version']).cell_standardized_mae.mean().unstack('signature_version')
      delta=(d['original']-d['expanded']).dropna() if {'original','expanded'}<=set(d.columns) else pd.Series(dtype=float)
      mean,lo,hi=_delta_ci(delta,rng_ci,2000) if biological_unit_key else (float(delta.mean()) if len(delta) else np.nan,np.nan,np.nan)
      interactions.append(dict(feature_id=pid,retention=retention,dimension=dimension,stratum=label,n_independent_units=len(delta) if biological_unit_key else 0,delta_original_minus_expanded=mean,ci_low=lo,ci_high=hi,scope='Version-by-input-depth/detection contrast; no cutoff fitted; technical stability only'))
    # Actual version × depth/detection contrast: high-minus-low version difference
    # within each unit, after averaging technical seeds. Not a fitted UMI cutoff.
    for (pid,retention,dimension),z in strata[strata.dimension.isin(['UMI','Original_detection'])].groupby(['feature_id','retention','dimension']):
      d=z.groupby(['unit','stratum','signature_version']).cell_standardized_mae.mean().unstack('signature_version')
      if not {'original','expanded'}<=set(d.columns):continue
      differences=(d['original']-d['expanded']).unstack('stratum')
      if not {'lower','higher'}<=set(differences.columns):continue
      contrast=(differences['higher']-differences['lower']).dropna()
      mean,lo,hi=_delta_ci(contrast,rng_ci,2000) if biological_unit_key else (float(contrast.mean()) if len(contrast) else np.nan,np.nan,np.nan)
      interactions.append(dict(feature_id=pid,retention=retention,dimension=dimension,stratum='higher_minus_lower_interaction',n_independent_units=len(contrast) if biological_unit_key else 0,delta_original_minus_expanded=mean,ci_low=lo,ci_high=hi,scope='Paired within-unit difference-in-differences: (Original-Expanded)higher minus (Original-Expanded)lower; own-scale technical error, no biological superiority claim'))
    interaction=pd.DataFrame(interactions)
    tolerance=ADVISOR_POLICY['practical_standardized_difference']
    for pid in REFERENCE_ORDER:
      ok=all(eligibility.get((pid,v),False) for v in versions);z=paired_metrics[paired_metrics.feature_id.eq(pid)];finite=z.delta_original_minus_expanded.notna();reason=''
      if not ok:state='NOT_EVALUABLE';reason='One or both signatures fail frozen coverage/mapping gates'
      elif not finite.all():state='NOT_EVALUABLE';reason='Baseline variance near zero or no valid paired cells/sample scores'
      elif not biological_unit_key or z.n_independent_units.min()<3:state='INCONCLUSIVE';reason='Technical perturbation measured, but fewer than 3 declared biological units; no reliable paired uncertainty'
      else:
        preferred=np.where(z.ci_low>tolerance,1,np.where(z.ci_high < -tolerance,-1,0))
        context=interaction[(interaction.feature_id==pid)&interaction.ci_low.notna() & ~interaction.stratum.eq('higher_minus_lower_interaction')]
        signs=set(np.where(context.ci_low>tolerance,1,np.where(context.ci_high < -tolerance,-1,0)))-{0}
        if (1 in preferred and -1 in preferred) or len(signs)>1:state='CONTEXT_DEPENDENT';reason='Depth/detection/identity or cell/sample metrics give opposing paired directions'
        elif np.all(preferred==1):state='EXPANDED_TECHNICALLY_PREFERRED';reason='Paired biological-unit CIs favor Expanded across all tested depths and both metrics'
        elif np.all(preferred==-1):state='ORIGINAL_TECHNICALLY_PREFERRED';reason='Paired biological-unit CIs favor Original across all tested depths and both metrics'
        elif np.all((z.ci_low>=-tolerance)&(z.ci_high<=tolerance)):state='COMPARABLE';reason='All paired CIs fit the frozen practical-equivalence tolerance'
        else:state='INCONCLUSIVE';reason='Measured differences do not consistently resolve paired uncertainty across all depths/metrics'
      phenotype=bool(condition_key and biological_unit_key and reference and case)
      recommendations.append(dict(feature_id=pid,recommendation=state,measurement_eligibility='ELIGIBLE_BOTH' if ok else 'INCOMPLETE_COVERAGE_OR_MAPPING',technical_stability_evidence='TESTED_CURRENT_INPUT',phenotype_response_evidence='DESCRIPTIVE_RNA_RESPONSE_ONLY' if phenotype else 'NOT_TESTED',biological_validation_status='NOT_ESTABLISHED_BY_THINNING_OR_CONDITION_LABEL',reason=reason,missing_evidence='Independent functional readouts and transfer replication',automatic_version_switch=False,reference_scope='Original default; Expanded L/A developmental even when technically preferred'))
    return {'metrics':covered,'robustness':robustness,'unit_metrics':u,'paired_metrics':paired_metrics,'stratum_metrics':strata,'interaction':interaction,'recommendations':pd.DataFrame(recommendations),'policy':dict(ADVISOR_POLICY),'test_status':'TESTED_CURRENT_INPUT','depth_boundaries':boundaries,'input_fingerprint':m['input_fingerprint'],'limitations':'Molecule retention relative to observed full-depth baseline is technical evidence, not biological truth or new independent validation'}
