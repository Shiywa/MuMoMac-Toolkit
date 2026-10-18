"""Biological units, not cells or replicate libraries, determine inference."""
import numpy as np,pandas as pd
from scipy import stats
from .statistics import welch_effect,adjust_bh
from .adapters import known
STRATA=['study_id','source_block_id','tissue','coarse_identity','species','modality']
METRICS=['A_fraction','B_fraction','joint_fraction','logOR','spearman','adjusted_association']
def biological_summary(samples):
 rows=[]
 for keys,g in samples.groupby(STRATA+['condition','pair_id'],observed=True,dropna=False,sort=False):
  r=dict(zip(STRATA+['condition','pair_id'],keys));valid=g[g.continuous_estimable];bio=g[g.biological_unit_confirmed]
  r.update(evaluable_sample_count=int(len(valid)),sample_record_count=int(g.sample_id[g.sample_id.map(known)].nunique()),supported_biological_unit_count=int(bio.biological_unit_id.nunique()),confirmed_biological_replicate_count=int(bio.biological_unit_id.nunique()) if len(bio)==len(g) and len(g) else None,study_count=int(g.study_id[g.study_id.map(known)].nunique()),source_count=int(g.source_block_id[g.source_block_id.map(known)].nunique()),sample_effects=g[['sample_id','lift','logOR','spearman','adjusted_association','estimability_status']].to_json(orient='records'),direction_agreement=float(max((valid.spearman>0).mean(),(valid.spearman<0).mean())) if len(valid) else np.nan,effect_heterogeneity=float(valid.spearman.std(ddof=1)) if len(valid)>1 else np.nan,heterogeneity_definition='SD of record-level Spearman, not between-study tau²',uncertainty='Unknown animal totals; no biological CI' if len(bio)<len(g) else 'Declared biological IDs; no automatic independent-source certification',evidence_level='DESCRIPTIVE_RECORD_LEVEL' if len(bio)<len(g) else 'BIOLOGICAL_UNIT_DESCRIPTIVE',limitations='Condition/context, gene overlap, technical dependence and source selection remain; no molecular synergy claim')
  rows.append(r)
 return pd.DataFrame(rows)

def condition_effects(samples,contrasts,n_bootstrap=2000,seed=20261016):
 out=[];rng=np.random.default_rng(seed)
 for contrast in contrasts:
  control=contrast['reference'];case=contrast['case'];paired=bool(contrast.get('paired',False))
  if control==case:raise ValueError('Conditions must differ.')
  for keys,g in samples.groupby(STRATA+['pair_id'],observed=True,dropna=False,sort=False):
   base=dict(zip(STRATA+['pair_id'],keys));g=g[g.condition.isin([control,case])];a=g[g.condition==case];b=g[g.condition==control]
   design='OK'
   if g.empty or a.empty or b.empty:design='NOT_ESTIMABLE_MISSING_CONDITION_WITHIN_SOURCE'
   elif not g.biological_unit_confirmed.all():design='NOT_ESTIMABLE_UNVERIFIED_BIOLOGICAL_UNITS'
   elif not known(base['study_id']) or not known(base['coarse_identity']):design='NOT_ESTIMABLE_MISSING_SOURCE_OR_IDENTITY'
   elif g.batch.eq('Unknown').any():design='NOT_ESTIMABLE_BATCH_UNVERIFIED'
   elif set(a.batch).isdisjoint(set(b.batch)):design='NOT_ESTIMABLE_COMPLETE_BATCH_CONFOUNDING'
   elif g.batch.nunique()>1 or g.batch.eq('MIXED_WITHIN_RECORD').any():design='NOT_ESTIMABLE_MULTIBATCH_REQUIRES_STRATIFIED_MODEL'
   elif not paired and set(a.biological_unit_id)&set(b.biological_unit_id):design='NOT_ESTIMABLE_SHARED_ANIMALS_REQUIRES_PAIRED'
   for metric in METRICS:
    z=g.copy()
    if metric in ['logOR']:z=z[z.binary_estimable]
    elif metric in ['spearman','adjusted_association']:z=z[z.continuous_estimable]
    # Equal library estimates within an actual animal; duplicate records do not increase n.
    descriptive_a=z.loc[z.condition==case,metric].mean();descriptive_b=z.loc[z.condition==control,metric].mean()
    z=z[z.biological_unit_confirmed & z.biological_unit_id.map(known)]
    units=z.groupby(['condition','biological_unit_id'],observed=True)[metric].mean().unstack('condition') if len(z) else pd.DataFrame()
    av=units[case].dropna().to_numpy() if case in units else np.array([]);bv=units[control].dropna().to_numpy() if control in units else np.array([])
    row=dict(**base,reference=control,case=case,paired=paired,metric=metric,statistical_unit='Declared animal/donor within one source; equal library estimate',descriptive_record_difference=descriptive_a-descriptive_b,n_case_confirmed_units_available=int(a.loc[a.biological_unit_confirmed,'biological_unit_id'].nunique()),n_control_confirmed_units_available=int(b.loc[b.biological_unit_confirmed,'biological_unit_id'].nunique()),p_adjusted=np.nan,bootstrap_ci_low=np.nan,bootstrap_ci_high=np.nan)
    result=welch_effect(av,bv)
    if paired and case in units and control in units:
     delta=(units[case]-units[control]).dropna().to_numpy();n=len(delta);result.update(n_case_units=n,n_control_units=n,effect=float(delta.mean()) if n else np.nan,status='NOT_ESTIMABLE_FEWER_THAN3_COMPLETE_PAIRS',p_value=np.nan,ci_low=np.nan,ci_high=np.nan)
     if n>=3 and delta.std(ddof=1)>1e-12:
      se=delta.std(ddof=1)/np.sqrt(n);crit=stats.t.ppf(.975,n-1);result.update(ci_low=delta.mean()-crit*se,ci_high=delta.mean()+crit*se,p_value=float(2*stats.t.sf(abs(delta.mean()/se),n-1)),df=n-1,status='BIOLOGICAL_UNIT_PAIRED_T; symmetry_and_small_n_limits_apply')
    if design!='OK':result.update(status=design,p_value=np.nan,ci_low=np.nan,ci_high=np.nan)
    elif metric in ['A_fraction','B_fraction','joint_fraction','logOR'] and ('threshold_inference_status' not in g or not g.threshold_inference_status.eq('INDEPENDENT_FROZEN_RULE').all()):result.update(status='NOT_ESTIMABLE_THRESHOLD_FIT_OVERLAPS_EVALUATION_OR_UNVERIFIED',p_value=np.nan,ci_low=np.nan,ci_high=np.nan)
    elif np.isfinite(result['p_value']) and n_bootstrap>0:
     if paired:boot=delta[rng.integers(0,len(delta),(n_bootstrap,len(delta)))].mean(1)
     else:boot=av[rng.integers(0,len(av),(n_bootstrap,len(av)))].mean(1)-bv[rng.integers(0,len(bv),(n_bootstrap,len(bv)))].mean(1)
     row['bootstrap_ci_low'],row['bootstrap_ci_high']=np.quantile(boot,[.025,.975])
    if not g.empty and not g.biological_unit_confirmed.all():result.update(n_case_units=None,n_control_units=None,effect=np.nan)
    row.update(result);out.append(row)
 result=pd.DataFrame(out)
 if len(result):
  for _,g in result.groupby(['reference','case','metric'],observed=True,dropna=False):result.loc[g.index,'p_adjusted']=adjust_bh(g.p_value)
  result['multiple_testing']='BH across all scheduled pairs and strata per contrast/endpoint; inestimable p=1 internally, q=NA externally'
 else:result=pd.DataFrame(columns=STRATA+['pair_id','reference','case','metric','effect','p_value','p_adjusted','status'])
 return result
