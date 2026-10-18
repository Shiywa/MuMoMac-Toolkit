import json,copy
from pathlib import Path
import numpy as np,pandas as pd
from mumomac_coupling.registry import ORDER,digest,canonicalize,pairs
from mumomac_coupling.engine import sample_statistics
from mumomac_coupling.statistics import four_table,spearman,rank_residual
from mumomac_coupling.states import pair_codes,PAIR_LABELS
from mumomac_apsc.calibration import select,key,independence
from mumomac_apsc.engine import applicable
from mumomac_apsc.distribution import posterior
from .config import load,inherited,merge
from .evidence import assess,codes
from .ecdf import freeze,percentile,quantile,tail_thresholds
def continuous_applicable(meta,p,c):
 # Reuse C2 lineage/modality/identity rules; continuous L features need no binary selection.
 return applicable(meta,'S01' if p=='A01' else p,{**inherited(c),'optional_lifecycle':list(ORDER[16:21])})
def fit(bundle,config=None):
 c=load(config);cal=select(bundle,inherited(c));rows=[];keys=c['calibration']['group_by']
 for vals,g in cal.groupby(keys,sort=False,observed=True,dropna=False):
  vals=vals if isinstance(vals,tuple) else (vals,);group=key(vals)
  for p in ORDER:
   ix=g.index[continuous_applicable(g,p,c)];x=bundle.scores.loc[ix,p].to_numpy(float);det=None if bundle.detected_members is None else bundle.detected_members.loc[ix,p].to_numpy(float);active=p in ORDER[:16] or p in c['optional_lifecycle']
   evidence=assess(x,cal.loc[ix],c,float(bundle.coverage.loc[p,'coverage']),bool(bundle.coverage.loc[p,'scoring_eligible']),det,active)
   row=dict(calibration_group=group,program_id=p,**evidence,ecdf=freeze(x,cal.loc[ix]),coverage=float(bundle.coverage.loc[p,'coverage']),calibration_scope=c['calibration']['scope'],calibration_interpretation='QUERY_RELATIVE' if c['calibration']['scope']=='all_query' else 'CONTROL_RELATIVE')
   rows.append(row)
 model=canonicalize(dict(version=c['version'],frozen_config=c,config_hash=digest(c),score_spec_hash=digest(bundle.spec.to_dict()),score_spec=bundle.spec.to_dict(),namespace=bundle.spec.namespace,keys=keys,calibration_cell_ids=list(cal.index),calibration_scores_hash=digest(bundle.scores.loc[cal.index].to_dict('split')),calibration_units=cal[['source_block_id','sample_id','biological_unit_id','biological_unit_confirmed','condition']].drop_duplicates().to_dict('records'),rows=rows,normal_reference=False,biological_validation=False))
 model['model_hash']=digest(model);return model
def validate(bundle,model):
 if digest({k:v for k,v in model.items() if k!='model_hash'})!=model.get('model_hash'):raise ValueError('Frozen evidence/calibration hash modified.')
 if digest(bundle.spec.to_dict())!=model['score_spec_hash']:raise ValueError('Score specification/namespace/gene universe differs from frozen calibration.')
def configure_mode(model,mode=None,overrides=None,relative=None):
 c=copy.deepcopy(model['frozen_config'])
 if mode is not None:c['default_mode']=mode
 if overrides is not None:c['overrides']=overrides
 if relative is not None:c['relative']=merge(c['relative'],relative)
 return load(c)
def call(bundle,model,mode=None,overrides=None,relative=None):
 validate(bundle,model);c=configure_mode(model,mode,overrides,relative);n=len(bundle.scores);state_codes=np.full((n,22),3,np.int8);percentiles=np.full((n,22),np.nan);records=[];lookup={(r['calibration_group'],r['program_id']):r for r in model['rows']};groups=bundle.metadata[model['keys']].astype(str).apply(lambda r:key(r),axis=1);bio=independence(bundle.metadata,model);request_hash=digest({'model_hash':model['model_hash'],'mode':c['default_mode'],'overrides':c['overrides'],'relative':c['relative']});samples=bundle.metadata.sample_id.to_numpy();animals=bundle.metadata.biological_unit_id.to_numpy();confirmed=bundle.metadata.biological_unit_confirmed.to_numpy()
 for j,p in enumerate(ORDER):
  x=bundle.scores[p].to_numpy(float);app=continuous_applicable(bundle.metadata,p,c);elig=bool(bundle.coverage.loc[p,'scoring_eligible']) and bundle.coverage.loc[p,'coverage']>=inherited(c)['qc']['min_gene_coverage'];valid=np.isfinite(x)&app&elig;det=np.ones(n,bool) if bundle.detected_members is None else bundle.detected_members[p].to_numpy()>=inherited(c)['qc']['min_detected_members']
  for group,indices in groups.groupby(groups,sort=False).indices.items():
   ii=np.asarray(indices);r=lookup.get((group,p));o=c['overrides'].get(p,{});requested=o.get('mode',c['default_mode']);active=p in ORDER[:16] or p in c['optional_lifecycle'];actual='continuous';rep='UNASSESSABLE' if not r else r['auto_representation'];evidence='CONTINUOUS_BASE_ONLY';lo=hi=None;pr=np.full(len(ii),np.nan);why='MISSING_CALIBRATION_STRATUM' if not r else ';'.join(r['decision_reasons']);binary_ok=valid[ii]&det[ii];vi=ii[binary_ok]
   if r and r['ecdf'].get('valid'):
    mapped=percentile(x[ii],r['ecdf']);mapped[~(valid[ii]&det[ii])]=np.nan;percentiles[ii,j]=mapped
   if r and active and p!='A01' and r['basic_QC_pass']:
    if requested=='auto' and rep=='DISCRETE':actual='evidence_supported_discrete';evidence=r['partition_evidence']
    elif requested=='force_gmm' and r['gmm'].get('estimable'):actual='evidence_supported_discrete' if rep=='DISCRETE' else 'forced_GMM_exploratory';evidence=r['partition_evidence'] if rep=='DISCRETE' else 'USER_FORCED_EXPLORATORY; NOT_SUPPORTED_DISCRETE'
    elif requested in ['force_tail','relative_high']:actual='relative_high';evidence='EXPLORATORY_RELATIVE_HIGH'
    if actual in ['evidence_supported_discrete','forced_GMM_exploratory']:
     gm=r['gmm'];po=posterior(x[vi],gm);pr[binary_ok]=po;state_codes[vi,j]=codes(x[vi],gm);lo=None;hi=gm.get('threshold_high')
    elif actual=='relative_high':
     fixed=o.get('fixed_threshold',c['relative']['fixed_thresholds'].get(p));q=o.get('quantile',c['relative']['quantile'])
     if fixed is not None:lo=hi=float(fixed)
     elif requested=='force_tail':lo,hi=tail_thresholds(r['ecdf'],merge(c['relative']['tail'],{**o.get('tail',{}),**{k:v for k,v in o.items() if k in c['relative']['tail']}}))
     else:lo=hi=quantile(r['ecdf'],q)
     v=x[vi];sc=np.full(len(vi),2,np.int8)
     if requested=='force_tail' and fixed is None:sc[v<lo]=0;sc[v>=hi]=1
     else:sc[v<=lo]=0;sc[v>hi]=1
     state_codes[vi,j]=sc
   if requested=='force_gmm' and actual=='continuous' and active:why+='; MODEL_NOT_ESTIMABLE_OR_BASIC_QC_FAILURE'
   for k,i in enumerate(ii):
    sc=int(state_codes[i,j]);localrep=rep if valid[i] else 'UNASSESSABLE';label=None
    if actual=='evidence_supported_discrete' and sc<3:label=['Model-supported low','Model-supported high','Uncertain'][sc]
    elif actual=='forced_GMM_exploratory' and sc<3:label=['Exploratory GMM low','Exploratory GMM high','Uncertain'][sc]
    elif actual=='relative_high' and sc<3:label=['Relative not-high','Relative high','Uncertain'][sc]
    records.append(dict(cell_id=bundle.scores.index[i],sample_id=samples[i],biological_unit_id=animals[i],biological_unit_confirmed=bool(confirmed[i]),program_id=p,raw_program_score=x[i],scoring_method=bundle.spec.scoring_method,score_namespace=bundle.spec.namespace,gene_set_version=bundle.spec.reference_version+'/'+bundle.spec.signature_version,gene_set_hash=bundle.spec.program_hashes[p],gene_coverage=float(bundle.coverage.loc[p,'coverage']),program_applicable=bool(app[i]),continuous_score_eligible=bool(valid[i]),representation_type=localrep,requested_mode=requested,actual_method=actual,discrete_or_relative_state=label,state_code=sc,posterior_high=pr[k],calibration_percentile=percentiles[i,j],threshold_low=lo,threshold_high=hi,state_definition='High p>=.9 and score>=midpoint; Low p<=.1 and score<midpoint; otherwise uncertain' if 'GMM' in actual or actual=='evidence_supported_discrete' else 'strict > quantile/fixed threshold; conservative Tail uses inherited dual MAD rule' if actual=='relative_high' else 'No binary assignment',calibration_scope=model['frozen_config']['calibration']['scope'],calibration_id=group,evidence_level=evidence,user_override=requested!='auto',biological_evidence=bio[i]+'; NO_FUNCTIONAL_VALIDATION',QC=why if valid[i] and det[i] else 'CELL_COVERAGE_APPLICABILITY_SCORE_OR_DETECTION_QC; raw score retained',method_version=c['version'],model_hash=model['model_hash'],parameter_hash=request_hash))
 return dict(table=pd.DataFrame(records),codes=state_codes,percentiles=percentiles,eligible_programs=bundle.coverage.scoring_eligible.to_numpy(bool)&(bundle.coverage.coverage.to_numpy(float)>=inherited(c)['qc']['min_gene_coverage']),config=c,parameter_hash=request_hash)
def continuous_coupling(bundle,calls,model):
 ps=pairs(model['frozen_config']['optional_lifecycle']);s=sample_statistics(bundle,calls,model['frozen_config']['coupling'],ps);m=model['frozen_config'];lookup={(r['calibration_group'],r['program_id']):r for r in model['rows']};typ=[]
 for row in s.to_dict('records'):
  group=key([row.get(k,'Unknown') for k in model['keys']]);a=lookup.get((group,row['program_A']),{});b=lookup.get((group,row['program_B']),{});typ.append((a.get('auto_representation','UNASSESSABLE'),b.get('auto_representation','UNASSESSABLE')))
 s['program_A_representation']=[a for a,b in typ];s['program_B_representation']=[b for a,b in typ];s['pair_representation']=[a+' × '+b for a,b in typ]
 # A table exists only when both program state rules were actually requested/supported.
 off=s.n_evaluable.eq(0)
 for col in ['n11','n10','n01','n00','n_evaluable','A_fraction','B_fraction','joint_fraction','lift','OR','logOR','OR_raw','logOR_raw']:s.loc[off,col]=np.nan
 s['binary_status']=np.where(off,'NOT_APPLICABLE_NO_JOINT_BINARY_RULE','DESCRIPTIVE_BINARY; evidence in per-program table');s['calibration_scope']=m['calibration']['scope'];s['model_hash']=model['model_hash'];s['parameter_hash']=calls['parameter_hash'];s['method_version']=m['version'];s['biological_evidence']='WITHIN_RECORD_DESCRIPTIVE; not cell-level biological replication';return s
def tail_profiles(bundle,model):
 validate(bundle,model);c=model['frozen_config'];lookup={(r['calibration_group'],r['program_id']):r for r in model['rows']};rows=[];ps=pairs(c['optional_lifecycle'],only=c['tail_profile']['pairs'])
 from mumomac_coupling.engine import GROUPS
 for vals,g in bundle.metadata.groupby(GROUPS,sort=False,observed=True,dropna=False):
  info=dict(zip(GROUPS,vals));group=key([g.iloc[0][k] for k in model['keys']]);ix=bundle.scores.index.get_indexer(g.index)
  for p in ps:
   a=p['program_A'];b=p['program_B'];ra=lookup.get((group,a));rb=lookup.get((group,b));xa=bundle.scores.loc[g.index,a].to_numpy(float);xb=bundle.scores.loc[g.index,b].to_numpy(float);valid=np.isfinite(xa)&np.isfinite(xb)&continuous_applicable(g,a,c)&continuous_applicable(g,b,c);ready=bool(ra and rb and ra['basic_QC_pass'] and rb['basic_QC_pass'] and ra['ecdf']['valid'] and rb['ecdf']['valid'])
   if bundle.detected_members is not None:valid &=(bundle.detected_members.loc[g.index,a].to_numpy()>=1)&(bundle.detected_members.loc[g.index,b].to_numpy()>=1)
   for q in c['tail_profile']['quantiles']:
    ta=quantile(ra['ecdf'],q) if ready else None;tb=quantile(rb['ecdf'],q) if ready else None
    if ready:
     ah=xa[valid]>ta;bh=xb[valid]>tb;ct=[int(np.sum(ah&bh)),int(np.sum(ah&~bh)),int(np.sum(~ah&bh)),int(np.sum(~ah&~bh))];f=four_table(*ct,min_cells=c['coupling']['min_cells'],min_table_cell=c['coupling']['min_table_cell'])
    else:f={k:None for k in ['n11','n10','n01','n00','n_evaluable','A_fraction','B_fraction','joint_fraction','lift','OR','logOR','OR_raw','logOR_raw','binary_estimable','estimability_status']}
    rows.append(dict(**info,**p,quantile_policy=q,threshold_A=ta,threshold_B=tb,**f,n_total=len(g),n_excluded=int(len(g)-valid.sum()),sparse_cell_warning=not bool(f['binary_estimable']),biological_unit_id=g.iloc[0].biological_unit_id,biological_unit_confirmed=bool(g.biological_unit_confirmed.all()),calibration_scope=c['calibration']['scope'],evidence_level='EXPLORATORY_TAIL_PROFILE; no threshold search or hypothesis testing',status='PASS_DESCRIPTIVE' if ready else 'UNASSESSABLE_CALIBRATION_OR_COVERAGE',threshold_rule=c['tail_profile']['high_rule'],model_hash=model['model_hash'],method_version=c['version']))
 return pd.DataFrame(rows)
def joint_features(bundle,calls,model):
 c=model['frozen_config'];ps=pairs(c['optional_lifecycle'],only=c['joint_index']['pairs']);rows=[]
 for p in ps:
  a=ORDER.index(p['program_A']);b=ORDER.index(p['program_B']);u=calls['percentiles'][:,a];v=calls['percentiles'][:,b];js=np.minimum(u,v);pc=pair_codes(calls['codes'][:,a],calls['codes'][:,b]);z=bundle.metadata[['sample_id','biological_unit_id','biological_unit_confirmed','condition','source_block_id']].copy();z.insert(0,'cell_id',z.index);z=z.reset_index(drop=True);z['pair_id']=p['pair_id'];z['program_A']=p['program_A'];z['program_B']=p['program_B'];z['raw_A']=bundle.scores.iloc[:,a].to_numpy();z['raw_B']=bundle.scores.iloc[:,b].to_numpy();z['u_A']=u;z['u_B']=v;z['J_soft']=js;z['binary_pair_state']=PAIR_LABELS[pc];z.loc[(calls['codes'][:,a]==3)|(calls['codes'][:,b]==3),'binary_pair_state']='Not applicable: no complete binary rule';z['J_soft_interpretation']='DESCRIPTIVE_RELATIVE_CO_HIGH; not activation probability or statistical coupling';z['calibration_scope']=c['calibration']['scope'];z['model_hash']=model['model_hash'];z['parameter_hash']=calls['parameter_hash'];z['method_version']=c['version'];z['coarse_identity']=bundle.metadata.coarse_identity.to_numpy();z['QC']=np.where(np.isfinite(js),'DESCRIPTIVE; see per-program distribution and provenance QC','NON_EVALUABLE_SCORE_COVERAGE_DETECTION_OR_CALIBRATION');rows.append(z)
 return pd.concat(rows,ignore_index=True)
def diagnostics(model):
 rows=[];comparisons=[]
 for r in model['rows']:
  row={k:v for k,v in r.items() if k not in ['models','gmm','validation_folds','initializations','ecdf']};row.update(gmm_reliable=r['gmm'].get('reliable'),gmm_reasons=r['gmm'].get('reasons'),gmm_separation=r['gmm'].get('separation'),gmm_kde_peak_counts=r['gmm'].get('kde_peak_counts'),gmm_weights=r['gmm'].get('weights'),model_hash=model['model_hash'],method_version=model['version']);rows.append({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()})
  for m in r['models']:
   held=[f[m['name']+'_heldout_mean_loglik'] for f in r['validation_folds'] if f.get(m['name']+'_heldout_mean_loglik') is not None];comparisons.append(dict(program_id=r['program_id'],calibration_group=r['calibration_group'],model_hash=model['model_hash'],blocked_mean_heldout_loglik=float(np.mean(held)) if held else None,n_valid_heldout_blocks=len(held),validation_unit=r['unit_level'],**{k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in m.items()}))
 return pd.DataFrame(rows),pd.DataFrame(comparisons)
def run(bundle,config=None,output=None,model=None,mode=None,overrides=None,relative=None,include_joint=None,include_tail=None):
 if model is not None and config is not None and digest(load(config))!=model['config_hash']:raise ValueError('Frozen calibration configuration differs; request only an explicitly audited mode/relative override or refit separately.')
 model=model or fit(bundle,config);calls=call(bundle,model,mode,overrides,relative);continuous=continuous_coupling(bundle,calls,model)
 if include_joint is None:include_joint=model['frozen_config']['joint_index']['enabled'] and model['frozen_config']['joint_index']['default_output']
 if include_tail is None:include_tail=model['frozen_config']['tail_profile']['enabled']
 tail=tail_profiles(bundle,model) if include_tail else pd.DataFrame(columns=['sample_id','pair_id','quantile_policy','status'])
 selected=pairs(model['frozen_config']['optional_lifecycle'],only=model['frozen_config']['joint_index']['pairs']);has_binary=any(np.any((calls['codes'][:,ORDER.index(p['program_A'])]<3)&(calls['codes'][:,ORDER.index(p['program_B'])]<3)) for p in selected)
 joint=joint_features(bundle,calls,model) if include_joint or has_binary else pd.DataFrame(columns=['cell_id','sample_id','pair_id','J_soft','method_version'])
 if has_binary and not include_joint:
  joint[['u_A','u_B','J_soft']]=np.nan;joint['J_soft_interpretation']='NOT_REQUESTED; binary pair states only'
 dec,cmp=diagnostics(model)
 if output:
  out=Path(output);out.mkdir(parents=True,exist_ok=True)
  if any(out.iterdir()):raise FileExistsError('Use a new empty output directory; historical results must not be overwritten.')
  for name,frame in [('cell_program_representation',calls['table']),('selected_pair_joint_features',joint),('sample_coupling_summary',continuous),('continuous_coupling_results',continuous),('tail_coenrichment_profiles',tail)]:frame.to_parquet(out/(name+'.parquet'),index=False)
  dec.to_csv(out/'program_representation_decisions.tsv',sep='\t',index=False,na_rep='NA');cmp.to_csv(out/'distribution_model_comparison.tsv',sep='\t',index=False,na_rep='NA');(out/'frozen_representation.json').write_text(json.dumps(model,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
  (out/'representation_request.json').write_text(json.dumps(canonicalize(dict(requested_config=calls['config'],parameter_hash=calls['parameter_hash'],model_hash=model['model_hash'],J_soft_requested=bool(include_joint),tail_profile_requested=bool(include_tail))),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
 return dict(model=model,calls=calls,continuous=continuous,tail=tail,joint=joint,decisions=dec,comparison=cmp)
