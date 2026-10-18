"""State-only extension. All coupling mathematics is imported from frozen C1."""
import json
from pathlib import Path
import numpy as np,pandas as pd
from mumomac_coupling.registry import ORDER,digest,canonicalize,pairs
from mumomac_coupling.engine import sample_statistics
from mumomac_coupling.states import pair_codes,PAIR_LABELS
from .config import load_config,for_program,DEFAULT
from .calibration import select,key,independence
from .distribution import basic_qc,fit_gmm,fit_tail,posterior
def applicable(meta,p,c):
 ok=meta.species.isin(['mouse','human']).to_numpy() & meta.modality.str.lower().isin(['scrna','scrna-seq','snrna','snrna-seq']).to_numpy()
 # Current scope is macrophage/monocyte single-cell scores; spots, other lineages fail binary applicability.
 ok &= meta.coarse_identity.str.contains('macrophage|monocyte|kupffer|microglia|osteoclast',case=False,regex=True).to_numpy()
 if p.startswith('L'):ok &= ~meta.coarse_identity.str.contains('monocyte',case=False).to_numpy() & (p in c['optional_lifecycle'])
 if p=='A01':ok[:]=False
 return ok
def fit(bundle,config=None):
 c=load_config(config);meta=select(bundle,c);keys=c['calibration']['group_by'];rows=[]
 for vals,g in meta.groupby(keys,sort=False,observed=True,dropna=False):
  vals=vals if isinstance(vals,tuple) else (vals,);group=key(vals)
  for p in ORDER:
   mode,tc,reason=for_program(c,p);active=p in ORDER[:16] or p in c['optional_lifecycle'];ix=g.index[applicable(g,p,c)] if active else g.index[:0]
   x=bundle.scores.loc[ix,p].to_numpy(float);finite=np.isfinite(x);ix=ix[finite];x=x[finite];gm={};tail={}
   det=None if bundle.detected_members is None else bundle.detected_members.loc[ix,p].to_numpy(float)
   qc=basic_qc(x,float(bundle.coverage.loc[p,'coverage']),bool(bundle.coverage.loc[p,'scoring_eligible']),det,c)
   r=dict(calibration_group=group,program_id=p,mode=mode,user_override=mode!='auto',override_reason=reason,**qc,auto_recommended_method='continuous_only',actual_method='continuous_only',evidence_level='UNASSESSABLE' if active and not qc['basic_QC_pass'] else 'CONTINUOUS_ONLY',QC_status=';'.join(qc['basic_QC_reasons']) if active else 'A01_CONTINUOUS_COVARIATE' if p=='A01' else 'LIFECYCLE_NOT_SELECTED',gmm={},tail={},calibration_scope=c['calibration']['scope'],calibration_id=group)
   if active and qc['basic_QC_pass']:
    use=meta.loc[ix];gm=fit_gmm(x,use,c);tail=fit_tail(x,use,tc)
    auto='gmm' if gm['reliable'] else 'tail' if c['auto']['allow_tail_fallback'] else 'continuous_only'
    actual={'auto':auto,'force_gmm':'gmm','force_tail':'tail','continuous_only':'continuous_only'}[mode]
    quality='PASS' if gm['reliable'] else ';'.join(gm['reasons'])
    if actual=='gmm' and not gm['estimable']:actual='continuous_only';quality='MODEL_NOT_ESTIMABLE;'+quality
    evidence='MODEL_SUPPORTED_HIGH' if actual=='gmm' and gm['reliable'] else 'USER_FORCED_EXPLORATORY' if actual=='gmm' else 'RELATIVE_HIGH' if actual=='tail' else 'CONTINUOUS_ONLY'
    if mode=='force_tail' and tc!=DEFAULT['tail']:evidence='USER_FORCED_EXPLORATORY'
    # Fitting success and replication evidence are distinct; sample support cannot establish biological validity.
    if actual=='gmm' and (not use.biological_unit_confirmed.all() or use[['source_block_id','biological_unit_id']].drop_duplicates().shape[0]<2):
     evidence='USER_FORCED_EXPLORATORY' if mode=='force_gmm' else 'MODEL_SUPPORTED_HIGH_TECHNICAL_ONLY'
    if det is None:quality+=';NOT_DONE_SCORE_ONLY'
    r.update(auto_recommended_method=auto,actual_method=actual,evidence_level=evidence,QC_status=quality,gmm=gm,tail=tail)
   rows.append(r)
 model=dict(version=c['version'],config=c,config_hash=digest(c),score_spec_hash=digest(bundle.spec.to_dict()),namespace=bundle.spec.namespace,keys=keys,calibration_cell_ids=list(meta.index),calibration_score_hash=digest(bundle.scores.loc[meta.index].to_dict('split')),calibration_units=meta[['source_block_id','sample_id','biological_unit_id','biological_unit_confirmed','condition']].drop_duplicates().to_dict('records'),rows=rows,absolute_functional_activation=False)
 model=canonicalize(model);model['parameter_hash']=digest(model);return model
def call(bundle,model):
 if model['parameter_hash']!=digest({k:v for k,v in model.items() if k!='parameter_hash'}):raise ValueError('Frozen calibration/config hash changed.')
 if model['score_spec_hash']!=digest(bundle.spec.to_dict()):raise ValueError('Scoring/gene namespace differs from frozen calibration.')
 c=model['config'];n=len(bundle.scores);codes=np.full((n,22),3,np.int8);rows=[];lookup={(r['calibration_group'],r['program_id']):r for r in model['rows']}
 groups=bundle.metadata[model['keys']].astype(str).apply(lambda r:key(r),axis=1)
 independent=independence(bundle.metadata,model)
 for j,p in enumerate(ORDER):
  x=bundle.scores[p].to_numpy(float);eligible=bool(bundle.coverage.loc[p,'scoring_eligible']) and float(bundle.coverage.loc[p,'coverage'])>=c['qc']['min_gene_coverage'];ok=np.isfinite(x)&applicable(bundle.metadata,p,c)&eligible
  if bundle.detected_members is not None:ok &= bundle.detected_members[p].to_numpy()>=c['qc']['min_detected_members']
  for group,indices in groups.groupby(groups,sort=False).indices.items():
   ii=np.asarray(indices);r=lookup.get((group,p));pr=np.full(len(ii),np.nan);lo=hi=None;method='continuous_only';evidence='UNASSESSABLE';qc='NO_CALIBRATION_STRATUM';auto='continuous_only';override=False
   if r:
    method=r['actual_method'];evidence=r['evidence_level'];qc=r['QC_status'];auto=r['auto_recommended_method'];override=r['user_override']
    valid=ii[ok[ii]];xx=x[valid];out=np.full(len(valid),2,np.int8)
    if method=='tail':
     lo=r['tail']['threshold_low'];hi=r['tail']['threshold_high'];out[xx<lo]=0;out[xx>=hi]=1;codes[valid,j]=out
    elif method=='gmm':
     g=r['gmm'];po=posterior(xx,g);pr[np.isin(ii,valid)]=po;lo=g['threshold_low'];hi=g['threshold_high'];out[(po<c['gmm']['posterior_borderline'])|(xx<g['midpoint'])]=0;out[(po>=c['gmm']['posterior_high'])&(xx>=g['midpoint'])]=1;codes[valid,j]=out
   for k,i in enumerate(ii):
    code=int(codes[i,j]);state=['Not-high','High','Uncertain','Unassessable'][code]
    if code==1:state='Model-supported High' if evidence.startswith('MODEL_SUPPORTED') else 'User-forced exploratory High' if evidence=='USER_FORCED_EXPLORATORY' else 'Relative High'
    if method=='continuous_only' and evidence=='CONTINUOUS_ONLY':state='Continuous-only'
    rows.append(dict(cell_id=bundle.scores.index[i],program_id=p,score=x[i],state=state,state_code=code,state_call_method=method,auto_recommended_method=auto,actual_method=method,user_override=override,posterior_high=pr[k],threshold_low=lo,threshold_high=hi,calibration_scope=c['calibration']['scope'],calibration_id=group,evidence_level=evidence if code in [0,1] else 'UNCERTAIN' if code==2 else 'CONTINUOUS_ONLY' if state=='Continuous-only' else 'UNASSESSABLE',methodological_evidence=evidence,biological_evidence=independent[i]+'; BIOLOGICAL_VALIDATION_NOT_ESTABLISHED',program_coverage=float(bundle.coverage.loc[p,'coverage']),QC_status=qc if ok[i] or method=='continuous_only' else 'CELL_SCORE_COVERAGE_APPLICABILITY_OR_DETECTION_FAILED',parameter_hash=model['parameter_hash']))
 # C1 continuous eligibility only reflects coverage. Mode/shape never masks valid scores.
 eligible=(bundle.coverage.scoring_eligible.to_numpy(bool)&(bundle.coverage.coverage.to_numpy(float)>=c['qc']['min_gene_coverage']))
 return dict(codes=codes,eligible_programs=eligible,table=pd.DataFrame(rows))
def run(bundle,config=None,output=None,model=None):
 if model is not None and config is not None and digest(load_config(config))!=model['config_hash']:raise ValueError('Requested config/override differs from the frozen calibration; refitting must be explicit and separately audited.')
 model=model or fit(bundle,config);calls=call(bundle,model);c=model['config'];ps=pairs(c['optional_lifecycle']);stats=sample_statistics(bundle,calls,c['coupling'],ps)
 # C1 computes the same underlying counts. Public continuous-only metrics are explicitly not applicable.
 inactive={p for p in ORDER if all(r['actual_method']=='continuous_only' for r in model['rows'] if r['program_id']==p)}
 off=stats.program_A.isin(inactive)|stats.program_B.isin(inactive)|stats.n_evaluable.eq(0)
 metrics=['n11','n10','n01','n00','n_evaluable','joint_fraction','A_fraction','B_fraction','lift','OR','logOR','OR_raw','logOR_raw']
 for col in metrics:
  if col in stats:stats.loc[off,col]=np.nan
 stats['binary_output_status']=np.where(off,'NOT_APPLICABLE_CONTINUOUS_OR_UNASSESSABLE','DESCRIPTIVE_RELATIVE_STATE')
 stats.loc[off,'binary_estimable']=False
 stats['threshold_inference_status']='DESCRIPTIVE_APSC; no independent biological validation claimed'
 diag=[]
 for r in model['rows']:
  row={k:v for k,v in r.items() if k not in ['gmm','tail']};row.update({'gmm_'+k:v for k,v in r['gmm'].items()});row.update({'tail_'+k:v for k,v in r['tail'].items()});row.update(parameter_hash=model['parameter_hash'],config_hash=model['config_hash'],score_namespace=model['namespace'],auto_decision_reason='DISTRIBUTION_MODEL_PASSED_FROZEN_RULES' if r['auto_recommended_method']=='gmm' else r['QC_status'],calibration_interpretation='EXPLORATORY_RELATIVE' if r['calibration_scope']=='all_query' else 'CONTROL_RELATIVE; independent biology not established');diag.append({k:json.dumps(canonicalize(v),ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in row.items()})
 diag=pd.DataFrame(diag);audit=diag[diag.user_override].copy()
 if output:
  out=Path(output);out.mkdir(parents=True,exist_ok=True)
  if any(out.iterdir()):raise FileExistsError('APSC output must be a new directory.')
  calls['table'].to_parquet(out/'cell_program_states.parquet',index=False);stats.to_parquet(out/'sample_coupling_statistics.parquet',index=False);diag.to_csv(out/'program_distribution_diagnostics.tsv',sep='\t',index=False,na_rep='NA');audit.to_csv(out/'manual_override_audit.tsv',sep='\t',index=False,na_rep='NA');(out/'frozen_calibration.json').write_text(json.dumps(model,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
  selected=[p for p in ps if p['pair_id'] in ['C01','C02','C03','C04']];comb=[]
  for p in selected:
   a=ORDER.index(p['program_A']);b=ORDER.index(p['program_B']);z=pd.DataFrame({'cell_id':bundle.scores.index,'pair_id':p['pair_id'],'pair_state':PAIR_LABELS[pair_codes(calls['codes'][:,a],calls['codes'][:,b])]})
   methods=calls['table'].pivot(index='cell_id',columns='program_id',values='actual_method').loc[bundle.scores.index]
   continuous=methods[p['program_A']].eq('continuous_only')|methods[p['program_B']].eq('continuous_only')
   z.loc[continuous.to_numpy(),'pair_state']='Not applicable (continuous-only)';comb.append(z)
  pd.concat(comb,ignore_index=True).to_parquet(out/'cell_combination_states.parquet',index=False)
 return dict(model=model,calls=calls,diagnostics=diag,audit=audit,samples=stats)
