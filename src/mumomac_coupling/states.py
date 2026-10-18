"""Explicit, frozen relative-high calling. No fitted normal reference is available."""
import json
import numpy as np,pandas as pd
from .registry import ORDER,digest
from .adapters import known
NOT_HIGH,HIGH,UNCERTAIN,UNASSESSABLE=0,1,2,3
LABELS=np.array(['Not-high','High','Uncertain','Unassessable'])
KEYS=['coarse_identity','tissue','species','modality']
def group_key(row,keys):return json.dumps([str(row[k]) for k in keys],ensure_ascii=False)
def weighted_quantile(values,weights,q):
 o=np.argsort(values,kind='stable');v=np.asarray(values)[o];w=np.asarray(weights)[o];cdf=np.cumsum(w)/np.sum(w)
 return float(v[min(np.searchsorted(cdf,q,side='left'),len(v)-1)])
def fit_thresholds(bundle,config):
 mode=config.get('threshold_mode','none')
 if mode=='reference':raise ValueError('Reference-calibrated Current UCell thresholds have not been independently calibrated.')
 if mode not in ['none','control','exploratory']:raise ValueError('Unknown threshold mode.')
 keys=config.get('threshold_group_keys',KEYS)
 if any(k in keys for k in ['condition','sample_id','biological_unit_id','study_id','source_block_id']):raise ValueError('Group/source-specific percentile fitting is prohibited in this prototype; use one frozen rule across conditions and sources.')
 if set(keys)!=set(KEYS):raise ValueError('Freeze identity/tissue/species/modality calibration strata explicitly.')
 q=float(config.get('quantile',.8));band=float(config.get('quantile_band',.02))
 if not 0<q-band<=q+band<1:raise ValueError('Quantile and uncertainty band must lie inside(0,1).')
 fitmeta=bundle.metadata.copy()
 if mode=='control':
  control=config.get('control_condition')
  if not known(control):raise ValueError('Control calibration requires an explicit control condition.')
  fitmeta=fitmeta[fitmeta.condition==control]
 if config.get('calibration_sample_ids') is not None:fitmeta=fitmeta[fitmeta.sample_id.isin(config['calibration_sample_ids'])]
 if mode!='none' and fitmeta.empty:raise ValueError('No declared calibration observations.')
 rows=[];active=set(ORDER[:16])|set(config.get('optional_lifecycle',[]))
 for vals,g in fitmeta.groupby(keys,observed=True,dropna=False,sort=False):
  ids=g.index;block=g.source_block_id.where(g.source_block_id.map(known),g.study_id);block=block.where(block.map(known),'UNVERIFIED_SOURCE')
  # Confirmed animals combine technical records; otherwise retain explicit library weighting.
  unit=np.where(g.biological_unit_confirmed,g.biological_unit_id,'RECORD::'+g.sample_id)
  tmp=pd.DataFrame({'block':block,'unit':unit},index=ids);nb=tmp.block.nunique();nu=tmp.groupby('block').unit.transform('nunique');nc=tmp.groupby(['block','unit']).unit.transform('size');weights=1/(nb*nu*nc)
  key=json.dumps(list(vals),ensure_ascii=False)
  for pid in ORDER:
   x=bundle.scores.loc[ids,pid].to_numpy(float);w=weights.to_numpy();eligible=bool(bundle.coverage.loc[pid,'scoring_eligible']) and bundle.coverage.loc[pid,'coverage']>=config.get('min_gene_coverage',.8);valid=np.isfinite(x)&eligible
   if bundle.detected_members is not None:valid &= bundle.detected_members.loc[ids,pid].to_numpy()>=config.get('min_detected_members',0)
   elif config.get('min_detected_members',0)>0:valid &= False
   row=dict(calibration_group=key,program_id=pid,mode=mode,quantile=q,quantile_band=band,threshold=np.nan,lower=np.nan,upper=np.nan,status='NO_VALID_THRESHOLD',n_calibration_cells=int(valid.sum()),sample_ids=sorted(set(g.sample_id)),source_blocks=sorted(set(block)),confirmed_animals_only=bool(g.biological_unit_confirmed.all()),weighting='source equal; confirmed animal or explicit record equal; cell equal within unit',interpretation='CONTROL_RELATIVE_HIGH' if mode=='control' else 'EXPLORATORY_RELATIVE_HIGH' if mode=='exploratory' else 'CONTINUOUS_ONLY')
   if mode!='none' and pid in active and valid.sum()>=3 and np.ptp(x[valid])>1e-12:
    row.update(threshold=weighted_quantile(x[valid],w[valid],q),lower=weighted_quantile(x[valid],w[valid],q-band),upper=weighted_quantile(x[valid],w[valid],q+band),status='FROZEN_RELATIVE_THRESHOLD')
   if pid=='A01':row['status']='CONTINUOUS_COVARIATE_NOT_BINARIZED'
   rows.append(row)
 model=dict(mode=mode,namespace=bundle.spec.namespace,species=bundle.spec.species,signature_version=bundle.spec.signature_version,program_hashes=bundle.spec.program_hashes,score_spec_hash=digest(bundle.spec.to_dict()),keys=keys,calibration_input_fingerprint=bundle.fingerprint() if mode=='exploratory' else digest({'cells':list(fitmeta.index),'scores':bundle.scores.loc[fitmeta.index].to_json()}),calibration_sample_unit_mapping=fitmeta[['sample_id','study_id','source_block_id','biological_unit_id','biological_unit_confirmed']].drop_duplicates().to_dict('records'),rows=rows,absolute_activation=False)
 model['threshold_definition_hash']=digest(model);return model

def call_states(bundle,model,config):
 if digest({k:v for k,v in model.items() if k!='threshold_definition_hash'})!=model.get('threshold_definition_hash'):raise ValueError('Frozen threshold definition was modified or its hash is invalid.')
 if model['namespace']!=bundle.spec.namespace or model['program_hashes']!=bundle.spec.program_hashes or model['score_spec_hash']!=digest(bundle.spec.to_dict()):raise ValueError('Frozen threshold/score definitions or parameters mismatch.')
 if 'threshold_mode' in config and config['threshold_mode']!=model['mode']:raise ValueError('Requested threshold mode conflicts with the frozen model.')
 for key in ['quantile','quantile_band']:
  if key in config and model['rows'] and config[key]!=model['rows'][0][key]:raise ValueError('Requested threshold parameter conflicts with the frozen model: '+key)
 n=len(bundle.scores);state=np.full((n,22),UNASSESSABLE,np.int8);threshold=np.full((n,22),np.nan,np.float32);reason=np.full((n,22),'COVERAGE_OR_SCORE_UNASSESSABLE',object)
 active=set(ORDER[:16])|set(config.get('optional_lifecycle',[]));look={(r['calibration_group'],r['program_id']):r for r in model['rows']}
 eligibility=(bundle.coverage.coverage.to_numpy(float)>=config.get('min_gene_coverage',.8))&bundle.coverage.scoring_eligible.to_numpy(bool)
 groups=bundle.metadata.apply(lambda r:group_key(r,model['keys']),axis=1)
 for j,pid in enumerate(ORDER):
  if pid not in active:reason[:,j]='CONTINUOUS_COVARIATE' if pid=='A01' else 'LIFECYCLE_NOT_SELECTED_FOR_BINARY_CALL';continue
  x=bundle.scores[pid].to_numpy();ok=np.isfinite(x)&eligibility[j]
  if pid.startswith('L') and not config.get('allow_lifecycle_in_monocytes',False):ok &= ~bundle.metadata.coarse_identity.str.contains('monocyte',case=False).to_numpy()
  if bundle.detected_members is not None:ok &= bundle.detected_members[pid].to_numpy()>=config.get('min_detected_members',0)
  elif config.get('min_detected_members',0)>0:ok &= False;reason[:,j]='MISSING_REQUESTED_DETECTION_QC'
  for key,indices in groups.groupby(groups,sort=False).indices.items():
   ix=np.asarray(indices);ix=ix[ok[ix]];r=look.get((key,pid));state[ix,j]=UNCERTAIN;reason[ix,j]='NO_VALID_THRESHOLD_CONTINUOUS_ONLY'
   if not r or r['status']!='FROZEN_RELATIVE_THRESHOLD':continue
   threshold[ix,j]=r['threshold'];reason[ix,j]='NEAR_FROZEN_THRESHOLD';state[ix[x[ix]<r['lower']],j]=NOT_HIGH;state[ix[x[ix]>r['upper']],j]=HIGH
   reason[ix[state[ix,j]==NOT_HIGH],j]='BELOW_LOWER_THRESHOLD';reason[ix[state[ix,j]==HIGH],j]='ABOVE_UPPER_THRESHOLD'
 return dict(codes=state,thresholds=threshold,reasons=reason,eligible_programs=eligibility)

def pair_codes(a,b):
 out=np.where((a==HIGH)&(b==HIGH),0,np.where((a==HIGH)&(b==NOT_HIGH),1,np.where((a==NOT_HIGH)&(b==HIGH),2,3))).astype(np.int8)
 out[(a==UNCERTAIN)|(b==UNCERTAIN)]=4;out[(a==UNASSESSABLE)|(b==UNASSESSABLE)]=5
 return out
PAIR_LABELS=np.array(['Joint','A-only','B-only','Neither','Uncertain','Unassessable'])
