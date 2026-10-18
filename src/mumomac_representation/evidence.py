"""Working-density evidence plus source/animal-blocked stability; no modality p-value."""
import copy,warnings
import numpy as np,pandas as pd
from scipy.stats import norm,skewnorm
from sklearn.mixture import GaussianMixture
from mumomac_apsc.distribution import fit_gmm,basic_qc,balance,posterior
from mumomac_coupling.registry import canonicalize
from .config import inherited
def logpdf(x,model):
 x=np.asarray(x,float)
 if model['name']=='Gaussian':return norm.logpdf(x,model['mean'],model['sd'])
 if model['name']=='SkewNormal':return skewnorm.logpdf(x,*model['parameters'])
 logs=np.column_stack([np.log(model['weights'][i])+norm.logpdf(x,model['means'][i],model['sds'][i]) for i in range(2)])
 return np.logaddexp(logs[:,0],logs[:,1])
def models(x,meta,c2,gm=None):
 z=balance(x,meta,c2['gmm']['max_fit_cells']);gm=gm if gm is not None else fit_gmm(x,meta,c2,stability=False);out=[]
 mu=float(np.mean(z));sd=float(np.sqrt(np.var(z)+c2['gmm']['reg_covar']));ga=dict(name='Gaussian',mean=mu,sd=sd,k=2,valid=True,outside_mass=float(norm.cdf(0,mu,sd)+norm.sf(1,mu,sd)));out.append(ga)
 try:
  with warnings.catch_warnings(record=True) as caught:params=skewnorm.fit(z)
  valid=bool(np.isfinite(params).all() and params[2]>1e-8)
  sk=dict(name='SkewNormal',parameters=list(params),k=3,valid=valid,fit_warnings=[str(w.message) for w in caught],outside_mass=float(skewnorm.cdf(0,*params)+skewnorm.sf(1,*params)));out.append(sk)
 except (ValueError,FloatingPointError,RuntimeError) as e:out.append(dict(name='SkewNormal',valid=False,k=3,error=str(e)))
 if gm.get('estimable'):
  g=dict(name='GMM2',k=5,valid=True,means=gm['means'],sds=gm['sds'],weights=gm['weights']);g['outside_mass']=float(sum(w*(norm.cdf(0,m,s)+norm.sf(1,m,s)) for m,s,w in zip(g['means'],g['sds'],g['weights'])));out.append(g)
 else:out.append(dict(name='GMM2',valid=False,k=5,error='MODEL_NOT_ESTIMABLE'))
 for model in out:
  if model['valid']:
   ll=float(logpdf(z,model).sum());model.update(loglik=ll,n_balanced=len(z),aic=2*model['k']-2*ll,bic=np.log(len(z))*model['k']-2*ll)
 return out,gm
def codes(x,gm,confidence=.9):
 p=posterior(x,gm);out=np.full(len(x),2,np.int8);mid=gm['midpoint'];out[(p<=1-confidence)&(x<mid)]=0;out[(p>=confidence)&(x>=mid)]=1;return out
def blocks(meta):
 if meta.biological_unit_confirmed.all():return (meta.source_block_id+'::'+meta.biological_unit_id).to_numpy(),'CONFIRMED_BIOLOGICAL_UNIT'
 if meta.source_block_id.nunique()>1:return meta.source_block_id.to_numpy(),'TECHNICAL_SOURCE_BLOCK; ANIMAL_RELATIONSHIPS_UNVERIFIED'
 return (meta.source_block_id+'::RECORD::'+meta.sample_id).to_numpy(),'RECORD_ONLY; ANIMAL_RELATIONSHIPS_UNVERIFIED'
def assess(x,meta,cfg,coverage=1.,eligible=True,detected=None,active=True):
 x=np.asarray(x,float);ok=np.isfinite(x);meta=meta.iloc[np.flatnonzero(ok)].copy();detected=None if detected is None else np.asarray(detected)[ok];x=x[ok];c2=inherited(cfg);d=cfg['distribution'];qc=basic_qc(x,coverage,eligible,detected,c2)
 result=dict(**qc,auto_representation='CONTINUOUS_PREFERRED',decision_reasons=[],gmm={},models=[],validation_folds=[],initializations=[],partition_evidence='NOT_ESTABLISHED',unit_level='NOT_ASSESSED',n_calibration_blocks=0)
 hard={'GENE_COVERAGE_FAILED','NEAR_CONSTANT','MEMBER_DETECTION_FAILED'}
 if not len(x) or set(qc['basic_QC_reasons'])&hard:
  result.update(auto_representation='UNASSESSABLE',decision_reasons=qc['basic_QC_reasons'] or ['MISSING_SCORE']);return canonicalize(result)
 if not active:result['decision_reasons']=['A01_CONTINUOUS_COVARIATE_OR_L_PARTITION_NOT_SELECTED'];return canonicalize(result)
 if not qc['basic_QC_pass']:
  result.update(auto_representation='INCONCLUSIVE',decision_reasons=qc['basic_QC_reasons']);return canonicalize(result)
 mm,gm=models(x,meta,c2);result.update(gmm=gm,models=mm)
 valid=[m for m in mm if m['valid'] and m['name']!='GMM2'];unimodal=min(valid,key=lambda m:m['bic']);result['best_unimodal']=unimodal['name'];result['bic_gain_over_best_unimodal']=unimodal['bic']-gm.get('bic2',np.inf)
 bb,level=blocks(meta);unique=sorted(set(bb));result.update(unit_level=level,n_calibration_blocks=len(unique));folds=[]
 # Every model is compared on the same held-out block; no disease labels enter selection.
 for block in unique[:d['max_validation_folds']]:
  train=bb!=block;test=~train
  if min(train.sum(),test.sum())<d['min_validation_cells'] or train.sum()<c2['qc']['min_cells']:folds.append(dict(block=block,status='INSUFFICIENT_BLOCK_CELLS'));continue
  tt,gg=models(x[train],meta.iloc[np.flatnonzero(train)],c2);vv=[m for m in tt if m['valid'] and m['name']!='GMM2'];uu=min(vv,key=lambda m:m['bic']);mix=next(m for m in tt if m['name']=='GMM2');row=dict(block=block,status='PASS' if mix['valid'] else 'GMM_NOT_ESTIMABLE',n_train=int(train.sum()),n_test=int(test.sum()),best_training_unimodal=uu['name'],refit_gmm_reliable=gg.get('reliable',False),evaluation_unit='HELD_OUT_CALIBRATION_BLOCK; not disease holdout tuning')
  for candidate in tt:
   row[candidate['name']+'_heldout_mean_loglik']=float(np.mean(logpdf(x[test],candidate))) if candidate['valid'] else None
  if mix['valid']:
   row['gmm_minus_best_unimodal_loglik']=float(np.mean(logpdf(x[test],mix)-logpdf(x[test],uu)))
   row['fixed_cells_label_agreement']=float(np.mean(codes(x,gg)==codes(x,gm))) if gm.get('estimable') else None
  folds.append(row)
 result['validation_folds']=folds;ff=[f for f in folds if f['status']=='PASS'];result['n_valid_folds']=len(ff)
 gains=[f['gmm_minus_best_unimodal_loglik'] for f in ff];result['mean_block_loglik_gain']=float(np.mean(gains)) if gains else None;result['positive_block_fraction']=float(np.mean(np.asarray(gains)>0)) if gains else None;result['LOO_shape_fraction']=float(np.mean([f['refit_gmm_reliable'] for f in ff])) if ff else None;result['LOO_min_label_agreement']=min([f['fixed_cells_label_agreement'] for f in ff if f.get('fixed_cells_label_agreement') is not None],default=None)
 # Seed sensitivity is distinct from sampling biological units.
 if gm.get('estimable'):
  z=balance(x,meta,c2['gmm']['max_fit_cells'])
  for seed in d['initialization_seeds']:
   with warnings.catch_warnings():
    warnings.simplefilter('ignore');g=GaussianMixture(2,n_init=1,reg_covar=c2['gmm']['reg_covar'],random_state=seed).fit(z[:,None])
   order=np.argsort(g.means_.ravel());mu=g.means_.ravel()[order];sd=np.sqrt(g.covariances_.reshape(-1))[order];gg=dict(means=mu.tolist(),sds=sd.tolist(),weights=g.weights_[order].tolist(),midpoint=float(mu.mean()));result['initializations'].append(dict(seed=seed,converged=bool(g.converged_),label_agreement=float(np.mean(codes(x,gg)==codes(x,gm)))))
  result['initialization_min_agreement']=min(r['label_agreement'] for r in result['initializations'])
 reasons=[];confounded=any(any(token in r for token in ['CONFOUNDING','INSTABILITY','UNSTABLE','SMALL_OR_DOMINANT']) for r in gm.get('reasons',[]))
 if confounded:result.update(auto_representation='INCONCLUSIVE',decision_reasons=gm['reasons'])
 elif not gm.get('estimable'):result.update(auto_representation='INCONCLUSIVE',decision_reasons=['NUMERIC_GMM_NOT_ESTIMABLE'])
 elif not gm['reliable']:result.update(auto_representation='CONTINUOUS_PREFERRED',decision_reasons=['INSUFFICIENT_STABLE_PARTITION_EVIDENCE; not proof of unimodality']+gm['reasons'])
 else:
  if result['bic_gain_over_best_unimodal']<d['min_bic_gain_over_best_unimodal']:reasons.append('REASONABLE_UNIMODAL_MODEL_NOT_CLEARLY_INFERIOR')
  if len(unique)<d['min_calibration_blocks'] or len(ff)<d['min_calibration_blocks']:reasons.append('INSUFFICIENT_CALIBRATION_BLOCKS')
  if not gains or result['mean_block_loglik_gain']<d['min_mean_block_loglik_gain'] or result['positive_block_fraction']<d['min_positive_block_fraction']:reasons.append('NO_CONSISTENT_BLOCK_HELDOUT_GAIN')
  if result['LOO_shape_fraction'] is None or result['LOO_shape_fraction']<d['min_LOO_shape_fraction']:reasons.append('LEAVE_BLOCK_OUT_SHAPE_UNSTABLE')
  if result['LOO_min_label_agreement'] is None or result['LOO_min_label_agreement']<d['min_LOO_label_agreement']:reasons.append('LEAVE_BLOCK_OUT_LABELS_UNSTABLE')
  if result['initialization_min_agreement']<d['min_initialization_agreement']:reasons.append('INITIALIZATION_UNSTABLE')
  if any(m.get('outside_mass',1)>d['max_outside_support_mass'] for m in mm if m['name'] in [unimodal['name'],'GMM2']):reasons.append('WORKING_DENSITY_SUPPORT_MISMATCH')
  result['decision_reasons']=reasons or ['SEPARATION_SHAPE_UNIMODAL_COMPARISON_AND_BLOCK_STABILITY_SUPPORTED']
  result['auto_representation']='DISCRETE' if not reasons else 'CONTINUOUS_PREFERRED' if reasons==['REASONABLE_UNIMODAL_MODEL_NOT_CLEARLY_INFERIOR'] else 'INCONCLUSIVE'
  if not reasons:result['partition_evidence']='MODEL_SUPPORTED; '+level+'; NOT_FUNCTIONALLY_VALIDATED'
 return canonicalize(result)
