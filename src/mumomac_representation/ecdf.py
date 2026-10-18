import numpy as np
from mumomac_coupling.states import weighted_quantile
from mumomac_apsc.distribution import weights
def freeze(x,meta):
 x=np.asarray(x,float);ok=np.isfinite(x);x=x[ok];meta=meta.iloc[np.flatnonzero(ok)]
 if len(x)<3 or np.ptp(x)<1e-12:return {'valid':False,'reason':'NO_CALIBRATION_VARIATION','n':len(x)}
 w=weights(meta);w=w/sum(w);order=np.argsort(x,kind='stable');v=x[order];sorted_weights=w[order];unique,start=np.unique(v,return_index=True);mass=np.add.reduceat(sorted_weights,start);cumulative=np.cumsum(mass)
 return dict(valid=True,n=len(x),values=unique.tolist(),mass=mass.tolist(),cumulative=cumulative.tolist(),median=weighted_quantile(x,w,.5),scaled_MAD=1.4826*weighted_quantile(abs(x-weighted_quantile(x,w,.5)),w,.5),weighting='source equal; confirmed animal or explicit record equal; cell equal within unit')
def percentile(x,model):
 x=np.asarray(x,float);out=np.full(x.shape,np.nan)
 if not model.get('valid'):return out
 v=np.asarray(model['values']);cum=np.asarray(model['cumulative']);mass=np.asarray(model['mass']);ok=np.isfinite(x);xx=x[ok];left=np.searchsorted(v,xx,side='left');right=np.searchsorted(v,xx,side='right');before=np.where(left>0,cum[np.maximum(left-1,0)],0.);same=right>left
 result=before.copy();result[same]+=.5*mass[left[same]];out[ok]=np.clip(result,0,1);return out
def quantile(model,q):
 if not 0<q<1:raise ValueError('Quantile must be inside (0,1).')
 if not model.get('valid'):return None
 return float(model['values'][min(int(np.searchsorted(model['cumulative'],q,side='left')),len(model['values'])-1)])
def tail_thresholds(model,parameters):
 a=max(quantile(model,parameters['borderline_quantile']),model['median']+parameters['borderline_mad_multiplier']*model['scaled_MAD']);b=max(quantile(model,parameters['high_quantile']),model['median']+parameters['high_mad_multiplier']*model['scaled_MAD']);return a,b
