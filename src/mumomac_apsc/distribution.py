"""Shape diagnostics describe a score distribution, not biological activation."""
import warnings
import numpy as np,pandas as pd
from scipy.stats import norm,spearmanr
from scipy.signal import find_peaks
from sklearn.mixture import GaussianMixture
from mumomac_coupling.states import weighted_quantile
def weights(meta):
 block=meta.source_block_id.astype(str);unit=np.where(meta.biological_unit_confirmed,meta.biological_unit_id,'RECORD::'+meta.sample_id)
 t=pd.DataFrame({'block':block,'unit':unit},index=meta.index)
 return (1/(t.block.nunique()*t.groupby('block').unit.transform('nunique')*t.groupby(['block','unit']).unit.transform('size'))).to_numpy(float)
def balance(x,m,max_n):
 w=weights(m);n=min(len(x),max_n)
 order=np.argsort(x,kind='stable');cdf=np.cumsum(w[order])/sum(w)
 return x[order[np.minimum(np.searchsorted(cdf,(np.arange(n)+.5)/n),len(x)-1)]]
def basic_qc(x,coverage,eligible,detected,cfg):
 q=cfg['qc'];v=x[np.isfinite(x)];reasons=[]
 if not eligible or not np.isfinite(coverage) or coverage<q['min_gene_coverage']:reasons.append('GENE_COVERAGE_FAILED')
 if len(v)<q['min_cells']:reasons.append('INSUFFICIENT_FINITE_CELLS')
 u,ct=np.unique(v,return_counts=True);tied=float(ct.max()/len(v)) if len(v) else 1.
 boundary=float(np.mean((v<=1e-8)|(v>=1-1e-8))) if len(v) else 1.
 if len(v) and np.std(v)<q['min_sd']:reasons.append('NEAR_CONSTANT')
 if len(u)<q['min_unique']:reasons.append('INSUFFICIENT_UNIQUE_RANK_SCORES')
 if tied>q['max_tied_mass']:reasons.append('TIES_DEGENERATE')
 if boundary>q['max_boundary_mass']:reasons.append('BOUNDARY_PILEUP')
 detfrac=None
 if detected is not None:
  detfrac=float(np.mean(detected>=q['min_detected_members'])) if len(detected) else None
  if detfrac is not None and detfrac<q['min_detection_fraction']:reasons.append('MEMBER_DETECTION_FAILED')
 return dict(n_finite=len(v),n_unique=len(u),max_tied_mass=tied,boundary_mass=boundary,score_min=float(v.min()) if len(v) else None,score_max=float(v.max()) if len(v) else None,score_sd=float(np.std(v)) if len(v) else None,detection_QC='NOT_DONE_SCORE_ONLY' if detected is None else 'PER_CELL_MEMBERS_CHECKED',detection_fraction=detfrac,basic_QC_pass=not reasons,basic_QC_reasons=reasons)
def posterior(x,m):
 x=np.asarray(x,float);mu=np.asarray(m['means']);sd=np.asarray(m['sds']);w=np.asarray(m['weights'])
 logs=np.column_stack([np.log(w[k])+norm.logpdf(x,mu[k],sd[k]) for k in range(2)])
 return np.exp(logs[:,1]-np.logaddexp(logs[:,0],logs[:,1]))
def threshold(m,p):
 # First admissible interval; posterior and midpoint are both checked at calling.
 grid=np.linspace(m['midpoint'],1.,10001);pr=posterior(grid,m);ix=np.flatnonzero(pr>=p)
 return float(grid[ix[0]]) if len(ix) else None
def eta2(x,labels):
 labels=np.asarray(labels).astype(str);total=np.sum((x-x.mean())**2)
 if len(np.unique(labels))<2 or total<=0:return 0.
 return float(sum(np.sum(labels==k)*(x[labels==k].mean()-x.mean())**2 for k in np.unique(labels))/total)
def fit_gmm(x,meta,cfg,stability=True):
 c=cfg['gmm'];z=balance(x,meta,c['max_fit_cells']);out={'estimable':False,'reliable':False,'reasons':[],'n_fit':len(z),'weighting':'source equal / confirmed unit or record equal; inverse-CDF deterministic balanced fit'}
 try:
  with warnings.catch_warnings(record=True) as caught:
   one=GaussianMixture(1,n_init=c['n_init'],reg_covar=c['reg_covar'],random_state=c['random_state']).fit(z[:,None])
   two=GaussianMixture(2,n_init=c['n_init'],reg_covar=c['reg_covar'],random_state=c['random_state']).fit(z[:,None])
  order=np.argsort(two.means_.ravel());mu=two.means_.ravel()[order];sd=np.sqrt(two.covariances_.reshape(-1))[order];w=two.weights_[order]
  out.update(means=mu.tolist(),sds=sd.tolist(),weights=w.tolist(),midpoint=float(mu.mean()),separation=float(np.diff(mu)[0]/np.sqrt(np.mean(sd**2))),converged=bool(one.converged_ and two.converged_),aic1=one.aic(z[:,None]),aic2=two.aic(z[:,None]),bic1=one.bic(z[:,None]),bic2=two.bic(z[:,None]),fit_warnings=[str(w.message) for w in caught])
  out['bic_gain']=out['bic1']-out['bic2']
  if not out['converged'] or np.min(w)*len(z)<2 or np.diff(mu)[0]<1e-5:out['reasons'].append('MODEL_NOT_ESTIMABLE');return out
  out['estimable']=True
  if out['separation']<c['min_separation']:out['reasons'].append('LOW_SEPARATION')
  if not c['min_high_weight']<=w[1]<=c['max_high_weight'] or np.min(w)*len(z)<c['min_component_cells']:out['reasons'].append('SMALL_OR_DOMINANT_COMPONENT')
  if out['bic_gain']<c['min_bic_gain']:out['reasons'].append('INSUFFICIENT_BIC_GAIN_AUXILIARY')
  # Fixed Gaussian KDE on histogram grid; Scott bandwidth and three frozen factors.
  grid=np.linspace(max(0.,z.min()-np.std(z)*.5),min(1.,z.max()+np.std(z)*.5),256)
  hist,edges=np.histogram(z,bins=256,range=(grid[0],grid[-1]),density=True);cent=(edges[1:]+edges[:-1])/2
  bw=max(np.std(z,ddof=1)*len(z)**(-.2),1e-5);peaks=[];valleys=[];modes=[]
  for factor in c['kde_bandwidth_factors']:
   density=norm.pdf((grid[:,None]-cent[None,:])/(bw*factor))@hist*(edges[1]-edges[0])/(bw*factor)
   pp,_=find_peaks(density,prominence=density.max()*c['kde_min_prominence'])
   left=[i for i in pp if grid[i]<out['midpoint']];right=[i for i in pp if grid[i]>=out['midpoint']]
   valid=False;vr=1.
   if left and right:
    a=max(left,key=lambda i:density[i]);b=max(right,key=lambda i:density[i]);vr=float(density[a:b+1].min()/min(density[a],density[b]));valid=vr<=c['kde_max_valley_ratio']
   peaks.append(bool(valid));valleys.append(vr);modes.append(len(pp))
  out.update(kde_peak_counts=modes,kde_valley_ratios=valleys,kde_bimodal_recurrence=sum(peaks),kde_scott_bandwidth=bw)
  if sum(peaks)<c['kde_min_peak_recurrence']:out['reasons'].append('NO_BANDWIDTH_STABLE_TWO_PEAKS')
  conf={}
  for k in ['study_id','source_block_id','coarse_identity','platform','batch']:
   if k in meta:conf[k]=eta2(x,meta[k].to_numpy())
  for k,v in conf.items():
   if cfg['auto']['check_technical_confounding'] and v>c['max_group_explained_variance']:out['reasons'].append('GROUP_CONFOUNDING:'+k)
  out['group_explained_variance']=conf;depth={}
  for k in ['nUMI','nGenes']:
   if k in meta:
    v=pd.to_numeric(meta[k],errors='coerce').to_numpy();ok=np.isfinite(v)&np.isfinite(x)
    depth[k]=float(spearmanr(x[ok],v[ok]).statistic) if ok.sum()>=30 and np.std(v[ok])>0 else None
    if cfg['auto']['check_technical_confounding'] and depth[k] is not None and abs(depth[k])>c['max_depth_rho']:out['reasons'].append('DEPTH_CONFOUNDING:'+k)
  out['depth_rho']=depth
  prob=posterior(x,out);group_fr=[]
  for _,ix in meta.groupby(['source_block_id','sample_id'],sort=False).indices.items():
   if len(ix)>=c['min_sample_cells']:group_fr.append(float(np.mean(prob[np.asarray(ix)]>=.5)))
  mixed=sum(.05<v<.95 for v in group_fr);exclusive=float(np.mean([v<=.05 or v>=.95 for v in group_fr])) if group_fr else None
  out.update(sample_high_fractions=group_fr,n_mixed_samples=mixed,sample_exclusive_fraction=exclusive,sample_stability='NOT_VERIFIABLE' if len(group_fr)<2 else 'ASSESSED')
  if cfg['auto']['check_sample_stability'] and len(group_fr)>=2 and (mixed<c['min_mixed_samples'] or exclusive>c['max_sample_exclusive_fraction']):out['reasons'].append('SAMPLE_COMPONENT_INSTABILITY')
  shifts=[];bad=0
  if stability and c['stability_resamples']:
   rng=np.random.default_rng(c['random_state']+100)
   for _ in range(c['stability_resamples']):
    zz=rng.choice(z,len(z),replace=True)
    with warnings.catch_warnings():
     warnings.simplefilter('ignore');g=GaussianMixture(2,n_init=1,reg_covar=c['reg_covar'],random_state=c['random_state']).fit(zz[:,None])
    mm=np.sort(g.means_.ravel());shifts.append(float(np.max(np.abs(mm-mu)/sd)))
    if not g.converged_:bad+=1
  out['resample_mean_shift_sd']=shifts
  if bad or (shifts and max(shifts)>c['max_mean_shift_sd']):out['reasons'].append('RESAMPLE_UNSTABLE')
  out['reliable']=not out['reasons'];out['threshold_low']=threshold(out,c['posterior_borderline']);out['threshold_high']=threshold(out,c['posterior_high'])
 except (ValueError,FloatingPointError,np.linalg.LinAlgError) as e:out['reasons'].append('MODEL_NOT_ESTIMABLE:'+str(e))
 return out
def fit_tail(x,meta,c):
 w=weights(meta);med=weighted_quantile(x,w,.5);mad=1.4826*weighted_quantile(abs(x-med),w,.5)
 lo=max(weighted_quantile(x,w,c['borderline_quantile']),med+c['borderline_mad_multiplier']*mad)
 hi=max(weighted_quantile(x,w,c['high_quantile']),med+c['high_mad_multiplier']*mad)
 return dict(threshold_low=lo,threshold_high=hi,median=med,scaled_MAD=mad,parameters=c,rule='max(weighted_quantile, median + multiplier * scaled_MAD)',no_observed_high=bool(hi>np.max(x)))
