"""Descriptive within-record estimands; no automatic cell-level significance."""
import numpy as np,pandas as pd
from scipy import stats
def four_table(n11,n10,n01,n00,*,min_cells=30,min_table_cell=5):
 c=np.asarray([n11,n10,n01,n00],dtype=float)
 if (c<0).any() or not np.isfinite(c).all() or not np.equal(c,np.rint(c)).all():raise ValueError('Four-table counts must be finite nonnegative integers.')
 n=int(c.sum());a=n11+n10;b=n11+n01
 row=dict(n11=int(n11),n10=int(n10),n01=int(n01),n00=int(n00),n_evaluable=n,A_fraction=a/n if n else np.nan,B_fraction=b/n if n else np.nan,joint_fraction=n11/n if n else np.nan,lift=(n11*n/(a*b)) if n and a and b else np.nan,OR_raw=np.nan,logOR_raw=np.nan,OR=np.nan,logOR=np.nan,continuity_correction=0.,raw_OR_status='UNDEFINED',binary_estimable=False,estimability_status='NO_EVALUABLE_CELLS')
 if not n:return row
 if a==0 or b==0 or a==n or b==n:row['estimability_status']='DEGENERATE_MARGIN';return row
 num=float(n11)*n00;den=float(n10)*n01
 if den>0:
  row['OR_raw']=num/den;row['logOR_raw']=float(np.log(num/den)) if num>0 else -np.inf;row['raw_OR_status']='FINITE_POSITIVE' if num>0 else 'ZERO_OR'
 elif num>0:row.update(OR_raw=np.inf,logOR_raw=np.inf,raw_OR_status='INFINITE_ZERO_DENOMINATOR')
 if (c==0).any():
  q=c+.5;row.update(OR=float(q[0]*q[3]/(q[1]*q[2])),logOR=float(np.log(q[0]*q[3]/(q[1]*q[2]))),continuity_correction=.5,estimability_status='SPARSE_HALDANE_DESCRIPTIVE_ONLY')
 else:row.update(OR=row['OR_raw'],logOR=row['logOR_raw'],estimability_status='DESCRIPTIVE_SPARSE' if c.min()<min_table_cell else 'DESCRIPTIVE_SUPPORTED')
 row['binary_estimable']=bool(n>=min_cells and c.min()>=min_table_cell)
 if n<min_cells:row['estimability_status']='TOO_FEW_CERTAIN_CELLS'
 return row

def spearman(x,y,min_cells=3):
 x=np.asarray(x,float);y=np.asarray(y,float);ok=np.isfinite(x)&np.isfinite(y)
 if ok.sum()<min_cells or np.ptp(x[ok])<1e-12 or np.ptp(y[ok])<1e-12:return np.nan
 return float(stats.spearmanr(x[ok],y[ok]).statistic)

def rank_residual(x,y,metadata,covariates,min_cells=30):
 x=np.asarray(x,float);y=np.asarray(y,float);parts=[];used=[];dropped=[]
 base=np.isfinite(x)&np.isfinite(y);meta=metadata.copy()
 for cov in covariates:
  if cov not in meta:return dict(adjusted_association=np.nan,status='MISSING_COVARIATE:'+cov,n_adjusted=0,covariates_used='',dropped='')
  s=meta[cov]
  if pd.api.types.is_numeric_dtype(s):
   v=s.to_numpy(float);base &= np.isfinite(v)
  else:
   valid=s.notna()&~s.astype(str).str.lower().isin(['','unknown','nan','na']);base &= valid.to_numpy()
 if base.sum()<min_cells:return dict(adjusted_association=np.nan,status='INSUFFICIENT_COMPLETE_COVARIATES',n_adjusted=int(base.sum()),covariates_used='',dropped='')
 for cov in covariates:
  s=meta.loc[base,cov]
  if s.nunique()<2:dropped.append(cov+':constant');continue
  if pd.api.types.is_numeric_dtype(s):new=stats.rankdata(s.to_numpy(float))[:,None]
  else:new=pd.get_dummies(s.astype(str),drop_first=True,dtype=float).to_numpy()
  for k in range(new.shape[1]):parts.append(new[:,k]);used.append(cov+(f':{k}' if new.shape[1]>1 else ''))
 n=int(base.sum());design=np.column_stack([np.ones(n)]+parts)
 if len(parts)==0:return dict(adjusted_association=np.nan,status='NO_VARIABLE_COVARIATE',n_adjusted=n,covariates_used='',dropped=';'.join(dropped))
 rank=int(np.linalg.matrix_rank(design))
 if rank!=design.shape[1] or n<max(min_cells,5*rank) or n-rank<10:return dict(adjusted_association=np.nan,status='RANK_DEFICIENT_OR_TOO_COMPLEX',n_adjusted=n,covariates_used=';'.join(used),dropped=';'.join(dropped))
 xr=stats.rankdata(x[base]);yr=stats.rankdata(y[base]);rx=xr-design@np.linalg.lstsq(design,xr,rcond=None)[0];ry=yr-design@np.linalg.lstsq(design,yr,rcond=None)[0]
 rho=float(np.corrcoef(rx,ry)[0,1]) if np.std(rx)>1e-10 and np.std(ry)>1e-10 else np.nan
 return dict(adjusted_association=rho,status='DESCRIPTIVE_PARTIAL_RANK' if np.isfinite(rho) else 'NO_RESIDUAL_VARIATION',n_adjusted=n,covariates_used=';'.join(used),dropped=';'.join(dropped))

def welch_effect(case,control):
 a=np.asarray(case,float);b=np.asarray(control,float);a=a[np.isfinite(a)];b=b[np.isfinite(b)];na=len(a);nb=len(b)
 out=dict(n_case_units=na,n_control_units=nb,effect=float(a.mean()-b.mean()) if na and nb else np.nan,ci_low=np.nan,ci_high=np.nan,p_value=np.nan,df=np.nan,status='NOT_ESTIMABLE_INSUFFICIENT_BIOLOGICAL_UNITS')
 if min(na,nb)<3:return out
 va=a.var(ddof=1)/na;vb=b.var(ddof=1)/nb;v=va+vb
 if not v>1e-15:out['status']='NOT_ESTIMABLE_ZERO_BIOLOGICAL_VARIANCE';return out
 df=v*v/(va*va/(na-1)+vb*vb/(nb-1));se=np.sqrt(v);crit=stats.t.ppf(.975,df)
 out.update(ci_low=out['effect']-crit*se,ci_high=out['effect']+crit*se,p_value=float(2*stats.t.sf(abs(out['effect']/se),df)),df=float(df),status='BIOLOGICAL_UNIT_WELCH; small_n_and_design_limits_apply')
 return out

def adjust_bh(pvalues):
 p=np.asarray(pvalues,float);valid=np.isfinite(p)
 if not len(p):return p
 x=np.where(valid,p,1.);order=np.argsort(x,kind='stable');q=np.empty(len(x));q[order]=np.minimum(1.,np.minimum.accumulate((x[order]*len(x)/np.arange(1,len(x)+1))[::-1])[::-1]);q[~valid]=np.nan;return q
