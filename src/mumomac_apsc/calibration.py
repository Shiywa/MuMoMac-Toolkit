import numpy as np,pandas as pd
from mumomac_coupling.registry import digest
def select(bundle,cfg):
 c=cfg['calibration'];m=bundle.metadata
 if c['scope']=='frozen_reference':raise ValueError('No qualified frozen Current UCell normal reference exists; historical Mean/q20 thresholds are incompatible.')
 if any(k not in m for k in c['group_by']):raise ValueError('Missing declared calibration stratum.')
 keep=np.ones(len(m),bool)
 if c['scope']=='control':
  if c['condition_column'] not in m or not c['control_labels']:raise ValueError('Explicit control labels required.')
  keep &= m[c['condition_column']].isin(c['control_labels']).to_numpy()
 if c['sample_ids'] is not None:keep &= m.sample_id.isin(c['sample_ids']).to_numpy()
 fit=m.iloc[np.flatnonzero(keep)].copy()
 if fit.empty:raise ValueError('No calibration observations selected.')
 if c['scope']=='control':
  if not fit.biological_unit_confirmed.all():raise ValueError('Control calibration requires confirmed biological-unit mapping; sample IDs are not animals.')
  unit=c['biological_unit_column']
  if unit not in fit or fit[unit].astype(str).isin(['Unknown','','nan']).any():raise ValueError('Invalid declared biological-unit key.')
  # An alternate key cannot silently split a verified animal into technical libraries.
  if fit.groupby(['source_block_id','biological_unit_id'])[unit].nunique().max()>1:raise ValueError('Declared biological-unit key splits a confirmed animal.')
 return fit
def key(values):return '|'.join(str(v).replace('|','%7C') for v in values)
def independence(meta,model):
 cal={(r['source_block_id'],r['biological_unit_id']) for r in model['calibration_units'] if r['biological_unit_confirmed']}
 out=[]
 for r in meta.itertuples():out.append('INDEPENDENT_OF_THRESHOLD_FIT' if r.biological_unit_confirmed and cal and (r.source_block_id,r.biological_unit_id) not in cal else 'CALIBRATION_MEMBER_OR_UNIT_UNVERIFIED')
 return np.array(out)
