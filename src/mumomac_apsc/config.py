"""Frozen experimental rules; no outcome-directed selection."""
import copy
import yaml
from mumomac_coupling.registry import ORDER,digest
DEFAULT={
 'version':'C2-experimental-0.1.0','default_mode':'auto','optional_lifecycle':[],
 'calibration':{'scope':'all_query','condition_column':'condition','control_labels':['Control'],'biological_unit_column':'biological_unit_id','sample_ids':None,'group_by':['tissue','coarse_identity','species','modality'],'freeze_before_evaluation':True},
 'qc':{'min_cells':100,'min_unique':20,'min_sd':1e-5,'max_tied_mass':.5,'max_boundary_mass':.5,'min_gene_coverage':.8,'min_detected_members':1,'min_detection_fraction':.5},
 'auto':{'assess_distribution':True,'check_sample_stability':True,'check_technical_confounding':True,'allow_tail_fallback':True,'allow_continuous_fallback':True},
 'gmm':{'n_components':2,'n_init':3,'reg_covar':1e-5,'random_state':17,'min_separation':1.5,'min_high_weight':.005,'max_high_weight':.95,'posterior_high':.90,'posterior_borderline':.50,'require_model_diagnostics':True,'min_component_cells':30,'min_bic_gain':10.,'kde_bandwidth_factors':[.8,1.,1.2],'kde_min_peak_recurrence':2,'kde_max_valley_ratio':.8,'kde_min_prominence':.05,'max_fit_cells':4000,'max_group_explained_variance':.7,'max_depth_rho':.65,'min_sample_cells':30,'min_mixed_samples':2,'max_sample_exclusive_fraction':.8,'stability_resamples':3,'max_mean_shift_sd':.75},
 'tail':{'high_quantile':.975,'borderline_quantile':.90,'high_mad_multiplier':3.,'borderline_mad_multiplier':1.5},
 'overrides':{},'coupling':{'continuous_always':True,'use_C1_engine':True,'min_cells':30,'min_table_cell':5,'covariates':['nUMI','nGenes'],'adjust_activity_secondary':True},
 'sensitivity':{'enabled':True,'tail_quantiles':[.8,.9,.975]}}
MODES={'auto','force_gmm','force_tail','continuous_only'}
def merge(a,b):
 z=copy.deepcopy(a)
 for k,v in b.items():z[k]=merge(z[k],v) if k in z and isinstance(z[k],dict) and isinstance(v,dict) else copy.deepcopy(v)
 return z
def tail_validate(t):
 if not 0<t['borderline_quantile']<=t['high_quantile']<1:raise ValueError('Tail quantiles must satisfy 0 < borderline <= high < 1.')
 if not 0<=t['borderline_mad_multiplier']<=t['high_mad_multiplier']:raise ValueError('MAD multipliers must be nonnegative and ordered.')
def load_config(value=None):
 if isinstance(value,(str,bytes)):
  with open(value) as f:value=yaml.safe_load(f)
 value=value or {};migration=None
 if 'apsc' in value:value=value['apsc']
 elif 'threshold_mode' in value:
  old=value;scope={'exploratory':'all_query','control':'control','none':'all_query','reference':'frozen_reference'}[old['threshold_mode']]
  # This migrates interface intent, not historical q80 +/- .02 numeric labels.
  value={'default_mode':'continuous_only' if old['threshold_mode']=='none' else 'auto','calibration':{'scope':scope,'control_labels':[old.get('control_condition','Control')]},'optional_lifecycle':old.get('optional_lifecycle',[]),'coupling':{k:old[k] for k in ['min_cells','min_table_cell','covariates','adjust_activity_secondary'] if k in old}}
  migration='C1 interface migrated; original q80 labels retained as separate baseline, not silently reproduced by APSC'
 z=merge(DEFAULT,value);z['migration_note']=migration
 if z['default_mode'] not in MODES:raise ValueError('Unknown APSC mode.')
 if z['calibration']['scope'] not in ['all_query','control','frozen_reference']:raise ValueError('Unknown calibration scope.')
 if not z['calibration']['freeze_before_evaluation']:raise ValueError('Calibration must freeze before evaluation.')
 forbidden={'condition',z['calibration']['condition_column'],'sample_id','biological_unit_id',z['calibration']['biological_unit_column']}
 keys=z['calibration']['group_by']
 if not keys or set(keys)&forbidden:raise ValueError('No condition/sample/unit-specific threshold refitting.')
 if not z['coupling']['continuous_always'] or not z['coupling']['use_C1_engine']:raise ValueError('Frozen C1 continuous engine is mandatory.')
 if not z['auto']['assess_distribution'] or not z['gmm']['require_model_diagnostics']:raise ValueError('Distribution diagnostics cannot be disabled in this experimental release.')
 if not z['auto']['check_sample_stability'] or not z['auto']['check_technical_confounding']:raise ValueError('Sample/source QC cannot be bypassed in this experimental release.')
 if set(z['optional_lifecycle'])-set(ORDER[16:21]):raise ValueError('A01 cannot be binary; lifecycle requires explicit applicability.')
 tail_validate(z['tail'])
 for p,o in z['overrides'].items():
  if p not in ORDER or o.get('mode','auto') not in MODES:raise ValueError('Invalid program override.')
  if p=='A01' and o.get('mode')!='continuous_only':raise ValueError('A01 must remain continuous.')
  tail_validate(merge(z['tail'],{k:v for k,v in o.items() if k in z['tail']}))
 g=z['gmm']
 if g['n_components']!=2 or g['n_init']<1 or g['reg_covar']<=0 or not 0<g['posterior_borderline']<g['posterior_high']<1:raise ValueError('Invalid GMM parameters.')
 if not 0<g['min_high_weight']<g['max_high_weight']<1 or g['min_separation']<=0:raise ValueError('Invalid component gates.')
 if z['qc']['min_cells']<3 or not 0<=z['qc']['min_gene_coverage']<=1:raise ValueError('Invalid basic QC.')
 for name in ['max_tied_mass','max_boundary_mass','min_detection_fraction']:
  if not 0<z['qc'][name]<=1:raise ValueError('Invalid fraction QC.')
 return z
def for_program(c,p):
 o=c['overrides'].get(p,{})
 return o.get('mode',c['default_mode']),merge(c['tail'],{k:v for k,v in o.items() if k in c['tail']}),o.get('reason','User selected '+o.get('mode',c['default_mode']))
