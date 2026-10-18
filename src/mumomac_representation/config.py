import copy,yaml
from mumomac_apsc.config import load_config as c2_config,merge,tail_validate
from mumomac_coupling.registry import ORDER
DEFAULT={'version':'C2.1-experimental-0.1.0','policy':'continuous_first','default_mode':'auto','optional_lifecycle':['L01','L03'],'calibration':{'scope':'all_query','condition_column':'condition','control_labels':['Control'],'biological_unit_column':'biological_unit_id','sample_ids':None,'group_by':['tissue','coarse_identity','species','modality'],'freeze_before_evaluation':True},'distribution':{'min_calibration_blocks':3,'min_validation_cells':30,'max_validation_folds':8,'min_bic_gain_over_best_unimodal':10.,'min_mean_block_loglik_gain':.02,'min_positive_block_fraction':.75,'min_LOO_shape_fraction':.75,'min_LOO_label_agreement':.90,'initialization_seeds':[7,17,37],'min_initialization_agreement':.95,'max_outside_support_mass':.05,'candidate_models':['Gaussian','SkewNormal','GMM2'],'optional_dip_test':'NOT_INSTALLED_NOT_RUN','support_policy':'unbounded_working_density_with_outside_mass_QC; raw UCell unchanged'},'relative':{'quantile':.8,'fixed_thresholds':{},'tail':{'high_quantile':.975,'borderline_quantile':.90,'high_mad_multiplier':3.,'borderline_mad_multiplier':1.5}},'tail_profile':{'enabled':True,'quantiles':[.70,.80,.85,.90,.95],'pairs':['C01','C02','C03','C04'],'high_rule':'strict_score_gt_weighted_quantile; never split ties'},'joint_index':{'enabled':True,'default_output':False,'percentile_rule':'weighted_mid_CDF','pairs':['C01','C02','C03','C04']},'overrides':{},'coupling':{'min_cells':30,'min_table_cell':5,'covariates':['nUMI','nGenes'],'adjust_activity_secondary':True}}
MODES={'auto','force_gmm','force_tail','continuous_only','relative_high'}
def load(value=None):
 if isinstance(value,(str,bytes)):
  with open(value) as f:value=yaml.safe_load(f)
 value=value or {};value=value.get('representation',value)
 z=merge(DEFAULT,value)
 if z['policy']!='continuous_first' or z['default_mode'] not in MODES:raise ValueError('Invalid continuous-first policy/mode.')
 if set(z['optional_lifecycle'])-set(ORDER[16:21]):raise ValueError('Invalid lifecycle selection; A01 must stay continuous.')
 if not 0<z['relative']['quantile']<1:raise ValueError('Relative quantile must lie in (0,1).')
 if z['tail_profile']['quantiles']!=sorted(set(z['tail_profile']['quantiles'])) or any(not 0<q<1 for q in z['tail_profile']['quantiles']):raise ValueError('Tail quantiles must be ordered unique fractions.')
 for p,v in z['relative']['fixed_thresholds'].items():
  if p not in ORDER or p=='A01' or not 0<=v<=1:raise ValueError('Invalid fixed UCell threshold.')
 for p,o in z['overrides'].items():
  if p not in ORDER or o.get('mode') not in MODES:raise ValueError('Invalid per-program mode.')
  if p=='A01' and o['mode']!='continuous_only':raise ValueError('A01 cannot be binarized.')
  if 'quantile' in o and not 0<o['quantile']<1:raise ValueError('Invalid override quantile.')
  if 'fixed_threshold' in o and not 0<=o['fixed_threshold']<=1:raise ValueError('Invalid fixed threshold.')
  tail_validate(merge(z['relative']['tail'],{**o.get('tail',{}),**{k:v for k,v in o.items() if k in z['relative']['tail']}}))
 # Validate the inherited scope, ordered MAD/quantiles and non-bypassable C2 diagnostics.
 c2_config({'calibration':z['calibration'],'optional_lifecycle':z['optional_lifecycle'],'tail':z['relative']['tail']})
 d=z['distribution']
 if d['min_calibration_blocks']<3 or d['max_validation_folds']<3 or d['min_validation_cells']<3:raise ValueError('At least three calibration blocks and valid held-out cells are required for supported partitioning.')
 for k in ['min_positive_block_fraction','min_LOO_shape_fraction','min_LOO_label_agreement','min_initialization_agreement','max_outside_support_mass']:
  if not 0<d[k]<=1:raise ValueError('Invalid distribution quality fraction.')
 return z
def inherited(c):return c2_config({'calibration':c['calibration'],'optional_lifecycle':c['optional_lifecycle'],'tail':c['relative']['tail']})
