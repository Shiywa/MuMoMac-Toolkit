"""One orchestrator, with exact pinned engines and a stable export contract."""
from pathlib import Path
import json,hashlib,time,traceback,resource,importlib.metadata as im
import numpy as np,pandas as pd,yaml
from scipy import sparse
from mumomac_coupling.registry import ORDER,canonicalize,digest,pairs
from mumomac_coupling.adapters import ScoreBundle,ScoreSpec,from_official
from mumomac_representation import engine as representation
from mumomac_representation.config import load as load_representation_config,DEFAULT as REPRESENTATION_DEFAULT
from .input import load_input,validate_input,select_momac,harmonize_genes,read_table,MODALITIES

def write_json(obj,path):Path(path).write_text(json.dumps(canonicalize(obj),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def file_hash(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def load_config(config):
 c=yaml.safe_load(Path(config).read_text()) if isinstance(config,(str,Path)) else dict(config)
 if not isinstance(c,dict):raise ValueError('Configuration must be a mapping.')
 return c
def score_programs(a,c,scorer_key=None):
 import mumomac
 return mumomac.score_signature_versions(a,species=c['species'],signature_version=c.get('signature_version','expanded'),layer=c['layer'],expression_type=c['expression_type'],gene_symbol_key=scorer_key,duplicate_policy='error',coverage_policy='reference',max_rank=1500,chunk_size=int(c.get('scoring_chunk_size',500)),copy=True)
def normalize_metadata(m,c):
 m=m.copy();m.index=m.index.astype(str);m.index.name='cell_id'
 for dest,opt in [('sample_id','sample_column'),('coarse_identity','identity_column'),('study_id','study_column'),('source_block_id','source_column'),('condition','condition_column'),('tissue','tissue_column'),('biological_unit_id','biological_unit_column'),('biological_unit_confirmed','biological_unit_confirmed_column')]:
  if c.get(opt):
   if c[opt] not in m:raise ValueError('Missing configured metadata column: '+c[opt])
   m[dest]=m[c[opt]]
 if 'coarse_identity' not in m:m['coarse_identity']=m[c['cell_type_column']]
 for k in ['species','modality']:
  if k in m and not m[k].eq(c[k]).all():raise ValueError('Configuration conflicts with metadata '+k)
  m[k]=c[k]
 # Biological IDs are only confirmed by an explicit input column; never from sample ID.
 if 'biological_unit_confirmed' not in m:m['biological_unit_confirmed']=False
 return m
def representation_config(c):
 r=dict(c.get('representation',{}));r.setdefault('optional_lifecycle',[])
 selected=c.get('selected_pairs',[])
 for key in ['joint_index','tail_profile']:r[key]={**r.get(key,{}),'pairs':selected}
 r['tail_profile']['enabled']=bool(c.get('include_tail',False));r['joint_index']['default_output']=bool(c.get('include_joint_index',False))
 cfg=load_representation_config(r)
 if cfg['distribution']!=REPRESENTATION_DEFAULT['distribution']:raise ValueError('C3 freezes C2.1 evidence criteria; change only explicit state modes/manual rules, not distribution acceptance thresholds.')
 pairs(cfg['optional_lifecycle'],only=selected) # Unknown or inapplicable L pairs fail early.
 return cfg
def _new_output(output):
 out=Path(output)
 if out.exists() and any(out.iterdir()):raise FileExistsError('Use a new empty output directory. Failed runs are preserved; restart in a new directory.')
 out.mkdir(parents=True,exist_ok=True);return out
def _failure(out,error):
 write_json(dict(status='FAILED',error_type=type(error).__name__,error=str(error),partial_results_not_certified=True),out/'run_qc.json')
 (out/'run_log.txt').write_text(traceback.format_exc())
def _input_manifest(paths):return [{'path':str(Path(p).resolve()),'bytes':Path(p).stat().st_size,'sha256':file_hash(p)} for p in paths if p]

def _legacy_analyze_expanded(input_path,config,output,metadata=None):
 c=load_config(config);out=_new_output(output);start=time.perf_counter()
 try:
  a=load_input(input_path,metadata);input_qc=validate_input(a,c);selected,audit=select_momac(a,c)
  input_qc['original_observation_axis_name']=selected.obs_names.name
  selected.obs_names=selected.obs_names.rename('cell_id') # Axis label only; native IDs unchanged.
  selected,genes,mapping,key=harmonize_genes(selected,c);m=normalize_metadata(selected.obs,c)
  x=selected.X if c['layer']=='X' else selected.layers[c['layer']]
  m['input_expression_total']=np.asarray(x.sum(1)).ravel();m['input_detected_genes']=np.asarray((x>0).sum(1)).ravel()
  if c['expression_type']=='counts':
   m['nUMI']=m.input_expression_total;m['nGenes']=m.input_detected_genes
  # Formal scorer receives every native gene including opaque background.
  scored=score_programs(selected,c,key);bundle=from_official(scored,m,'expanded')
  symbols=selected.var_names.astype(str).tolist() if key is None else selected.var[key].tolist();gene_detected=np.asarray((x>0).sum(0)).ravel()>0;symbol_index={str(v):j for j,v in enumerate(symbols) if not pd.isna(v)}
  bundle.coverage['n_detected_genes']=[sum(bool(gene_detected[symbol_index[g]]) for g in str(v).split(';') if g in symbol_index) for v in bundle.coverage.present_genes]
  bundle.coverage['detection_coverage']=bundle.coverage.n_detected_genes/bundle.coverage.n_mapped_genes
  coords=None
  if c['modality']=='cell_resolved_spatial':
   xy=np.asarray(selected.obsm[c.get('coordinate_key','spatial')]);coords=pd.DataFrame(xy,index=selected.obs_names,columns=[f'coordinate_{j}' for j in range(xy.shape[1])]);coords.index.name='cell_id'
  audit.to_parquet(out/'selection_audit.parquet',index=False);genes.to_parquet(out/'gene_harmonization.parquet',index=False)
  return _execute(bundle,c,out,start,_input_manifest([input_path,metadata]),input_qc,mapping,coords)
 except Exception as e:_failure(out,e);raise

def _analyze_scores_expanded(scores,metadata,spec,coverage,config,output,detected_members=None,frozen_model=None):
 c=load_config(config);out=_new_output(output);start=time.perf_counter()
 try:
  if c.get('species') not in {'mouse','human'} or c.get('modality') not in MODALITIES:raise ValueError('Explicit supported species/modality required.')
  s=read_table(scores);m=read_table(metadata)
  for frame in [s,m]:
   if 'cell_id' not in frame:raise ValueError('Score-only tables require explicit cell_id.')
  ss=ScoreSpec(**json.loads(Path(spec).read_text()))
  if ss.scoring_method!='UCell' or ss.signature_version!='expanded':raise ValueError('Unified Current Expanded22 entry rejects legacy/original/mixed scoring namespaces.')
  if set(ORDER)-set(s):raise ValueError('All22 score columns required; unassessable programs use explicit NA.')
  if not ss.parameters.get('gene_universe_hash'):raise ValueError('Score-only requires full native gene universe hash.')
  values=s[list(ORDER)].to_numpy(float)
  if np.isinf(values).any():raise ValueError('Infinite scores are invalid.')
  finite=values[np.isfinite(values)]
  if len(finite) and (finite.min() < -1e-6 or finite.max()>1+1e-6):raise ValueError('Declared UCell scores outside[0,1]; validate before float32 conversion.')
  m=normalize_metadata(m.set_index('cell_id'),c)
  label=c.get('cell_type_column');targets=c.get('target_labels')
  if not label or label not in m or not targets or not m[label].isin(targets).all():raise ValueError('Score-only input must contain explicitly provided target annotations for every observation.')
  if m.sample_id.isna().any() or m.sample_id.astype(str).isin(['Unknown','nan','']).any():raise ValueError('Complete sample IDs required.')
  if c['modality']=='mixed_spots':raise ValueError('Mixed-spot input refused by single-cell Mo/Mac workflow.')
  if c['modality']=='cell_resolved_spatial' and (not c.get('segmentation_confirmed') or not c.get('score_coordinate_columns')):raise ValueError('Spatial score-only requires segmentation confirmation and explicit coordinate columns.')
  det=read_table(detected_members).set_index('cell_id') if detected_members else None
  if det is None and c.get('representation',{}).get('default_mode','auto')!='continuous_only':raise ValueError('Member-detection table required for state evidence; otherwise explicitly choose continuous_only.')
  if det is not None:
   if det.index.has_duplicates or set(det.index)!=set(m.index) or set(ORDER)-set(det):raise ValueError('Detection table IDs/programs must match exactly.')
   dv=det[list(ORDER)].to_numpy(float)
   if not np.isfinite(dv).all() or (dv<0).any() or (dv!=np.rint(dv)).any():raise ValueError('Detection member counts must be finite nonnegative integers.')
  bundle=ScoreBundle(s,m,ss,read_table(coverage),det)
  coords=m[c['score_coordinate_columns']].rename(columns={v:f'coordinate_{i}' for i,v in enumerate(c['score_coordinate_columns'])}) if c.get('score_coordinate_columns') else None
  if coords is not None and not np.isfinite(coords.to_numpy(float)).all():raise ValueError('Invalid spatial coordinates.')
  model=json.loads(Path(frozen_model).read_text()) if frozen_model else None
  return _execute(bundle,c,out,start,_input_manifest([scores,metadata,spec,coverage,detected_members,frozen_model]),{'status':'PASS','n_cells':len(m),'n_native_genes':'NO_EXPRESSION_LOADED','mode':'scores_only; no re-scoring'},dict(species=ss.species,mapping_version=ss.mapping_version,mapping_hash=ss.mapping_hash,source='validated input ScoreSpec'),coords,model)
 except Exception as e:_failure(out,e);raise

def _analysis_bundle(bundle,c):
 if c['modality']!='cell_resolved_spatial':return bundle
 # C2 has an explicit RNA-only modality guard. The integration presents segmented
 # cells as compatible RNA observations for engineering applicability only, retaining
 # real modality in input/export and declaring cross-domain evidence unvalidated.
 meta=bundle.metadata.copy();meta['original_modality']=meta.modality;meta['modality']='scRNA'
 return ScoreBundle(bundle.scores,meta,bundle.spec,bundle.coverage,bundle.detected_members)

def _execute(bundle,c,out,start,manifest,input_qc,mapping,coords=None,model=None):
 import pyarrow as pa,pyarrow.parquet as pq
 from .versions import validate_pinned_cores
 core_manifest=validate_pinned_cores()
 cfg=representation_config(c);analysis=_analysis_bundle(bundle,c)
 if model is not None and model['config_hash']!=digest(cfg):raise ValueError('Frozen representation config differs; provide exactly the frozen analysis config.')
 model=model or representation.fit(analysis,cfg)
 representation.validate(analysis,model)
 n=len(bundle.scores);chunk=int(c.get('export_chunk_size',2048))
 if chunk<1:raise ValueError('export_chunk_size must be positive.')
 codes=np.full((n,22),3,np.int8);percentiles=np.full((n,22),np.nan);writer=None;jwriter=None;param=None;elig=None;state_counts={};include_j=bool(c.get('include_joint_index',False))
 try:
  for lo in range(0,n,chunk):
   hi=min(lo+chunk,n);part=ScoreBundle(analysis.scores.iloc[lo:hi],analysis.metadata.iloc[lo:hi],analysis.spec,analysis.coverage,None if analysis.detected_members is None else analysis.detected_members.iloc[lo:hi])
   calls=representation.call(part,model);codes[lo:hi]=calls['codes'];percentiles[lo:hi]=calls['percentiles'];param=calls['parameter_hash'];elig=calls['eligible_programs']
   frame=calls['table'];frame['input_modality']=c['modality'];frame['analysis_unit']='segmented_cell' if c['modality']=='cell_resolved_spatial' else 'cell_or_nucleus'
   if c['modality']=='cell_resolved_spatial':frame['domain_evidence']='TARGETED_SPATIAL_ENGINEERING_EXTENSION; NO_CROSS_MODALITY_BIOLOGICAL_CERTIFICATION'
   else:frame['domain_evidence']='RNA_EXPLORATORY'
   for k,v in frame.representation_type.value_counts().items():state_counts[k]=state_counts.get(k,0)+int(v)
   # All chunks share explicit schema; None-only state columns must not become Arrow null.
   for k in frame:
    if frame[k].dtype=='object':frame[k]=frame[k].astype('string')
   table=pa.Table.from_pandas(frame,preserve_index=False)
   if writer is None:writer=pq.ParquetWriter(out/'cell_program_representation.parquet',table.schema)
   writer.write_table(table.cast(writer.schema))
   has_binary=any(np.any((calls['codes'][:,ORDER.index(p['program_A'])]<3)&(calls['codes'][:,ORDER.index(p['program_B'])]<3)) for p in pairs(cfg['optional_lifecycle'],only=c.get('selected_pairs',[])))
   if c.get('selected_pairs'):
    joint=representation.joint_features(part,calls,model)
    if not include_j:joint[['u_A','u_B','J_soft']]=np.nan;joint['J_soft_interpretation']='NOT_REQUESTED; selected raw scores and optional binary state only'
    joint['input_modality']=c['modality']
    for k in joint:
     if joint[k].dtype=='object':joint[k]=joint[k].astype('string')
    jt=pa.Table.from_pandas(joint,preserve_index=False)
    if jwriter is None:jwriter=pq.ParquetWriter(out/'selected_pair_joint_features.parquet',jt.schema)
    jwriter.write_table(jt.cast(jwriter.schema))
 finally:
  if writer:writer.close()
  if jwriter:jwriter.close()
 if jwriter is None:pd.DataFrame(columns=['cell_id','sample_id','pair_id','raw_A','raw_B','binary_pair_state','J_soft','QC']).to_parquet(out/'selected_pair_joint_features.parquet',index=False)
 calls=dict(codes=codes,percentiles=percentiles,eligible_programs=elig,parameter_hash=param)
 continuous=representation.continuous_coupling(analysis,calls,model)
 # Prevent optional lifecycle pairs in monocytes from escaping applicability guards.
 bad=continuous.coarse_identity.str.lower().str.contains('monocyte')&(continuous.program_A.str.startswith('L')|continuous.program_B.str.startswith('L'))
 if bad.any():
  for col in ['spearman','adjusted_association','activity_adjusted_association','joint_fraction','lift','OR','logOR']:continuous.loc[bad,col]=np.nan
  continuous.loc[bad,'continuous_estimable']=False;continuous.loc[bad,'technical_QC']='LIFECYCLE_PAIR_NOT_APPLICABLE_IN_MONOCYTES'
 continuous['input_modality']=c['modality'];continuous['analysis_unit']='segmented_cell' if c['modality']=='cell_resolved_spatial' else 'cell_or_nucleus'
 continuous['QC_reference']='run_qc.json; program_coverage.tsv; pairwise_gene_overlap.tsv'
 continuous['evidence_level']=continuous.biological_evidence
 continuous.to_parquet(out/'sample_coupling_statistics.parquet',index=False)
 if c.get('include_tail'):representation.tail_profiles(analysis,model).assign(input_modality=c['modality']).to_parquet(out/'tail_coenrichment_profiles.parquet',index=False)
 dec,cmp=representation.diagnostics(model);dec.to_csv(out/'program_representation_decisions.tsv',sep='\t',index=False);cmp.to_csv(out/'distribution_model_comparison.tsv',sep='\t',index=False)
 wide=bundle.scores.copy();wide.insert(0,'identity',bundle.metadata.coarse_identity);wide.insert(0,'sample_id',bundle.metadata.sample_id);wide=wide.reset_index(names='cell_id')
 for k,v in [('scoring_method',bundle.spec.scoring_method),('scoring_version',bundle.spec.scoring_version),('program_version',bundle.spec.reference_version+'/'+bundle.spec.signature_version),('QC_reference','run_qc.json; program_coverage.tsv')]:wide[k]=v
 wide.to_parquet(out/'cell_program_scores.parquet',index=False)
 bundle.metadata.reset_index(names='cell_id').to_parquet(out/'cell_metadata.parquet',index=False)
 if bundle.detected_members is not None:bundle.detected_members.reset_index(names='cell_id').to_parquet(out/'detected_members.parquet',index=False)
 if coords is not None:coords.reset_index(names='cell_id').to_parquet(out/'cell_coordinates.parquet',index=False)
 cov=bundle.coverage.reset_index(names='feature_id');cov.to_csv(out/'program_coverage.tsv',sep='\t',index=False)
 overlap=continuous[['pair_id','program_A','program_B','gene_overlap_n','gene_overlap_jaccard','gene_overlap_fraction_A','gene_overlap_fraction_B','shared_genes']].drop_duplicates();overlap.to_csv(out/'pairwise_gene_overlap.tsv',sep='\t',index=False)
 from mumomac_coupling.engine import GROUPS
 from mumomac_coupling.statistics import spearman
 technical=[]
 for vals,g in bundle.metadata.groupby(GROUPS,sort=False,observed=True,dropna=False):
  for p in ORDER:
   for covariate in ['nUMI','nGenes']:
    rho=spearman(bundle.scores.loc[g.index,p].to_numpy(),g[covariate].to_numpy())
    technical.append(dict(**dict(zip(GROUPS,vals)),program_id=p,covariate=covariate,n_cells=len(g),spearman=rho,status='DESCRIPTIVE_HIGH_DEPTH_DEPENDENCY' if np.isfinite(rho) and abs(rho)>.65 else 'DESCRIPTIVE' if np.isfinite(rho) else 'NOT_ESTIMABLE',interpretation='C1 Spearman; .65 inherited C2 technical warning, not biological certification'))
 pd.DataFrame(technical).to_csv(out/'technical_dependencies.tsv',sep='\t',index=False,na_rep='NA')
 write_json(bundle.spec.to_dict(),out/'score_spec.json');write_json(model,out/'frozen_representation.json')
 write_json({'parameter_hash':param,'model_hash':model['model_hash'],'config':cfg,'J_soft_requested':include_j,'tail_profile_requested':bool(c.get('include_tail'))},out/'representation_request.json')
 (out/'run_config.yaml').write_text(yaml.safe_dump(c,allow_unicode=True,sort_keys=False))
 write_json(manifest,out/'input_manifest.json');write_json(mapping,out/'ortholog_mapping_manifest.json')
 write_json({'reference':'1.0.0','signature':bundle.spec.signature_version,'feature_order':list(ORDER),'program_hashes':bundle.spec.program_hashes,'score_parameters':bundle.spec.parameters,'feature_classes':{'S':'16 functional programs','L':'5 lifecycle/residence features','A01':'continuous activity covariate'}},out/'gene_set_manifest.json')
 write_json({'python':__import__('sys').version,'dependencies':{p:im.version(p) for p in ['mumomac-toolkit','numpy','pandas','scipy','scikit-learn','anndata','pyucell','pyarrow','PyYAML']},'experimental_version':'0.1.1rc1','pinned_core_versions':core_manifest},out/'software_versions.json')
 qc=dict(status='PASS_WITH_LIMITATIONS',input=input_qc,n_input=input_qc['n_cells'],n_selected=n,n_sample_records=int(bundle.metadata[['study_id','source_block_id','sample_id']].drop_duplicates().shape[0]),n_confirmed_biological_units=int(bundle.metadata.loc[bundle.metadata.biological_unit_confirmed,['source_block_id','biological_unit_id']].drop_duplicates().shape[0]),n_programs_scoring_eligible=int(bundle.coverage.scoring_eligible.sum()),representation_rows=n*22,representation_counts=state_counts,model_hash=model['model_hash'],parameter_hash=param,score_spec_hash=digest(bundle.spec.to_dict()),gene_universe_hash=bundle.spec.parameters['gene_universe_hash'],binary_rule_not_available_means_NA_not_zero=True,p_values_generated=False,auto_classifier_used=False,reference_lookup_required=False,scoring_not_recomputed=input_qc.get('mode','').startswith('scores_only'),coordinates_preserved=coords is not None,spatial_domain_validated=False if c['modality']=='cell_resolved_spatial' else None,warning='RNA proxies; no absolute activation/normality thresholds; sample records are not automatically animals; historical exposure retained',runtime_seconds=time.perf_counter()-start,peak_rss_GiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,output_bytes_before_report=sum(p.stat().st_size for p in out.iterdir() if p.is_file()))
 write_json(qc,out/'run_qc.json')
 from .report import export_report
 export_report(out,bundle,dec,continuous,qc,c)
 (out/'run_log.txt').write_text('Completed pinned official scorer → C2.1 representation → C1 statistics.\nNo source data modified; no biological disease testing.\n'+json.dumps(canonicalize(qc),ensure_ascii=False,indent=2)+'\n')
 return {'output':str(out.resolve()),'qc':qc,'model':model,'decisions':dec,'continuous':continuous,'bundle':bundle}


def analyze(input_path,config,output,metadata=None,signature_version=None):
 from .integration import analyze as go
 return go(input_path,config,output,metadata,signature_version)

def analyze_scores(scores,metadata,spec,coverage,config,output,detected_members=None,frozen_model=None,seed_qc_run=None):
 from .integration import analyze_scores as go
 return go(scores,metadata,spec,coverage,config,output,detected_members,frozen_model,seed_qc_run)
