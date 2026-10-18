"""Sample/identity-stratified coupling and bounded, on-demand cell exports."""
from pathlib import Path
import json,time,resource
import numpy as np,pandas as pd
from .registry import ORDER,VERSION,pairs,signatures,digest,canonicalize
from .adapters import known
from .states import fit_thresholds,call_states,pair_codes,LABELS,PAIR_LABELS
from .statistics import four_table,spearman,rank_residual
from .inference import biological_summary,condition_effects
GROUPS=['study_id','source_block_id','sample_id','tissue','coarse_identity','condition','species','modality']
def sample_statistics(bundle,calls,config,pair_list=None):
 pair_list=pair_list if pair_list is not None else pairs(config.get('optional_lifecycle',[]));rows=[]
 sig=signatures(bundle.spec.species,bundle.spec.signature_version);x=bundle.scores.to_numpy();m=bundle.metadata;idx={p:i for i,p in enumerate(ORDER)}
 for keys,indices in m.groupby(GROUPS,observed=True,dropna=False,sort=False).indices.items():
  indices=np.asarray(indices);base=dict(zip(GROUPS,keys));meta=m.iloc[indices];r=meta.iloc[0]
  for p in pair_list:
   a=idx[p['program_A']];b=idx[p['program_B']];pc=pair_codes(calls['codes'][indices,a],calls['codes'][indices,b]);ct=np.bincount(pc,minlength=6)
   binary=four_table(*ct[:4],min_cells=config.get('min_cells',30),min_table_cell=config.get('min_table_cell',5))
   xa=x[indices,a].astype(float).copy();xb=x[indices,b].astype(float).copy()
   if not calls['eligible_programs'][a] or not calls['eligible_programs'][b]:xa[:]=np.nan;xb[:]=np.nan
   n_cont=int((np.isfinite(xa)&np.isfinite(xb)).sum());rho=spearman(xa,xb);covmeta=meta.copy();covmeta['A01']=bundle.scores.iloc[indices].A01.to_numpy()
   adjusted=rank_residual(xa,xb,covmeta,config.get('covariates',['nUMI','nGenes']),config.get('min_cells',30))
   activity=rank_residual(xa,xb,covmeta,config.get('covariates',['nUMI','nGenes'])+['A01'],config.get('min_cells',30)) if config.get('adjust_activity_secondary',True) else {}
   ga=set(sig[p['program_A']]['genes']);gb=set(sig[p['program_B']]['genes']);overlap=sorted(ga&gb)
   rows.append(dict(**base,**p,**binary,n_total=len(indices),n_uncertain=int(ct[4]),n_unassessable=int(ct[5]),n_continuous=n_cont,biological_unit_id=r.biological_unit_id,biological_unit_confirmed=bool(meta.biological_unit_confirmed.all()),batch=r.batch if meta.batch.nunique()==1 else 'MIXED_WITHIN_RECORD',spearman=rho,continuous_estimable=bool(n_cont>=config.get('min_cells',30) and np.isfinite(rho)),adjusted_association=adjusted['adjusted_association'],adjustment_status=adjusted['status'],n_adjusted=adjusted['n_adjusted'],covariates_used=adjusted['covariates_used'],dropped_covariates=adjusted['dropped'],activity_adjusted_association=activity.get('adjusted_association',np.nan),activity_adjustment_status=activity.get('status','NOT_REQUESTED'),gene_overlap_n=len(overlap),gene_overlap_jaccard=len(overlap)/len(ga|gb),gene_overlap_fraction_A=len(overlap)/len(ga),gene_overlap_fraction_B=len(overlap)/len(gb),shared_genes=';'.join(overlap),coverage_A=float(bundle.coverage.loc[p['program_A'],'coverage']),coverage_B=float(bundle.coverage.loc[p['program_B'],'coverage']),threshold_method=config.get('threshold_mode','none'),score_namespace=bundle.spec.namespace,technical_QC='DESCRIPTIVE_WITHIN_RECORD; overlap/depth/threshold sensitivity required',p_value=np.nan))
 result=pd.DataFrame(rows)
 if len(result):assert (result[['n11','n10','n01','n00','n_uncertain','n_unassessable']].sum(axis=1)==result.n_total).all()
 return result

def export_outputs(bundle,calls,model,samples,rep,effects,out,config,pair_list):
 import pyarrow as pa,pyarrow.parquet as pq
 out=Path(out);out.mkdir(parents=True,exist_ok=True)
 if any(out.iterdir()):raise FileExistsError('Output directory is not empty; old results are not overwritten.')
 writers={};n=len(bundle.scores);chunk=int(config.get('export_chunk_size',1024));idx={p:i for i,p in enumerate(ORDER)}
 wanted=config.get('export_pairs',['C01','C03','C04']);selected=[p for p in pair_list if p['pair_id'] in wanted]
 if set(wanted)-{p['pair_id'] for p in selected}:raise ValueError('Export pair missing; lifecycle pairs require explicit applicability.')
 if len(selected)>8 and not config.get('allow_large_pair_export',False):raise ValueError('More than8 cell-pair exports requires explicit allow_large_pair_export; sample statistics remain all-pair.')
 def append(name,df):
  table=pa.Table.from_pandas(df,preserve_index=False)
  if name not in writers:writers[name]=pq.ParquetWriter(out/(name+'.parquet'),table.schema,compression='zstd')
  writers[name].write_table(table)
 try:
  for lo in range(0,n,chunk):
   hi=min(n,lo+chunk);meta=bundle.metadata.iloc[lo:hi].copy();meta.insert(0,'cell_id',meta.index);meta=meta.reset_index(drop=True)
   scores=pd.concat([meta,bundle.scores.iloc[lo:hi].reset_index(drop=True)],axis=1);scores['scoring_method']=bundle.spec.scoring_method;scores['scoring_version']=bundle.spec.scoring_version;scores['score_namespace']=bundle.spec.namespace;scores['program_version']=bundle.spec.reference_version+'/'+bundle.spec.signature_version;scores['QC_references']='coverage_and_state_QC.tsv; run_manifest.json'
   append('cell_program_scores',scores)
   state_columns={}
   for j,pid in enumerate(ORDER):
    state_columns.update({pid+'__score':bundle.scores.iloc[lo:hi,j].to_numpy(),pid+'__state':LABELS[calls['codes'][lo:hi,j]],pid+'__threshold':calls['thresholds'][lo:hi,j],pid+'__coverage':float(bundle.coverage.iloc[j].coverage),pid+'__uncertainty_status':calls['reasons'][lo:hi,j]})
   states=pd.concat([meta[['cell_id','sample_id','coarse_identity','species','modality']],pd.DataFrame(state_columns,index=meta.index)],axis=1)
   states['threshold_method']=model['mode'];states['threshold_definition_hash']=model['threshold_definition_hash'];states['score_namespace']=bundle.spec.namespace;append('cell_program_states',states)
   for p in selected:
    a=idx[p['program_A']];b=idx[p['program_B']];z=meta.copy();z['pair_id']=p['pair_id'];z['program_A']=p['program_A'];z['program_B']=p['program_B'];z['A_score']=bundle.scores.iloc[lo:hi,a].to_numpy();z['B_score']=bundle.scores.iloc[lo:hi,b].to_numpy();z['A_state']=LABELS[calls['codes'][lo:hi,a]];z['B_state']=LABELS[calls['codes'][lo:hi,b]];z['pair_state']=PAIR_LABELS[pair_codes(calls['codes'][lo:hi,a],calls['codes'][lo:hi,b])];z['threshold_A']=calls['thresholds'][lo:hi,a];z['threshold_B']=calls['thresholds'][lo:hi,b];z['coverage_A']=float(bundle.coverage.iloc[a].coverage);z['coverage_B']=float(bundle.coverage.iloc[b].coverage);z['score_namespace']=bundle.spec.namespace;z['scoring_method']=bundle.spec.scoring_method;z['state_definition_version']=VERSION;z['threshold_definition_hash']=model['threshold_definition_hash'];z['QC']='RELATIVE_STATE_ONLY; not function or molecular synergy'
    spots=z.modality.str.contains('visium|spot',case=False)&~z.modality.str.contains('hd.*cell',case=False)
    z['observation_level']=np.where(spots,'MIXED_SPOT_CO_SIGNAL','CELL_OR_SEGMENTED_CELL');z.loc[spots,'pair_state']=z.loc[spots,'pair_state'].map(lambda s:'Spot-co-signal' if s=='Joint' else 'Spot-'+s);z['spatial_autocorrelation_controlled']=False;append('cell_combination_states',z)
 finally:
  for w in writers.values():w.close()
 samples.to_parquet(out/'sample_coupling_statistics.parquet',index=False);rep.to_csv(out/'coupling_reproducibility.tsv',sep='\t',index=False,na_rep='NA');effects.to_csv(out/'condition_effects.tsv',sep='\t',index=False,na_rep='NA')
 qc=[]
 for j,pid in enumerate(ORDER):
  counts=np.bincount(calls['codes'][:,j],minlength=4);qc.append(dict(program_id=pid,**dict(zip(LABELS,counts)),coverage=bundle.coverage.iloc[j].coverage,scoring_eligible=bool(calls['eligible_programs'][j]),missing_score_column=pid in bundle.missing_programs))
 pd.DataFrame(qc).to_csv(out/'coverage_and_state_QC.tsv',sep='\t',index=False)
 def write(name,obj):(out/name).write_text(json.dumps(canonicalize(obj),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
 write('frozen_thresholds.json',model);write('score_spec.json',bundle.spec.to_dict());write('run_config.json',config)
 write('run_manifest.json',dict(engine_version=VERSION,input_fingerprint=bundle.fingerprint(),n_cells=n,n_programs=22,n_primary_S_pairs=120,n_total_pairs=len(pair_list),n_sample_pair_rows=len(samples),n_cell_export_pairs=len(selected),score_namespace=bundle.spec.namespace,threshold_definition_hash=model['threshold_definition_hash'],no_cell_level_pvalues=bool(samples.p_value.isna().all()),biological_unit_known_cells=int(bundle.metadata.biological_unit_confirmed.sum()),gene_overlap_rescoring='NOT_PERFORMED_BY_CORE; optional expression diagnostic is separate',spatial_permutation_inference=False,independent_validation=False))

def run_engine(bundle,config=None,output=None,frozen_thresholds=None):
 config=dict(config or {});start=time.perf_counter()
 if config.get('min_cells',30)<3 or not 0<=config.get('min_gene_coverage',.8)<=1:raise ValueError('Invalid cell-count or coverage gate.')
 plist=pairs(config.get('optional_lifecycle',[]));model=frozen_thresholds or fit_thresholds(bundle,config);config.setdefault('threshold_mode',model['mode'])
 if model['rows']:
  config.setdefault('quantile',model['rows'][0]['quantile']);config.setdefault('quantile_band',model['rows'][0]['quantile_band'])
 calls=call_states(bundle,model,config);sample=sample_statistics(bundle,calls,config,plist)
 cal=model.get('calibration_sample_unit_mapping',[])
 calkeys={(r['study_id'],r['source_block_id'],r['sample_id']) for r in cal};evalkeys=set(bundle.metadata[['study_id','source_block_id','sample_id']].itertuples(index=False,name=None))
 calunits={(r['source_block_id'],r['biological_unit_id']) for r in cal if r['biological_unit_confirmed']};evalunits={(r.source_block_id,r.biological_unit_id) for r in bundle.metadata.itertuples() if r.biological_unit_confirmed}
 disjoint=bool(calkeys and not calkeys&evalkeys and calunits and evalunits and not calunits&evalunits and all(r['biological_unit_confirmed'] and known(r['source_block_id']) for r in cal) and bundle.metadata.biological_unit_confirmed.all())
 sample['threshold_inference_status']='INDEPENDENT_FROZEN_RULE' if disjoint else 'CALIBRATION_OVERLAP_OR_UNVERIFIED; binary condition inference disabled'
 rep=biological_summary(sample);effect=condition_effects(sample,config.get('condition_contrasts',[]),config.get('bootstrap_replicates',2000),config.get('seed',20261016))
 if output is not None:export_outputs(bundle,calls,model,sample,rep,effect,output,config,plist)
 return dict(samples=sample,reproducibility=rep,effects=effect,thresholds=model,calls=calls,pairs=plist,resource=dict(wall_seconds=time.perf_counter()-start,peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2))
