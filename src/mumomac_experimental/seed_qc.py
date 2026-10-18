"""Original-seed evidence adapter. No UCell mathematics or frozen sets changed."""
from pathlib import Path
from importlib import resources
import json,hashlib
import numpy as np,pandas as pd
from mumomac_coupling.registry import ORDER,digest
from mumomac_coupling.adapters import ScoreSpec,ScoreBundle,from_official

class SeedScoreSpec(ScoreSpec):
 def validate(self):
  if self.namespace!='CURRENT_UCELL_V1_ORIGINAL_SEED_QC_V1' or self.signature_version!='original':raise ValueError('Invalid Original SeedQC namespace.')
  p=self.parameters
  if p.get('original_coverage_policy')!='partial_seed3_v1' or p.get('original_min_present')!=3 or p.get('original_min_detected_per_cell')!=3:raise ValueError('Frozen Seed3 technical policy required.')
  if not p.get('effective_program_hashes') or not p.get('seed_qc_fingerprint'):raise ValueError('Effective-subset and detection provenance required.')
  if set(p['effective_program_members'])!=set(ORDER):raise ValueError('Effective memberships need all22 frozen programs.')
  import mumomac
  frozen=mumomac.reference_signatures('original',self.species)
  for k,v in p['effective_program_members'].items():
   if digest(sorted(v))!=p['effective_program_hashes'][k]:raise ValueError('Effective subset hash changed.')
   if len(v)!=len(set(v)) or not set(v)<=set(frozen[k]['genes']):raise ValueError('Effective members must be a subset of frozen Original members.')
  ScoreSpec('CURRENT_UCELL_V1_ORIGINAL','original',self.species,self.scoring_method,self.scoring_version,self.program_hashes,self.reference_version,p,self.mapping_version,self.mapping_hash).validate()
  return True

def spec_from_dict(x):
 if x.get('parameters',{}).get('original_coverage_policy')=='partial_seed3_v1' and x.get('namespace')!='CURRENT_UCELL_V1_ORIGINAL_SEED_QC_V1':raise ValueError('Partial SeedQC cannot masquerade as full Original.')
 return SeedScoreSpec(**x) if x.get('namespace')=='CURRENT_UCELL_V1_ORIGINAL_SEED_QC_V1' else ScoreSpec(**x)
def validate_policy(c):
 if c.get('signature_version','expanded') not in {'original','expanded','both'}:raise ValueError('signature_version must be original,expanded or both.')
 for k in ['original_min_present','original_min_detected_per_cell']:
  if c.get(k,3)!=3:raise ValueError('partial_seed3_v1 freezes '+k+'=3; compare other cutoffs only as QC sensitivity.')
 if c.get('original_coverage_policy','partial_seed3_v1')!='partial_seed3_v1':raise ValueError('Unknown Original coverage policy.')
 if c.get('scoring_qc_mode','strict') not in {'strict','audit'}:raise ValueError('scoring_qc_mode must be strict or audit.')

def seed_evidence(a,m,c,symbol_key=None):
 import mumomac
 symbols=a.var_names.tolist() if symbol_key is None else a.var[symbol_key].tolist();ix={str(v):j for j,v in enumerate(symbols) if not pd.isna(v)}
 root=resources.files('mumomac.resources')
 with root.joinpath('human_one_to_one.tsv').open() as f:smap=pd.read_csv(f,sep='\t')
 with root.joinpath('human_la_one_to_one_v1.tsv').open() as f:lamap=pd.read_csv(f,sep='\t')
 smap=dict(zip(smap.mouse_gene,smap.human_gene));lamap=dict(zip(lamap.mouse_gene,lamap.human_gene))
 native=mumomac.reference_signatures('original','mouse');mapped=mumomac.reference_signatures('original',c['species']);genes=sorted({g for v in mapped.values() for g in v['genes'] if g in ix})
 detection_layer=c.get('seed_detection_layer','counts' if 'counts' in a.layers else c['layer']);x=a.X if detection_layer=='X' else a.layers.get(detection_layer)
 if x is None:raise ValueError('Declared seed detection layer absent.')
 counts_layer=detection_layer=='counts' or detection_layer==c['layer'] and c['expression_type']=='counts'
 if not counts_layer and c.get('detection_zero_preserving') is not True:raise ValueError('Original detection without counts requires explicit detection_zero_preserving=true and transformation provenance.')
 for lo in range(0,a.n_obs,2048):
  b=x[lo:lo+2048];v=b.data if hasattr(b,'tocsr') else np.asarray(b)
  if not np.isfinite(v).all() or (v<0).any() or counts_layer and (v!=np.rint(v)).any():raise ValueError('Seed detection expression is not valid finite nonnegative declared counts/expression.')
 detect=np.asarray((x[:,[ix[g] for g in genes]]>0).toarray() if hasattr(x,'tocsr') else (x[:,[ix[g] for g in genes]]>0),dtype=bool);lookup={g:j for j,g in enumerate(genes)};counts=np.zeros((a.n_obs,22),np.int16);catalog=[];cover=[];subsets={}
 for j,p in enumerate(ORDER):
  frozen=native[p]['genes'];mapper={g:g for g in frozen} if c['species']=='mouse' else smap if p.startswith('S') else lamap;members=[mapper[g] for g in frozen if g in mapper and mapper[g] in ix];subsets[p]=members
  counts[:,j]=detect[:,[lookup[g] for g in members]].sum(1) if members else 0
  unmapped=[g for g in frozen if g not in mapper];absent=[g for g in frozen if g in mapper and mapper[g] not in ix]
  cover.append(dict(program_id=p,feature_class='Functional Program' if p.startswith('S') else 'Activity Covariate' if p=='A01' else 'Lifecycle/Residence',original_seed_count=len(frozen),mapped_seed_count=len(frozen)-len(unmapped),measured_seed_count=len(members),absent_seed_count=len(absent),mapping_loss_count=len(unmapped),coverage_fraction=len(members)/len(frozen),mapped_coverage_fraction=len(members)/(len(frozen)-len(unmapped)) if len(frozen)>len(unmapped) else 0.,exact_measured_members=';'.join(members),exact_missing_members=';'.join(absent),exact_unmapped_members=';'.join(unmapped),eligibility='TECHNICAL_SEED3_ELIGIBLE' if p.startswith('S') and len(members)>=3 else 'UNASSESSABLE_LT3_MEASURED' if p.startswith('S') else 'FROZEN_LA_POLICY_UNCHANGED',completeness_status='FULL_ORIGINAL_SEED' if len(members)==len(frozen) else 'PARTIAL_ORIGINAL_SEED',definition_namespace='CURRENT_ORIGINAL_FULL' if len(members)==len(frozen) else 'CURRENT_ORIGINAL_PARTIAL_SEED3' if len(members)>=3 else 'ORIGINAL_UNASSESSABLE',original_signature_hash=native[p]['gene_set_hash'],measured_subset_hash=digest(sorted(members))))
  for g in frozen:
   mg=mapper.get(g);measured=mg in ix if mg else False;d=detect[:,lookup[mg]] if measured else None;by=[]
   if measured:
    for sid,indices in m.groupby(['study_id','source_block_id','sample_id'],sort=False,dropna=False).indices.items():by.append(dict(study_id=str(sid[0]),source_block_id=str(sid[1]),sample_id=str(sid[2]),n_cells=len(indices),n_detected=int(d[indices].sum()),fraction=float(d[indices].mean())))
   catalog.append(dict(program_id=p,original_gene=g,mapped_gene=mg,mapping_status='NATIVE_MOUSE' if c['species']=='mouse' else 'FROZEN_ONE_TO_ONE' if mg else 'UNMAPPED',measured_in_dataset=measured,measurement_status='NOT_MAPPED' if not mg else 'NOT_MEASURED' if not measured else 'MEASURED_ZERO_IN_ALL_SELECTED_CELLS' if not d.any() else 'MEASURED',n_cells_detected=int(d.sum()) if measured else None,fraction_cells_detected=float(d.mean()) if measured else None,relevant_sample_detection=json.dumps(by,ensure_ascii=False),exclusion_reason='NO_ACCEPTABLE_ORTHOLOG' if not mg else 'ABSENT_FROM_MEASURED_AXIS' if not measured else 'MEASURED; ZERO_IS_NOT_ASSAY_ABSENCE' if not d.any() else ''))
 packed=np.packbits(detect,axis=1,bitorder='little');fingerprint=digest(dict(cell_ids=list(a.obs_names),gene_universe=symbols,subsets=subsets,detected_counts=counts,packed_sha256=hashlib.sha256(packed.tobytes()).hexdigest(),detection_layer=detection_layer,expression_type='counts' if counts_layer else c['expression_type']))
 return dict(catalog=pd.DataFrame(catalog),coverage=pd.DataFrame(cover),counts=pd.DataFrame(counts,index=a.obs_names,columns=ORDER),subsets=subsets,packed=packed,genes=genes,fingerprint=fingerprint,cell_ids=list(a.obs_names),detection_layer=detection_layer,expression_type='counts' if counts_layer else c['expression_type'],detection_evidence='EXPRESSION_NONZERO_AUDITED' if counts_layer else 'DECLARED_ZERO_PRESERVING_EXPRESSION; TRANSFORMATION_ASSERTED_BY_CALLER')

def original_bundle(reference_scored,available_scored,m,evidence,c):
 ref=from_official(reference_scored,m,'original');avail=from_official(available_scored,m,'original');raw=avail.scores.copy();raw.loc[:,list(ORDER[16:])]=ref.scores.loc[:,list(ORDER[16:])]
 cover=evidence['coverage'].set_index('program_id');cov=ref.coverage.copy();valid=raw.copy()
 for p in ORDER[:16]:
  ok=int(cover.loc[p,'measured_seed_count'])>=3;valid.loc[evidence['counts'][p]<3,p]=np.nan
  if not ok:valid[p]=np.nan
  cov.loc[p,'full_original_coverage']=cover.loc[p,'coverage_fraction'];cov.loc[p,'measured_subset_size']=cover.loc[p,'measured_seed_count'];cov.loc[p,'coverage']=1. if ok else 0.;cov.loc[p,'scoring_eligible']=ok;cov.loc[p,'coverage_scope']='COMPLETE_OF_EFFECTIVE_MEASURED_SUBSET; NOT_FULL_ORIGINAL';cov.loc[p,'present_genes']=';'.join(evidence['subsets'][p])
 cov['seed3_policy_scope']='S_ORIGINAL_ONLY; LA_UNCHANGED'
 params={**ref.spec.parameters,'coverage_policy':'partial_seed3_v1; LA_reference','original_coverage_policy':'partial_seed3_v1','original_min_present':3,'original_min_detected_per_cell':3,'seed_qc_fingerprint':evidence['fingerprint'],'effective_program_members':evidence['subsets'],'effective_program_hashes':{p:digest(sorted(v)) for p,v in evidence['subsets'].items()},'raw_scorer_parameters':avail.spec.parameters,'full_reference_coverage':cover.coverage_fraction.to_dict(),'valid_score_policy':'S Dataset K>=3 and cell detected>=3; NA otherwise; LA formal reference unchanged'}
 spec=SeedScoreSpec('CURRENT_UCELL_V1_ORIGINAL_SEED_QC_V1','original',ref.spec.species,'UCell',ref.spec.scoring_version,ref.spec.program_hashes,ref.spec.reference_version,params,ref.spec.mapping_version,ref.spec.mapping_hash)
 return ScoreBundle(valid,m,spec,cov,evidence['counts']),raw

def detection_table(evidence,m,raw,valid):
 rows=[];n=len(m);cover=evidence['coverage'].set_index('program_id')
 for p in ORDER:
  cov=cover.loc[p];ok=np.isfinite(valid[p].to_numpy());det=evidence['counts'][p].to_numpy();reason=np.where(ok,'PASS_TECHNICAL_NOT_FUNCTIONAL_POSITIVE','FROZEN_LA_UNASSESSABLE_OR_ZERO_CELL')
  if p.startswith('S'):reason=np.where(cov.measured_seed_count<3,'LT3_MEASURED_SEEDS',np.where(det<3,'LT3_DETECTED_SEEDS','RAW_SCORE_NOT_FINITE'));reason=np.where(ok,'PASS_TECHNICAL_NOT_FUNCTIONAL_POSITIVE',reason)
  rows.append(pd.DataFrame(dict(cell_id=m.index,sample_id=m.sample_id.to_numpy(),program_id=p,n_seed_total=cov.original_seed_count,n_seed_measured=cov.measured_seed_count,n_seed_detected=det,coverage_fraction=cov.coverage_fraction,eligible_for_scoring=ok,reason_if_invalid=reason,raw_score_for_qc=raw[p].to_numpy(),valid_program_score=valid[p].to_numpy(),definition_namespace=cov.definition_namespace,original_signature_hash=cov.original_signature_hash,measured_subset_hash=cov.measured_subset_hash)))
 return pd.concat(rows,ignore_index=True)

def export_seed_evidence(e,out,m,raw,valid):
 out=Path(out);e['catalog'].to_csv(out/'program_seed_gene_qc.tsv',sep='\t',index=False,na_rep='NA');e['coverage'].to_csv(out/'program_seed_coverage.tsv',sep='\t',index=False,na_rep='NA');dt=detection_table(e,m,raw,valid);dt.to_parquet(out/'cell_seed_detection.parquet',index=False)
 np.savez_compressed(out/'seed_detection_bits.npz',packed=e['packed']);m.reset_index(names='cell_id').to_parquet(out/'seed_cell_index.parquet',index=False)
 manifest=dict(schema='SeedQC1',fingerprint=e['fingerprint'],n_cells=len(m),genes=e['genes'],cell_id_order_hash=digest(list(m.index)),gene_universe_hash=e.get('gene_universe_hash'),effective_program_hashes={p:digest(sorted(v)) for p,v in e['subsets'].items()},raw_score_hash=digest(raw.to_dict('split')),valid_score_hash=digest(valid.to_dict('split')),bitorder='little',detection_layer=e['detection_layer'],expression_type=e['expression_type'],detection_evidence=e['detection_evidence'],file_sha256={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['program_seed_gene_qc.tsv','program_seed_coverage.tsv','cell_seed_detection.parquet','seed_detection_bits.npz','seed_cell_index.parquet']})
 (out/'seed_qc_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');return dt,manifest

def verify_seed_evidence(run):
 run=Path(run)
 if not (run/'seed_qc_manifest.json').exists() and (run.parent/'seed_qc_manifest.json').exists():run=run.parent
 manifest=json.loads((run/'seed_qc_manifest.json').read_text())
 for name,h in manifest['file_sha256'].items():
  if hashlib.sha256((run/name).read_bytes()).hexdigest()!=h:raise ValueError('SeedQC artifact changed: '+name)
 return run,manifest

def seed_details(run,cell_id,program_id):
 run,manifest=verify_seed_evidence(run);cells=pd.read_parquet(run/'seed_cell_index.parquet',columns=['cell_id']);ix=pd.Index(cells.cell_id.astype(str));pos=ix.get_indexer([str(cell_id)])[0]
 if pos<0 or ix.has_duplicates:raise ValueError('Unknown or non-unique cell ID.')
 cat=pd.read_csv(run/'program_seed_gene_qc.tsv',sep='\t');cat=cat[cat.program_id.eq(program_id)]
 if not len(cat):raise ValueError('Unknown program.')
 bits=np.unpackbits(np.load(run/'seed_detection_bits.npz')['packed'][pos],bitorder=manifest['bitorder'])[:len(manifest['genes'])];hit={g:bool(v) for g,v in zip(manifest['genes'],bits)}
 result=dict(cell_id=str(cell_id),program_id=program_id,detected_genes=[],measured_but_zero_genes=[],not_measured_genes=[],unmapped_genes=[],detection_evidence=manifest['detection_evidence'])
 for r in cat.itertuples():
  key='unmapped_genes' if pd.isna(r.mapped_gene) else 'not_measured_genes' if not r.measured_in_dataset else 'detected_genes' if hit[r.mapped_gene] else 'measured_but_zero_genes';result[key].append(dict(original_gene=r.original_gene,mapped_gene=None if pd.isna(r.mapped_gene) else r.mapped_gene))
 return result
