from pathlib import Path
import json,time,shutil,html
import numpy as np,pandas as pd,yaml
from mumomac_coupling.registry import ORDER,digest
from mumomac_coupling.adapters import ScoreBundle,ScoreSpec,from_official
from mumomac_coupling.engine import GROUPS
from mumomac_representation.engine import continuous_applicable
from .input import load_input,validate_input,select_momac,harmonize_genes,read_table
from .seed_qc import seed_evidence,original_bundle,export_seed_evidence,validate_policy,verify_seed_evidence,spec_from_dict,SeedScoreSpec

def run_versions(a,m,c,key):
 """Official reference Both once; available Original only when some members absent."""
 from .workflow import score_programs
 import mumomac
 ref=score_programs(a,{**c,'signature_version':'both'},key)
 r=from_official(ref,m,'original');m=r.metadata
 e=seed_evidence(a,m,c,key);e['gene_universe_hash']=r.spec.parameters['gene_universe_hash']
 missing=r.coverage.loc[list(ORDER[:16])].n_present_genes.ne(r.coverage.loc[list(ORDER[:16])].n_mapped_genes).any()
 available=mumomac.score_signature_versions(a,species=c['species'],signature_version='original',layer=c['layer'],expression_type=c['expression_type'],gene_symbol_key=key,duplicate_policy='error',coverage_policy='available',max_rank=1500,chunk_size=int(c.get('scoring_chunk_size',500)),copy=True) if missing else ref
 original,raw=original_bundle(ref,available,m,e,c);expanded=from_official(ref,m,'expanded')
 # Expanded detection/eligibility stays exactly the official full-signature policy.
 return {'original':original,'expanded':expanded},raw,e,ref

def analyze(input_path,config,output,metadata=None,signature_version=None):
 from .workflow import load_config,_new_output,_failure,normalize_metadata,_input_manifest,_execute
 c=load_config(config)
 if signature_version is not None:c['signature_version']=signature_version
 c.setdefault('signature_version','expanded');validate_policy(c);out=_new_output(output);start=time.perf_counter()
 try:
  a=load_input(input_path,metadata);iq=validate_input(a,c);a,selection=select_momac(a,c);iq['original_observation_axis_name']=a.obs_names.name;a.obs_names=a.obs_names.rename('cell_id');a,genes,mapping,key=harmonize_genes(a,c);m=normalize_metadata(a.obs,c)
  x=a.X if c['layer']=='X' else a.layers[c['layer']];m['input_expression_total']=np.asarray(x.sum(1)).ravel();m['input_detected_genes']=np.asarray((x>0).sum(1)).ravel()
  if c['expression_type']=='counts':m['nUMI']=m.input_expression_total;m['nGenes']=m.input_detected_genes
  bundles,raw,e,ref=run_versions(a,m,c,key);m=bundles['original'].metadata;dt,manifest=export_seed_evidence(e,out,m,raw,bundles['original'].scores)
  coords=None
  if c['modality']=='cell_resolved_spatial':
   xy=np.asarray(a.obsm[c.get('coordinate_key','spatial')]);coords=pd.DataFrame(xy,index=a.obs_names,columns=[f'coordinate_{j}' for j in range(xy.shape[1])]);coords.index.name='cell_id';coords.reset_index().to_parquet(out/'cell_coordinates.parquet',index=False)
  selection.to_parquet(out/'selection_audit.parquet',index=False);genes.to_parquet(out/'gene_harmonization.parquet',index=False)
  requested=['original','expanded'] if c['signature_version']=='both' else [c['signature_version']];results={};inputs=_input_manifest([input_path,metadata])
  for v in requested:
   target=out/v if len(requested)>1 else out;target.mkdir(exist_ok=True);cfg={**c,'signature_version':v};result=_execute(bundles[v],cfg,target,time.perf_counter(),inputs,iq,mapping,coords)
   _annotate(result,v,c,raw,e,dt,manifest,out);results[v]=result
  _finish_run(out,c,results,bundles,raw,e,dt,manifest,iq,inputs,start)
  result=results[requested[0]] if len(requested)==1 else dict(output=str(out.resolve()),qc=json.loads((out/'run_qc.json').read_text()),versions=results)
  result['seed_evidence']=e;return result
 except Exception as ex:_failure(out,ex);raise

def annotate_denominators(result,version,c,seed_dt=None):
 b=result['bundle'];s=result['continuous'].copy();m=b.metadata;scope_m=m.copy()
 if c['modality']=='cell_resolved_spatial':scope_m['modality']='scRNA' # same transparent frozen RNA adapter as C3
 cfg=__import__('mumomac_experimental.workflow',fromlist=['representation_config']).representation_config(c);groups={tuple(map(str,k)):v for k,v in scope_m.groupby(GROUPS,sort=False,observed=True,dropna=False).indices.items()}
 for j,r in s.iterrows():
  ix=np.asarray(groups[tuple(str(r[k]) for k in GROUPS)]);g=scope_m.iloc[ix];a=r.program_A;bb=r.program_B;va=np.isfinite(b.scores.iloc[ix][a].to_numpy())&continuous_applicable(g,a,cfg);vb=np.isfinite(b.scores.iloc[ix][bb].to_numpy())&continuous_applicable(g,bb,cfg)
  s.loc[j,'n_target_cells']=len(ix);s.loc[j,'n_A_valid']=int(va.sum());s.loc[j,'n_B_valid']=int(vb.sum());s.loc[j,'n_pair_valid']=int((va&vb).sum());s.loc[j,'fraction_pair_valid']=float((va&vb).mean());s.loc[j,'failure_reason_distribution']=json.dumps(dict(A_invalid=int((~va).sum()),B_invalid=int((~vb).sum()),both_invalid=int((~va&~vb).sum()),score_or_scope_failure=True))
  if not (va&vb).any():
   s.loc[j,['spearman','adjusted_association','activity_adjusted_association']]=np.nan;s.loc[j,'continuous_estimable']=False;s.loc[j,'n_continuous']=0
  if version=='original':
   member=b.spec.parameters.get('effective_program_members',{p:str(b.coverage.loc[p,'present_genes']).split(';') for p in ORDER});ga=set(member[a]);gb=set(member[bb]);shared=sorted(ga&gb)
   s.loc[j,'gene_overlap_n_full_definition']=r.gene_overlap_n;s.loc[j,'gene_overlap_n']=len(shared);s.loc[j,'gene_overlap_jaccard']=len(shared)/len(ga|gb) if ga|gb else np.nan;s.loc[j,'gene_overlap_fraction_A']=len(shared)/len(ga) if ga else np.nan;s.loc[j,'gene_overlap_fraction_B']=len(shared)/len(gb) if gb else np.nan;s.loc[j,'shared_genes']=';'.join(shared);s.loc[j,'overlap_scope']='ACTUAL_MEASURED_MEMBERS; full frozen overlap separately retained'
 s['signature_version']=version;s['n_binary_evaluable']=s.n_evaluable;s['fraction_denominator_policy']='Continuous n_pair_valid; binary n_evaluable excludes uncertain/NA; n_target_cells never substituted'
 for k in ['n_target_cells','n_A_valid','n_B_valid','n_pair_valid']:s[k]=s[k].astype(int)
 s.to_parquet(Path(result['output'])/'sample_coupling_statistics.parquet',index=False)
 overlap_cols=['pair_id','program_A','program_B','gene_overlap_n','gene_overlap_jaccard','gene_overlap_fraction_A','gene_overlap_fraction_B','shared_genes']
 if version=='original':overlap_cols+=['gene_overlap_n_full_definition','overlap_scope']
 s[overlap_cols].drop_duplicates().to_csv(Path(result['output'])/'pairwise_gene_overlap.tsv',sep='\t',index=False)
 result['continuous']=s;return s

def _annotate(result,version,c,raw,e,dt,manifest,root):
 from .workflow import write_json
 p=Path(result['output']);b=result['bundle'];annotate_denominators(result,version,c,dt);frame=pd.read_parquet(p/'cell_program_representation.parquet');frame['signature_version']=version
 if version=='original':
  lookup=dt.set_index(['cell_id','program_id']);keys=pd.MultiIndex.from_frame(frame[['cell_id','program_id']]);z=lookup.loc[keys]
  for col in ['raw_score_for_qc','n_seed_detected','n_seed_measured','n_seed_total','definition_namespace','measured_subset_hash','reason_if_invalid']:frame[col]=z[col].to_numpy()
  frame['full_original_coverage_fraction']=z.coverage_fraction.to_numpy();frame['gene_coverage_interpretation']='Core coverage over measured effective subset; full original fraction is a separate field'
  raw.reset_index(names='cell_id').to_parquet(p/'raw_cell_program_scores_for_qc.parquet',index=False)
 frame.to_parquet(p/'cell_program_representation.parquet',index=False)
 pair=pd.read_parquet(p/'selected_pair_joint_features.parquet')
 if len(pair):
  bad=~np.isfinite(pair.raw_A.to_numpy())|~np.isfinite(pair.raw_B.to_numpy());pair.loc[bad,'binary_pair_state']='Unassessable: invalid program score';pair['signature_version']=version;pair.to_parquet(p/'selected_pair_joint_features.parquet',index=False)
 if p!=root:write_json({'relative_seed_qc_root':'..','seed_qc_fingerprint':manifest['fingerprint']},p/'seed_qc_location.json')
 wide=pd.read_parquet(p/'cell_program_scores.parquet');wide['signature_version']=version
 if version=='original':
  for prog in ORDER:wide['raw_score_for_qc__'+prog]=raw.loc[wide.cell_id,prog].to_numpy()
 wide.to_parquet(p/'cell_program_scores.parquet',index=False)
 qc=result['qc'];qc.update(signature_version=version,original_qc_mode=c.get('scoring_qc_mode','strict'),seed_qc_fingerprint=manifest['fingerprint'],original_strict_min3_applied=version=='original',original_technical_qualification_not_activation=True,version_namespace=b.spec.namespace)
 write_json(qc,p/'run_qc.json');_seed_report(p,e,dt,version,c)

def _seed_report(out,e,dt,version,c):
 from .report import _table
 cov=e['coverage'];summary=cov[['program_id','original_seed_count','mapped_seed_count','measured_seed_count','coverage_fraction','completeness_status','exact_missing_members','exact_unmapped_members']]
 mode_text={'original':'本运行版本Original：cell_program_scores中的S为有效分数；raw_score_for_qc列仅审计，未达标不作为Low/Neither。program_coverage.coverage是有效已测子集的覆盖，不能解释成完整Seed覆盖；真实完整定义比例见full_original_coverage_fraction和program_seed_coverage.tsv。','expanded':'本运行版本Expanded：SeedQC仅提供Original技术审计，不对Expanded追加三基因门控。','both':'本运行版本Both：original__列为三基因门控后的Original有效分数；expanded__列遵循冻结Expanded规则。两个版本分别判定状态及Coupling，不能相互填补缺失或直接相减作生物学效应。'}[version]
 text='\n\n## Seed Gene Coverage & Reliability\n\n'+_table(summary)+'\n\nOriginal S需要至少3个可测Seed及每细胞至少3个非零Seed。达标只是最低技术资格，不是High阈值或功能准确性认证。K<N或同源映射损失均为Partial；不能与完整Original无标签混用。Expanded及L/A沿用冻结规则，不受S的三基因门槛重定义。\n\n'+mode_text+'\n明确区分Unmapped、Not measured与Measured but zero；Not measured不意味着生物学不表达。精确基因列表可用seed-details按需查询，位图默认不导出细胞×全基因文本。Human旧S映射release仍缺失；空间低覆盖不强行补值。\n'
 sample=dt[dt.program_id.str.startswith('S')].groupby(['sample_id','program_id'],sort=False).agg(n_target_cells=('cell_id','size'),n_eligible=('eligible_for_scoring','sum')).reset_index();sample['n_failed_seed_QC']=sample.n_target_cells-sample.n_eligible;sample['failure_fraction']=sample.n_failed_seed_QC/sample.n_target_cells
 text+='\n逐样本Original S资格（Expanded不追加本门槛）：\n\n'+_table(sample)+'\n\n完整成员列表、映射失败及未测量原因见program_seed_gene_qc.tsv；逐细胞状态见cell_seed_detection.parquet，Both子目录的SeedQC见seed_qc_location.json所指根目录。Original减少Expanded新增成员的覆盖依赖，但缩减与三基因门槛可能改变数值、排名及纳入细胞。Expanded保留冻结加权定义，不能由Original补齐NA。新增Partial只作探索性技术表示。\n'
 catalog=e.get('catalog')
 if catalog is not None:
  measured=catalog[catalog.measured_in_dataset].sort_values('fraction_cells_detected').head(12)
  text+='\n本数据中检出率最低的已测量Seed（零值不等于未测量）：\n\n'+_table(measured[['program_id','original_gene','mapped_gene','n_cells_detected','fraction_cells_detected']])+'\n'
 md=out/'MuMoMac_Report_ZH.md';md.write_text(md.read_text()+text);h=out/'MuMoMac_Report.html';htmltext=h.read_text();section='<h2>Seed Gene Coverage &amp; Reliability</h2>'+summary.to_html(index=False)+'<pre>'+html.escape(text)+'</pre>'
 h.write_text(htmltext.replace('</body>',section+'</body>') if '</body>' in htmltext else htmltext.replace('</html>',section+'</html>'))

def _finish_run(out,c,results,bundles,raw,e,dt,manifest,iq,inputs,start):
 from .workflow import write_json
 from .report import _table
 rows=[]
 for version,r in results.items():
  for vals,g in r['bundle'].metadata.groupby(['study_id','source_block_id','sample_id','coarse_identity','condition','biological_unit_id','biological_unit_confirmed'],dropna=False,sort=False):
   for prog in ORDER:
    n=int(np.isfinite(r['bundle'].scores.loc[g.index,prog]).sum());rows.append(dict(zip(['study_id','source_block_id','sample_id','coarse_identity','condition','biological_unit_id','biological_unit_confirmed'],vals))|dict(signature_version=version,program_id=prog,n_target_cells=len(g),n_valid=n,fraction_valid=n/len(g)))
 eligibility=pd.DataFrame(rows);eligibility.to_csv(out/'score_eligibility_by_sample.tsv',sep='\t',index=False)
 comparison=[]
 for p in ORDER:
  ro=bundles['original'].scores[p];ex=bundles['expanded'].scores[p];valid=np.isfinite(ro)&np.isfinite(ex)
  from mumomac_coupling.statistics import spearman
  comparison.append(dict(program_id=p,n_original_valid=int(np.isfinite(ro).sum()),n_expanded_valid=int(np.isfinite(ex).sum()),n_both_valid=int(valid.sum()),rank_association_within_query=spearman(ro[valid].to_numpy(),ex[valid].to_numpy()),raw_score_difference_as_biological_effect='NOT_COMPUTED; different definitions'))
 pd.DataFrame(comparison).to_csv(out/'scoring_version_comparison.tsv',sep='\t',index=False)
 if len(results)>1:
  wide=bundles['original'].metadata[['sample_id','coarse_identity']].reset_index(names='cell_id')
  for v,r in results.items():
   for p in ORDER:wide[v+'__'+p]=r['bundle'].scores[p].to_numpy()
  for p in ORDER:wide['raw_original_for_qc__'+p]=raw[p].to_numpy()
  wide.to_parquet(out/'cell_program_scores.parquet',index=False);bundles['original'].metadata.reset_index(names='cell_id').to_parquet(out/'cell_metadata.parquet',index=False)
  for name in ['cell_program_representation.parquet','sample_coupling_statistics.parquet','selected_pair_joint_features.parquet']:
   pd.concat([pd.read_parquet(Path(r['output'])/name).assign(signature_version=v) for v,r in results.items()],ignore_index=True).to_parquet(out/name,index=False)
  write_json({'versions':{v:r['bundle'].spec.to_dict() for v,r in results.items()},'mixed_version_pairs_allowed':False},out/'score_spec.json');write_json({'signature_version':'both','versions':{v:json.loads((Path(r['output'])/'gene_set_manifest.json').read_text()) for v,r in results.items()},'same_gene_universe':True},out/'gene_set_manifest.json')
  for name in ['software_versions.json','ortholog_mapping_manifest.json']:shutil.copyfile(Path(results['expanded']['output'])/name,out/name)
  write_json(inputs,out/'input_manifest.json');(out/'run_config.yaml').write_text(yaml.safe_dump(c,allow_unicode=True,sort_keys=False))
  text='# MuMoMac Original / Expanded / Both\n\n两个版本分别进行状态表示和Coupling；禁止跨版本默认程序对，缺失不能由另一版本自动补齐。原始分数不直接相减作为生物学效应。三基因只是Original技术QC，不是激活阈值。\n\n'+_table(pd.DataFrame(comparison))+'\n\n详情分别见 original/MuMoMac_Report.html 与 expanded/MuMoMac_Report.html。\n'
  (out/'MuMoMac_Report_ZH.md').write_text(text);(out/'MuMoMac_Report.html').write_text('<!doctype html><meta charset="utf-8"><body style="background:white;font-family:Arial"><h1>MuMoMac Both</h1><p><a href="original/MuMoMac_Report.html">Original</a> | <a href="expanded/MuMoMac_Report.html">Expanded</a></p><pre>'+html.escape(text)+'</pre></body></html>');_seed_report(out,e,dt,'both',c)
  q=dict(status='PASS_WITH_LIMITATIONS',signature_version='both',n_selected=len(raw),n_input=iq['n_cells'],versions={v:r['qc'] for v,r in results.items()},same_cell_ids=True,same_gene_universe=True,mixed_version_pairs=False,seed_qc_fingerprint=manifest['fingerprint'],runtime_seconds=time.perf_counter()-start)
  write_json(q,out/'run_qc.json');(out/'run_log.txt').write_text('Both mode: two separately frozen scoring/representation/coupling namespaces.\n'+json.dumps(q,ensure_ascii=False,indent=2)+'\n')
 else:
  q=next(iter(results.values()))['qc'];q['runtime_seconds_including_seed_QC']=time.perf_counter()-start;write_json(q,out/'run_qc.json')
 if c.get('scoring_qc_mode','strict')=='audit':dt.groupby('program_id').agg(n_cells=('cell_id','size'),n_valid=('eligible_for_scoring','sum'),raw_mean=('raw_score_for_qc','mean'),valid_mean=('valid_program_score','mean')).reset_index().to_csv(out/'audit_raw_vs_valid.tsv',sep='\t',index=False)

def analyze_scores(scores,metadata,spec,coverage,config,output,detected_members=None,frozen_model=None,seed_qc_run=None):
 from .workflow import _analyze_scores_expanded,load_config,_new_output,_failure,_execute,normalize_metadata,_input_manifest,write_json
 c=load_config(config);sd=json.loads(Path(spec).read_text())
 if 'versions' in sd:raise ValueError('Both score-only requires score-run --input-run; do not mix version-prefixed columns into one ScoreBundle.')
 if sd.get('signature_version')=='expanded':
  if c.get('signature_version','expanded') not in {'expanded'}:raise ValueError('Requested scoring version differs from input ScoreSpec.')
  r=_analyze_scores_expanded(scores,metadata,spec,coverage,c,output,detected_members,frozen_model);annotate_denominators(r,'expanded',c)
  from .missing_qc import export_unknown_seed_qc
  export_unknown_seed_qc(r['bundle'],Path(output),'expanded')
  for name in ['cell_program_scores.parquet','cell_program_representation.parquet','selected_pair_joint_features.parquet']:
   path=Path(output)/name;frame=pd.read_parquet(path);frame['signature_version']='expanded'
   if name=='selected_pair_joint_features.parquet' and len(frame):
    bad=~np.isfinite(frame.raw_A.to_numpy())|~np.isfinite(frame.raw_B.to_numpy());frame.loc[bad,'binary_pair_state']='Unassessable: invalid program score'
   frame.to_parquet(path,index=False)
  _score_only_summaries(r,Path(output),'expanded');return r
 out=_new_output(output);start=time.perf_counter()
 try:
  c.setdefault('signature_version','original');validate_policy(c)
  if c['signature_version']!='original':raise ValueError('Score-only request and Original input version differ.')
  ss=spec_from_dict(sd);ss.validate()
  if ss.signature_version!='original' or ss.scoring_method!='UCell':raise ValueError('Current Original UCell required,never Mean/q20.')
  s=read_table(scores);m=normalize_metadata(read_table(metadata).set_index('cell_id'),c);values=s[list(ORDER)].to_numpy(float);f=values[np.isfinite(values)]
  if 'sample_id' not in m or m.sample_id.isna().any() or m.sample_id.astype(str).isin(['','Unknown','nan']).any():raise ValueError('Explicit complete sample IDs required.')
  if np.isinf(values).any() or len(f) and (f.min()<-1e-6 or f.max()>1+1e-6):raise ValueError('Invalid UCell score range.')
  if not ss.parameters.get('gene_universe_hash'):raise ValueError('Full native gene universe provenance required.')
  if not m[c['cell_type_column']].isin(c['target_labels']).all():raise ValueError('Provided target annotations required.')
  cov=read_table(coverage);det=read_table(detected_members).set_index('cell_id') if detected_members else None
  if det is not None:
   dv=det[list(ORDER)].to_numpy(float)
   if det.index.has_duplicates or set(det.index)!=set(m.index) or not np.isfinite(dv).all() or (dv<0).any() or (dv!=np.rint(dv)).any():raise ValueError('Invalid member detection IDs/counts.')
  base=ScoreBundle(s,m,ss,cov,det);raw=base.scores.copy();valid=raw.copy()
  if seed_qc_run:
   qcroot,man=verify_seed_evidence(seed_qc_run);idx=pd.read_parquet(qcroot/'seed_cell_index.parquet').cell_id.astype(str).tolist()
   if set(idx)!=set(base.scores.index) or man['gene_universe_hash']!=ss.parameters['gene_universe_hash']:raise ValueError('SeedQC cells/native gene universe differ.')
   dt=pd.read_parquet(qcroot/'cell_seed_detection.parquet');vals=dt.pivot(index='cell_id',columns='program_id',values='valid_program_score').reindex(index=base.scores.index,columns=ORDER);rr=dt.pivot(index='cell_id',columns='program_id',values='raw_score_for_qc').reindex(index=base.scores.index,columns=ORDER)
   if not (np.array_equal(raw.to_numpy(),vals.to_numpy(),equal_nan=True) or np.array_equal(raw.to_numpy(),rr.to_numpy(),equal_nan=True)):raise ValueError('Scores do not match the audited raw or valid SeedQC artifact.')
   raw=rr;valid=vals;counts=dt.pivot(index='cell_id',columns='program_id',values='n_seed_detected').reindex(index=base.scores.index,columns=ORDER)
   if isinstance(ss,SeedScoreSpec) and man['fingerprint']!=ss.parameters['seed_qc_fingerprint']:raise ValueError('SeedQC policy fingerprint differs.')
   if not isinstance(ss,SeedScoreSpec):raise ValueError('Audited SeedQC score-only requires its corresponding C3.1 SeedScoreSpec,not a relabeled formal score.')
   b=ScoreBundle(valid,base.metadata,ss,base.coverage,counts)
   for name in man['file_sha256']:shutil.copyfile(qcroot/name,out/name)
   shutil.copyfile(qcroot/'seed_qc_manifest.json',out/'seed_qc_manifest.json')
  else:
   # Missing evidence never becomes a fake zero detection count or a certified score.
   valid.loc[:,list(ORDER[:16])]=np.nan
   if det is None:det=pd.DataFrame(np.nan,index=base.scores.index,columns=ORDER)
   b=ScoreBundle(valid,base.metadata,ss,base.coverage,det);write_json({'status':'MISSING_DETECTION_EVIDENCE','Original_S_valid_scores':'ALL_NA','raw_scores_retained_for_qc':True},out/'seed_qc_missing_evidence.json')
  coords=None
  if c['modality']=='mixed_spots':raise ValueError('Mixed spots are not single-cell Mo/Mac.')
  if c['modality']=='cell_resolved_spatial':
   keys=c.get('score_coordinate_columns')
   if not c.get('segmentation_confirmed') or not keys:raise ValueError('Spatial score-only requires explicit coordinates and segmentation confirmation.')
   coords=m[keys].rename(columns={v:f'coordinate_{j}' for j,v in enumerate(keys)})
  model=json.loads(Path(frozen_model).read_text()) if frozen_model else None
  result=_execute(b,c,out,start,_input_manifest([scores,metadata,spec,coverage,detected_members,frozen_model]),dict(status='PASS',n_cells=len(m),mode='scores_only; no re-scoring'),{'source':'input ScoreSpec','mapping_version':ss.mapping_version,'mapping_hash':ss.mapping_hash},coords,model);annotate_denominators(result,'original',c)
  raw.reset_index(names='cell_id').to_parquet(out/'raw_cell_program_scores_for_qc.parquet',index=False)
  if seed_qc_run:
   e=dict(coverage=pd.read_csv(out/'program_seed_coverage.tsv',sep='\t'));_annotate(result,'original',c,raw,e,dt,man,out)
  else:
   frame=pd.read_parquet(out/'cell_program_representation.parquet');frame['seed_detection_evidence']='MISSING; Original S masked toNA';frame.to_parquet(out/'cell_program_representation.parquet',index=False)
   result['qc']['seed_detection_evidence']='MISSING; Original S not certified';write_json(result['qc'],out/'run_qc.json')
   from .missing_qc import export_unknown_seed_qc
   export_unknown_seed_qc(b,out,'original',raw)
  _score_only_summaries(result,out,'original');return result
 except Exception as ex:_failure(out,ex);raise

def _score_only_summaries(result,out,version):
 b=result['bundle'];rows=[]
 for (sample,identity),g in b.metadata.groupby(['sample_id','coarse_identity'],sort=False,dropna=False):
  for p in ORDER:
   n=int(b.scores.loc[g.index,p].notna().sum());rows.append(dict(sample_id=sample,coarse_identity=identity,signature_version=version,program_id=p,n_target_cells=len(g),n_valid=n,fraction_valid=n/len(g)))
 pd.DataFrame(rows).to_csv(out/'score_eligibility_by_sample.tsv',sep='\t',index=False)
 pd.DataFrame([dict(program_id=p,input_signature_version=version,n_valid=int(b.scores[p].notna().sum()),cross_version_comparison='NOT_EVALUABLE_SINGLE_VERSION_SCORE_ONLY') for p in ORDER]).to_csv(out/'scoring_version_comparison.tsv',sep='\t',index=False)

def analyze_score_run(input_run,output,config=None):
 """Both/single score-only replay with exact previously exported SeedQC evidence."""
 from .workflow import _new_output,load_config,_input_manifest
 src=Path(input_run);c=load_config(config or src/'run_config.yaml');root=_new_output(output);start=time.perf_counter();qcroot,manifest=verify_seed_evidence(src)
 for name in [*manifest['file_sha256'],'seed_qc_manifest.json']:shutil.copyfile(qcroot/name,root/name)
 versions=['original','expanded'] if c.get('signature_version')=='both' else [c.get('signature_version','expanded')];results={}
 for v in versions:
  source=src/v if len(versions)>1 else src;cfg={**c,'signature_version':v};meta=source/'cell_metadata.parquet'
  if c['modality']=='cell_resolved_spatial':
   m=pd.read_parquet(meta).set_index('cell_id');xy=pd.read_parquet(source/'cell_coordinates.parquet').set_index('cell_id');m=m.join(xy);cfg['score_coordinate_columns']=list(xy.columns);meta=root/(v+'_score_input_metadata.parquet');m.reset_index().to_parquet(meta,index=False)
  dest=root/v if len(versions)>1 else root/'scored_version'
  results[v]=analyze_scores(source/'cell_program_scores.parquet',meta,source/'score_spec.json',source/'program_coverage.tsv',cfg,dest,source/'detected_members.parquet',source/'frozen_representation.json',qcroot if v=='original' else None)
 e=dict(coverage=pd.read_csv(root/'program_seed_coverage.tsv',sep='\t'));dt=pd.read_parquet(root/'cell_seed_detection.parquet');raw=dt.pivot(index='cell_id',columns='program_id',values='raw_score_for_qc');cells=pd.read_parquet(root/'seed_cell_index.parquet').cell_id;raw=raw.reindex(index=cells,columns=ORDER);raw.index.name='cell_id'
 if len(versions)>1:
  _finish_run(root,c,results,{v:r['bundle'] for v,r in results.items()},raw,e,dt,manifest,dict(n_cells=len(raw)),_input_manifest([src/'cell_program_scores.parquet',src/'score_spec.json',src/'seed_qc_manifest.json']),start)
 else:
  for p in Path(results[versions[0]]['output']).iterdir():
   if p.is_file() and not (root/p.name).exists():shutil.copyfile(p,root/p.name)
 return dict(output=str(root.resolve()),qc=json.loads((root/'run_qc.json').read_text()),versions=results)
