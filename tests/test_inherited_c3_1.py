import json,copy

import numpy as np,pandas as pd,pytest,anndata as ad

from scipy import sparse

from mumomac_coupling.registry import ORDER,digest

from mumomac_experimental import analyze,seed_details

from mumomac_experimental.integration import run_versions

from mumomac_experimental.seed_qc import seed_evidence,SeedScoreSpec,spec_from_dict,validate_policy,verify_seed_evidence

import mumomac

def fixture(k=3,n=8):
 genes=mumomac.reference_signatures('original','mouse')['S03']['genes'];g=genes[:k]+['background_'+str(i) for i in range(30)];x=np.ones((n,len(g)),dtype=np.int32)*5
 if k>=3:x[0,:k]=0;x[0,:2]=5
 a=ad.AnnData(sparse.csr_matrix(x),obs=pd.DataFrame(dict(sample_id=['sample']*n,coarse_identity=['Macrophage']*n,study_id='study',source_block_id='source',tissue='Unknown',condition='Unknown',biological_unit_id='Unknown',biological_unit_confirmed=False),index=pd.Index(['cell'+str(i) for i in range(n)],name='cell_id')),var=pd.DataFrame(index=g));a.layers['counts']=a.X.copy()
 c=dict(species='mouse',modality='scRNA-seq',layer='counts',expression_type='counts',sample_column='sample_id',cell_type_column='coarse_identity',target_labels=['Macrophage'],selected_pairs=['C01'],representation={'default_mode':'continuous_only'},signature_version='both');return a,c

def test_three_measured_two_detected_is_NA_three_detected_valid():
 a,c=fixture();b,raw,e,_=run_versions(a,a.obs,c,None);v=b['original'].scores.S03;assert pd.isna(v.iloc[0]);assert v.iloc[1:].notna().all();assert np.isfinite(raw.S03.iloc[0]);assert (e['counts'].S03.iloc[:2].to_numpy()==[2,3]).all()

def test_two_measured_is_not_formally_scorable():
 a,c=fixture(2);b,raw,e,_=run_versions(a,a.obs,c,None);assert b['original'].scores.S03.isna().all();assert raw.S03.notna().all();assert not b['original'].coverage.loc['S03','scoring_eligible']

def test_one_missing_does_not_destroy_whole_original():
 genes=mumomac.reference_signatures('original','mouse')['S03']['genes'];a,c=fixture(len(genes)-1);b,raw,e,ref=run_versions(a,a.obs,c,None);assert np.isnan(ref.obsm['X_mumomac_reference_original'][:,2]).all();assert b['original'].scores.S03.iloc[1:].notna().all();assert e['coverage'].set_index('program_id').loc['S03','definition_namespace']=='CURRENT_ORIGINAL_PARTIAL_SEED3'

def test_full_original_raw_exact_with_official():
 n=len(mumomac.reference_signatures('original','mouse')['S03']['genes']);a,c=fixture(n);b,raw,e,ref=run_versions(a,a.obs,c,None);assert np.array_equal(raw.S03.to_numpy(),ref.obsm['X_mumomac_reference_original'][:,2],equal_nan=True);ok=b['original'].scores.S03.notna();assert np.array_equal(b['original'].scores.loc[ok,'S03'],raw.loc[ok,'S03'])

def test_expanded_formal_scores_unmodified():
 a,c=fixture();b,raw,e,ref=run_versions(a,a.obs,c,None);expected=mumomac.score_signature_versions(a,species='mouse',signature_version='expanded',layer='counts',expression_type='counts',coverage_policy='reference',copy=True);assert np.array_equal(b['expanded'].scores.to_numpy(),expected.obsm['X_mumomac_reference_expanded'],equal_nan=True)

def test_not_measured_and_measured_zero_distinct(tmp_path):
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);r=analyze(p,c,tmp_path/'out');detail=seed_details(tmp_path/'out','cell0','S03');assert len(detail['detected_genes'])==2;assert len(detail['measured_but_zero_genes'])==1;assert len(detail['not_measured_genes'])>0;assert not detail['unmapped_genes']

def test_query_unknown_cell_or_program(tmp_path):
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);analyze(p,c,tmp_path/'out')
 with pytest.raises(ValueError):seed_details(tmp_path/'out','not_a_cell','S03')
 with pytest.raises(ValueError):seed_details(tmp_path/'out','cell0','S99')

def test_human_unmapped_is_not_measured_zero():
 a,c=fixture();c['species']='human';a.var_names=mumomac.reference_signatures('original','human')['S03']['genes'][:3]+['human_fixture_background_'+str(i) for i in range(a.n_vars-3)]
 e=seed_evidence(a,a.obs,c);cat=e['catalog'];un=cat[cat.mapping_status=='UNMAPPED'];assert len(un)>0;assert un.n_cells_detected.isna().all();measured=cat[cat.measured_in_dataset];assert measured.n_cells_detected.notna().all();cov=e['coverage'];assert (cov.original_seed_count==cov.mapped_seed_count+cov.mapping_loss_count).all();assert (cov.mapped_seed_count==cov.measured_seed_count+cov.absent_seed_count).all()

def test_counts_nonzero_priority_over_normalized_X():
 a,c=fixture();c.update(layer='X',expression_type='log1p');a.X=a.X.astype(float)+sparse.csr_matrix(np.ones(a.shape))
 e=seed_evidence(a,a.obs,c);assert e['detection_layer']=='counts';assert e['counts'].S03.iloc[0]==2

def test_normalized_only_requires_zero_semantics():
 a,c=fixture();del a.layers['counts'];c.update(layer='X',expression_type='log1p')
 with pytest.raises(ValueError):seed_evidence(a,a.obs,c)
 c['detection_zero_preserving']=True;assert seed_evidence(a,a.obs,c)['detection_evidence'].startswith('DECLARED')

def test_valid_zero_score_is_not_called_high():
 a,c=fixture();# fill background with high expression; the minimum gate is not an activation call
 a.layers['counts']=a.layers['counts'].copy();a.layers['counts'][:,3:]=100
 b,raw,e,_=run_versions(a,a.obs,c,None);assert b['original'].scores.S03.iloc[1:].notna().all()

def test_seed_partial_namespace_and_hash_cannot_be_forged():
 a,c=fixture();b,*_=run_versions(a,a.obs,c,None);sp=b['original'].spec.to_dict();sp['namespace']='CURRENT_UCELL_V1_ORIGINAL'
 with pytest.raises(ValueError):spec_from_dict(sp)
 sp=b['original'].spec.to_dict();sp['parameters']['effective_program_members']['S03'][0]='NOT_A_FROZEN_SEED'
 with pytest.raises(ValueError):spec_from_dict(sp).validate()

@pytest.mark.parametrize('key,value',[('signature_version','legacy'),('original_min_present',2),('original_min_present',4),('original_min_detected_per_cell',2),('original_min_detected_per_cell',4),('scoring_qc_mode','ignore'),('original_coverage_policy','available')])
def test_invalid_seed_policy(key,value):
 with pytest.raises(ValueError):validate_policy({key:value})

def test_both_equals_independent_versions(tmp_path):
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);both=analyze(p,c,tmp_path/'both')
 for version in ['original','expanded']:
  one=analyze(p,{**c,'signature_version':version},tmp_path/version);assert np.array_equal(both['versions'][version]['bundle'].scores.to_numpy(),one['bundle'].scores.to_numpy(),equal_nan=True)

def test_no_NA_enters_neither_and_denominator_exact(tmp_path):
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);r=analyze(p,c,tmp_path/'out');s=r['versions']['original']['continuous'];b=r['versions']['original']['bundle'];pair=pd.read_parquet(tmp_path/'out/original/selected_pair_joint_features.parquet');bad=~np.isfinite(pair.raw_A)|~np.isfinite(pair.raw_B);assert pair.loc[bad,'binary_pair_state'].str.contains('Unassessable').all();assert not pair.loc[bad,'binary_pair_state'].eq('Neither').any()
 for row in s.itertuples():
  n=int((b.scores[row.program_A].notna()&b.scores[row.program_B].notna()).sum());assert row.n_pair_valid==n;assert row.n_target_cells==len(a)

def test_audit_never_upgrades_two_detected(tmp_path):
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);r=analyze(p,{**c,'signature_version':'original','scoring_qc_mode':'audit'},tmp_path/'out');assert pd.isna(r['bundle'].scores.S03.iloc[0]);assert (tmp_path/'out/audit_raw_vs_valid.tsv').exists()

@pytest.mark.parametrize('mode',['auto','force_gmm','force_tail','relative_high','continuous_only'])
def test_manual_modes_do_not_rescue_NA(tmp_path,mode):
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);r=analyze(p,{**c,'signature_version':'original','representation':{'default_mode':mode}},tmp_path/'out');f=pd.read_parquet(tmp_path/'out/cell_program_representation.parquet');row=f[(f.cell_id=='cell0')&(f.program_id=='S03')].iloc[0];assert row.state_code==3 and pd.isna(row.discrete_or_relative_state)

def test_bitmaps_have_no_all_cell_gene_text(tmp_path):
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);analyze(p,c,tmp_path/'out');bits=np.load(tmp_path/'out/seed_detection_bits.npz')['packed'];assert bits.shape[0]==len(a);assert bits.shape[1]<=a.n_vars;assert (tmp_path/'out/program_seed_gene_qc.tsv').stat().st_size<1000000

def test_seed_qc_tamper_rejected(tmp_path):
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);analyze(p,c,tmp_path/'out');f=tmp_path/'out/program_seed_gene_qc.tsv';f.write_text(f.read_text()+'\n')
 with pytest.raises(ValueError):verify_seed_evidence(tmp_path/'out')

def test_score_only_missing_detection_never_certified(tmp_path):
 from mumomac_experimental import analyze_scores
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);src=tmp_path/'expr';analyze(p,c,src);src=src/'original';out=tmp_path/'scores'
 analyze_scores(src/'raw_cell_program_scores_for_qc.parquet',src/'cell_metadata.parquet',src/'score_spec.json',src/'program_coverage.tsv',{**c,'signature_version':'original'},out)
 valid=pd.read_parquet(out/'cell_program_scores.parquet');assert valid[list(ORDER[:16])].isna().all().all();dt=pd.read_parquet(out/'cell_seed_detection.parquet');assert dt.n_seed_detected.isna().all();raw=pd.read_parquet(out/'raw_cell_program_scores_for_qc.parquet');assert raw.S03.notna().all()

def test_trusted_both_score_only_replay_exact(tmp_path):
 from mumomac_experimental import analyze_score_run
 a,c=fixture();p=tmp_path/'a.h5ad';a.write_h5ad(p);src=tmp_path/'expr';analyze(p,c,src);out=tmp_path/'replay';analyze_score_run(src,out)
 for v in ['original','expanded']:
  for name in ['cell_program_scores.parquet','sample_coupling_statistics.parquet']:
   pd.testing.assert_frame_equal(pd.read_parquet(src/v/name),pd.read_parquet(out/v/name))

def test_native_background_axis_retained():
 a,c=fixture();a.layers['counts'][:,3:]=100;b,raw,e,_=run_versions(a,a.obs,c,None);small,raw_small,e_small,_=run_versions(a[:,:3].copy(),a.obs,c,None);assert e['gene_universe_hash']==b['expanded'].spec.parameters['gene_universe_hash'];assert e['gene_universe_hash']!=e_small['gene_universe_hash'];assert not np.array_equal(raw.S03.to_numpy(),raw_small.S03.to_numpy(),equal_nan=True)

def test_L_A_do_not_inherit_S_per_cell_gate():
 a,c=fixture();b,raw,e,ref=run_versions(a,a.obs,c,None);assert np.array_equal(b['original'].scores[list(ORDER[16:])].to_numpy(),ref.obsm['X_mumomac_reference_original'][:,16:],equal_nan=True)

