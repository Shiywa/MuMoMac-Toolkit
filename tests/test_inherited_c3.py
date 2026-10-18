import copy,json

import numpy as np,pandas as pd,pytest,anndata as ad

from scipy import sparse

from mumomac_experimental.input import validate_input,select_momac,harmonize_genes,load_input

from mumomac_experimental.workflow import normalize_metadata,representation_config,analyze_scores,_analysis_bundle

from mumomac_coupling.adapters import ScoreBundle,ScoreSpec

from mumomac_coupling.registry import ORDER

def fixture():
 a=ad.AnnData(sparse.csr_matrix([[1,2,0],[0,3,4],[2,0,0]]),obs=pd.DataFrame({'label':['Macrophage','T','Monocyte'],'sample':['a','a','b']},index=['c1','c2','c3']),var=pd.DataFrame(index=['Cd14','Csf1r','Lyz2']));a.layers['counts']=a.X.copy()
 c=dict(species='mouse',modality='scRNA',layer='counts',expression_type='counts',sample_column='sample',cell_type_column='label',target_labels=['Macrophage','Monocyte'])
 return a,c

def test_provided_exact_selection_universe_and_annotation():
 a,c=fixture();s,t=select_momac(a,c);assert list(s.obs_names)==['c1','c3'];assert s.n_vars==a.n_vars;assert s.obs.label.tolist()==['Macrophage','Monocyte'];assert t.selected.sum()==2;assert a.n_obs==3

@pytest.mark.parametrize('field,value',[('species','rat'),('modality','VisiumUnknown'),('layer',None),('expression_type',None),('sample_column',None),('layer','corrected_absent'),('expression_type','scaled')])
def test_invalid_input_contract(field,value):
 a,c=fixture();c[field]=value
 with pytest.raises(ValueError):validate_input(a,c)

@pytest.mark.parametrize('which',['duplicate_cell','duplicate_gene','negative','noninteger','nan','inf','sample_missing','conflicting_cell_id'])
def test_invalid_matrix(which):
 a,c=fixture()
 if which=='duplicate_cell':a.obs_names=['c1','c1','c3']
 elif which=='duplicate_gene':a.var_names=['A','A','B']
 elif which=='sample_missing':a.obs.loc['c1','sample']=None
 elif which=='conflicting_cell_id':a.obs['cell_id']=['x','y','z']
 else:
  x=a.layers['counts'].astype(float);x.data[0]={'negative':-1,'noninteger':1.5,'nan':np.nan,'inf':np.inf}[which];a.layers['counts']=x
 with pytest.raises(ValueError):validate_input(a,c)

@pytest.mark.parametrize('which',['missing_column','missing_targets','empty_targets','auto'])
def test_selection_refusals(which):
 a,c=fixture()
 if which=='missing_column':c['cell_type_column']='absent'
 elif which=='missing_targets':c.pop('target_labels')
 elif which=='empty_targets':c['target_labels']=['never_seen']
 else:c['selection_mode']='auto'
 with pytest.raises(ValueError):select_momac(a,c)

def test_metadata_exact_join_reorders_not_loses(tmp_path):
 a,c=fixture();p=tmp_path/'x.h5ad';a.write_h5ad(p);m=a.obs.iloc[::-1].reset_index(names='cell_id');m.to_csv(tmp_path/'meta.csv',index=False);b=load_input(p,tmp_path/'meta.csv');assert list(b.obs_names)==list(a.obs_names);assert b.obs.label.equals(a.obs.label)

def test_metadata_missing_id_refused(tmp_path):
 a,c=fixture();p=tmp_path/'x.h5ad';a.write_h5ad(p);pd.DataFrame({'cell_id':['c1','c2']}).to_csv(tmp_path/'m.csv',index=False)
 with pytest.raises(ValueError):load_input(p,tmp_path/'m.csv')

def test_text_ids_preserve_leading_zeros(tmp_path):
 from mumomac_experimental.input import read_table
 p=tmp_path/'m.tsv';p.write_text('cell_id\tsample_id\n001\t007\n');m=read_table(p);assert m.cell_id.iloc[0]=='001' and m.sample_id.iloc[0]=='007'

def test_no_sample_to_animal_inference():
 a,c=fixture();m=normalize_metadata(a.obs,c);assert not m.biological_unit_confirmed.any();assert 'biological_unit_id' not in m

def test_native_symbols_not_case_converted():
 a,c=fixture();b,genes,mp,key=harmonize_genes(a,c);assert genes.scoring_symbol.tolist()==['Cd14','Csf1r','Lyz2'];assert mp['case_conversion'] is False;assert key is None

def test_ensembl_subset_retains_all_background():
 a,c=fixture();a.var_names=['ENSMUSG00000051439.1','ENSMUSG99999999999','opaque'];c['gene_id_type']='ensembl';b,t,mp,key=harmonize_genes(a,c);assert b.n_vars==3;assert t.scoring_symbol.iloc[0]=='Cd14';assert t.scoring_symbol.iloc[1:].isna().all();assert mp['unmapped_native_genes']==2;assert key=='mumomac_scoring_symbol'

def test_ensembl_not_guessed():
 a,c=fixture();a.var_names=['ENSMUSG00000051439','other','third']
 with pytest.raises(ValueError):harmonize_genes(a,c)

def test_duplicate_mapped_symbols_refused():
 a,c=fixture();a.var['symbol']=['Cd14','Cd14','Lyz2'];c['gene_symbol_key']='symbol'
 with pytest.raises(ValueError):harmonize_genes(a,c)

def test_unknown_symbols_opaque_not_dropped():
 a,c=fixture();a.var['symbol']=['Cd14',None,'Lyz2'];c['gene_symbol_key']='symbol';b,t,mp,key=harmonize_genes(a,c);assert b.n_vars==3;assert mp['unmapped_native_genes']==1

def test_mixed_spot_is_refused_before_output():
 a,c=fixture();c['modality']='mixed_spots'
 with pytest.raises(ValueError,match='Mixed spots'):validate_input(a,c)

@pytest.mark.parametrize('defect',['no_segmentation','no_coordinates','bad_coordinates'])
def test_spatial_guards(defect):
 a,c=fixture();c.update(modality='cell_resolved_spatial',segmentation_confirmed=True);a.obsm['spatial']=np.ones((3,2))
 if defect=='no_segmentation':c['segmentation_confirmed']=False
 elif defect=='no_coordinates':del a.obsm['spatial']
 else:a.obsm['spatial'][0,0]=np.nan
 with pytest.raises(ValueError):validate_input(a,c)

def test_spatial_coordinates_retained():
 a,c=fixture();c.update(modality='cell_resolved_spatial',segmentation_confirmed=True);a.obsm['spatial']=np.array([[1,2],[3,4],[5,6]]);validate_input(a,c);s,_=select_momac(a,c);assert np.array_equal(s.obsm['spatial'],[[1,2],[5,6]])

def test_defaults_continuous_S_only_no_soft_no_tail():
 r=representation_config({});assert r['optional_lifecycle']==[];assert r['default_mode']=='auto';assert r['joint_index']['default_output'] is False;assert r['tail_profile']['enabled'] is False

@pytest.mark.parametrize('pair',['C02','A01__S01','unknown'])
def test_pairs_need_valid_scope(pair):
 with pytest.raises(ValueError):representation_config({'selected_pairs':[pair]})

@pytest.mark.parametrize('minimum_blocks',[1,4])
def test_distribution_gates_not_tunable_in_integration(minimum_blocks):
 with pytest.raises(ValueError):representation_config({'representation':{'distribution':{'min_calibration_blocks':minimum_blocks}}})

def test_explicit_non_symbol_background_preserved():
 a,c=fixture();a.var_names=['Cd14','MHC I','Lyz2'];c['unmapped_gene_ids']=['MHC I'];b,t,mp,key=harmonize_genes(a,c);assert b.n_vars==3;assert pd.isna(t.scoring_symbol.iloc[1]);assert key=='mumomac_scoring_symbol'

