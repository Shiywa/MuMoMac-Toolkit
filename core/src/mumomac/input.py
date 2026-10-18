"""Explicit matrix import and count/QC characterization."""
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,io
from .validation import prepare_input
from .io import read_adata

def from_matrix(matrix,*,genes,cells,obs=None,var=None,expression_type,umi_counts=False):
    if expression_type not in ('counts','log1p','nonnegative_expression'):raise ValueError('Declare expression_type')
    if len(set(cells))!=len(cells):raise ValueError('Cell IDs must be unique')
    if matrix.shape!=(len(cells),len(genes)):raise ValueError('Expected cells by genes matrix')
    if umi_counts and expression_type!='counts':raise ValueError('UMI declaration needs integer counts')
    a=ad.AnnData(sparse.csr_matrix(matrix),obs=pd.DataFrame(index=pd.Index(cells)) if obs is None else obs.copy(),var=pd.DataFrame(index=pd.Index(genes)) if var is None else var.copy())
    if list(a.obs_names)!=list(cells) or list(a.var_names)!=list(genes):raise ValueError('Metadata IDs/order differ from matrix IDs')
    if expression_type=='counts':a.layers['counts']=a.X.copy()
    a.uns['input_contract']={'expression_type':expression_type,'umi_counts_verified':bool(umi_counts),'source':'explicit caller declaration'}
    prepare_input(a,layer='counts' if expression_type=='counts' else 'X',expression_type=expression_type)
    return a

def read_matrix_market(matrix_path,genes_path,cells_path,*,orientation,expression_type,umi_counts=False):
    if orientation not in ('genes_by_cells','cells_by_genes'):raise ValueError('Declare Matrix Market orientation')
    genes=pd.read_csv(genes_path,sep='\t',index_col=0);cells=pd.read_csv(cells_path,sep='\t',index_col=0)
    x=io.mmread(matrix_path).tocsr()
    if orientation=='genes_by_cells':x=x.T.tocsr()
    return from_matrix(x,genes=genes.index.astype(str).tolist(),cells=cells.index.astype(str).tolist(),obs=cells,var=genes,expression_type=expression_type,umi_counts=umi_counts)

def profile_input(adata,*,layer='counts',expression_type=None,gene_symbol_key=None,duplicate_policy='error',sample_key=None,study_key=None,chunk_size=500):
    p=prepare_input(adata,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key,duplicate_policy=duplicate_policy,chunk_size=chunk_size)
    totals=np.asarray(p.matrix.sum(1)).ravel();det=np.asarray((p.matrix>0).sum(1)).ravel()
    symbols=np.array(p.symbols);mt=np.array([str(g).lower().startswith('mt-') for g in symbols]);ribo=np.array([str(g).lower().startswith(('rpl','rps')) for g in symbols])
    def fraction(mask):
        if not mask.any():return np.full(p.adata.n_obs,np.nan)
        return np.divide(np.asarray(p.matrix[:,mask].sum(1)).ravel(),totals,out=np.full(p.adata.n_obs,np.nan),where=totals>0)
    cells=pd.DataFrame({'cell_id':p.adata.obs_names,'input_expression_total':totals,'detected_genes':det,'mitochondrial_fraction':fraction(mt),'ribosomal_fraction':fraction(ribo)})
    if sample_key:
        if sample_key not in p.adata.obs:raise KeyError(sample_key)
        cells[sample_key]=p.adata.obs[sample_key].to_numpy()
        grouping=[sample_key]
        if study_key and study_key!=sample_key:
            if study_key not in p.adata.obs:raise KeyError(study_key)
            cells[study_key]=p.adata.obs[study_key].to_numpy();grouping=[study_key,sample_key]
        samples=cells.groupby(grouping,dropna=False).agg(n_cells=('cell_id','size'),median_expression_total=('input_expression_total','median'),median_detected_genes=('detected_genes','median'),mitochondrial_fraction=('mitochondrial_fraction','median'),ribosomal_fraction=('ribosomal_fraction','median')).reset_index()
    else:samples=pd.DataFrame()
    summary=dict(p.provenance,expression_total_unit='UMI only when explicitly verified; otherwise declared expression units',median_expression_total=float(np.median(totals)),median_detected_genes=float(np.median(det)),mitochondrial_genes_measured=int(mt.sum()),ribosomal_genes_measured=int(ribo.sum()))
    return {'summary':summary,'cell_qc':cells,'sample_qc':samples}

def validate_input_metadata(adata,*,sample_key=None,biological_unit_key=None,condition_key=None,identity_key=None,study_key=None):
    a=read_adata(adata);rows=[]
    for role,key in [('sample',sample_key),('biological_unit',biological_unit_key),('condition',condition_key),('identity',identity_key),('study',study_key)]:
        present=key is not None and key in a.obs
        rows.append(dict(role=role,key=key or '',status='PRESENT' if present else 'MISSING',missing_values=int(a.obs[key].isna().sum()) if present else a.n_obs,blocked_analysis='formal inference' if role=='biological_unit' and not present else ''))
    if sample_key and sample_key in a.obs:
        group=[sample_key]+([study_key] if study_key and study_key!=sample_key else [])
        for key in [condition_key,biological_unit_key]:
            if key and key in a.obs and a.obs.groupby(group,observed=True,dropna=False)[key].nunique(dropna=False).gt(1).any():raise ValueError(f'Sample maps to multiple {key}; pooled mapping must be resolved')
    return pd.DataFrame(rows)
