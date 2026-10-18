from pathlib import Path
import re,json,hashlib
import numpy as np,pandas as pd
from scipy import sparse
import anndata as ad

MODALITIES={'scRNA','scRNA-seq','snRNA','snRNA-seq','cell_resolved_spatial','mixed_spots'}
def read_table(p):
 p=Path(p)
 # Text metadata identifiers are lexical (e.g.001), never guessed numeric IDs.
 return pd.read_parquet(p) if p.suffix=='.parquet' else pd.read_csv(p,sep=',' if p.suffix=='.csv' else '\t',dtype=object)
def load_input(path,metadata=None):
 if Path(path).suffix!='.h5ad':raise ValueError('Native expression entry currently requires H5AD; use the documented conversion route.')
 a=ad.read_h5ad(path)
 if metadata:
  m=read_table(metadata)
  if 'cell_id' not in m:raise ValueError('Metadata needs explicit cell_id.')
  m=m.set_index('cell_id');m.index=m.index.astype(str)
  if m.index.has_duplicates or set(m.index)!=set(a.obs_names):raise ValueError('Exact unique cell_id match required; no inner join or silent row loss.')
  # Explicit supplied metadata takes precedence, while original annotations remain available.
  for c in a.obs:
   if c in m and not a.obs[c].astype(str).equals(m.loc[a.obs_names,c].astype(str)):a.obs['input_obs__'+c]=a.obs[c]
  for c in m:a.obs[c]=m.loc[a.obs_names,c].to_numpy()
 return a
def validate_input(a,c):
 if c.get('species') not in {'mouse','human'}:raise ValueError('Explicit mouse/human species required.')
 if c.get('modality') not in MODALITIES:raise ValueError('Declare supported modality explicitly.')
 if not c.get('layer') or not c.get('expression_type'):raise ValueError('Explicit layer and expression_type required; corrected/normalized layers are never guessed.')
 if a.n_obs<1 or a.n_vars<1 or not a.obs_names.is_unique or not a.var_names.is_unique:raise ValueError('Nonempty cell-by-gene matrix with unique native cell/gene IDs required.')
 if 'cell_id' in a.obs and not np.array_equal(a.obs.cell_id.astype(str),a.obs_names.astype(str)):raise ValueError('obs cell_id conflicts with matrix observation IDs.')
 layer=c['layer'];x=a.X if layer=='X' else a.layers.get(layer)
 if x is None:raise ValueError('Selected expression layer absent.')
 if c['expression_type'] not in {'counts','log1p','nonnegative_expression'}:raise ValueError('Unsupported expression type; scaled/corrected data cannot be assumed suitable.')
 if x.shape!=(a.n_obs,a.n_vars) or np.dtype(x.dtype).kind not in 'biuf':raise ValueError('Matrix orientation/dimensions/type invalid.')
 for start in range(0,a.n_obs,2048):
  b=x[start:start+2048];v=b.data if sparse.issparse(b) else np.asarray(b)
  if not np.isfinite(v).all() or (v<0).any():raise ValueError('Expression must be finite and nonnegative.')
  if c['expression_type']=='counts' and (v!=np.rint(v)).any():raise ValueError('Declared counts must be integer-valued.')
 sample=c.get('sample_column')
 if not sample or sample not in a.obs or a.obs[sample].isna().any() or a.obs[sample].astype(str).isin(['','Unknown','nan']).any():raise ValueError('An explicit complete sample column is required; sample is not inferred as animal.')
 if c['modality']=='cell_resolved_spatial':
  if c.get('segmentation_confirmed') is not True:raise ValueError('Cell-resolved spatial requires explicit segmentation provenance confirmation.')
  xy=a.obsm.get(c.get('coordinate_key','spatial'))
  if xy is None or xy.shape[0]!=a.n_obs or xy.shape[1]<2 or not np.isfinite(xy).all():raise ValueError('Finite per-cell coordinates required for cell-resolved spatial.')
 if c['modality']=='mixed_spots':raise ValueError('Mixed spots are not single-cell Mo/Mac: this target-specific workflow refuses them. Use an explicitly spot-level scoring analysis; no cell Joint is emitted.')
 return dict(status='PASS',n_cells=a.n_obs,n_native_genes=a.n_vars,sparse=bool(sparse.issparse(x)),orientation='cells_by_genes',layer=layer,expression_type=c['expression_type'],expression_provenance=c.get('expression_provenance','CALLER_DECLARED_NOT_INDEPENDENTLY_CERTIFIED'))
def select_momac(a,c):
 if c.get('selection_mode','provided')!='provided':raise ValueError('Auto-assisted remains a separate frozen upstream workflow; provide its audited annotation explicitly.')
 col=c.get('cell_type_column');labels=c.get('target_labels')
 if not col or col not in a.obs or not labels:raise ValueError('Provided selection requires explicit cell_type_column and nonempty target_labels.')
 keep=a.obs[col].isin(labels).to_numpy()
 audit=pd.DataFrame({'cell_id':a.obs_names,'original_annotation':a.obs[col].astype(str).to_numpy(),'selected':keep,'selection_reason':np.where(keep,'USER_PROVIDED_TARGET_LABEL','NON_TARGET_OR_MISSING_LABEL')})
 if not keep.any():raise ValueError('Provided target labels select no cells.')
 return a[keep,:].copy(),audit
def harmonize_genes(a,c):
 """Native universe retained; only explicit, unambiguous symbol assignments used."""
 from importlib import resources
 root=resources.files('mumomac.resources');key=c.get('gene_symbol_key');mode=c.get('gene_id_type','symbol')
 if key and key not in a.var:raise ValueError('Gene symbol column not found.')
 native=a.var_names.astype(str);raw=a.var[key] if key else pd.Series(native,index=a.var_names)
 symbols=raw.copy();status=np.full(a.n_vars,'EXPLICIT_NATIVE_SYMBOL',object)
 if mode=='ensembl':
  if key:raise ValueError('Ensembl mode uses native gene IDs; use symbol mode for an explicit author symbol column.')
  with root.joinpath('human_one_to_one.tsv').open() as f:mp=pd.read_csv(f,sep='\t')
  idcol='mouse_ensembl_id' if c['species']=='mouse' else 'human_ensembl_id';symcol='mouse_gene' if c['species']=='mouse' else 'human_gene'
  pairs=mp[[idcol,symcol]].dropna().drop_duplicates();bad=pairs.groupby(idcol)[symcol].nunique();bad=set(bad[bad>1].index)
  pairs=pairs[~pairs[idcol].isin(bad)];d=dict(zip(pairs[idcol],pairs[symcol]));symbols=pd.Series([d.get(re.sub(r'\.\d+$','',v),np.nan) for v in native],index=a.var_names);status[:]='FROZEN_ENSEMBL_SUBSET_MAPPING'
 elif mode!='symbol':raise ValueError('gene_id_type must be symbol or ensembl.')
 elif key is None and any(re.fullmatch(r'ENS[A-Z]*G\d+(\.\d+)?',v) for v in native):raise ValueError('Ensembl IDs need explicit gene_id_type=ensembl or an author symbol column.')
 s=pd.Series(symbols.to_numpy(),index=a.var_names);bad=s.notna()&s.astype(str).duplicated(keep=False)
 explicitly_unmapped=c.get('unmapped_gene_ids',[])
 if set(explicitly_unmapped)-set(native):raise ValueError('Explicit unmapped gene IDs absent from native axis.')
 if explicitly_unmapped:s.loc[explicitly_unmapped]=np.nan
 if bad.any():raise ValueError('Duplicate mapped symbols: resolve explicitly upstream; no arbitrary ortholog edge choice.')
 missing=s.isna()|s.astype(str).eq('');status[missing]='UNMAPPED_OPAQUE_BACKGROUND_RETAINED'
 audit=pd.DataFrame({'native_gene_id':native,'scoring_symbol':s.to_numpy(),'status':status})
 a.var['mumomac_scoring_symbol']=s.to_numpy()
 manifest={'species':c['species'],'native_gene_axis_retained':True,'native_gene_count':a.n_vars,'symbol_source':key or mode,'unmapped_native_genes':int(missing.sum()),'conflicts_selected':0,'case_conversion':False,'database_release_status':'MISSING_LEGACY_S_MAPPING' if c['species']=='human' else 'MOUSE_NATIVE; optional partial Ensembl map frozen snapshot','human_support':'EXPERIMENTAL' if c['species']=='human' else 'MOUSE','resources_sha256':{n:hashlib.sha256(root.joinpath(n).read_bytes()).hexdigest() for n in ['human_one_to_one.tsv','human_la_one_to_one_v1.tsv','signatures_manifest.json']}}
 # Preserve exact official symbol parameter for native-symbol inputs (parity hashes).
 scorer_key=None if key is None and mode=='symbol' and not missing.any() else 'mumomac_scoring_symbol'
 if key is not None and mode=='symbol' and not explicitly_unmapped:scorer_key=key
 return a,audit,manifest,scorer_key
