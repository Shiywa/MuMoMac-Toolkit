"""AnnData/TSV IO with explicit no-overwrite behavior."""
from pathlib import Path
import json,os,tempfile
import anndata as ad
import pandas as pd

def read_adata(value):
    """Return AnnData; objects are returned unchanged and paths are loaded in memory."""
    if isinstance(value,ad.AnnData):return value
    if isinstance(value,(str,Path)):
        p=Path(value)
        if p.suffix!='.h5ad':raise ValueError('Input file must be an AnnData .h5ad file.')
        return ad.read_h5ad(p)
    raise TypeError('Expected AnnData or a .h5ad path.')

def _commit(path,writer,overwrite=False):
    p=Path(path);p.parent.mkdir(parents=True,exist_ok=True)
    if p.exists() and not overwrite:raise FileExistsError(f'Output exists: {p}; explicit overwrite is required.')
    fd,name=tempfile.mkstemp(prefix='.mumomac-',suffix=p.suffix,dir=p.parent);os.close(fd)
    t=Path(name)
    try:
        writer(t)
        if overwrite:os.replace(t,p)
        else:os.link(t,p);t.unlink() # atomic exclusive destination creation
    finally:
        if t.exists():t.unlink()

def write_adata(adata,path,overwrite=False):
    _commit(path,lambda p:adata.write_h5ad(p,compression='gzip'),overwrite)

def write_table(table,path,overwrite=False):
    meta=Path(str(path)+'.metadata.json')
    if table.attrs and meta.exists() and not overwrite:raise FileExistsError(f'Output exists: {meta}')
    _commit(path,lambda p:table.to_csv(p,sep='\t',index=False,na_rep='Missing'),overwrite)
    if table.attrs:_commit(meta,lambda p:p.write_text(json.dumps(table.attrs,ensure_ascii=False,indent=2)+'\n'),overwrite)

def read_table(path):
    # Empty gene lists mean no missing members; they are distinct from unknown.
    # Only the writer's explicit Missing sentinel is parsed as a missing value.
    d=pd.read_csv(path,sep='\t',na_values=['Missing'],keep_default_na=False);meta=Path(str(path)+'.metadata.json')
    if meta.exists():d.attrs.update(json.loads(meta.read_text()))
    return d
