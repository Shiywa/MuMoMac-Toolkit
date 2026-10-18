"""Engineering input checks. Computability is not biological reliability."""
from dataclasses import dataclass
from pathlib import Path
import re,warnings
import numpy as np,pandas as pd
from scipy import sparse
from .io import read_adata
from .signatures import load_signatures,PROGRAM_ORDER

@dataclass
class PreparedInput:
    adata: object
    matrix: object
    symbols: tuple
    expression_type: str
    provenance: dict
    zero_cells: np.ndarray
    gene_detection: np.ndarray
    gene_sum: np.ndarray

def prepare_input(adata,*,layer='counts',expression_type=None,gene_symbol_key=None,
                  duplicate_policy='error',chunk_size=1000):
    origin=str(Path(adata).resolve()) if isinstance(adata,(str,Path)) else 'AnnData object'
    a=read_adata(adata)
    if a.n_obs==0 or a.n_vars==0:raise ValueError('Input must contain at least one cell and one feature.')
    if not isinstance(chunk_size,int) or chunk_size<1:raise ValueError('chunk_size must be a positive integer.')
    if layer=='X':layer=None
    if expression_type is None:
        if layer=='counts':expression_type='counts'
        else:raise ValueError('Declare expression_type for X or a non-counts layer; normalization is not guessed.')
    if expression_type not in ('counts','log1p','nonnegative_expression'):
        raise ValueError('expression_type must be counts, log1p, or nonnegative_expression; scaled data is unsupported.')
    if layer is not None and layer not in a.layers:raise ValueError(f'Expression layer {layer!r} is absent. For X explicitly declare expression_type.')
    X=a.X if layer is None else a.layers[layer]
    if getattr(X,'ndim',2)!=2 or X.shape!=(a.n_obs,a.n_vars):raise ValueError('Expression matrix dimensions do not match AnnData.')
    if np.dtype(X.dtype).kind not in 'biuf':raise ValueError('Expression matrix must contain real numeric values.')
    if gene_symbol_key is not None:
        if gene_symbol_key not in a.var:raise KeyError(f'Gene-symbol column {gene_symbol_key!r} is absent.')
        raw=a.var[gene_symbol_key].tolist()
    else:raw=a.var_names.tolist()
    symbols=[];unknown=[]
    for j,value in enumerate(raw):
        if pd.isna(value) or str(value)=='':symbols.append(f'__mumomac_unmapped_{j}');unknown.append(j);continue
        name=str(value)
        if name!=name.strip() or any(c.isspace() for c in name):raise ValueError('Gene identifiers contain whitespace; provide a cleaned, documented symbol column.')
        if gene_symbol_key is None and (re.fullmatch(r'ENS[A-Z]*G\d+(\.\d+)?',name) or name.isnumeric()):
            raise ValueError('var_names contain Ensembl/numeric identifiers; specify gene_symbol_key instead of guessing symbols.')
        symbols.append(name)
    if unknown:warnings.warn(f'{len(unknown)} features have no symbols; they remain opaque background features.',UserWarning,stacklevel=2)
    idx=pd.Index(symbols);duplicate_symbols=idx[idx.duplicated(keep=False)].unique().tolist()
    if duplicate_policy not in ('error','sum','exclude'):raise ValueError('duplicate_policy must be error, sum, or exclude.')
    if duplicate_symbols and duplicate_policy=='error':raise ValueError(f'Duplicated gene symbols ({len(duplicate_symbols)}): {duplicate_symbols[:8]}; select an explicit duplicate_policy.')
    if duplicate_symbols and duplicate_policy=='sum' and expression_type!='counts':raise ValueError('Summing duplicate features is only supported for declared counts, not normalized/log expression.')
    # Validate native values before any grouping or rank computation; no matrix mutation.
    for start in range(0,a.n_obs,chunk_size):
        block=X[start:start+chunk_size,:];values=block.data if sparse.issparse(block) else np.asarray(block)
        if not np.isfinite(values).all():raise ValueError('Expression contains NaN/Inf; scoring is refused.')
        if (values<0).any():raise ValueError('Negative expression is unsupported (e.g. z-score/scaled input).')
        if expression_type=='counts' and not np.equal(values,np.rint(values)).all():raise ValueError('Declared counts contain noninteger values; supply correct expression_type/source.')
    original_nvars=a.n_vars
    if duplicate_symbols and duplicate_policy=='exclude':
        keep=~idx.duplicated(keep=False);X=X[:,keep];symbols=list(idx[keep])
    elif duplicate_symbols and duplicate_policy=='sum':
        codes,unique=pd.factorize(idx,sort=False);G=sparse.csr_matrix((np.ones(len(codes)),(np.arange(len(codes)),codes)),shape=(len(codes),len(unique)))
        # Sparse multiplication; dense input converts to sparse rather than densifying a sparse input.
        X=sparse.csr_matrix(X)@G;symbols=list(unique)
    if len(symbols)==0:raise ValueError('No feature universe remains after duplicate exclusion.')
    detected=np.zeros(len(symbols),dtype=np.int64);gsum=np.zeros(len(symbols),dtype=np.float64);zero=np.zeros(a.n_obs,bool)
    for start in range(0,a.n_obs,chunk_size):
        b=X[start:start+chunk_size,:]
        detected+=np.asarray((b>0).sum(axis=0)).ravel().astype(np.int64);gsum+=np.asarray(b.sum(axis=0)).ravel()
        zero[start:start+b.shape[0]]=np.asarray((b>0).sum(axis=1)).ravel()==0
    prov={'source':origin,'layer':'X' if layer is None else str(layer),'expression_type':expression_type,
          'gene_symbol_key':'var_names' if gene_symbol_key is None else str(gene_symbol_key),
          'duplicate_policy':duplicate_policy,'duplicate_symbols':duplicate_symbols,
          'unknown_symbol_feature_count':len(unknown),'input_features':original_nvars,
          'effective_features':len(symbols),'input_cells':a.n_obs,'zero_expression_cells':int(zero.sum()),
          'rank_universe_policy':'Every prepared native feature including opaque unmapped background; no HVG selection.'}
    return PreparedInput(a,X,tuple(symbols),expression_type,prov,zero,detected/a.n_obs,gsum/a.n_obs)

def _program_plan(prepared,species,signature_version,coverage_policy,min_genes):
    if coverage_policy not in ('complete','available'):raise ValueError('coverage_policy must be complete or available.')
    if not isinstance(min_genes,int) or min_genes<1:raise ValueError('min_genes must be a positive integer engineering parameter.')
    registry=load_signatures(signature_version,species);index={g:i for i,g in enumerate(prepared.symbols)};rows=[];plans={}
    for pid in PROGRAM_ORDER:
        p=registry[pid];present=[g for g in p.genes if g in index];missing=[g for g in p.genes if g not in index]
        original=[g for g in p.original_genes if g in index];added=[g for g in p.added_genes if g in index]
        has_added=signature_version=='expanded' and len(p.mouse_genes)>len(load_signatures('original',species)[pid].mouse_genes)
        structural=len(present)>=min_genes and len(original)>=min_genes and (not has_added or len(added)>=min_genes)
        if coverage_policy=='complete':structural=structural and not missing
        structural=structural and not prepared.zero_cells.all()
        fr=prepared.gene_detection[[index[g] for g in present]] if present else np.array([])
        messages=[]
        if p.unmapped_genes:messages.append('Mouse members lost in one-to-one mapping; Experimental human extension.')
        if missing:messages.append('Declared mapped/signature genes absent; no zero imputation.')
        if coverage_policy=='available' and missing:messages.append('Explicit engineering reduced signature; no biological reliability certification.')
        if present and not np.any(fr>0):messages.append('Present members have no observed detection; measured zero is not functional absence.')
        if not structural:messages.append('Not Evaluable under declared engineering/input conditions.')
        if prepared.zero_cells.any():messages.append('Wholly zero-expression cells return missing scores.')
        rows.append({'program_id':pid,'program_name':p.name,'species':species,'signature_version':signature_version,
                     'definition_version':p.definition_version,'n_signature_genes':len(p.mouse_genes),
                     'n_mapped_genes':len(p.genes),'n_present_genes':len(present),
                     'mapping_retention':len(p.genes)/len(p.mouse_genes),
                     'gene_coverage_fraction':len(present)/len(p.genes) if p.genes else 0.,
                     'n_detected_genes':int(np.sum(fr>0)),
                     'mean_gene_detection_fraction':float(np.mean(fr)) if len(fr) else np.nan,
                     'median_gene_detection_fraction':float(np.median(fr)) if len(fr) else np.nan,
                     'genes_detected_in_lt1pct_cells':int(np.sum(fr<.01)),
                     'unmapped_genes':';'.join(p.unmapped_genes),'missing_genes':';'.join(missing),
                     'ambiguous_genes':';'.join(g for g in p.genes if g in prepared.provenance['duplicate_symbols']),
                     'present_genes':';'.join(present),'scoring_eligibility':bool(structural),
                     'coverage_policy':coverage_policy,'reliability_status':'UNVERIFIED: engineering computability is not biological reliability',
                     'scope':'Experimental cross-species extension' if species=='human' else 'Exploratory RNA signature',
                     'warning_messages':' '.join(messages)})
        plans[pid]={'original':original,'added':added,'has_added':has_added,'eligible':bool(structural),'present':present}
    return pd.DataFrame(rows),plans

def validate_programs(adata,*,species='mouse',signature_version='original',layer='counts',
                      expression_type=None,gene_symbol_key=None,duplicate_policy='error',
                      coverage_policy='complete',min_genes=1,chunk_size=1000):
    """Return presence, detection and declared computability; never certify reliability."""
    p=prepare_input(adata,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key,
                    duplicate_policy=duplicate_policy,chunk_size=chunk_size)
    report,_=_program_plan(p,species,signature_version,coverage_policy,min_genes)
    report.attrs['mumomac_validation']={'input_provenance':p.provenance,'biological_reliability':'Unverified'}
    return report
