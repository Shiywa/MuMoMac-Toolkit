"""Descriptive aggregation; sample IDs are never automatically biological units."""
import numpy as np,pandas as pd
from .io import read_adata
from .signatures import PROGRAM_ORDER

def summarize_programs(adata,*,sample_key,identity_key=None,tissue_key=None,group_key=None,
                       biological_unit_key=None,min_cells=1,missing_metadata='error'):
    """Long table of cell summaries per sample and optional identity/tissue.

    The biological key is carried as a user assertion, not inferred from UID.
    Statistical inference is deferred to compare_programs.
    """
    a=read_adata(adata)
    if 'X_mumomac_s16' not in a.obsm or 'mumomac' not in a.uns:raise ValueError('AnnData has no MuMoMac scores/provenance; run score_programs first.')
    if tuple(a.uns['mumomac']['program_order'])!=PROGRAM_ORDER:raise ValueError('Score program order does not match frozen S01–S16.')
    values=np.asarray(a.obsm['X_mumomac_s16'])
    if values.shape!=(a.n_obs,16) or np.isinf(values).any():raise ValueError('Invalid score shape or infinite scores.')
    if not isinstance(min_cells,int) or min_cells<1:raise ValueError('min_cells must be a positive descriptive aggregation threshold.')
    if missing_metadata not in ('error','label'):raise ValueError('missing_metadata must be error or label.')
    required=list(dict.fromkeys(k for k in [sample_key,identity_key,tissue_key,group_key,biological_unit_key] if k is not None))
    if len({sample_key,identity_key,tissue_key,group_key}-{None})!=len([k for k in [sample_key,identity_key,tissue_key,group_key] if k is not None]):raise ValueError('Sample/identity/tissue/group keys must name distinct concepts.')
    for key in required:
        if key not in a.obs:raise KeyError(f'Metadata column {key!r} is missing.')
    obs=a.obs[required].copy();missing=obs.isna().any(axis=1)
    if missing.any() and missing_metadata=='error':raise ValueError('Required aggregation metadata contain missing values; fix them or explicitly use missing_metadata=label for descriptive output.')
    biovalid=np.ones(a.n_obs,bool) if biological_unit_key else np.zeros(a.n_obs,bool)
    if biological_unit_key:biovalid=obs[biological_unit_key].notna().to_numpy()
    for key in required:
        if obs[key].isna().any():obs[key]=obs[key].astype(object).where(obs[key].notna(),'(Missing)')
    # Group and donor assignments must be constant within a technical sample.
    for key in [group_key,biological_unit_key]:
        if key and key!=sample_key and obs.groupby(sample_key,observed=True)[key].nunique().gt(1).any():raise ValueError(f'One sample maps to multiple {key} values; pooled/paired designs require an explicit corrected sample manifest.')
    strata=[k for k in [identity_key,tissue_key] if k is not None];keys=[sample_key]+strata
    rows=[]
    for key,ids in obs.groupby(keys,sort=False,observed=True,dropna=False).indices.items():
        key=key if isinstance(key,tuple) else (key,);base=dict(zip(keys,key));ids=np.asarray(ids);cellobs=obs.iloc[ids]
        for name in [group_key,biological_unit_key]:
            if name:base[name]=cellobs[name].iloc[0]
        for j,pid in enumerate(PROGRAM_ORDER):
            v=values[ids,j];valid=np.isfinite(v);nv=int(valid.sum());flags=[]
            if len(ids)<min_cells:flags.append('Below declared descriptive minimum cells')
            if not nv:flags.append('Not Evaluable')
            elif nv<len(ids):flags.append('Some cells have missing scores')
            if missing.iloc[ids].any():flags.append('Missing metadata explicitly labelled')
            rows.append(dict(**base,program_id=pid,n_cells=len(ids),n_valid_cells=nv,
                             mean_score=float(v[valid].mean()) if nv else np.nan,
                             median_score=float(np.median(v[valid])) if nv else np.nan,
                             score_sd=float(v[valid].std(ddof=1)) if nv>1 else np.nan,
                             available_program=bool(nv and len(ids)>=min_cells),
                             biological_unit_valid=bool(biovalid[ids].all()) if biological_unit_key else False,
                             quality_flags='; '.join(flags) or 'Descriptive score available; biological reliability unverified',
                             score_config_hash=a.uns['mumomac']['score_config_hash']))
    out=pd.DataFrame(rows)
    out.attrs['mumomac_summary']={'schema':'sample_summary_v1','sample_key':sample_key,'strata_keys':strata,
                                  'group_key':group_key or '', 'biological_unit_key':biological_unit_key or '',
                                  'biological_unit_assertion':'User-declared; not externally certified' if biological_unit_key else 'Not provided; no independent unit assumed',
                                  'minimum_cells_engineering':min_cells,'aggregation':'Mean/median/SD of observed cells within sample and declared strata; no cell-level inference',
                                  'score_config_hash':a.uns['mumomac']['score_config_hash']}
    return out
