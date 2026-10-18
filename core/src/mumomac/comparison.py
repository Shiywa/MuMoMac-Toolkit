"""Simple independent two-group unit comparison, without cell pseudoreplication."""
import itertools,math
import numpy as np,pandas as pd
from .io import read_table

def _adjust(p,method):
    n=len(p);order=np.argsort(p);sortedp=p[order]
    if method=='holm':q=np.maximum.accumulate(sortedp*(n-np.arange(n)))
    elif method=='bh':q=np.minimum.accumulate((sortedp*n/(np.arange(n)+1))[::-1])[::-1]
    else:raise ValueError('correction must be holm or bh.')
    result=np.empty(n);result[order]=np.minimum(q,1.);return result

def compare_programs(summary,*,group_key,reference,case,biological_unit_key=None,
                     strata_keys=None,n_bootstrap=2000,n_permutations=9999,
                     seed=20261011,correction='holm'):
    """Compare independent biological units; repeated samples are equally averaged.

    Paired/crossover/mixed designs are rejected. Passing a biological key is a
    user assertion of biological independence; no UID inference is performed.
    Without it only sample-level descriptive effects are returned.
    """
    if not isinstance(summary,pd.DataFrame):summary=read_table(summary)
    d=summary.copy();meta=summary.attrs.get('mumomac_summary',{})
    if reference==case:raise ValueError('case and reference must differ.')
    if n_bootstrap<1 or n_permutations<1:raise ValueError('Resampling counts must be positive.')
    if correction not in ('holm','bh'):raise ValueError('correction must be holm or bh.')
    needed=[group_key,'program_id','mean_score','available_program','score_config_hash']
    for k in needed:
        if k not in d:raise KeyError(f'Required summary column {k!r} is absent; use summarize_programs.')
    if d.score_config_hash.nunique()!=1:raise ValueError('Summary mixes incompatible score definitions/configurations/gene universes; absolute calibration is not supported.')
    declared=list(meta.get('strata_keys',[]))
    if strata_keys is None:
        if not meta:raise ValueError('Summary metadata are absent; provide explicit strata_keys (including [] for an unstratified table).')
        strata=list(declared)
    else:
        strata=list(strata_keys)
        if any(k not in strata and d[k].nunique()>1 for k in declared):raise ValueError('Cannot pool distinct declared identity/tissue strata as independent observations.')
    for k in strata:
        if k not in d:raise KeyError(f'Stratum {k!r} is absent.')
    sample_key=meta.get('sample_key','')
    if not sample_key or sample_key not in d:raise ValueError('A traceable sample_key in summary metadata is required.')
    if d.duplicated([sample_key]+declared+['program_id']).any():raise ValueError('Duplicate sample/program summary records; do not concatenate repeated copies as new samples.')
    d=d[d[group_key].isin([reference,case])]
    if set(d[group_key])!={reference,case}:raise ValueError('Both requested groups must be present.')
    infer=biological_unit_key is not None
    if infer:
        if meta.get('biological_unit_key') and meta['biological_unit_key']!=biological_unit_key:
            raise ValueError('Requested biological key differs from the declared summary mapping; technical sample IDs cannot replace donor/animal IDs.')
        if biological_unit_key not in d:raise KeyError(f'Biological unit column {biological_unit_key!r} is absent; no independent animal/donor can be inferred.')
        if d[biological_unit_key].isna().any() or d[biological_unit_key].eq('(Missing)').any():raise ValueError('Biological unit IDs are missing; formal inference is unsupported.')
        if 'biological_unit_valid' in d and not d.biological_unit_valid.all():raise ValueError('Summary did not preserve a valid declared biological-unit mapping; regenerate it with biological_unit_key.')
        if d.groupby(biological_unit_key,observed=True)[group_key].nunique().gt(1).any():raise ValueError('A biological unit appears in both groups (paired/crossover design); v0.1 supports independent groups only.')
    unit=biological_unit_key if infer else sample_key;rng=np.random.default_rng(seed);out=[]
    by=strata+['program_id']
    for keys,z in d.groupby(by,sort=False,observed=True,dropna=False):
        keys=keys if isinstance(keys,tuple) else (keys,);base=dict(zip(by,keys))
        usable=z[z.available_program&z.mean_score.notna()]
        # One vote per declared biological unit, independent of cell/library count.
        scores=usable.groupby([unit,group_key],observed=True,sort=False).mean_score.mean().reset_index()
        a=scores[scores[group_key].eq(case)].mean_score.to_numpy();b=scores[scores[group_key].eq(reference)].mean_score.to_numpy()
        effect=float(a.mean()-b.mean()) if len(a) and len(b) else np.nan
        row=dict(**base,reference=str(reference),case=str(case),n_case_units=len(a) if infer else np.nan,
                 n_reference_units=len(b) if infer else np.nan,n_case_samples=z[z[group_key].eq(case)][sample_key].nunique(),
                 n_reference_samples=z[z[group_key].eq(reference)][sample_key].nunique(),
                 mean_difference=effect,median_difference=float(np.median(a)-np.median(b)) if len(a) and len(b) else np.nan,
                 effect_direction='Increase' if np.isfinite(effect) and effect>0 else 'Decrease' if np.isfinite(effect) and effect<0 else 'Zero' if np.isfinite(effect) else 'Not Evaluable',
                 standardized_effect=np.nan,ci_low=np.nan,ci_high=np.nan,p_value=np.nan,p_adjusted=np.nan,
                 statistical_unit='User-declared biological unit' if infer else 'Technical sample descriptive only',
                 repeated_observation_policy='Equal sample means within unit and declared stratum; never cell-weighted inference',
                 inference_status='No declared biological unit: descriptive only' if not infer else 'Inference Limited: fewer than 3 evaluable units per group',
                 biological_reliability='Unverified RNA association; not measured function')
        if infer and min(len(a),len(b))>=2:
            df=len(a)+len(b)-2;sp=np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/df)
            if sp>0:row['standardized_effect']=(1-3/(4*(len(a)+len(b))-9))*effect/sp
        if infer and min(len(a),len(b))>=3:
            boots=a[rng.integers(0,len(a),(n_bootstrap,len(a)))].mean(1)-b[rng.integers(0,len(b),(n_bootstrap,len(b)))].mean(1)
            row['ci_low'],row['ci_high']=np.quantile(boots,[.025,.975]);y=np.r_[a,b];k=len(a);N=len(y);alloc=math.comb(N,k)
            if alloc<=100000:
                ix=np.array(list(itertools.combinations(range(N),k)));v=y[ix].mean(1)-(y.sum()-y[ix].sum(1))/(N-k)
                row['p_value']=min(1.,2*min((v<=effect+1e-12).mean(),(v>=effect-1e-12).mean()));row['p_method']='Exhaustive permutation, twice smaller tail';row['permutation_count']=alloc
            else:
                lt=gt=0
                for start in range(0,n_permutations,2000):
                    m=min(2000,n_permutations-start);ix=np.argpartition(rng.random((m,N)),k-1,axis=1)[:,:k];v=y[ix].mean(1)-(y.sum()-y[ix].sum(1))/(N-k);lt+=int((v<=effect+1e-12).sum());gt+=int((v>=effect-1e-12).sum())
                row['p_value']=min(1.,2*min((lt+1)/(n_permutations+1),(gt+1)/(n_permutations+1)));row['p_method']='Monte Carlo permutation, double-tail +1';row['permutation_count']=n_permutations
            row['inference_status']='Simple independent-group RNA association; metadata/confounding not certified'
            if alloc<=100000 and 2/alloc>.05:row['inference_status']+='; limited exact p resolution'
        out.append(row)
    result=pd.DataFrame(out)
    if len(result):
        p=result.p_value.to_numpy(dtype=float);q=_adjust(np.where(np.isfinite(p),p,1.),correction);q[~np.isfinite(p)]=np.nan;result['p_adjusted']=q
    result.attrs['mumomac_comparison']={'schema':'independent_two_group_v1','family_size':len(result),'correction':correction,
                                      'stopped_tests':'Displayed missing; conservative1 padding only for family adjustment',
                                      'biological_unit_key':biological_unit_key or '', 'strata_keys':strata,
                                      'seed':seed,'n_bootstrap':n_bootstrap,'n_permutations':n_permutations,
                                      'complex_designs':'Paired/crossover/mixed/adjusted inference is not supported in v0.1'}
    return result
