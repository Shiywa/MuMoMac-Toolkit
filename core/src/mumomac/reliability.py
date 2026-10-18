"""Context-bound measurement quality, measured robustness and historical evidence."""
from dataclasses import dataclass
from pathlib import Path
import hashlib,json,datetime
import numpy as np,pandas as pd
from scipy import sparse
from .io import read_adata,write_table,_commit
from .signatures import PROGRAM_ORDER,signature_manifest
from .validation import validate_programs

def input_fingerprint(adata,layer='counts'):
    """Canonical ordered cell/gene/expression digest; no inferred biological IDs."""
    a=read_adata(adata);h=hashlib.sha256(json.dumps({'cells':a.obs_names.tolist(),'genes':a.var_names.tolist(),'layer':layer},ensure_ascii=False).encode())
    X=a.X if layer in (None,'X') else a.layers[layer]
    for lo in range(0,a.n_obs,500):
        b=sparse.csr_matrix(X[lo:lo+500]).astype(np.float64);b.sum_duplicates();b.eliminate_zeros();b.sort_indices()
        for z in [b.indptr.astype(np.int64),b.indices.astype(np.int64),b.data]:h.update(z.tobytes())
    return h.hexdigest()

def symbol_fingerprint(a,key=None):
    values=a.var_names.tolist() if key is None or key=='var_names' else [None if pd.isna(x) else str(x) for x in a.var[key]]
    return hashlib.sha256(json.dumps(values,ensure_ascii=False).encode()).hexdigest()

def _score_digest(scores):
    x=np.asarray(scores,np.float32).copy();x[np.isnan(x)]=np.float32(-1)
    return hashlib.sha256(x.tobytes()).hexdigest()

def _clean(value):
    if isinstance(value,dict):return {str(k):_clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,np.ndarray)):return [_clean(v) for v in value]
    if isinstance(value,(float,np.floating)):return float(value) if np.isfinite(value) else None
    if isinstance(value,(bool,np.bool_)):return bool(value)
    if isinstance(value,(int,np.integer)):return int(value)
    if pd.isna(value):return None
    return value

@dataclass
class RobustnessBenchmark:
    directory: str
    context: dict
    metrics: pd.DataFrame
    configurations: dict
    score_digests: dict
    metadata_fields: dict
    resource_hashes: dict

def load_robustness_benchmark(directory):
    """Read an actual local frozen benchmark, verifying its complete references.

    Expected files: execution_summary.json, program_metrics.tsv,
    baseline_original.h5ad and baseline_expanded.h5ad. This loader grants no
    biological validation and does not turn a reference into current-data tests.
    """
    root=Path(directory);context=json.loads((root/'execution_summary.json').read_text())
    if context.get('all_conditions_status')!='PASS':raise ValueError('Benchmark execution is incomplete/failed; it cannot attest current-input robustness.')
    metrics=pd.read_csv(root/'program_metrics.tsv',sep='\t',na_values=['Missing'])
    required={'program_id','signature_version','condition_id','experiment','scheme','level','seed','cell_MAE','cell_spearman','n_paired_cells','evaluation_status'}
    if not required<=set(metrics):raise ValueError('Benchmark metrics schema is incomplete.')
    if metrics.duplicated(['signature_version','program_id','condition_id']).any():raise ValueError('Duplicate benchmark metric rows; technical repeats must be explicit.')
    if set(metrics.program_id)!=set(PROGRAM_ORDER):raise ValueError('Benchmark must preserve all16programs including stopped rows.')
    configs={};digests={};meta={};hashes={}
    for v in ['original','expanded']:
        a=read_adata(root/f'baseline_{v}.h5ad');u=a.uns['mumomac'];cfg=dict(u['configuration'])
        if input_fingerprint(a,cfg['input_layer'])!=context['input_fingerprint']:raise ValueError('Benchmark complete-reference expression/cell/gene fingerprint differs from execution record.')
        if cfg['signature_version']!=v or cfg['coverage_policy']!='complete':raise ValueError('Benchmark baseline must be the declared complete definition.')
        if tuple(u['program_order'])!=PROGRAM_ORDER:raise ValueError('Benchmark program order differs from frozen resources.')
        configs[v]=cfg;digests[v]=_score_digest(a.obsm['X_mumomac_s16']);hashes[v]=dict(u['definition_resource_hashes'])
        meta[v]={k:a.obs[k].astype(str).tolist() for k in ['sample_uid','source_block','identity','tissue','depth_stratum','study','animal_key','condition'] if k in a.obs}
    context['files_sha256']={f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in ['execution_summary.json','program_metrics.tsv']}
    return RobustnessBenchmark(str(root.resolve()),context,metrics,configs,digests,meta,hashes)

_MATH_KEYS=('scorer_version','signature_version','definition_version','species','scoring_method','max_rank','ties_method','coverage_policy','min_genes_engineering','effective_gene_universe_sha256','input_layer','expression_type','duplicate_policy','expanded_weights')

def assess_reliability(adata,*,benchmark=None,reference_benchmark=None):
    """Return16rows separating input quality, local robustness and historical validity.

    `benchmark` must match actual counts/cells/genes, score math, complete baseline
    score vectors and source/sample metadata. Mismatch raises an error. A
    `reference_benchmark` provides contextual risk notes only: current robustness
    remains NOT TESTED. No empirical stability threshold or probability is used.
    """
    a=read_adata(adata)
    if 'mumomac' not in a.uns or 'X_mumomac_s16' not in a.obsm:raise ValueError('Run score_programs before assess_reliability; existing scores/provenance are required.')
    u=a.uns['mumomac'];cfg=dict(u['configuration']);v=cfg['signature_version'];species=cfg['species'];layer=cfg['input_layer'];key=u['input_provenance'].get('gene_symbol_key','var_names');key=None if key=='var_names' else key
    current=input_fingerprint(a,layer);warnings=[]
    if 'input_fingerprint' in u:
        if current!=u['input_fingerprint'] or symbol_fingerprint(a,key)!=u.get('symbol_metadata_fingerprint'):raise ValueError('Expression/cell/gene/symbol input changed or was subset after scoring; recompute scores for the current input before quality assessment.')
    else:warnings.append('Legacy scored object lacks input attestation; measurement quality uses current data. Recompute withv0.2 before attaching direct local benchmark.')
    if tuple(u['program_order'])!=PROGRAM_ORDER:raise ValueError('Scored program order differs from frozen16programs.')
    matrix=np.asarray(a.obsm['X_mumomac_s16'])
    if matrix.shape!=(a.n_obs,16) or np.isinf(matrix).any():raise ValueError('Invalid score matrix shape or infinite values.')
    if ((matrix[np.isfinite(matrix)]<0)|(matrix[np.isfinite(matrix)]>1)).any():raise ValueError('UCell scores must be in [0,1]; corrupted scores are refused.')
    coverage=validate_programs(a,species=species,signature_version=v,layer=layer,expression_type=cfg['expression_type'],gene_symbol_key=key,duplicate_policy=cfg['duplicate_policy'],coverage_policy=cfg['coverage_policy'],min_genes=cfg['min_genes_engineering']).set_index('program_id')
    local=load_robustness_benchmark(benchmark) if isinstance(benchmark,(str,Path)) else benchmark
    reference=load_robustness_benchmark(reference_benchmark) if isinstance(reference_benchmark,(str,Path)) else reference_benchmark
    for x in [local,reference]:
        if x is not None and not isinstance(x,RobustnessBenchmark):raise TypeError('Benchmark must be a validated RobustnessBenchmark or its directory.')
    manifest=signature_manifest();limits=[]
    if local is not None:
        if 'input_fingerprint' not in u:raise ValueError('Direct local robustness attachment requires v0.2 input attestation; recompute scores.')
        if current!=local.context['input_fingerprint']:raise ValueError('Benchmark belongs to different counts/cells/genes; use reference_benchmark for risk context, not current-data TESTED.')
        if any(cfg.get(k)!=local.configurations[v].get(k) for k in _MATH_KEYS):raise ValueError('Benchmark scorer/definition/gene universe/coverage configuration does not match current scores.')
        if dict(u['definition_resource_hashes'])!=local.resource_hashes[v]:raise ValueError('Benchmark definition/mapping resource hashes differ.')
        if _score_digest(matrix)!=local.score_digests[v]:raise ValueError('Current baseline score vector differs from benchmark reference; numerical inconsistency cannot certify robustness.')
        for name,values in local.metadata_fields[v].items():
            if name not in a.obs or a.obs[name].astype(str).tolist()!=values:raise ValueError(f'Benchmark source/sample/identity mapping differs for {name}; biological units cannot be reassigned silently.')
    if reference is not None:
        if reference.resource_hashes[v]!=dict(u['definition_resource_hashes']):warnings.append('Reference benchmark definition/mapping differs; risk context is not directly applicable.')
        warnings.append('Reference benchmark is contextual evidence only; current data robustness remains NOT TESTED unless a matching direct local benchmark is provided.')
    rows=[]
    for j,pid in enumerate(PROGRAM_ORDER):
        c=coverage.loc[pid];definition=manifest['programs'][pid];score=matrix[:,j];finite=score[np.isfinite(score)]
        quality='NOT EVALUABLE under declared engineering input conditions' if not c.scoring_eligibility else 'COMPUTABLE; no observed signature detection' if c.n_detected_genes==0 else 'COMPUTABLE; biological reliability unverified'
        status='NOT TESTED';metrics={};n=0
        if local is not None:
            z=local.metrics[(local.metrics.signature_version==v)&(local.metrics.program_id==pid)&~local.metrics.experiment.eq('baseline')].copy()
            z=z[~(z.experiment.eq('missing_genes')&z.level.eq(0))&~(z.experiment.eq('depth')&z.level.eq(1))];n=len(z)
            status='TESTED: matching current input; inspect continuous conditions' if z.n_paired_cells.gt(0).any() else 'NOT EVALUABLE against complete same-input baseline'
            records=z.astype(object).where(pd.notna(z),None).to_dict('records');metrics={'condition_rows':_clean(records),'no_pass_threshold':True,'missing_tests_retained':True}
        msg=' '.join([str(c.warning_messages),*warnings,'Technical stability does not establish biological function or cross-study absolute calibration.'])
        rows.append({'program_id':pid,'program_name':definition['program_name'],'signature_version':v,'species':species,'definition_hash':manifest['source_files'][v]['sha256'],'definition_version':cfg['definition_version'],'scorer_version':cfg['scorer_version'],'gene_universe_size':u['input_provenance']['effective_features'],'genes_expected_mouse':int(c.n_signature_genes),'genes_expected_mapped':int(c.n_mapped_genes),'genes_present':int(c.n_present_genes),'genes_detected_any':int(c.n_detected_genes),'genes_missing':c.missing_genes,'genes_unmapped':c.unmapped_genes,'genes_ambiguous':c.ambiguous_genes,'mapping_retention':float(c.mapping_retention),'coverage_fraction':float(c.gene_coverage_fraction),'mean_present_gene_detection_fraction':float(c.mean_gene_detection_fraction),'median_present_gene_detection_fraction':float(c.median_gene_detection_fraction),'input_quality_status':quality,'n_cells':a.n_obs,'n_finite_scores':len(finite),'score_mean':float(finite.mean()) if len(finite) else np.nan,'score_median':float(np.median(finite)) if len(finite) else np.nan,'robustness_evaluation_status':status,'n_perturbation_rows':n,'robustness_metrics':json.dumps(_clean(metrics),ensure_ascii=False,allow_nan=False),'historical_evidence_status':definition['evidence_status'],'historical_external_evidence':definition['external_evidence'],'biological_validity_status':'NOT CERTIFIED by technical robustness','limitations':definition['known_limitations'],'warning_messages':msg})
    result=pd.DataFrame(rows);result.attrs['mumomac_reliability']={'schema':'context_bound_reliability_v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_fingerprint':current,'score_config_hash':u['score_config_hash'],'software_versions':dict(u['software_versions']),'resource_hashes':manifest['resource_sha256'],'local_benchmark_context':_clean(local.context) if local else None,'reference_benchmark_context':_clean(reference.context) if reference else None,'measurement_quality':'Current gene presence/detection and declared engineering computability; not biological certification','robustness':'Only actual matching current-input tests; no benchmark meansNOT TESTED','biological_validity':'Historical evidence only; never upgraded by numerical stability','warnings':warnings}
    return result

def write_reliability_report(report,prefix,*,overwrite=False):
    """Write protected TSV+metadata, strict JSON and a readable Markdown report."""
    if not isinstance(report,pd.DataFrame) or 'mumomac_reliability' not in report.attrs:raise ValueError('Expected a report from assess_reliability.')
    prefix=Path(prefix);paths={k:Path(str(prefix)+'.'+k) for k in ['tsv','json','md']}
    if not overwrite and any(p.exists() for p in [*paths.values(),Path(str(paths['tsv'])+'.metadata.json')]):raise FileExistsError('Reliability report output exists; no file was overwritten.')
    rows=report.astype(object).where(pd.notna(report),None).to_dict('records');metadata=report.attrs['mumomac_reliability'];payload={'metadata':_clean(metadata),'programs':_clean(rows)}
    _commit(paths['json'],lambda p:p.write_text(json.dumps(payload,ensure_ascii=False,allow_nan=False,indent=2)+'\n'),overwrite)
    write_table(report,paths['tsv'],overwrite)
    lines=['# MuMoMac program reliability report','', 'Measurement quality, measured technical robustness and historical biological evidence are separate. Continuous metrics are conditional on this input; no permanent Stable label,80%biological rule or95%Reliable probability is assigned.','',f"Input fingerprint: `{metadata['input_fingerprint']}`. Score configuration: `{metadata['score_config_hash']}`.",'','| Program | Definition | Presence | Detection across cells | Input quality | Robustness |','|---|---|---:|---:|---|---|']
    for _,r in report.iterrows():lines.append(f"| {r.program_id} | {r.signature_version} | {r.coverage_fraction:.3f} | {r.mean_present_gene_detection_fraction:.3f} | {r.input_quality_status} | {r.robustness_evaluation_status} |")
    lines += ['', 'Presence is mapped-member presence. Detection is mean per-present-gene positive-expression fraction across current cells. Both differ from program biological reliability. Missing scores are never presented as0.', '', 'Local benchmark context: '+json.dumps(_clean(metadata['local_benchmark_context']),ensure_ascii=False), '', 'Reference-only context: '+json.dumps(_clean(metadata['reference_benchmark_context']),ensure_ascii=False), '', 'Program-specific frozen evidence, adverse results, mapping losses, limitations and every measured condition are preserved in the TSV/JSON. RNA scores are not measured protein activity, function or flux; human extension remainsExperimental when mapping release is missing.','']
    _commit(paths['md'],lambda p:p.write_text('\n'.join(lines)),overwrite)
    return {k:str(p) for k,p in paths.items()}
