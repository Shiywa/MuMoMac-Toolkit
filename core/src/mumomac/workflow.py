"""Portable end-to-end reference workflow; no hidden Atlas dependency."""
from pathlib import Path
from importlib import metadata
import html,json,hashlib,sys,time
import numpy as np,pandas as pd
from scipy import stats
from .io import read_adata,write_adata,write_table
from .reference import load_reference,REFERENCE_ORDER
from .reference_scoring import score_signature_versions
from .input import profile_input,validate_input_metadata
from .analysis import aggregate_program_scores,biological_unit_scores,compare_conditions,analyze_combinations,assess_identity,summarize_composition
from .advisor import advise_signature_version

def _jsonable(x):
    if isinstance(x,dict):return {k:_jsonable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple,np.ndarray)):return [_jsonable(v) for v in x]
    if isinstance(x,np.generic):x=x.item()
    if isinstance(x,float) and not np.isfinite(x):return None
    return x

def write_analysis_report(result,output_dir):
    out=Path(output_dir);sections=[]
    def section(title,text,frame=None):sections.append((title,text,frame))
    section('1. Input Data & QC',json.dumps(result['profile']['summary'],ensure_ascii=False,indent=2),result['profile']['sample_qc'])
    section('2. Reference Mapping','Reference 1.0.0; user-provided identity. Unknown categories are retained. Human projection is experimental one-to-one.',result['identity'])
    cv=result['coverage']
    section('3. Functional Programs','Sixteen functional RNA proxies; continuous scores do not directly measure function.',cv[cv.feature_id.str.startswith('S')])
    section('4. Lifecycle / Residence','L04 is transition-associated and refine_required; L05 is an auxiliary identity feature. L scores in monocytes are descriptive; no historical L positive calls are transferred.',cv[cv.feature_id.str.startswith('L')])
    section('5. Global Activity','A01 is continuous. QC correlations below describe this exact input; they do not establish protein translation flux.',result['activity_qc'])
    section('6. Seed vs Expanded','Different score definitions/scales are not ranked by their absolute mean. Full Expanded L/A remain developmental. Advisor never changes the default.',result['advisor']['recommendations'])
    section('7. Condition Comparison',result['comparison_status']+' Composition and within-identity scores are exported separately; biological n comes from declared donor/animal IDs.',result['comparison'])
    section('8. Program Combinations','C01-C04 continuous exploratory associations. C02 describes residence co-variation; no confirmed interaction advantage. A01 is only a sensitivity covariate.',result['combinations'])
    section('9. Evidence & Limitations','Software operation, technical perturbation, RNA-response association and biological function validation are distinct. '+result['evidence_scope']+' No new gene selection or threshold calibration was performed. Missing columns and non-evaluable features are retained.',result['metadata_audit'])
    section('10. Reproducibility',json.dumps(result['provenance'],ensure_ascii=False,indent=2))
    md=['# MuMoMac reference analysis'];ht=['<!doctype html><html><head><meta charset="utf-8"><title>MuMoMac analysis</title><style>body{font-family:Arial,sans-serif;max-width:1200px;margin:2rem auto;background:white;color:#202124}table{border-collapse:collapse;font-size:12px;display:block;overflow:auto}td,th{padding:6px;border:1px solid #ddd}pre{white-space:pre-wrap}img{max-width:100%}</style></head><body><h1>MuMoMac reference analysis</h1>']
    for title,text,frame in sections:
      md.extend(['',f'## {title}','',text]);ht.extend([f'<h2>{html.escape(title)}</h2>',f'<pre>{html.escape(text)}</pre>'])
      if frame is not None and len(frame):
        # Full tables are TSV; report tables retain every feature, not a significant subset.
        cols=[c for c in frame.columns if c not in ('present_genes','absent_genes','unmapped_genes','biological_limits')]
        view=frame[cols];md.append(view.to_csv(sep='\t',index=False,na_rep='Missing'));ht.append(view.to_html(index=False,escape=True,na_rep='Missing'))
    figures=sorted(out.glob('*.png'))
    for p in figures:md.extend(['',f'![{p.stem}]({p.name})']);ht.append(f'<img src="{html.escape(p.name)}" alt="{html.escape(p.stem)}">')
    ht.append('</body></html>');(out/'REPORT.md').write_text('\n'.join(md));(out/'REPORT.html').write_text('\n'.join(ht))
    return out/'REPORT.html'

def _figures(result,out):
    try:
      import matplotlib
      matplotlib.use('Agg');import matplotlib.pyplot as plt
    except ImportError:return 'NOT_INSTALLED: install report extra for plots'
    plt.rcParams.update({'figure.facecolor':'white','axes.facecolor':'white','pdf.fonttype':42,'font.family':'Arial'})
    c=result['coverage'];labels=[p+' '+v for p,v in zip(c.feature_id,c.signature_version)]
    fig,ax=plt.subplots(figsize=(12,5));ax.bar(range(len(c)),c.coverage,color=np.where(c.scoring_eligible,'#2878a1','#bbbbbb'));ax.set_xticks(range(len(c)),labels,rotation=90,fontsize=6);ax.set_ylabel('Gene-axis coverage');fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(out/f'coverage.{ext}',dpi=180)
    plt.close(fig)
    comp=result['comparison']
    if len(comp):
      for version,z in comp.groupby('signature_version'):
        fig,ax=plt.subplots(figsize=(10,5));ax.axhline(0,color='gray',lw=1)
        for i,row in enumerate(z.itertuples()):
          if np.isfinite(row.mean_difference):
            ax.scatter(i,row.mean_difference,c='#2878a1')
            if np.isfinite(row.ci_low) and np.isfinite(row.ci_high):ax.plot([i,i],[row.ci_low,row.ci_high],color='#2878a1')
        ax.set_xticks(range(len(z)),z.feature_id,rotation=90,fontsize=7);ax.set_ylabel('Within-study condition mean difference');fig.tight_layout()
        for ext in ('png','pdf'):fig.savefig(out/f'condition_{version}.{ext}',dpi=180)
        plt.close(fig)
    return 'GENERATED'

def run_workflow(adata,output_dir,*,species='mouse',layer='counts',expression_type=None,gene_symbol_key=None,sample_key=None,identity_key=None,condition_key=None,biological_unit_key=None,study_key=None,tissue_key=None,reference=None,case=None,paired=False,advisor_mode='coverage',umi_counts_verified=False,coverage_policy='reference',chunk_size=500,evidence_scope='Application/technical analysis; independent biological validation not established',export_h5ad=True,identity_marker_sets=None):
    output=Path(output_dir)
    if output.exists() and any(output.iterdir()):raise FileExistsError('Use an empty/new output directory; original outputs are protected')
    output.mkdir(parents=True,exist_ok=True);start=time.perf_counter();a=read_adata(adata)
    data_sha=hashlib.sha256(Path(adata).read_bytes()).hexdigest() if isinstance(adata,(str,Path)) else 'Object; canonical matrix fingerprint recorded'
    audit=validate_input_metadata(a,sample_key=sample_key,biological_unit_key=biological_unit_key,condition_key=condition_key,identity_key=identity_key,study_key=study_key)
    profile=profile_input(a,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key,sample_key=sample_key,study_key=study_key)
    scored=score_signature_versions(a,species=species,signature_version='both',layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key,coverage_policy=coverage_policy,chunk_size=chunk_size,copy=False)
    identity=assess_identity(scored,identity_key=identity_key,marker_sets=identity_marker_sets,layer=layer,expression_type=expression_type,gene_symbol_key=gene_symbol_key) if identity_key else pd.DataFrame([{'status':'MISSING_USER_IDENTITY','analysis_scope':'Unstratified descriptive RNA only; no classifier run'}])
    summary=pd.DataFrame();comparison=pd.DataFrame();combinations=pd.DataFrame();composition=pd.DataFrame();unit_summary=pd.DataFrame();comparison_status='NOT_EVALUABLE: sample/condition metadata absent'
    if sample_key:
      summary=aggregate_program_scores(scored,sample_key=sample_key,identity_key=identity_key,condition_key=condition_key,biological_unit_key=biological_unit_key,study_key=study_key,tissue_key=tissue_key)
      qcgroup=list(dict.fromkeys(k for k in [study_key,sample_key] if k))
      sq=profile['sample_qc'];qcols=['median_expression_total','median_detected_genes','ribosomal_fraction']
      sq=sq[qcgroup+qcols].rename(columns={k:'qc_'+k for k in qcols});summary_attrs=dict(summary.attrs)
      summary=summary.merge(sq,on=qcgroup,how='left',validate='many_to_one');summary.attrs.update(summary_attrs)
      combinations=analyze_combinations(summary,qc_covariates=['qc_median_expression_total','qc_median_detected_genes'])
      if biological_unit_key:unit_summary=biological_unit_scores(summary)
      if identity_key:composition=summarize_composition(scored,sample_key=sample_key,identity_key=identity_key,condition_key=condition_key,biological_unit_key=biological_unit_key,study_key=study_key)
      if reference and case and condition_key:
        comparison=compare_conditions(summary,reference=reference,case=case,paired=paired)
        comparison_status='Biological-unit comparison' if biological_unit_key else 'DESCRIPTIVE_ONLY: actual independent biological-unit metadata absent'
    activity=[];qc=profile['cell_qc']
    for version in ['original','expanded']:
      values=scored.obsm[f'X_mumomac_reference_{version}'][:,-1]
      contexts=scored.obs.groupby(sample_key,observed=True,dropna=False).indices if sample_key else {'input':np.arange(scored.n_obs)}
      for context,idx in contexts.items():
       for cov in ['input_expression_total','detected_genes','mitochondrial_fraction','ribosomal_fraction']:
        x=qc[cov].to_numpy()[idx];y=values[idx];valid=np.isfinite(x)&np.isfinite(y)
        rho=float(stats.spearmanr(x[valid],y[valid]).statistic) if valid.sum()>=3 and np.std(x[valid])>1e-8 and np.std(y[valid])>1e-8 else np.nan
        activity.append(dict(signature_version=version,context=str(context),covariate=cov,n_cells=int(valid.sum()),spearman_rho=rho,scope='Cell descriptive within library; no cell-level p value'))
    advisor=advise_signature_version(scored,mode=advisor_mode,sample_key=sample_key,biological_unit_key=biological_unit_key,condition_key=condition_key,study_key=study_key,identity_key=identity_key,reference=reference,case=case,umi_counts_verified=umi_counts_verified,chunk_size=chunk_size)
    from . import __version__
    versions={p:metadata.version(p) for p in ['pyucell','numpy','pandas','scipy','anndata']};versions['mumomac-toolkit']=__version__
    provenance={'reference_version':'1.0.0','toolkit_version':__version__,'input_sha256':data_sha,'input_fingerprint':scored.uns['mumomac_reference']['input_fingerprint'],'score_configuration':scored.uns['mumomac_reference']['configuration'],'software_versions':versions,'argv':sys.argv,'advisor_policy':advisor['policy'],'biological_unit_key':biological_unit_key or 'Missing','evidence_scope':evidence_scope,'elapsed_seconds':time.perf_counter()-start}
    result={'scored':scored,'coverage':scored.uns['mumomac_reference']['gene_coverage'],'profile':profile,'identity':identity,'sample_scores':summary,'biological_unit_scores':unit_summary,'comparison':comparison,'comparison_status':comparison_status,'combinations':combinations,'composition':composition,'activity_qc':pd.DataFrame(activity),'advisor':advisor,'metadata_audit':audit,'provenance':provenance,'evidence_scope':evidence_scope}
    tables={'PROGRAM_COVERAGE':result['coverage'],'INPUT_CELL_QC':qc,'INPUT_SAMPLE_QC':profile['sample_qc'],'METADATA_AUDIT':audit,'IDENTITY_AUDIT':identity,'SAMPLE_PROGRAM_SCORES':summary,'BIOLOGICAL_UNIT_SCORES':unit_summary,'CONDITION_COMPARISON':comparison,'PROGRAM_COMBINATIONS':combinations,'CELL_COMPOSITION':composition,'A01_QC_ASSOCIATIONS':result['activity_qc'],'SIGNATURE_SUITABILITY_METRICS':advisor['metrics'],'SEED_EXPANDED_PAIRED_ROBUSTNESS':advisor['robustness'],'DEPTH_DETECTION_INTERACTION':advisor['interaction'],'SIGNATURE_RECOMMENDATIONS':advisor['recommendations']}
    for key in ['unit_metrics','paired_metrics','stratum_metrics']:
      if key in advisor:tables['ADVISOR_'+key.upper()]=advisor[key]
    for name,frame in tables.items():write_table(frame,output/f'{name}.tsv')
    for version in ['original','expanded']:
      values=pd.DataFrame(scored.obsm[f'X_mumomac_reference_{version}'],index=scored.obs_names,columns=REFERENCE_ORDER)
      values.index.name='cell_id';values.reset_index().to_csv(output/f'CELL_SCORES_{version}.tsv.gz',sep='\t',index=False)
    (output/'REPRODUCIBILITY.json').write_text(json.dumps(_jsonable(provenance),ensure_ascii=False,indent=2)+'\n')
    (output/'REFERENCE_REGISTRY.json').write_text(json.dumps(load_reference(),ensure_ascii=False,indent=2)+'\n')
    if export_h5ad:write_adata(scored,output/'scored.h5ad')
    result['plot_status']=_figures(result,output);write_analysis_report(result,output)
    files=[]
    for p in sorted(output.iterdir()):
      if p.is_file():files.append({'path':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    (output/'OUTPUT_MANIFEST.json').write_text(json.dumps(files,indent=2)+'\n')
    return result
