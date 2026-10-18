import argparse,json,sys
from pathlib import Path
def main(argv=None):
 p=argparse.ArgumentParser(prog='mumomac-experimental',description='Experimental Original / Expanded / Both → Seed QC → continuous-first representation → sample coupling.')
 p.add_argument('--version',action='version',version='0.1.1rc1')
 subs=p.add_subparsers(dest='command',required=True)
 a=subs.add_parser('analyze',help='Provided mouse/human H5AD workflow');a.add_argument('--input',required=True);a.add_argument('--metadata');a.add_argument('--config',required=True);a.add_argument('--output',required=True)
 a.add_argument('--signature-version',choices=['original','expanded','both']);a.add_argument('--original-min-present',type=int);a.add_argument('--original-min-detected-per-cell',type=int);a.add_argument('--original-coverage-policy');a.add_argument('--scoring-qc-mode',choices=['strict','audit'])
 s=subs.add_parser('coupling',help='Current Original or Expanded score-only workflow; Original S requires trusted SeedQC evidence')
 for k in ['scores','metadata','spec','coverage','config','output']:s.add_argument('--'+k,required=True)
 s.add_argument('--detected-members');s.add_argument('--frozen-model');s.add_argument('--seed-qc-run')
 d=subs.add_parser('seed-details',help='Exact detected/zero/not-measured/unmapped genes for one cell/program');d.add_argument('--run',required=True);d.add_argument('--cell-id',required=True);d.add_argument('--program-id',required=True);d.add_argument('--output')
 b=subs.add_parser('score-run',help='Replay one exported Original/Expanded/Both run without expression or re-scoring');b.add_argument('--input-run',required=True);b.add_argument('--output',required=True);b.add_argument('--config')
 r=subs.add_parser('reference-lookup',help='Optional descriptive metadata index lookup, never a normality test');r.add_argument('--index',required=True);r.add_argument('--output',required=True)
 for k in ['tissue','identity','condition']:r.add_argument('--'+k)
 x=p.parse_args(argv)
 try:
  if x.command=='analyze':
   from .workflow import analyze,load_config
   cfg=load_config(x.config)
   for k in ['signature_version','original_min_present','original_min_detected_per_cell','original_coverage_policy','scoring_qc_mode']:
    if getattr(x,k) is not None:cfg[k]=getattr(x,k)
   result=analyze(x.input,cfg,x.output,x.metadata)
  elif x.command=='coupling':
   from .workflow import analyze_scores
   result=analyze_scores(x.scores,x.metadata,x.spec,x.coverage,x.config,x.output,x.detected_members,x.frozen_model,x.seed_qc_run)
  elif x.command=='seed-details':
   from .seed_qc import seed_details
   result=seed_details(x.run,x.cell_id,x.program_id);text=json.dumps(result,ensure_ascii=False,indent=2)
   if x.output:
    if Path(x.output).exists():raise FileExistsError('Seed detail output exists.')
    Path(x.output).write_text(text+'\n')
   else:print(text)
   return 0
  elif x.command=='score-run':
   from .integration import analyze_score_run
   result=analyze_score_run(x.input_run,x.output,x.config)
  else:
   import pandas as pd
   from .input import read_table
   out=Path(x.output)
   if out.exists():raise FileExistsError('Lookup output must be new.')
   t=read_table(x.index)
   for value,alternatives in [(x.tissue,['tissue','Tissue']),(x.identity,['coarse_identity','identity','Identity']),(x.condition,['condition','Condition'])]:
    if value:
     key=next((v for v in alternatives if v in t),None)
     if key is None:raise ValueError('Requested lookup field unavailable.')
     t=t[t[key].astype(str).eq(value)]
   t['lookup_interpretation']='DESCRIPTIVE_HISTORICAL_METADATA; not calibrated Current UCell normality or independent biological validation';t.to_csv(out,sep='\t',index=False)
   print(json.dumps({'rows':len(t),'output':str(out.resolve())}));return 0
  print(json.dumps({'status':result['qc']['status'],'output':result['output'],'n_selected':result['qc']['n_selected']}));return 0
 except Exception as e:print(type(e).__name__+': '+str(e),file=sys.stderr);return 2
if __name__=='__main__':raise SystemExit(main())
