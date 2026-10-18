"""Additive reference CLI; previous v0.2 commands keep their defaults."""
import argparse,json
from pathlib import Path
from .reference import load_reference,resolve_feature,validate_reference
from .reference_scoring import get_program_coverage
from .io import read_adata,read_table,write_table
from .analysis import compare_conditions,analyze_combinations
from .advisor import advise_signature_version
from .workflow import run_workflow,_jsonable
from .input import profile_input

def metadata_options(p):
    for name in ('sample','identity','condition','biological-unit','study','tissue'):p.add_argument('--'+name+'-col')
    p.add_argument('--reference');p.add_argument('--case');p.add_argument('--paired',action='store_true')

def main(argv=None):
    parser=argparse.ArgumentParser(prog='mumomac');sub=parser.add_subparsers(dest='command',required=True)
    r=sub.add_parser('reference');r.add_argument('--feature');r.add_argument('--validate',action='store_true')
    for cmd in ('run','report','profile'):
      p=sub.add_parser(cmd);p.add_argument('input');p.add_argument('--output-dir',required=True);p.add_argument('--species',choices=['mouse','human'],default='mouse');p.add_argument('--layer',default='counts');p.add_argument('--expression-type',choices=['counts','log1p','nonnegative_expression']);p.add_argument('--gene-symbol-col');p.add_argument('--coverage',choices=['reference','complete','available'],default='reference');p.add_argument('--chunk-size',type=int,default=500);p.add_argument('--advisor-mode',choices=['coverage','umi'],default='coverage');p.add_argument('--raw-umi',action='store_true',help='Explicit assay provenance assertion; integers alone do not certify UMI');metadata_options(p)
    a=sub.add_parser('advise');a.add_argument('input');a.add_argument('--output-dir',required=True);a.add_argument('--mode',choices=['coverage','umi'],default='coverage');a.add_argument('--raw-umi',action='store_true');metadata_options(a)
    c=sub.add_parser('compare-conditions');c.add_argument('input');c.add_argument('--output',required=True);c.add_argument('--reference',required=True);c.add_argument('--case',required=True);c.add_argument('--paired',action='store_true')
    b=sub.add_parser('combination');b.add_argument('input');b.add_argument('--output',required=True)
    args=parser.parse_args(argv)
    try:
      if args.command=='reference':
        value=validate_reference() if args.validate else resolve_feature(args.feature) if args.feature else load_reference()
        print(json.dumps(value,ensure_ascii=False,indent=2));return 0
      if args.command in ('compare-conditions','combination'):
        d=read_table(args.input)
        out=compare_conditions(d,reference=args.reference,case=args.case,paired=args.paired) if args.command=='compare-conditions' else analyze_combinations(d)
        write_table(out,args.output);return 0
      out=Path(args.output_dir)
      if out.exists() and any(out.iterdir()):raise FileExistsError('Output directory must be empty/new')
      if args.command in ('run','report'):
        run_workflow(args.input,args.output_dir,species=args.species,layer=args.layer,expression_type=args.expression_type,gene_symbol_key=args.gene_symbol_col,sample_key=args.sample_col,identity_key=args.identity_col,condition_key=args.condition_col,biological_unit_key=args.biological_unit_col,study_key=args.study_col,tissue_key=args.tissue_col,reference=args.reference,case=args.case,paired=args.paired,advisor_mode=args.advisor_mode,umi_counts_verified=args.raw_umi,coverage_policy=args.coverage,chunk_size=args.chunk_size)
      else:
        out.mkdir(parents=True,exist_ok=True)
        if args.command=='profile':
          result=profile_input(args.input,layer=args.layer,expression_type=args.expression_type,gene_symbol_key=args.gene_symbol_col,sample_key=args.sample_col)
          write_table(result['cell_qc'],out/'CELL_QC.tsv');write_table(result['sample_qc'],out/'SAMPLE_QC.tsv');(out/'PROFILE.json').write_text(json.dumps(_jsonable(result['summary']),indent=2))
        else:
          result=advise_signature_version(read_adata(args.input),mode=args.mode,sample_key=args.sample_col,biological_unit_key=args.biological_unit_col,condition_key=args.condition_col,study_key=args.study_col,identity_key=args.identity_col,reference=args.reference,case=args.case,umi_counts_verified=args.raw_umi)
          for name,value in result.items():
            if hasattr(value,'to_csv'):write_table(value,out/(name.upper()+'.tsv'))
          (out/'POLICY.json').write_text(json.dumps(_jsonable(result['policy']),indent=2))
      return 0
    except (ValueError,KeyError,FileExistsError,FileNotFoundError,TypeError,OSError) as error:parser.error(str(error))
