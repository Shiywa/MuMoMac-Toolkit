import argparse,sys
from pathlib import Path
from . import __version__,score_programs,validate_programs,summarize_programs,compare_programs
from .io import write_adata,write_table,read_table
from .reliability import assess_reliability,write_reliability_report

def _input_options(p):
    p.add_argument('input');p.add_argument('--species',choices=['mouse','human'],required=True)
    p.add_argument('--signature',choices=['original','expanded','both'],default='original')
    p.add_argument('--feature-scope',choices=['s16','reference'],default='s16')
    p.add_argument('--layer',default='counts',help='Named layer; X explicitly selects adata.X')
    p.add_argument('--expression-type',choices=['counts','log1p','nonnegative_expression'])
    p.add_argument('--gene-symbol-col');p.add_argument('--duplicates',choices=['error','sum','exclude'],default='error')
    p.add_argument('--coverage',choices=['complete','available'],default='complete')
    p.add_argument('--chunk-size',type=int,default=500)

def main(argv=None):
    values=sys.argv[1:] if argv is None else list(argv)
    if values and values[0] in ('reference','advise','run','report','compare-conditions','combination','profile'):
        from .cli_reference import main as reference_main
        return reference_main(values)
    parser=argparse.ArgumentParser(prog='mumomac',description='Exploratory MuMoMac RNA program definition transfer',epilog='Reference v1 commands: reference, run, report, profile, advise, compare-conditions, combination. Use mumomac COMMAND --help.')
    parser.add_argument('--version',action='version',version=__version__);sub=parser.add_subparsers(dest='command',required=True)
    v=sub.add_parser('validate');_input_options(v);v.add_argument('--output');v.add_argument('--overwrite',action='store_true')
    p=sub.add_parser('score');_input_options(p);p.add_argument('--output',required=True);p.add_argument('--max-rank',type=int,default=1500);p.add_argument('--n-jobs',type=int,default=1);p.add_argument('--overwrite',action='store_true');p.add_argument('--add-obs',action='store_true')
    s=sub.add_parser('summarize');s.add_argument('input');s.add_argument('--sample-col',required=True);s.add_argument('--identity-col');s.add_argument('--tissue-col');s.add_argument('--group-col');s.add_argument('--biological-unit-col');s.add_argument('--min-cells',type=int,default=1);s.add_argument('--output',required=True);s.add_argument('--overwrite',action='store_true')
    c=sub.add_parser('compare');c.add_argument('input');c.add_argument('--group-col',required=True);c.add_argument('--reference',required=True);c.add_argument('--case',required=True);c.add_argument('--biological-unit-col');c.add_argument('--correction',choices=['holm','bh'],default='holm');c.add_argument('--seed',type=int,default=20261011);c.add_argument('--output',required=True);c.add_argument('--overwrite',action='store_true')
    r=sub.add_parser('reliability');r.add_argument('input');r.add_argument('--benchmark-dir');r.add_argument('--reference-benchmark-dir');r.add_argument('--output-prefix',required=True);r.add_argument('--overwrite',action='store_true')
    args=parser.parse_args(argv)
    try:
        if getattr(args,'output',None):
            if Path(args.output).resolve()==Path(args.input).resolve():raise ValueError('Input and output paths must differ; original input is protected.')
            if Path(args.output).exists() and not args.overwrite:raise FileExistsError('Output already exists; no file was overwritten.')
        if args.command in ('score','validate'):
            kw=dict(species=args.species,signature_version=args.signature,layer=args.layer,expression_type=args.expression_type,gene_symbol_key=args.gene_symbol_col,duplicate_policy=args.duplicates,coverage_policy=args.coverage,chunk_size=args.chunk_size)
            if args.command=='score':
                out=score_programs(args.input,**kw,max_rank=args.max_rank,n_jobs=args.n_jobs,copy=False,add_obs=args.add_obs,feature_scope=args.feature_scope);write_adata(out,args.output,args.overwrite)
            else:
                if args.feature_scope=='reference' or args.signature=='both':
                    from .reference_scoring import get_program_coverage
                    kw['coverage_policy']='reference' if args.coverage=='complete' else args.coverage
                    out=get_program_coverage(args.input,**kw)
                else:out=validate_programs(args.input,**kw)
                if args.output:write_table(out,args.output,args.overwrite)
                else:print(out.to_csv(sep='\t',index=False,na_rep='Missing'),end='')
        elif args.command=='summarize':
            out=summarize_programs(args.input,sample_key=args.sample_col,identity_key=args.identity_col,tissue_key=args.tissue_col,group_key=args.group_col,biological_unit_key=args.biological_unit_col,min_cells=args.min_cells);write_table(out,args.output,args.overwrite)
        elif args.command=='reliability':
            out=assess_reliability(args.input,benchmark=args.benchmark_dir,reference_benchmark=args.reference_benchmark_dir);write_reliability_report(out,args.output_prefix,overwrite=args.overwrite)
        else:
            out=compare_programs(read_table(args.input),group_key=args.group_col,reference=args.reference,case=args.case,biological_unit_key=args.biological_unit_col,correction=args.correction,seed=args.seed);write_table(out,args.output,args.overwrite)
        return 0
    except (ValueError,KeyError,FileExistsError,FileNotFoundError,TypeError,OSError) as error:
        parser.error(str(error))
    return 2

if __name__=='__main__':sys.exit(main())
