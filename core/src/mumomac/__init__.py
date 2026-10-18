"""MuMoMac: exploratory frozen RNA signature transfer, not a classifier."""
__version__='0.3.0'
from .signatures import load_signatures,signature_manifest,PROGRAM_ORDER,Program
from .scoring import score_programs
from .validation import validate_programs
from .aggregation import summarize_programs
from .comparison import compare_programs
from .io import read_adata,write_adata,read_table,write_table
from .reliability import assess_reliability,load_robustness_benchmark,write_reliability_report,RobustnessBenchmark
__all__=['load_signatures','signature_manifest','PROGRAM_ORDER','Program','score_programs','validate_programs','summarize_programs','compare_programs','read_adata','write_adata','read_table','write_table','assess_reliability','load_robustness_benchmark','write_reliability_report','RobustnessBenchmark','__version__']
from .reference import load_reference,resolve_feature,reference_signatures,validate_reference,REFERENCE_ORDER
from .reference_scoring import score_signature_versions,get_program_coverage,score_historical_raw
from .input import from_matrix,read_matrix_market,profile_input,validate_input_metadata
from .analysis import aggregate_program_scores,biological_unit_scores,compare_conditions,compare_within_identity,assess_identity,summarize_composition,analyze_combinations,historical_combinations
from .advisor import advise_signature_version
from .workflow import run_workflow,write_analysis_report
__all__ += ['load_reference','resolve_feature','reference_signatures','validate_reference','REFERENCE_ORDER','score_signature_versions','get_program_coverage','score_historical_raw','from_matrix','read_matrix_market','profile_input','validate_input_metadata','aggregate_program_scores','biological_unit_scores','compare_conditions','compare_within_identity','assess_identity','summarize_composition','analyze_combinations','historical_combinations','advise_signature_version','run_workflow','write_analysis_report']
