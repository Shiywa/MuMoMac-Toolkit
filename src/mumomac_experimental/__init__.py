"""Experimental integration; frozen scientific engines remain unchanged."""
__version__='0.1.1rc1'
from .workflow import analyze, analyze_scores
from .input import load_input,validate_input,select_momac,harmonize_genes
from .integration import analyze_score_run
from .seed_qc import seed_details
