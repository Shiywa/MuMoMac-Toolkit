#!/usr/bin/env python3
"""Generate an artificial software fixture; contains no biological source data."""
import argparse,json
from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,yaml
from scipy import sparse
import mumomac

def make_example(output):
    out=Path(output)
    if out.exists() and any(out.iterdir()):raise FileExistsError('Use a new empty example directory.')
    out.mkdir(parents=True,exist_ok=True)
    ref=mumomac.reference_signatures('expanded','mouse')
    genes=sorted(set(g for x in ref.values() for g in x['genes']))+[f'SyntheticBackground{i:03d}' for i in range(80)]
    rng=np.random.default_rng(20261018);n=256
    counts=rng.poisson(rng.uniform(.8,5,(n,1)),(n,len(genes))).astype('int32')
    s03=mumomac.reference_signatures('original','mouse')['S03']['genes']
    ix=[genes.index(g) for g in s03];counts[0,ix]=0;counts[0,ix[:2]]=5
    labels=['Macrophage' if i%64<30 else 'Monocyte' if i%64<60 else 'Other' for i in range(n)]
    metadata=pd.DataFrame({'cell_id':[f'SYN{i:04d}' for i in range(n)],'sample_id':[f'sample{i//64+1:02d}' for i in range(n)],'provided_identity':labels,'study_id':'SYNTHETIC_ENGINEERING_FIXTURE','source_block_id':'synthetic_source','tissue':'Synthetic','condition':'Unknown','biological_unit_id':'Unknown','biological_unit_confirmed':False})
    # Labels are supplied externally; the matrix itself has only sample IDs.
    obs=metadata.set_index('cell_id')[['sample_id']].copy()
    a=ad.AnnData(sparse.csr_matrix(counts),obs=obs,var=pd.DataFrame(index=genes));a.layers['counts']=a.X.copy();a.uns['synthetic_fixture']=True
    a.write_h5ad(out/'synthetic_counts.h5ad',compression='gzip')
    metadata.iloc[::-1].to_csv(out/'provided_metadata.tsv',sep='\t',index=False)
    config=dict(species='mouse',modality='scRNA-seq',expression_type='counts',expression_provenance='Code-generated synthetic nonnegative integer counts; no biological validation',layer='counts',gene_id_type='symbol',selection_mode='provided',sample_column='sample_id',cell_type_column='provided_identity',identity_column='provided_identity',target_labels=['Macrophage','Monocyte'],study_column='study_id',source_column='source_block_id',condition_column='condition',tissue_column='tissue',biological_unit_column='biological_unit_id',biological_unit_confirmed_column='biological_unit_confirmed',signature_version='both',selected_pairs=['C01','C03','C04'],representation={'default_mode':'continuous_only'},scoring_qc_mode='strict',original_min_present=3,original_min_detected_per_cell=3,original_coverage_policy='partial_seed3_v1',include_joint_index=False,include_tail=False)
    (out/'config.yaml').write_text(yaml.safe_dump(config,sort_keys=False))
    manifest=dict(synthetic=True,seed=20261018,n_input_cells=n,n_target_cells=240,n_native_genes=len(genes),n_sample_records=4,n_independent_biological_replicates=0,purpose='Software execution, joins and QC only; no biological ground truth',controlled_failure={'cell_id':'SYN0000','program':'S03','measured_members':len(ix),'detected_members':2,'expected_original_strict_score':'NA'})
    (out/'example_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest))
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True);args=p.parse_args();make_example(args.output)
