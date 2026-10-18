"""Explicit unknown evidence for score-only input; never invent zero detection."""
import numpy as np,pandas as pd
from mumomac_coupling.registry import ORDER
def export_unknown_seed_qc(bundle,out,version,raw=None):
 import mumomac
 from importlib import resources
 raw=bundle.scores if raw is None else raw;m=bundle.metadata;original=mumomac.reference_signatures('original','mouse');mapped=mumomac.reference_signatures('original',bundle.spec.species);root=resources.files('mumomac.resources')
 with root.joinpath('human_one_to_one.tsv').open() as f:sm=pd.read_csv(f,sep='\t')
 with root.joinpath('human_la_one_to_one_v1.tsv').open() as f:la=pd.read_csv(f,sep='\t')
 sm=dict(zip(sm.mouse_gene,sm.human_gene));la=dict(zip(la.mouse_gene,la.human_gene));cat=[];cov=[];cells=[]
 for p in ORDER:
  genes=original[p]['genes'];mapper={g:g for g in genes} if bundle.spec.species=='mouse' else sm if p.startswith('S') else la;declared=str(bundle.coverage.loc[p,'present_genes']).split(';') if 'present_genes' in bundle.coverage else [];present=set(declared);members=[mapper[g] for g in genes if g in mapper and mapper[g] in present];un=[g for g in genes if g not in mapper];absent=[g for g in genes if g in mapper and mapper[g] not in present]
  for g in genes:cat.append(dict(program_id=p,original_gene=g,mapped_gene=mapper.get(g),mapping_status='UNMAPPED' if g not in mapper else 'FROZEN_MAPPING',measured_in_dataset=mapper.get(g) in present if g in mapper else False,measurement_status='UNMAPPED' if g not in mapper else 'DECLARED_MEASURED_DETECTION_UNKNOWN' if mapper[g] in present else 'DECLARED_NOT_MEASURED',n_cells_detected=np.nan,fraction_cells_detected=np.nan,relevant_sample_detection='MISSING',exclusion_reason='MEMBER_DETECTION_EVIDENCE_MISSING; availability inherited from supplied coverage,not independently audited'))
  cov.append(dict(program_id=p,original_seed_count=len(genes),mapped_seed_count=len(genes)-len(un),measured_seed_count=len(members),absent_seed_count=len(absent),mapping_loss_count=len(un),coverage_fraction=len(members)/len(genes),exact_measured_members=';'.join(members),exact_missing_members=';'.join(absent),exact_unmapped_members=';'.join(un),completeness_status='DECLARED_AVAILABILITY_ONLY',eligibility='DETECTION_UNVERIFIED',detection_evidence='MISSING'))
  cells.append(pd.DataFrame(dict(cell_id=m.index,sample_id=m.sample_id.to_numpy(),program_id=p,n_seed_total=len(genes),n_seed_measured=len(members),n_seed_detected=np.nan,coverage_fraction=len(members)/len(genes),eligible_for_scoring=False if version=='original' and p.startswith('S') else pd.NA,reason_if_invalid='DETECTION_EVIDENCE_MISSING; not proof of<3 detection',raw_score_for_qc=raw[p].to_numpy() if version=='original' else np.nan,valid_program_score=bundle.scores[p].to_numpy() if version=='original' else np.nan)))
 pd.DataFrame(cat).to_csv(out/'program_seed_gene_qc.tsv',sep='\t',index=False);pd.DataFrame(cov).to_csv(out/'program_seed_coverage.tsv',sep='\t',index=False);pd.concat(cells,ignore_index=True).to_parquet(out/'cell_seed_detection.parquet',index=False)
