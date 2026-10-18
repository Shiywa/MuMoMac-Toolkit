# Output schema and NA semantics

Single-definition runs place files in the requested output directory. Both mode places definition-specific files in `original/` and `expanded/`, and also writes combined summaries at the root. Do not read the root Both table as one shared score scale.

| File | Row/key and meaning |
| --- | --- |
| `selection_audit.parquet` | Input cell ID, original annotation, selection flag and reason, including non-target cells |
| `gene_harmonization.parquet` | Every native gene ID, scoring symbol and mapping status; background retained |
| `cell_program_scores.parquet` | Selected cell ID, sample/identity and score columns; single mode has S01–S16/L01–L05/A01, Both has version-prefixed score columns |
| `cell_metadata.parquet` | Selected cell ID and supplied/normalized sample, identity, source and biological-unit metadata |
| `raw_cell_program_scores_for_qc.parquet` | Raw audit values; availability does not imply Strict eligibility |
| `cell_seed_detection.parquet` | Cell × Original S feature detection evidence, measured/detected counts, eligibility and failure reason |
| `program_seed_coverage.tsv` | Original seeds, mapped/measured/absent members, coverage and definition namespace |
| `program_seed_gene_qc.tsv` | Program/member mapping, dataset measurement and detection evidence |
| `program_coverage.tsv` | Per-definition 22-feature qualification and measured-member information |
| `score_eligibility_by_sample.tsv` | Actual qualified counts/fractions per exported sample/identity/program/version grouping |
| `sample_coupling_statistics.parquet` | Per-sample/identity pair, raw and adjusted association, actual pair denominator and binary availability |
| `selected_pair_joint_features.parquet` | Cell-level selected pair values; state label is unavailable when no complete rule exists |
| `pairwise_gene_overlap.tsv` | Defining-pair shared genes and overlap measures |
| `technical_dependencies.tsv` | Descriptive score association with count/detection covariates |
| `cell_program_representation.parquet` | Cell/program continuous representation, optional explicitly qualified state, QC and applicability |
| `score_spec.json` | Definition namespace, program hashes, scoring parameters and input gene-universe hash |
| `gene_set_manifest.json` | Reference version, feature order, signature version and program hashes |
| `ortholog_mapping_manifest.json` | Mapping hashes/loss and known provenance gaps |
| `software_versions.json` | Workflow/core/dependency versions and frozen engine integrity |
| `run_qc.json` | Input/selected counts, actual units, qualification, state availability, performance and limitations |
| `seed_qc_manifest.json`, `seed_detection_bits.npz` | Auditable seed evidence and compact detection storage; use `seed-details` for exact lists |

Files are documented according to their actual creation scope; some single-definition files are not duplicated at the Both root. The smoke test records exact columns and dtypes in `release_audit/OBSERVED_OUTPUT_SCHEMAS.json`.

NA means the value cannot be assessed under the declared measurement, detection, applicability, or calibration rule. Measured zero, not measured, unmapped, invalid, and valid zero are different. In particular, absent detection evidence never certifies a zero detected count. No complete binary policy means Joint fraction/Lift/OR remain NA. Their missingness is not evidence that a joint population is absent. A01's discrete state is always unavailable.

No new scored H5AD writer is claimed by the wrapper. The supported wrapper exports are the actual TSV/Parquet/JSON/NPZ/report files above. Join scores back to the source AnnData by `cell_id`; unselected observations remain missing.
