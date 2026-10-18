# MuMoMac Toolkit

MuMoMac profiles 22 layered macrophage transcriptional features from user-supplied H5AD expression and explicit monocyte/macrophage annotations. The workflow exports versioned scores, measured-gene and per-cell detection QC, and continuous program associations within sample/identity groups.

**[Public prerelease 0.1.1rc1](https://github.com/Shiywa/MuMoMac-Toolkit/releases/tag/v0.1.1rc1), published 2026-10-18: scoring core 0.3.0 and Reference 1.0.0. MIT project software, copyright 2026 Shiyou Wang; lightweight resources retain separate scope and upstream notices. GitHub CI and anonymous fresh-install acceptance passed.**

The reference has 16 S functional transcriptional programs, five L lifecycle/residence-associated features, and continuous A01. Original contains 229 memberships/229 distinct genes. Current Expanded contains 820 memberships/715 distinct genes. Current Expanded22 is separate from the historical 601-gene expansion.

## Install the public prerelease

Use Python 3.11 on Ubuntu. In your Python environment, run one command:

```bash
python -m pip install -r https://raw.githubusercontent.com/Shiywa/MuMoMac-Toolkit/main/requirements-release.txt
```

Pip downloads and installs both required packages from the public `v0.1.1rc1` Release, verifies their published SHA256 hashes, and installs external dependencies from PyPI. No manual package download, GitHub login, Git executable, or Atlas download is required. For installation directly from a fixed GitHub tag, a local source checkout, or manual artifact verification, see [installation details](docs/installation.md).

To create a separate environment first, run `python3.11 -m venv .venv` and `source .venv/bin/activate`, then run the install command above. The project has not been published to PyPI; use the public requirements URL shown here.

## Quickstart

Download the synthetic example generator from the fixed release tag and run the workflow:

```bash
curl -fL https://raw.githubusercontent.com/Shiywa/MuMoMac-Toolkit/v0.1.1rc1/examples/make_synthetic_example.py -o make_synthetic_example.py
python make_synthetic_example.py --output example-data
mumomac-experimental analyze --input example-data/synthetic_counts.h5ad --metadata example-data/provided_metadata.tsv --config example-data/config.yaml --signature-version both --output example-results
```

Check the installation with `python -m pip check` or view available commands with `mumomac-experimental --help`. See [the complete quickstart](docs/quickstart.md) for seed details, replay, and separate Original/Expanded runs.

The synthetic input has 256 artificial profiles and 795 native genes; Provided labels select 240 targets. Scores and QC are in `example-results/original/` and `example-results/expanded/`; root tables preserve version-prefixed score namespaces. `SYN0000` has insufficient Original S03 detection and correctly retains NA. Missing Joint/Lift/OR values do not mean a joint population is absent.

The tested platforms are Python 3.11 on Ubuntu 20.04 and Ubuntu 24.04. Windows, macOS, and million-cell operation are not certified. The [public verification summary](release_audit/PUBLIC_RELEASE_VERIFIED.json) records the fixed software release and asset hashes; [GitHub Actions](https://github.com/Shiywa/MuMoMac-Toolkit/actions/workflows/ci.yml) checks builds, clean installs, and synthetic scoring.

## Inputs and results

H5AD matrices are cells × genes. Declare the expression layer, expression type, species, modality, sample column, provided cell-type column, and exact target labels in YAML. External TSV/CSV/Parquet metadata must contain unique `cell_id` values matching the full H5AD observation set; order may differ. Non-target cells are recorded in the selection audit. The native background gene axis is retained for scoring. Never pad unmeasured genes and present them as measured zeros.

Both mode produces separate `original/` and `expanded/` runs, plus a combined score table with `original__S01` / `expanded__S01` namespaces. Outputs include `cell_program_scores.parquet`, `cell_metadata.parquet`, `cell_seed_detection.parquet`, `program_seed_coverage.tsv`, `program_seed_gene_qc.tsv`, `score_eligibility_by_sample.tsv`, `sample_coupling_statistics.parquet`, `gene_set_manifest.json`, `ortholog_mapping_manifest.json`, `software_versions.json`, and `run_qc.json`. See [output schemas](docs/output_schema.md) for scopes, keys, and NA semantics. The inherited HTML/Markdown reports are auxiliary; structured numerical exports are the primary contract.

## Interpretation and limitations

Scores describe transcriptional features, rather than protein activity, functional positivity, antiviral protection, or a diagnosis. Original S qualification requires at least three measured seeds and three detected seeds per cell; this minimum is a technical gate, not an accuracy guarantee. Invalid scores remain NA. Expanded requires its frozen coverage rules and may be largely non-evaluable in targeted human spatial panels. Human mapping is experimental; the inherited S ortholog database release is unknown and remains explicitly recorded.

Provided annotations are the supported default. General all-cell automated macrophage identification is not certified by this package. Mixed spatial spots are refused. Cell-resolved spatial input requires declared segmentation and finite coordinates; the workflow does not perform neighborhood inference. Atlas lookup is descriptive and does not define healthy-normal z scores.

Continuous associations retain the actual finite-pair denominator, technical adjustments, and unavailable reasons within sample/identity groups. The synthetic quickstart explicitly selects `continuous_only`. No valid binary policy means Joint/Lift/OR remain unavailable; NA is never converted to zero or Neither. A01 is never binary. Sample IDs are not automatically animal or donor IDs, and the default workflow does not generate disease P values from cells.

The four historical engineering inputs were already used during development and are not new independent biological validation. Million-cell scalability, macOS, and Windows are not certified. See [limitations](docs/limitations.md), [QC](docs/coverage_and_qc.md), and [troubleshooting](docs/troubleshooting.md).

## Reference, citation, and rights

Scientific gene sets, mapping tables, score mathematics, and frozen inference engines are preserved. Only non-scientific source-path provenance and publication-rights fields in two manifests are updated for distribution. The canonical Reference registry SHA256 is `2c3ff81f43010caf11251053ec676d7d03ded3f941b9ac1660d836b5beeee515`. Project-authored software and documentation use the [MIT License](LICENSE), copyright 2026 Shiyou Wang. See [license scope](LICENSE_SCOPE.md), [Reference terms](REFERENCE_TERMS.md), and [third-party notices](THIRD_PARTY_NOTICES.md) for the separate resource permissions. Dependencies retain their own licenses and are installed from their original distributions. The Atlas is not included and its future hosting platform is undecided.

Human provenance compatibility requires attention: the inherited adapter includes the complete `signatures_manifest.json` and `reference_version_manifest.json` in its mapping digest. Source-path sanitization and current publication-rights metadata therefore change the human provenance digest, even though mapping tables, program members, scores, and NA masks are preserved. Historical human score-run specifications cannot be assumed hash-identical to this candidate. The legacy S-mapping database release remains unknown. See [third-party notices](THIRD_PARTY_NOTICES.md) for source attribution.

Cite Shiyou Wang, MuMoMac Toolkit 0.1.1rc1, [fixed release](https://github.com/Shiywa/MuMoMac-Toolkit/releases/tag/v0.1.1rc1) and [software source commit](https://github.com/Shiywa/MuMoMac-Toolkit/commit/61880256cd085296c38f04b05cdc34875fdbbde7). `CITATION.cff` supplies editable author/version metadata. No Zenodo DOI is assigned. The released archives retain that fixed source snapshot; subsequent documentation updates on main do not move the release tag or replace its assets. Earlier audit files describe their labeled historical staging checkpoint.
