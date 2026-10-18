# MuMoMac Toolkit

MuMoMac profiles 22 layered macrophage transcriptional features from user-supplied H5AD expression and explicit monocyte/macrophage annotations. The workflow exports versioned scores, measured-gene and per-cell detection QC, and continuous program associations within sample/identity groups.

**[Public prerelease 0.1.1rc1](https://github.com/Shiywa/MuMoMac-Toolkit/releases/tag/v0.1.1rc1), published 2026-10-18: scoring core 0.3.0 and Reference 1.0.0. MIT project software, copyright 2026 Shiyou Wang; lightweight resources retain separate scope and upstream notices. GitHub CI and anonymous fresh-install acceptance passed.**

The reference has 16 S functional transcriptional programs, five L lifecycle/residence-associated features, and continuous A01. Original contains 229 memberships/229 distinct genes. Current Expanded contains 820 memberships/715 distinct genes. Current Expanded22 is separate from the historical 601-gene expansion.

## Install the public prerelease

Use Python 3.11 on Ubuntu. Download both packages' wheels, their source distributions, checksums, and the fixed-tag synthetic example into a new directory outside any source checkout. The following downloader was executed anonymously and its outputs matched the published hashes:

```bash
mkdir mumomac-0.1.1rc1-work
cd mumomac-0.1.1rc1-work
python3.11 - <<'PY'
from urllib.request import urlretrieve

release = "https://github.com/Shiywa/MuMoMac-Toolkit/releases/download/v0.1.1rc1"
files = [
    "mumomac_toolkit-0.3.0-py3-none-any.whl",
    "mumomac_experimental-0.1.1rc1-py3-none-any.whl",
    "mumomac_toolkit-0.3.0.tar.gz",
    "mumomac_experimental-0.1.1rc1.tar.gz",
    "SHA256SUMS",
]
for name in files:
    urlretrieve(f"{release}/{name}", name)
urlretrieve(
    "https://raw.githubusercontent.com/Shiywa/MuMoMac-Toolkit/v0.1.1rc1/examples/make_synthetic_example.py",
    "make_synthetic_example.py",
)
PY
sha256sum --check SHA256SUMS
python3.11 -m venv venv
venv/bin/python -m pip --isolated install --index-url https://pypi.org/simple mumomac_toolkit-0.3.0-py3-none-any.whl mumomac_experimental-0.1.1rc1-py3-none-any.whl
venv/bin/python -m pip check
venv/bin/python make_synthetic_example.py --output example-data
venv/bin/mumomac-experimental --help
venv/bin/mumomac-experimental analyze --help
venv/bin/mumomac-experimental analyze --input example-data/synthetic_counts.h5ad --metadata example-data/provided_metadata.tsv --config example-data/config.yaml --signature-version both --output example-results
venv/bin/mumomac-experimental seed-details --run example-results --cell-id SYN0000 --program-id S03
venv/bin/mumomac-experimental score-run --input-run example-results --output replay-results
venv/bin/mumomac-experimental analyze --input example-data/synthetic_counts.h5ad --metadata example-data/provided_metadata.tsv --config example-data/config.yaml --signature-version original --output original-results
venv/bin/mumomac-experimental analyze --input example-data/synthetic_counts.h5ad --metadata example-data/provided_metadata.tsv --config example-data/config.yaml --signature-version expanded --output expanded-results
```

The two wheel filenames belong in the same pip call. MuMoMac itself is distributed through this GitHub Release; the index supplies external dependencies. No Atlas download is required. The synthetic input has 256 artificial profiles and 795 native genes; Provided labels select 240 targets. Structured scores and QC are in `example-results/original/` and `example-results/expanded/`; root tables preserve version-prefixed score namespaces. `SYN0000` has insufficient Original S03 detection and correctly retains NA. Missing Joint/Lift/OR values do not mean a joint population is absent.

Fresh wheel and sdist installs, the complete anonymous public-asset run, and [GitHub CI](https://github.com/Shiywa/MuMoMac-Toolkit/actions/runs/39130124445) passed. The Ubuntu 24.04/Python 3.11 CI and local Ubuntu 20.04/Python 3.11 checks do not certify Windows, macOS, or million-cell operation. The [public verification summary](release_audit/PUBLIC_RELEASE_VERIFIED.json) records the release-source commit and asset hashes.

## Build from source

The candidate consists of two distributions released together: `mumomac-toolkit==0.3.0` in `core/`, and `mumomac-experimental==0.1.1rc1` at the repository root. Both must be installed. The workflow's dependency on the core is exact; no local private wheel path is required. The tested environment and real install outcomes are recorded in `release_audit/INSTALLATION_STATUS.md`.

The tested platform is Python 3.11 on Ubuntu 20.04 x86_64; CI passed on Ubuntu 24.04/Python 3.11. For this local candidate, build both distributions, then install both wheels into a fresh environment. Start inside the candidate source directory:

```bash
MUMOMAC_SOURCE_DIR="$PWD"
python3.11 -m venv ../mumomac-build-env
../mumomac-build-env/bin/python -m pip install --index-url https://pypi.org/simple build 'setuptools>=77' wheel
../mumomac-build-env/bin/python -m build --outdir dist core
../mumomac-build-env/bin/python -m build --outdir dist .
python3.11 -m venv ../mumomac-run-env
../mumomac-run-env/bin/python -m pip install --index-url https://pypi.org/simple dist/mumomac_toolkit-0.3.0-py3-none-any.whl dist/mumomac_experimental-0.1.1rc1-py3-none-any.whl
MUMOMAC_PYTHON="$MUMOMAC_SOURCE_DIR/../mumomac-run-env/bin/python"
MUMOMAC_EXAMPLE_SCRIPT="$MUMOMAC_SOURCE_DIR/examples/make_synthetic_example.py"
mkdir ../mumomac-example-work
cd ../mumomac-example-work
```

The preceding commands set the script/environment paths and move outside the source tree. Continue in that working directory:

```bash
"$MUMOMAC_PYTHON" "$MUMOMAC_EXAMPLE_SCRIPT" --output example-data
"$MUMOMAC_PYTHON" -m mumomac_experimental.cli --help
"$MUMOMAC_PYTHON" -m mumomac_experimental.cli analyze --help
"$MUMOMAC_PYTHON" -m mumomac_experimental.cli analyze --input example-data/synthetic_counts.h5ad --metadata example-data/provided_metadata.tsv --config example-data/config.yaml --signature-version both --output example-results
"$MUMOMAC_PYTHON" -m mumomac_experimental.cli seed-details --run example-results --cell-id SYN0000 --program-id S03
"$MUMOMAC_PYTHON" -m mumomac_experimental.cli score-run --input-run example-results --output replay-results
```

The `mumomac-experimental` executable installed in the environment provides the same CLI. To run one definition independently, change `--signature-version` to `original` or `expanded` and use a new output directory. Existing nonempty output directories are refused; failed runs are retained.

The fixed public Release and anonymous installation are verified above. A source checkout or ZIP alone does not replace the wheel/sdist checks. PyPI, Conda, and BioConda publication are outside this release.

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
