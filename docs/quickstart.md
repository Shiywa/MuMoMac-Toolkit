# Quickstart

Use Python 3.11 on Ubuntu. Optionally create and activate a new environment with `python3.11 -m venv .venv` and `source .venv/bin/activate`.

## Install

```bash
python -m pip install -r https://raw.githubusercontent.com/Shiywa/MuMoMac-Toolkit/main/requirements-release.txt
```

This installs both required distributions from the fixed public release and checks their wheel hashes automatically. The full Atlas is not required. [Other installation routes](installation.md) include a single pip invocation using two fixed-tag GitHub source requirements.

## Generate input and score

```bash
curl -fL https://raw.githubusercontent.com/Shiywa/MuMoMac-Toolkit/v0.1.1rc1/examples/make_synthetic_example.py -o make_synthetic_example.py
python make_synthetic_example.py --output example-data
mumomac-experimental analyze --input example-data/synthetic_counts.h5ad --metadata example-data/provided_metadata.tsv --config example-data/config.yaml --signature-version both --output example-results
```

Each definition has its own `original/` or `expanded/` directory with `cell_program_scores.parquet`, coverage, detection QC, and sample/identity-level continuous coupling. The root also contains version-prefixed combined outputs. The example selects 240 target profiles; 22 features are retained, with invalid scores represented as NA. No valid binary policy is selected, so Joint/Lift/OR remain unavailable. A01 remains continuous.

## Inspect seeds and replay

```bash
mumomac-experimental seed-details --run example-results --cell-id SYN0000 --program-id S03
mumomac-experimental score-run --input-run example-results --output replay-results
```

`SYN0000` intentionally has insufficient Original S03 detection and retains NA. Replay uses the saved specification. Existing nonempty output directories are refused.

## Run each definition separately

```bash
mumomac-experimental analyze --input example-data/synthetic_counts.h5ad --metadata example-data/provided_metadata.tsv --config example-data/config.yaml --signature-version original --output original-results
mumomac-experimental analyze --input example-data/synthetic_counts.h5ad --metadata example-data/provided_metadata.tsv --config example-data/config.yaml --signature-version expanded --output expanded-results
```

Run `python -m pip check`, `mumomac --help`, or `mumomac-experimental --help` to inspect the installation. See [output schemas](output_schema.md), [input schema](input_schema.md), and [QC/limitations](qc_and_limitations.md) before using your own H5AD.

The released software source is `61880256cd085296c38f04b05cdc34875fdbbde7`, fixed by `v0.1.1rc1`. These installation and documentation updates preserve the software tag and existing release assets.
