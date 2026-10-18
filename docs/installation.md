# Installation details

The recommended install is one pip command in a Python 3.11 environment on Ubuntu:

```bash
python -m pip install -r https://raw.githubusercontent.com/Shiywa/MuMoMac-Toolkit/main/requirements-release.txt
```

The public manifest lists both required wheels from `v0.1.1rc1`, with their published SHA256 values in URL fragments. Pip checks these hashes automatically and resolves external dependencies from PyPI. This manifest pins the two MuMoMac distributions; external dependencies follow their package metadata. The full Atlas is not required. No GitHub login or Git executable is needed for this route. The `main` URL provides the current installation manifest; the wheel URLs and hashes inside it are fixed to this prerelease.

## Install from a fixed GitHub tag

If Git is available, pip can build both distributions directly from the fixed public tag in one invocation:

```bash
python -m pip install \
  "mumomac-toolkit @ git+https://github.com/Shiywa/MuMoMac-Toolkit.git@v0.1.1rc1#subdirectory=core" \
  "mumomac-experimental @ git+https://github.com/Shiywa/MuMoMac-Toolkit.git@v0.1.1rc1"
```

The core is located in `core/`; the workflow package is at the repository root. Supply both requirements in the same pip call so pip can satisfy the exact core dependency. Installing the root alone does not supply the core. Both routes install the same released scoring code and Reference resources.

The public wheel and fixed-tag routes are validated separately. The longer build-and-install procedure below is retained for artifact auditing and development.

The project has not been published to PyPI. A bare `pip install mumomac-toolkit` is not the supported installation command for this release. See [pip's VCS documentation](https://pip.pypa.io/en/stable/topics/vcs-support/) for the source URL and `subdirectory` syntax.

## Download and verify release artifacts

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

Fresh wheel and sdist installs, the complete anonymous public-asset run, and [GitHub CI](https://github.com/Shiywa/MuMoMac-Toolkit/actions/runs/39130124445) passed. The Ubuntu 24.04/Python 3.11 CI and local Ubuntu 20.04/Python 3.11 checks do not certify Windows, macOS, or million-cell operation. The [public verification summary](../release_audit/PUBLIC_RELEASE_VERIFIED.json) records the release-source commit and asset hashes.

## Build from source

The release consists of two distributions released together: `mumomac-toolkit==0.3.0` in `core/`, and `mumomac-experimental==0.1.1rc1` at the repository root. Both must be installed. The workflow's dependency on the core is exact; no local private wheel path is required. The [public verification summary](../release_audit/PUBLIC_RELEASE_VERIFIED.json) records online release acceptance; [historical staging results](../release_audit/INSTALLATION_STATUS.md) are retained separately.

The tested platform is Python 3.11 on Ubuntu 20.04 x86_64; CI passed on Ubuntu 24.04/Python 3.11. To build from source, build both distributions, then install both wheels into a fresh environment. Start inside a source checkout of tag `v0.1.1rc1`:

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
