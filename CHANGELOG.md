# Changelog

## Simplified public installation — October 18, 2026

Added a public requirements manifest that installs both fixed release wheels in one pip command and verifies their SHA256 hashes. The README now starts with this install command and a short synthetic quickstart. Fixed-tag GitHub installation and manual artifact auditing are documented separately. Scientific code, Reference resources, package versions, the release tag, and release assets are unchanged.

## Public verification and documentation — October 18, 2026

Published `v0.1.1rc1` from source commit `61880256cd085296c38f04b05cdc34875fdbbde7` after successful GitHub CI. Anonymous downloads, a new virtual environment, the full synthetic quickstart, exact outputs, and asset hashes passed. The verified public download instructions and citation metadata were then added on main without moving the software tag or replacing release assets.

## Publication packaging revision — October 18, 2026

MIT project software approved with copyright holder Shiyou Wang. Lightweight Reference scope and upstream mapping notices are separate. Code/license metadata and publication-rights statements are updated, CI is named `ci.yml`, and a combined QC guide is added. Frozen scientific code, members, thresholds, registry, and maps are unchanged. Online publication, CI, and anonymous public installation are still pending.

The staged packaging revision sanitizes private source-path provenance without changing scoring code or Reference members. Human mapping digests include the sanitized manifest and therefore differ from the private-source digest; numerical/NA parity is assessed separately. No human score-run hash compatibility is inferred.

## 0.1.1rc1 public-release candidate staging

The actual validated workflow version remains 0.1.1rc1, paired with scoring core 0.3.0 and Reference 1.0.0. No public tag has been created. The planned prerelease tag is `v0.1.1rc1` after the release gates pass.

Packaging changes add an isolated two-distribution source layout, synthetic example, portable regression tests, public-facing input/output documentation, and CI configuration. Private absolute source paths are removed only from two non-scientific provenance fields. Scientific member sets, maps, score mathematics, inference engines, and canonical Reference hashes are unchanged. The internal rights notice is retained pending author approval; no MIT or Apache-2.0 license has been selected by the packager.
