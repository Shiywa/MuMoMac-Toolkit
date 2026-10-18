# Changelog

The staged packaging revision sanitizes private source-path provenance without changing scoring code or Reference members. Human mapping digests include the sanitized manifest and therefore differ from the private-source digest; numerical/NA parity is assessed separately. No human score-run hash compatibility is inferred.

## 0.1.1rc1 public-release candidate staging

The actual validated workflow version remains 0.1.1rc1, paired with scoring core 0.3.0 and Reference 1.0.0. No public tag has been created. The planned prerelease tag is `v0.1.1rc1` after the release gates pass.

Packaging changes add an isolated two-distribution source layout, synthetic example, portable regression tests, public-facing input/output documentation, and CI configuration. Private absolute source paths are removed only from two non-scientific provenance fields. Scientific member sets, maps, score mathematics, inference engines, and canonical Reference hashes are unchanged. The internal rights notice is retained pending author approval; no MIT or Apache-2.0 license has been selected by the packager.
