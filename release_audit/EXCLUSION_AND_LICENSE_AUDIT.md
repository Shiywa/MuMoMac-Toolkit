# Exclusions and publication boundary

Included: the minimal two-package source closure; frozen TSV/JSON resources required at runtime; preserved inference engines; synthetic-example generator; portable synthetic tests; US English input/output/QC documentation; citation metadata; build/CI definitions; redacted audit summaries.

Excluded: the atlas, all real H5AD matrices, real sample metadata and engineering exports, manuscript/figures, private analysis history, developer environments, editable-install metadata, caches, credentials, downloaded dependency wheels and original private regression fixtures. Synthetic example data are generated outside the repository rather than committed.

The public test files retain assertions from portable C3/C3.1 tests. Private-fixture tests and original-tree pin-location checks are excluded only from the portable suite; the complete original C1/C2/C2.1/C3/C3.1 suites ran internally (189 passed). Installed resource hashes are tested separately. The new missing-resource harness test initially failed because it assumed an absent cache interface; that failure and its correction are retained in the internal evidence.

Only source-path provenance fields in two manifests were sanitized. Frozen member lists, maps, canonical registry, thresholds and executable scientific code are unchanged. Human mapping provenance digests necessarily change because the adapter hashes a sanitized manifest; exact human numeric/NA parity is separately recorded. No historical human score-run hash compatibility is asserted.

Code-license selection and curated Reference redistribution remain pending. The source trees retain their internal notice; no public LICENSE file, tag, Release, package-index publication, or public-data redistribution is authorized by this staging audit.
