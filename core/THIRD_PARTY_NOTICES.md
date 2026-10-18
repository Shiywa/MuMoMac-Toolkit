# Third-party resources and notices

Project-authored software is MIT-licensed, copyright 2026 Shiyou Wang. Project-authored lightweight Reference inclusion is approved for the Toolkit under `REFERENCE_TERMS.md`. This notice retains separate upstream terms and attribution; third-party material is not relicensed under MIT.

## Mouse Genome Informatics

The frozen L/A human ortholog resource `core/src/mumomac/resources/human_la_one_to_one_v1.tsv` is derived from Mouse Genome Informatics (MGI), The Jackson Laboratory. MGI data and annotations are distributed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/), according to the [MGI copyright policy](https://www.informatics.jax.org/mgihome/other/copyright.shtml). No MGI images are distributed here.

The retained table selects required L/A symbols and excludes all globally ambiguous mouse–human relationships. Every packaged pair was checked against the retained MGI homology report; the mapping is not regenerated. Packaged SHA256: `60db7f137335dd270de2b83455bd9f1d056d4e1f492a3f23680ea5870856fdb2`. The report's exact original access time is not reconstructed from inconsistent historical labels. No sequences, images, or real expression data accompany this table. Human mapped-signature JSON resources retain this attribution for their MGI-derived symbols.

## Ensembl

The inherited S-program human mapping `core/src/mumomac/resources/human_one_to_one.tsv` is described in the frozen provenance as an explicit Ensembl BioMart one-to-one mapping. Its database release is **unknown**. The original manifest records `Missing`; public documentation must not infer a release from the current Ensembl website or retrieval date.

The [Ensembl data disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html) states that Ensembl-generated data are available without restriction, while incorporated third-party data can have additional restrictions. This retained table contains Ensembl-generated homology relationships, symbols, and gene IDs, verified as a subset of the original BioMart audit export; no sequences, images, or third-party research matrices are included. Packaged SHA256: `5b87f7f3e66433c014e024d8dd5c89f2e294057486372e11ef999c028161b9e7`. Filtering and project curation do not convert upstream data to MIT. Ensembl code licensing is distinct from data terms. The unversioned historical endpoint and missing database release remain explicit; neither a release nor current mapping is substituted. Human mapped-signature JSON resources retain this provenance.

## Python dependencies

Dependencies are obtained from their original public package distributions. They are not vendored in this source repository and are not relicensed by a future project license. The installed-version and license-metadata inventory is in `release_audit/DEPENDENCY_LICENSE_AUDIT.md` and `DEPENDENCY_LICENSES.tsv`. Metadata are recorded as supplied by the distributions, with missing or ambiguous values identified for review.
