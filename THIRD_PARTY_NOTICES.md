# Third-party resources and notices

This notice records upstream terms and attribution. It does not grant a public license to project-authored code or curated Reference resources. Those rights remain pending author approval, as described in `release_audit/PUBLICATION_RIGHTS_REVIEW.md`.

## Mouse Genome Informatics

The frozen L/A human ortholog resource `core/src/mumomac/resources/human_la_one_to_one_v1.tsv` is derived from Mouse Genome Informatics (MGI), The Jackson Laboratory. MGI data and annotations are distributed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/), according to the [MGI copyright policy](https://www.informatics.jax.org/mgihome/other/copyright.shtml). No MGI images are distributed here.

The retained table is a filtered, globally unique mouse–human mapping used by the frozen scoring resource. Its exact inherited file and hash are retained; the staging process does not replace or regenerate it. Attribution and the general upstream terms do not resolve the separate authorization of the project's curated resources or all historical snapshot documentation.

## Ensembl

The inherited S-program human mapping `core/src/mumomac/resources/human_one_to_one.tsv` is described in the frozen provenance as an explicit Ensembl BioMart one-to-one mapping. Its database release is **unknown**. The original manifest records `Missing`; public documentation must not infer a release from the current Ensembl website or retrieval date.

The [Ensembl data disclaimer](https://www.ensembl.org/info/about/legal/disclaimer.html) states that Ensembl-generated data are available without restriction, while incorporated third-party data can have additional restrictions. Ensembl's code license is separate from its data terms. This candidate preserves the mapping and flags the historical provenance gap; it does not assert that all upstream rights for the retained snapshot have been independently reconstructed.

## Python dependencies

Dependencies are obtained from their original public package distributions. They are not vendored in this source repository and are not relicensed by a future project license. The installed-version and license-metadata inventory is in `release_audit/DEPENDENCY_LICENSE_AUDIT.md` and `DEPENDENCY_LICENSES.tsv`. Metadata are recorded as supplied by the distributions, with missing or ambiguous values identified for review.
