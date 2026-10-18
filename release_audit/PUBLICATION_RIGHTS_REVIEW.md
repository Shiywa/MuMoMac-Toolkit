# Publication rights review

Status: pending author/rights-holder decisions. No public redistribution license has been granted by this staging task. The inherited internal-research restriction remains in `INTERNAL_NOTICE.txt`; no public LICENSE is supplied.

| Asset | Status and required action |
| --- | --- |
| Project-authored scoring/integration/continuous engines | Owner must select and authorize a public code license, such as MIT or Apache-2.0. No license is chosen by the packager. |
| Original and Expanded gene-member registries | Owner must confirm public redistribution scope and any third-party contribution/curation restrictions. Scientific definitions remain frozen. |
| `human_one_to_one.tsv` and human-mapped S resources | Inherited explicit BioMart mapping with an unknown database release; confirm the retained mapping's provenance and redistribution terms. No new mapping is substituted. |
| `human_la_one_to_one_v1.tsv` and human-mapped L/A resources | Inherited MGI snapshot provenance; confirm redistribution terms and attribution for the actual retained snapshot. |
| Historical thresholds and identity/reference metadata | Internal evaluation rights do not establish public distribution rights. Confirm included records are cleared. |
| Dependencies | Installed from their public distributions, not relicensed or vendored under a future project license. License inventory is separate. |
| Synthetic example | Newly generated counts and artificial IDs, no real source data; author public license decision must also cover the example code. |
| Atlas matrices, private metadata, manuscript, real engineering inputs/exports and caches | Excluded from this candidate. No permission to redistribute them is inferred. |

Source availability, GEO accessibility and public paper access do not by themselves resolve all redistribution questions. The review identifies exact resources rather than inventing license approval. The planned public repository is a new independent `Shiywa/MuMoMac-Toolkit`; the existing private analysis repository is outside the publication scope.

General upstream terms have been checked: MGI data/annotations use CC BY 4.0, and Ensembl-generated data are unrestricted with a third-party-data caveat. Attribution and source links are recorded in [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md). These findings narrow the review but do not approve the project-authored code/curation or reconstruct the unknown historical S-mapping release.
