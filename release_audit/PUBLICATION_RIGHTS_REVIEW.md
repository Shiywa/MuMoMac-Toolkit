# Publication rights review

Status: project software MIT approved, copyright 2026 Shiyou Wang; lightweight runtime resources approved within the Toolkit scope. This revision supersedes the earlier pending-code-license staging status. The current scope is in LICENSE_SCOPE.md and REFERENCE_TERMS.md; upstream materials keep separate terms.

| Asset | Status and required action |
| --- | --- |
| Project-authored scoring/integration/continuous engines | MIT approved by the author; exact individual copyright name supplied as Shiyou Wang. No third-party scoring implementation is vendored. |
| Original and Expanded gene-member registries | Author-approved lightweight Toolkit inclusion; frozen member lists and canonical registry preserved. Separate Reference scope, not blanket MIT data licensing. |
| `human_one_to_one.tsv` and human-mapped S resources | 478 packaged mapping pairs verified against the original BioMart audit export. Ensembl data terms recorded; only generated homology facts/symbols/IDs. Database release remains unknown. No new mapping substituted. |
| `human_la_one_to_one_v1.tsv` and human-mapped L/A resources | 150 packaged pairs verified against the retained MGI homology report; CC BY 4.0, JAX/MGI attribution, filtering and global-ambiguity exclusions documented. |
| Historical thresholds and identity/reference metadata | Small aggregate legacy threshold table and 18-name identity vocabulary; no individual records or expression matrices. Approved only within lightweight Toolkit resource scope. |
| Dependencies | Installed from their public distributions, not relicensed or vendored under a future project license. License inventory is separate. |
| Synthetic example | Newly generated counts and artificial IDs, no real source data; project generator and synthetic tests use MIT. |
| Atlas matrices, private metadata, manuscript, real engineering inputs/exports and caches | Excluded from this candidate. No permission to redistribute them is inferred. |

Source availability, GEO accessibility and public paper access do not by themselves resolve all redistribution questions. The review identifies exact resources rather than inventing license approval. The planned public repository is a new independent `Shiywa/MuMoMac-Toolkit`; the existing private analysis repository is outside the publication scope.

General upstream terms and the actual retained mapping ancestry have been checked: MGI data/annotations use CC BY 4.0, and Ensembl-generated data are unrestricted with a third-party-data caveat. The packaged fields contain no copied sequences, images, or third-party research matrices. Attribution and source links are recorded in [THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md). Unknown historical S-mapping release provenance remains a documented limitation, not an invented version.
