# Program reference

Reference 1.0.0 contains S01–S16 functional transcriptional programs, L01–L05 lifecycle/residence-associated features, and continuous A01. Original has 229 memberships/229 unique genes. Current Expanded has 820 memberships/715 unique genes. The historical 601-gene expansion is a distinct namespace.

Original uses frozen Seed UCell scores. Expanded S uses `0.5 * Seed UCell + 0.5 * Added UCell`, with its frozen no-Added fallback. Expanded L/A01 use their later positive Full sets. UCell is official pyUCell 0.8.0 with average ties and maximum rank 1500. No member, coefficient, threshold, or calibrated historical call is retrained at runtime.

L04 is the Monocyte-to-macrophage transition-associated feature; it is not a monotonic maturation scale. L05 is an auxiliary mature macrophage feature. L scores have identity/applicability limitations. A01 is an RNA activity covariate, not a direct translation assay and never a binary state.

The frozen Reference registry and source-member hashes are shipped with the scoring core. A release audit checks their exact hashes. Original/Expanded absolute values, historical mean/q20 values, and historical binary thresholds must not be pooled or directly transferred.
