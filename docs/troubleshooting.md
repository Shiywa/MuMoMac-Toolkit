# Troubleshooting

**A program is NA.** Check measured members, mapping loss, per-cell nonzero Seeds, identity applicability, and the relevant score-version policy. Do not replace NA with zero or broaden a frozen signature to force a result.

**Human CosMx yields few Expanded scores.** A targeted panel can measure only a subset of the frozen signatures. Successful execution does not imply valid biological coverage. Retain the non-evaluable output and mapping provenance.

**Why no automatic macrophage discovery?** The supported contract requires reviewed Provided labels. Upstream auto-assisted work is an independent research procedure and must be reviewed before its labels are supplied here.

**Why are Joint/Lift/OR unavailable?** Continuous scores alone do not define a complete binary policy. A score threshold, reference population, applicability and table support are required. A01 cannot be binary.

**Metadata join fails.** Supply unique lexical `cell_id` values matching all H5AD observations exactly. Order can differ. Duplicate, missing or extra IDs are errors rather than silently dropped rows.

**Counts are rejected.** Verify nonnegative finite integers in the explicitly declared layer. Declare a documented non-count expression type only when appropriate; do not relabel invalid counts to bypass validation.

**Ensembl or symbols do not map.** Declare the ID type or an explicit author symbol column. The frozen Ensembl map is partial, and IDs are not guessed by case conversion. Resolve duplicate/ambiguous symbols upstream without choosing arbitrary ortholog edges.

**Import says the core is missing.** Install both built wheels in the same environment. The core is a separate pinned distribution. No private project directory or PYTHONPATH should be needed.

**The output directory is not empty.** Use a new directory. Failed/partial artifacts are retained and are not silently overwritten.

**Can I compare raw means across versions, batches, or studies?** No universal calibration is supplied. Keep definitions separate, retain source/sample structure, inspect technical dependence, and use a design appropriate to the biological unit.
