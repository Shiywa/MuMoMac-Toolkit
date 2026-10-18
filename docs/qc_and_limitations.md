# QC and limitations

Use the [coverage and QC guide](coverage_and_qc.md), [coupling/state rules](coupling_and_states.md), and [limitations](limitations.md) with every analysis. Provided annotations, declared measurement semantics, and explicit biological-unit metadata are required. The supported quickstart does not require the Atlas.

Original S features require at least three measured and three detected seeds per cell. This is technical eligibility, not functional activation or an accuracy guarantee. Keep valid scores separate from raw QC and NA. Measured zero, absent measurement, unmapped genes, and missing detection evidence have different meanings. Expanded follows its own frozen coverage rules; it cannot be filled from Original. A01 remains continuous.

Continuous associations use actual sample × identity finite-pair denominators and technical adjustment. They do not establish mechanisms or cell-level disease significance. Without a qualified binary policy, Joint/Lift/OR are unavailable. Human targeted panels can have little usable Expanded coverage; the legacy S ortholog database release remains unknown. Mixed spots are refused and cell-resolved spatial data require explicit segmentation and coordinate declarations.

The synthetic example tests execution and boundary behavior. It is not a biological validation cohort. General automatic identification, normal/abnormal z scores, protein activity, full 22-feature coverage on arbitrary panels, million-cell scalability, macOS, and Windows are not certified.
