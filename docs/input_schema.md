# Input schema

`analyze` accepts `.h5ad` expression only. Matrix orientation is observations × native genes. Observation and gene IDs must be unique and nonempty. `counts` must be finite nonnegative integers; `log1p` and `nonnegative_expression` must be finite and nonnegative. The layer is explicitly named (`counts` or `X`, for example). Scaled or corrected expression is not inferred to be suitable.

The YAML requires `species` (`mouse` or `human`), `modality` (`scRNA`, `scRNA-seq`, `snRNA`, `snRNA-seq`, or `cell_resolved_spatial`), `layer`, `expression_type`, `sample_column`, `cell_type_column`, and nonempty `target_labels`. `selection_mode` is `provided`. Optional `identity_column`, `study_column`, `source_column`, `condition_column`, and `tissue_column` map existing metadata to standard fields. `biological_unit_column` and `biological_unit_confirmed_column` are explicit; sample IDs are never treated as confirmed animals automatically.

External `.tsv`, `.csv`, or `.parquet` metadata must have one unique lexical `cell_id` for every observation in the H5AD. No silent inner join is used. Matching rows are reordered by ID, not position. External values take precedence, and conflicting original observation values are retained under `input_obs__...`.

Native gene symbols are used without case conversion. `gene_symbol_key` can identify an explicit symbol column. `gene_id_type: ensembl` uses the frozen partial Ensembl mapping only; it is not a genome-wide annotation service. Version suffixes are stripped in that explicit mode. Ambiguous or duplicate mapped symbols are refused. Unmapped genes remain as opaque background, rather than being deleted. `unmapped_gene_ids` must name actual native-axis IDs.

Cell-resolved spatial data additionally require `segmentation_confirmed: true` and finite coordinates in the configured `coordinate_key` (default `spatial`). `mixed_spots` is rejected. Coordinates are exported by cell ID without assigning physical units. A score-only workflow requires complete matching score/metadata tables, the exact ScoreSpec, coverage, and trustworthy detection evidence; the `score-run` replay route is preferred for existing package exports.
