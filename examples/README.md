# Synthetic example

`make_synthetic_example.py` creates 256 artificial count profiles, 240 provided target labels, four sample records, a shuffled metadata table, and a complete YAML configuration. The generator derives the measured gene axis from the installed frozen mouse reference and adds synthetic background genes. It includes a cell with only two detected S03 seeds so Strict Original must return NA for that cell/program.

The example tests software execution, exact cell-ID joins, measurement qualification, and exports. It contains no real cells, independent biological replicates, biological truth, or validation of a functional phenotype. No author dataset or atlas download is needed.

Run the exact commands in the root README from a writable directory outside the source tree. The generated H5AD and output tables are intentionally excluded from Git.
