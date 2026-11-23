# TissueTalk dataset guide

TissueTalk expects an AnnData `.h5ad` file.

## Required fields
- `adata.X`: gene-expression matrix (spots × genes).
- `adata.var_names`: gene identifiers.
- `adata.obsm["spatial"]`: numeric x/y coordinates, one row per spot.

## Before analysis
1. Read the dataset's original publication and data-use terms.
2. Check whether the expression matrix is raw counts, normalized counts or log-transformed values.
3. Inspect missing values, duplicated gene identifiers and coordinate validity.
4. Apply an appropriate quality-control and normalization workflow before comparisons.
5. Preserve the original object and document any filtering or transformation.

## Interpretation
The current starter app reports arithmetic mean expression and a region-minus-rest difference. It does not calculate p-values, adjust for multiple testing, correct for spatial autocorrelation or control for confounders. Treat the output as exploratory visualization, not formal differential-expression analysis.

Do not commit large or restricted datasets to this repository. Provide a source link and instructions for obtaining data instead.
