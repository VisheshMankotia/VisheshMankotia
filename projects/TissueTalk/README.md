# TissueTalk 🧬
### Explore spatial gene-expression maps through an interactive interface

TissueTalk is a research prototype for exploring spatial transcriptomics data. It lets a user select a region of tissue, compare gene expression in that region with the remaining tissue, and inspect genes that may be enriched locally.

> **Status:** Prototype / work in progress. The repository contains a starter interface and analysis outline. Dataset-specific validation, biological interpretation and model integration are ongoing.

## Goals
- Load a spatial transcriptomics dataset with spot coordinates and gene-expression values.
- Visualize the tissue spot map.
- Select a region using an interactive control.
- Calculate mean expression in the selected region and compare it with the rest of the tissue.
- Rank genes by regional expression difference.
- Add gene annotations and natural-language explanations as a future extension.

## Current analysis concept
For each gene:

`difference = mean(selected\_region) - mean(rest\_of\_tissue)`

A positive difference indicates higher average expression in the selected region. This is a descriptive comparison, not a statistical test or proof of biological significance.

## Project structure
```
TissueTalk/
├── app.py
├── requirements.txt
├── .gitignore
└── README.md
```

## Setup
1. Create and activate a Python virtual environment.
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Prepare a compatible AnnData `.h5ad` file with spatial coordinates in `adata.obsm["spatial"]`.
4. Run:
   ```bash
   streamlit run app.py
   ```
5. Upload the dataset in the app.

## Expected input
- AnnData object containing gene expression in `adata.X`.
- Gene identifiers in `adata.var_names`.
- Two-dimensional spatial coordinates in `adata.obsm["spatial"]`.

Check the dataset's documentation before assuming its expression values are normalized or log-transformed.

## Important interpretation notes
- Regional mean differences are exploratory and can be affected by spot density, library size, normalization, batch effects and tissue composition.
- Use appropriate normalization and statistical testing before drawing research conclusions.
- Gene-function annotations should be checked against authoritative biological databases.
- Do not treat generated natural-language summaries as validated biological findings.

## Roadmap
- [x] Define region-versus-rest expression comparison.
- [x] Create an interactive Streamlit starter app.
- [ ] Add dataset-specific preprocessing and quality-control checks.
- [ ] Add differential-expression statistics with appropriate assumptions.
- [ ] Integrate curated gene annotations.
- [ ] Add reproducible example data and screenshots.
- [ ] Evaluate explanations with domain expertise.

## Author
**Vishesh Mankotia** · [GitHub](https://github.com/VisheshMankotia) · [LinkedIn](https://www.linkedin.com/in/vishesh-mankotia-b8569b292/)
