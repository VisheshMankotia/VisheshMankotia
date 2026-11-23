# TissueTalk 🧬
### Interactive spatial gene-expression exploration

TissueTalk is a Streamlit research prototype for exploring spatial transcriptomics data. Select a region of tissue, compare its average gene expression with the rest of the sample, rank regional differences, and inspect the spatial distribution of an individual gene.

**Project status:** Working prototype. The app includes a synthetic demo mode and supports user-supplied AnnData `.h5ad` files. It is not a validated differential-expression pipeline and does not make diagnostic or clinical claims.

## Features
- Synthetic demo dataset to test the interface without downloading data.
- Upload support for AnnData `.h5ad` datasets.
- Spatial spot visualization and circular region selection.
- Region-versus-rest gene-wise mean comparison.
- Lists of genes with higher and lower mean expression in the selected region.
- Downloadable CSV results.
- Spatial expression map for an individual gene.
- Sparse-matrix-aware calculations.

## Quick start

```bash
git clone https://github.com/VisheshMankotia/VisheshMankotia.git
cd VisheshMankotia/projects/TissueTalk
python -m venv .venv
```

Activate the environment:

**Windows PowerShell**
```powershell
.venv\Scripts\Activate.ps1
```

**macOS / Linux**
```bash
source .venv/bin/activate
```

Install and run:

```bash
pip install -r requirements.txt
streamlit run app.py
```

The app starts in simulated demo mode. Switch off **Use simulated demo data** in the sidebar to upload your own dataset.

## Input format

The uploaded file must be an AnnData `.h5ad` object containing:
- `adata.X`: expression matrix (spots × genes).
- `adata.var_names`: gene identifiers.
- `adata.obsm["spatial"]`: two-dimensional spot coordinates.

See [DATA_GUIDE.md](DATA_GUIDE.md) for preprocessing and data-provenance guidance.

## Analysis

For each gene, TissueTalk calculates:

`difference = mean(selected region) - mean(remaining tissue)`

Sparse matrices are averaged without densifying the full expression matrix. Expression values are used as stored in the dataset; the app does not silently normalize or log-transform data.

## Project layout

```
TissueTalk/
├── app.py
├── requirements.txt
├── DATA_GUIDE.md
├── README.md
├── .gitignore
└── .streamlit/
    └── config.toml
```

## Limitations and responsible interpretation
- Regional mean differences are descriptive and are not equivalent to statistically tested differential expression.
- The app does not currently perform normalization, quality control, covariate adjustment, spatial statistical testing or multiple-testing correction.
- Results may be influenced by library size, spot density, batch effects, tissue composition and coordinate scaling.
- Gene names and annotations should be verified against authoritative sources.
- Synthetic demo data are not biological observations.
- Do not use the tool for diagnosis, treatment decisions or clinical interpretation.

## Roadmap
- [x] Interactive spatial spot map and region selection.
- [x] Synthetic demo dataset.
- [x] Region-versus-rest mean expression comparison.
- [x] CSV export and single-gene spatial visualization.
- [ ] Dataset-specific preprocessing and QC workflows.
- [ ] Curated gene annotations with source links.
- [ ] Statistically appropriate differential-expression methods.
- [ ] Reproducible screenshots, tests and example walkthrough.
- [ ] Optional natural-language assistant with transparent source attribution.

## Author
**Vishesh Mankotia** · [GitHub](https://github.com/VisheshMankotia) · [LinkedIn](https://www.linkedin.com/in/vishesh-mankotia-b8569b292/)
