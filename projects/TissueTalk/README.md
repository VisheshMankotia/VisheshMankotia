# TissueTalk 🧬
### AI-assisted spatial transcriptomics exploration

TissueTalk is a Streamlit research prototype that combines spatial gene-expression visualization, region-based expression summaries, gene-function annotations and a conversational AI assistant.

The original prototype was developed around a spatial transcriptomics dataset with **4,035 tissue spots and 36,601 genes**. This repository does not redistribute that dataset; upload your own compatible AnnData file or use the synthetic demo.

## Core capabilities
- Load AnnData `.h5ad` files with spatial coordinates.
- Explore tissue spots and select a circular region of interest.
- Compare gene-wise mean expression in the selected region against the remaining tissue.
- Rank genes by positive and negative mean-expression differences.
- Inspect a selected gene's spatial expression map.
- Retrieve gene symbols, names and summaries from MyGene.info.
- Ask a Groq-powered AI assistant questions grounded in the region-level results and retrieved annotations.
- Listen to answers using browser speech synthesis.
- Export expression comparisons and annotations as CSV files.
- Explore the interface with synthetic demo data.

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

The app opens in synthetic demo mode. Switch off **Use simulated demo data** to upload a real dataset.

## Enable the AI assistant

Set a Groq API key as an environment variable:

**macOS / Linux**
```bash
export GROQ_API_KEY="your_key_here"
```

**Windows PowerShell**
```powershell
$env:GROQ_API_KEY="your_key_here"
```

For Streamlit Community Cloud, add this under App settings → Secrets:

```toml
GROQ_API_KEY = "your_key_here"
```

The app uses Groq's `openai/gpt-oss-20b` model. Keep the key private; never commit it to GitHub. If no key is configured, spatial analysis and MyGene annotation remain available while AI Q&A displays a setup message.

## Dataset requirements

The uploaded AnnData file must contain:
- `adata.X`: expression matrix (spots × genes).
- `adata.var_names`: gene identifiers, preferably gene symbols for annotation.
- `adata.obsm["spatial"]`: two-dimensional spatial coordinates.

See [DATA_GUIDE.md](DATA_GUIDE.md) for preprocessing and provenance notes.

## How the AI assistant works

TissueTalk passes the assistant a compact evidence package containing:
- Dataset and region spot counts.
- Top positive and negative region-versus-rest mean differences.
- MyGene.info annotations retrieved for the most enriched genes.
- A prompt instructing the model to distinguish measured expression patterns from known gene functions and to avoid unsupported claims.

The model is not given the full expression matrix. Its explanations are limited to the supplied summary and annotations.

## Scientific interpretation and limitations
- Region-minus-rest values are descriptive arithmetic mean differences, not formal differential-expression statistics.
- The app does not automatically normalize expression, adjust for confounders, calculate p-values or perform multiple-testing correction.
- Expression values retain the scale stored in `adata.X`.
- MyGene.info annotations provide external gene knowledge and do not establish a causal explanation for the spatial pattern.
- AI-generated explanations can be wrong and must be checked against primary literature and biological expertise.
- Synthetic demo data are not biological observations.
- Do not use TissueTalk for diagnosis or clinical decisions.

## Project structure

```
TissueTalk/
├── app.py
├── requirements.txt
├── README.md
├── DATA_GUIDE.md
├── DEPLOYMENT.md
├── CHANGELOG.md
├── .gitignore
└── .streamlit/
    └── config.toml
```

## Roadmap
- [x] Spatial spot visualization and circular region selection.
- [x] Region-versus-rest expression comparison.
- [x] Gene ranking, CSV export and single-gene spatial visualization.
- [x] MyGene.info annotation integration.
- [x] Groq-based evidence-grounded Q&A.
- [x] Browser text-to-speech controls.
- [ ] Add reproducible formal differential-expression analysis.
- [ ] Add validation against known tissue markers and expert review.
- [ ] Add tests and a public deployment after end-to-end verification.

## Author
**Vishesh Mankotia** · [GitHub](https://github.com/VisheshMankotia) · [LinkedIn](https://www.linkedin.com/in/vishesh-mankotia-b8569b292/)
