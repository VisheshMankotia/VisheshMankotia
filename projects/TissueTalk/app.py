"""TissueTalk — AI-assisted spatial transcriptomics explorer.

Features:
- AnnData .h5ad upload and synthetic demo mode
- Spatial spot visualization and circular region selection
- Region-vs-rest expression summaries and CSV export
- Optional MyGene.info annotations
- Groq-powered, data-grounded Q&A
- Browser text-to-speech

The expression comparison is descriptive, not formal differential-expression testing.
"""
from __future__ import annotations

import html
import json
import os
from typing import Any

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

st.set_page_config(page_title="TissueTalk | Spatial AI Explorer", page_icon="🧬", layout="wide")

st.markdown(
    """
    <style>
    .block-container {padding-top: 1.8rem; padding-bottom: 3rem;}
    [data-testid="stMetric"] {background: #f4f7fb; padding: 14px 16px; border-radius: 12px;}
    .tissuetalk-subtitle {color: #64748b; font-size: 1.05rem;}
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🧬 TissueTalk")
st.markdown('<p class="tissuetalk-subtitle">Ask questions about the genes expressed in a selected region of your tissue.</p>', unsafe_allow_html=True)

with st.sidebar:
    st.header("1 · Dataset")
    demo_mode = st.toggle("Use simulated demo data", value=True)
    uploaded = None
    if not demo_mode:
        uploaded = st.file_uploader("Upload AnnData (.h5ad)", type=["h5ad"])
    st.caption("Demo values are synthetic. Use your own AnnData file for biological exploration.")

@st.cache_data(show_spinner=False)
def make_demo_data():
    from scipy import sparse
    from anndata import AnnData

    rng = np.random.default_rng(17)
    side = 45
    xx, yy = np.meshgrid(np.arange(side), np.arange(side))
    coords = np.column_stack([xx.ravel(), yy.ravel()]).astype(float)
    n = len(coords)
    n_genes = 120
    values = rng.poisson(1.2, size=(n, n_genes)).astype(np.float32)
    centre = np.hypot(coords[:, 0] - 22, coords[:, 1] - 22) < 10
    values[centre, 0] += rng.poisson(5, size=centre.sum())
    values[centre, 1] += rng.poisson(3, size=centre.sum())
    values[~centre, 2] += rng.poisson(4, size=(~centre).sum())
    genes = ["CXCL13", "MS4A1", "COL1A1"] + [f"DEMO_GENE_{i:03d}" for i in range(4, n_genes + 1)]
    adata = AnnData(X=sparse.csr_matrix(values))
    adata.var_names = genes
    adata.obsm["spatial"] = coords
    return adata

adata = None
if demo_mode:
    adata = make_demo_data()
elif uploaded is not None:
    try:
        import anndata as ad
        adata = ad.read_h5ad(uploaded)
    except Exception as exc:
        st.error(f"Could not read this file as AnnData: {exc}")
        st.stop()
else:
    st.info("Choose demo data or upload an AnnData .h5ad file.")
    st.stop()

if "spatial" not in adata.obsm:
    st.error("Missing adata.obsm['spatial']. The file must contain spatial coordinates.")
    st.stop()

coords = np.asarray(adata.obsm["spatial"])
if coords.ndim != 2 or coords.shape[0] != adata.n_obs or coords.shape[1] < 2:
    st.error("Spatial coordinates must have shape (number of spots, at least 2).")
    st.stop()
coords = coords[:, :2].astype(float)
valid = np.isfinite(coords).all(axis=1)
if not valid.all():
    st.warning(f"Excluded {(~valid).sum()} spots with invalid spatial coordinates.")
    adata = adata[valid].copy()
    coords = coords[valid]
if adata.n_obs < 2 or adata.n_vars < 1:
    st.error("The dataset needs at least two valid spots and one gene.")
    st.stop()

gene_names = pd.Index(adata.var_names.astype(str))
if gene_names.has_duplicates:
    st.warning("Duplicate gene identifiers found. Duplicate names have been made unique for display.")
    gene_names = pd.Index([f"{g}_{i}" for i, g in enumerate(gene_names)])

x_min, x_max = float(coords[:, 0].min()), float(coords[:, 0].max())
y_min, y_max = float(coords[:, 1].min()), float(coords[:, 1].max())
x_mid, y_mid = (x_min + x_max) / 2, (y_min + y_max) / 2
max_radius = max(float(np.hypot(x_max - x_min, y_max - y_min) / 2), 1.0)

with st.sidebar:
    st.divider()
    st.header("2 · Select tissue region")
    x_center = st.slider("Centre X", x_min, x_max, x_mid)
    y_center = st.slider("Centre Y", y_min, y_max, y_mid)
    radius = st.slider("Radius", 0.0, max_radius, max_radius * 0.25)

distance = np.hypot(coords[:, 0] - x_center, coords[:, 1] - y_center)
region_mask = distance <= radius
n_region = int(region_mask.sum())
n_rest = int((~region_mask).sum())

c1, c2, c3 = st.columns(3)
c1.metric("Tissue spots", f"{adata.n_obs:,}")
c2.metric("Selected spots", f"{n_region:,}")
c3.metric("Genes", f"{adata.n_vars:,}")

map_df = pd.DataFrame({
    "x": coords[:, 0], "y": coords[:, 1],
    "Region": np.where(region_mask, "Selected region", "Outside region"),
})
fig = px.scatter(
    map_df, x="x", y="y", color="Region",
    color_discrete_map={"Selected region": "#2563eb", "Outside region": "#cbd5e1"},
    title="Spatial tissue map",
)
fig.update_traces(marker={"size": 6, "opacity": 0.85})
fig.update_layout(height=460, margin=dict(l=8, r=8, t=55, b=8), legend_title_text="")
fig.update_yaxes(scaleanchor="x", scaleratio=1, autorange="reversed")
st.plotly_chart(fig, use_container_width=True)

if n_region == 0 or n_rest == 0:
    st.warning("Adjust the region centre or radius so both the selected region and the remaining tissue contain spots.")
    st.stop()

st.header("3 · Gene-expression analysis")
if st.button("Calculate region vs. rest", type="primary"):
    with st.spinner("Calculating gene-wise mean expression..."):
        matrix = adata.X
        region_mean = np.asarray(matrix[region_mask].mean(axis=0)).ravel()
        rest_mean = np.asarray(matrix[~region_mask].mean(axis=0)).ravel()
        result = pd.DataFrame({
            "gene": gene_names,
            "region_mean": region_mean,
            "rest_mean": rest_mean,
            "difference": region_mean - rest_mean,
        }).replace([np.inf, -np.inf], np.nan).dropna()
        result = result.sort_values("difference", ascending=False).reset_index(drop=True)
        st.session_state["tt_result"] = result
        st.session_state["tt_context"] = {
            "total_spots": int(adata.n_obs),
            "selected_spots": n_region,
            "remaining_spots": n_rest,
            "number_of_genes": int(adata.n_vars),
            "demo_mode": bool(demo_mode),
        }

result = st.session_state.get("tt_result")
if result is None:
    st.info("Select a region and calculate expression to activate gene annotations and the AI assistant.")
    st.stop()

top_n = st.slider("Number of genes to show", min_value=5, max_value=50, value=20)
left, right = st.columns(2)
with left:
    st.subheader("Higher in selected region")
    st.dataframe(result.head(top_n), use_container_width=True, hide_index=True)
with right:
    st.subheader("Lower in selected region")
    st.dataframe(result.tail(top_n).sort_values("difference"), use_container_width=True, hide_index=True)

st.download_button(
    "Download gene comparison (CSV)",
    result.to_csv(index=False).encode("utf-8"),
    file_name="tissuetalk_region_comparison.csv",
    mime="text/csv",
)

st.subheader("Explore an individual gene")
selected_gene = st.selectbox("Choose a gene", result["gene"].tolist())
gene_idx = list(gene_names).index(selected_gene)
gene_values = adata.X[:, gene_idx]
if hasattr(gene_values, "toarray"):
    gene_values = gene_values.toarray()
gene_values = np.asarray(gene_values).ravel()
gene_map = pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1], "expression": gene_values})
gene_fig = px.scatter(gene_map, x="x", y="y", color="expression", color_continuous_scale="Viridis", title=f"Spatial expression: {selected_gene}")
gene_fig.update_traces(marker={"size": 7})
gene_fig.update_layout(height=420, margin=dict(l=8, r=8, t=55, b=8))
gene_fig.update_yaxes(scaleanchor="x", scaleratio=1, autorange="reversed")
st.plotly_chart(gene_fig, use_container_width=True)

st.divider()
st.header("4 · Gene-function annotations")
st.write("Retrieve gene symbols, names and summaries from MyGene.info for the most enriched genes.")

if st.button("Annotate top genes with MyGene.info"):
    symbols = result.head(30)["gene"].astype(str).tolist()
    try:
        import mygene
        mg = mygene.MyGeneInfo()
        with st.spinner("Querying MyGene.info..."):
            annotations = mg.querymany(
                symbols,
                scopes="symbol",
                fields="symbol,name,summary",
                species="human",
                as_dataframe=False,
                verbose=False,
            )
        annotation_rows = []
        for item in annotations:
            if not isinstance(item, dict) or item.get("notfound"):
                continue
            annotation_rows.append({
                "gene": item.get("query", ""),
                "official_symbol": item.get("symbol", ""),
                "gene_name": item.get("name", ""),
                "summary": item.get("summary", ""),
            })
        annotation_df = pd.DataFrame(annotation_rows).drop_duplicates(subset=["gene"])
        st.session_state["tt_annotations"] = annotation_df
    except Exception as exc:
        st.error(f"Gene annotation lookup failed: {exc}")

annotations = st.session_state.get("tt_annotations")
if annotations is not None:
    if annotations.empty:
        st.warning("No annotations were returned for the selected gene symbols.")
    else:
        st.dataframe(annotations, use_container_width=True, hide_index=True)
        st.download_button(
            "Download gene annotations (CSV)",
            annotations.to_csv(index=False).encode("utf-8"),
            file_name="tissuetalk_gene_annotations.csv",
            mime="text/csv",
        )

st.divider()
st.header("5 · Ask TissueTalk")
st.caption("The assistant receives the selected-region summary and retrieved gene annotations. It should not invent findings beyond that evidence.")

question = st.text_area(
    "Ask a question about this region",
    key="tissuetalk_question",
    placeholder="Which genes have the largest positive expression differences? What is known about their functions?",
    height=90,
)

def build_evidence_context() -> str:
    top = result.head(20).to_dict(orient="records")
    bottom = result.tail(10).sort_values("difference").to_dict(orient="records")
    annotation_text: list[dict[str, Any]] = []
    ann = st.session_state.get("tt_annotations")
    if ann is not None and not ann.empty:
        annotation_text = ann.to_dict(orient="records")
    payload = {
        "dataset": st.session_state.get("tt_context", {}),
        "expression_scale": "Values are exactly as stored in adata.X; preprocessing is not inferred.",
        "top_positive_differences": top,
        "most_negative_differences": bottom,
        "gene_annotations": annotation_text,
    }
    return json.dumps(payload, ensure_ascii=False, default=str)

def get_groq_key() -> str | None:
    try:
        return st.secrets.get("GROQ_API_KEY") or os.environ.get("GROQ_API_KEY")
    except Exception:
        return os.environ.get("GROQ_API_KEY")

if st.button("Ask TissueTalk", type="primary", disabled=not bool(question.strip())):
    api_key = get_groq_key()
    if not api_key:
        st.error("Groq API key not found. Add GROQ_API_KEY to Streamlit Secrets or your environment variables. Never paste it into this source file.")
    else:
        try:
            from groq import Groq
            client = Groq(api_key=api_key)
            system_prompt = (
                "You are TissueTalk, a careful biomedical data assistant. "
                "Answer using only the supplied dataset summary and gene annotations. "
                "Clearly distinguish observed expression differences from known gene functions. "
                "Do not claim statistical significance, causality, diagnosis, treatment relevance, "
                "or biological validation because the supplied comparisons are descriptive means only. "
                "If evidence is missing, say so. Explain technical terms simply and mention that "
                "expression values retain the dataset's original scale."
            )
            with st.spinner("TissueTalk is preparing an evidence-grounded answer..."):
                completion = client.chat.completions.create(
                    model="openai/gpt-oss-20b",
                    messages=[
                        {"role": "system", "content": system_prompt},
                        {"role": "user", "content": "Question: " + question + "\n\nEvidence JSON:\n" + build_evidence_context()},
                    ],
                    temperature=0.2,
                    max_tokens=900,
                )
                answer = completion.choices[0].message.content or "The model returned an empty response."
            st.session_state["tt_answer"] = answer
        except Exception as exc:
            st.error(f"Could not generate an AI response: {exc}")

answer = st.session_state.get("tt_answer")
if answer:
    st.markdown("#### TissueTalk answer")
    st.markdown(answer)
    safe_answer = json.dumps(answer, ensure_ascii=False)
    tts_html = f"""
    <div style="display:flex;gap:10px;align-items:center;margin:8px 0 16px 0;">
      <button id="speak" style="background:#2563eb;color:white;border:0;border-radius:8px;padding:9px 15px;cursor:pointer;">🔊 Speak answer</button>
      <button id="stop" style="background:#e2e8f0;color:#172033;border:0;border-radius:8px;padding:9px 15px;cursor:pointer;">⏹ Stop</button>
      <span id="status" style="font:13px sans-serif;color:#64748b;"></span>
    </div>
    <script>
      const textToSpeak = {safe_answer};
      const status = document.getElementById("status");
      document.getElementById("speak").onclick = () => {{
        if (!("speechSynthesis" in window)) {{
          status.textContent = "Speech synthesis is not supported in this browser.";
          return;
        }}
        window.speechSynthesis.cancel();
        const utterance = new SpeechSynthesisUtterance(textToSpeak);
        utterance.onstart = () => status.textContent = "Speaking...";
        utterance.onend = () => status.textContent = "Finished.";
        utterance.onerror = () => status.textContent = "Speech could not be played.";
        window.speechSynthesis.speak(utterance);
      }};
      document.getElementById("stop").onclick = () => {{
        if ("speechSynthesis" in window) window.speechSynthesis.cancel();
        status.textContent = "Stopped.";
      }};
    </script>
    """
    st.components.v1.html(tts_html, height=65)

with st.expander("How to interpret TissueTalk results"):
    st.markdown(
        """
        - Difference = selected-region mean − remaining-tissue mean.
        - A positive difference indicates a higher arithmetic mean in the selected region.
        - This is exploratory descriptive analysis, not formal differential-expression testing.
        - The app does not automatically normalize data, adjust for covariates, calculate p-values or correct for multiple testing.
        - MyGene.info annotations describe known gene information; they do not prove that a gene causes the observed spatial pattern.
        """
    )

st.caption("TissueTalk · AI-assisted spatial transcriptomics prototype · Vishesh Mankotia")
