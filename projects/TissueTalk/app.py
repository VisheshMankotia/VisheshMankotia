"""TissueTalk — interactive spatial gene-expression explorer.

The app performs descriptive regional mean comparisons only. It does not
perform formal differential-expression testing.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

st.set_page_config(page_title="TissueTalk | Spatial Explorer", page_icon="🧬", layout="wide")

st.markdown(
    """
    <style>
    .block-container {padding-top: 2rem; padding-bottom: 3rem;}
    [data-testid="stMetric"] {background: #f4f7fb; padding: 14px 16px; border-radius: 12px;}
    </style>
    """,
    unsafe_allow_html=True,
)

st.title("🧬 TissueTalk")
st.markdown("### Spatial gene-expression exploration")
st.write("Select a region of tissue, compare its expression profile with the rest of the sample, and inspect spatial patterns gene by gene.")

with st.sidebar:
    st.header("1 · Choose a dataset")
    demo_mode = st.toggle("Use simulated demo data", value=True)
    uploaded = None
    if not demo_mode:
        uploaded = st.file_uploader("Upload AnnData (.h5ad)", type=["h5ad"])
    st.caption("Demo data are synthetic and intended only to test the interface.")

@st.cache_data(show_spinner=False)
def make_demo_data():
    from scipy import sparse
    from anndata import AnnData

    rng = np.random.default_rng(17)
    side = 45
    xx, yy = np.meshgrid(np.arange(side), np.arange(side))
    coords = np.column_stack([xx.ravel(), yy.ravel()]).astype(float)
    n = coords.shape[0]
    n_genes = 120
    base = rng.poisson(1.2, size=(n, n_genes)).astype(np.float32)
    # Synthetic spatial signal for a few clearly labelled demo features.
    centre = np.sqrt((coords[:, 0] - 22) ** 2 + (coords[:, 1] - 22) ** 2) < 10
    base[centre, 0] += rng.poisson(5, size=centre.sum())
    base[centre, 1] += rng.poisson(3, size=centre.sum())
    base[~centre, 2] += rng.poisson(4, size=(~centre).sum())
    genes = [f"DEMO_GENE_{i:03d}" for i in range(1, n_genes + 1)]
    adata = AnnData(X=sparse.csr_matrix(base))
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
        st.error(f"Unable to read this file as AnnData: {exc}")
        st.stop()
else:
    st.info("Choose the simulated demo dataset or upload an AnnData .h5ad file to begin.")
    st.stop()

if "spatial" not in adata.obsm:
    st.error("The dataset must contain spatial coordinates in adata.obsm['spatial'].")
    st.stop()

coords = np.asarray(adata.obsm["spatial"])
if coords.ndim != 2 or coords.shape[0] != adata.n_obs or coords.shape[1] < 2:
    st.error("Spatial coordinates must have shape (number of spots, at least 2).")
    st.stop()
coords = coords[:, :2].astype(float)
valid = np.isfinite(coords).all(axis=1)
if not valid.all():
    st.warning(f"Excluded {(~valid).sum()} spots with missing or invalid coordinates.")
    adata = adata[valid].copy()
    coords = coords[valid]
if adata.n_obs < 2:
    st.error("At least two valid spatial spots are required.")
    st.stop()

gene_names = pd.Index(adata.var_names.astype(str))
if gene_names.has_duplicates:
    st.warning("Duplicate gene identifiers detected. Display labels may not be unique.")
    gene_names = pd.Index([f"{g}_{i}" for i, g in enumerate(gene_names)])

x_min, x_max = float(coords[:, 0].min()), float(coords[:, 0].max())
y_min, y_max = float(coords[:, 1].min()), float(coords[:, 1].max())
x_mid, y_mid = (x_min + x_max) / 2, (y_min + y_max) / 2
max_radius = max(float(np.hypot(x_max - x_min, y_max - y_min)), 1.0) / 2

with st.sidebar:
    st.divider()
    st.header("2 · Select a tissue region")
    x_center = st.slider("Centre X", x_min, x_max, x_mid)
    y_center = st.slider("Centre Y", y_min, y_max, y_mid)
    radius = st.slider("Radius", 0.0, max_radius, max_radius * 0.25)
    st.caption("The circular selection is defined using Euclidean distance in the supplied coordinate units.")

dist = np.hypot(coords[:, 0] - x_center, coords[:, 1] - y_center)
region_mask = dist <= radius
n_region = int(region_mask.sum())
n_rest = int((~region_mask).sum())

m1, m2, m3 = st.columns(3)
m1.metric("Tissue spots", f"{adata.n_obs:,}")
m2.metric("Selected spots", f"{n_region:,}")
m3.metric("Genes", f"{adata.n_vars:,}")

if n_region == 0 or n_rest == 0:
    st.warning("Adjust the centre or radius so that both the selected region and the remaining tissue contain spots.")
    st.stop()

map_df = pd.DataFrame({
    "x": coords[:, 0], "y": coords[:, 1],
    "Selection": np.where(region_mask, "Selected region", "Outside region")
})
fig = px.scatter(
    map_df, x="x", y="y", color="Selection",
    color_discrete_map={"Selected region": "#2563eb", "Outside region": "#cbd5e1"},
    title="Spatial spot map",
)
fig.update_traces(marker={"size": 6, "opacity": 0.85})
fig.update_layout(height=480, margin=dict(l=10, r=10, t=55, b=10), legend_title_text="")
fig.update_yaxes(scaleanchor="x", scaleratio=1, autorange="reversed")
st.plotly_chart(fig, use_container_width=True)

st.divider()
st.header("3 · Explore regional expression")
if st.button("Calculate region vs. rest", type="primary", use_container_width=False):
    with st.spinner("Calculating gene-wise average expression..."):
        matrix = adata.X
        # Sparse matrices are averaged without converting the full matrix to dense.
        region_mean = np.asarray(matrix[region_mask].mean(axis=0)).ravel()
        rest_mean = np.asarray(matrix[~region_mask].mean(axis=0)).ravel()
        result = pd.DataFrame({
            "gene": gene_names,
            "region_mean": region_mean,
            "rest_mean": rest_mean,
            "difference": region_mean - rest_mean,
        }).replace([np.inf, -np.inf], np.nan).dropna()
        result = result.sort_values("difference", ascending=False).reset_index(drop=True)
        st.session_state["tissuetalk_result"] = result

result = st.session_state.get("tissuetalk_result")
if result is not None:
    pos, neg = st.columns(2)
    with pos:
        st.subheader("Higher in selected region")
        st.dataframe(result.head(15), use_container_width=True, hide_index=True)
    with neg:
        st.subheader("Lower in selected region")
        st.dataframe(result.tail(15).sort_values("difference"), use_container_width=True, hide_index=True)

    st.download_button(
        "Download complete comparison (CSV)",
        result.to_csv(index=False).encode("utf-8"),
        file_name="tissuetalk_region_comparison.csv",
        mime="text/csv",
    )

    st.subheader("Inspect one gene")
    selected_gene = st.selectbox("Gene", options=result["gene"].tolist())
    gene_idx = list(gene_names).index(selected_gene)
    expression = adata.X[:, gene_idx]
    if hasattr(expression, "toarray"):
        expression = expression.toarray()
    expression = np.asarray(expression).ravel()
    gene_df = pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1], "expression": expression})
    gene_fig = px.scatter(gene_df, x="x", y="y", color="expression", color_continuous_scale="Viridis", title=f"Spatial expression: {selected_gene}")
    gene_fig.update_traces(marker={"size": 7})
    gene_fig.update_layout(height=430, margin=dict(l=10, r=10, t=55, b=10))
    gene_fig.update_yaxes(scaleanchor="x", scaleratio=1, autorange="reversed")
    st.plotly_chart(gene_fig, use_container_width=True)

    st.caption("Expression values are shown as stored in the input matrix. Their scale depends on the dataset's preprocessing.")
else:
    st.info("Click **Calculate region vs. rest** to generate the descriptive gene comparison.")

with st.expander("How to interpret these results"):
    st.markdown(
        """
        - **Difference** = mean expression in selected spots − mean expression in remaining spots.
        - A positive value means a higher arithmetic mean in the selected region; a negative value means lower.
        - These are descriptive summaries, not formal differential-expression results.
        - Results depend on normalization, spot composition, tissue sampling and data quality.
        - Statistical testing, multiple-testing correction and biological validation are needed before making research claims.
        """
    )

st.caption("TissueTalk · Research prototype by Vishesh Mankotia")
