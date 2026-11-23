"""TissueTalk: exploratory spatial gene-expression viewer.

Input: AnnData .h5ad with adata.obsm["spatial"] containing x/y coordinates.
This starter app performs descriptive region-vs-rest comparisons only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

st.set_page_config(page_title="TissueTalk", page_icon="🧬", layout="wide")
st.title("🧬 TissueTalk")
st.caption("Explore spatial gene expression — select a tissue region and compare gene activity.")

uploaded = st.file_uploader("Upload an AnnData .h5ad file", type=["h5ad"])

if uploaded is None:
    st.info(
        "Upload a spatial transcriptomics AnnData file to begin. "
        "The file must contain spatial coordinates in adata.obsm['spatial']."
    )
    st.stop()

try:
    import anndata as ad
    adata = ad.read_h5ad(uploaded)
except Exception as exc:
    st.error(f"Could not read the AnnData file: {exc}")
    st.stop()

if "spatial" not in adata.obsm:
    st.error("Missing adata.obsm['spatial']. Please provide 2D spot coordinates.")
    st.stop()

coords = np.asarray(adata.obsm["spatial"])
if coords.ndim != 2 or coords.shape[0] != adata.n_obs or coords.shape[1] < 2:
    st.error("Spatial coordinates must have shape (number_of_spots, 2 or more).")
    st.stop()

coords = coords[:, :2]
x_min, x_max = float(np.nanmin(coords[:, 0])), float(np.nanmax(coords[:, 0]))
y_min, y_max = float(np.nanmin(coords[:, 1])), float(np.nanmax(coords[:, 1]))
x_center = st.sidebar.slider("Region centre — X", x_min, x_max, float((x_min + x_max) / 2))
y_center = st.sidebar.slider("Region centre — Y", y_min, y_max, float((y_min + y_max) / 2))
max_radius = max(float(np.hypot(x_max - x_min, y_max - y_min) / 2), 1e-6)
radius = st.sidebar.slider("Region radius", 0.0, max_radius, max_radius / 4)

distances = np.sqrt((coords[:, 0] - x_center) ** 2 + (coords[:, 1] - y_center) ** 2)
region_mask = distances <= radius
n_region = int(region_mask.sum())

left, right = st.columns(2)
left.metric("Total tissue spots", f"{adata.n_obs:,}")
right.metric("Selected region spots", f"{n_region:,}")

if n_region == 0 or n_region == adata.n_obs:
    st.warning("Adjust the centre or radius so the selected region contains some, but not all, spots.")
    st.stop()

st.subheader("Tissue region")
plot_df = pd.DataFrame({"x": coords[:, 0], "y": coords[:, 1],
                        "region": np.where(region_mask, "Selected region", "Outside region")})
st.scatter_chart(plot_df, x="x", y="y", color="region")

if st.button("Calculate regional expression", type="primary"):
    with st.spinner("Calculating gene-wise means..."):
        matrix = adata.X
        # Convert sparse matrices only for the selected calculation.
        region = matrix[region_mask]
        rest = matrix[~region_mask]
        if hasattr(region, "toarray"):
            region = region.toarray()
            rest = rest.toarray()

        region_mean = np.asarray(region.mean(axis=0)).ravel()
        rest_mean = np.asarray(rest.mean(axis=0)).ravel()
        result = pd.DataFrame({
            "gene": adata.var_names.astype(str),
            "region_mean": region_mean,
            "rest_mean": rest_mean,
            "difference": region_mean - rest_mean,
        }).sort_values("difference", ascending=False)

    st.subheader("Genes with higher mean expression in the selected region")
    st.dataframe(result.head(25), use_container_width=True, hide_index=True)
    st.download_button(
        "Download all gene comparisons (CSV)",
        result.to_csv(index=False).encode("utf-8"),
        file_name="tissuetalk_region_comparison.csv",
        mime="text/csv",
    )
    st.caption(
        "Exploratory descriptive results only. Confirm normalization, quality control, "
        "statistical assumptions and biological relevance before interpreting findings."
    )
