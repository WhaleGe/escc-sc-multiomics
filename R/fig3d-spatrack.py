# ── Import libraries ──────────────────────────────────────────────────────────────────
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import spaTrack as spt

# ── Global parameters ─────────────────────────────────────────────────────────────────
outdir     = "spatrack"       # Output directory for all results
inrds      = "neu.h5ad"       # Path to input AnnData (.h5ad) file
pre        = "neu"            # Filename prefix for saved outputs
dimpre     = "umap"           # Dimensionality reduction embedding to use
groupname  = "cluster"        # obs column name that stores cell-type / cluster labels
startgroup = "neu_start"      # Cell-type label used to define trajectory starting cells

# Color palette for cluster visualization in the streamplot panel
pal3 = [
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF", "#999999"
]

# ── Setup ─────────────────────────────────────────────────────────────────────────────
os.chdir(outdir)

# ── Step 1: Load data ─────────────────────────────────────────────────────────────────
adata = sc.read(inrds)

# Copy groupname column to a standardized "cluster" key for downstream use
adata.obs["cluster"] = adata.obs[groupname]

# Ensure gene names are unique (required by Scanpy)
adata.var_names_make_unique()

# ── Step 2: Preprocessing ─────────────────────────────────────────────────────────────
# Retain only genes expressed in at least 30 cells
sc.pp.filter_genes(adata, min_cells=30)

# Library-size normalization: scale each cell to a total count of 10,000
sc.pp.normalize_total(adata, target_sum=1e4)

# Log-transform normalized counts: log(x + 1)
sc.pp.log1p(adata)

# ── Step 3: Define starting cells for trajectory ──────────────────────────────────────
# Identify starting cells by matching the specified cell-type label
start_cells = spt.set_start_cells(
    adata,
    select_way="cell_type",
    cell_type=startgroup
)

# Visualize all cells (grey) and starting cells (orange) on the UMAP
fig, ax = plt.subplots(figsize=(6, 5))
ax.scatter(
    adata.obsm["X_umap"][:, 0],
    adata.obsm["X_umap"][:, 1],
    c="lightgrey", s=20, label="All cells"
)
ax.scatter(
    adata.obsm["X_umap"][start_cells, 0],
    adata.obsm["X_umap"][start_cells, 1],
    c="orange", s=25, label="Start cells"
)
ax.legend(frameon=False)
ax.set_title("Start cells on UMAP")
plt.tight_layout()
plt.savefig("startctumap.pdf")
plt.close()

# ── Step 4: Compute optimal-transport matrix and pseudotime ──────────────────────────
# Build cell-cell transport matrix using optimal transport (single-cell mode)
adata.obsp["trans"] = spt.get_ot_matrix(adata, data_type="single-cell")

# Assign pseudotime to each cell based on its distance from the starting cells
adata.obs["ptime"] = spt.get_ptime(adata, start_cells)

# ── Step 5: Compute velocity grid on UMAP embedding ──────────────────────────────────
# n_neigh_pos: number of positional neighbours; n_neigh_gene: gene-space neighbours
adata.uns["E_grid"], adata.uns["V_grid"] = spt.get_velocity(
    adata,
    basis="umap",
    n_neigh_pos=100,
    n_neigh_gene=0
)

# ── Step 6: Plot clusters and starting-cell overlay ──────────────────────────────────
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 5))

# Left panel: UMAP colored by cluster label
sc.pl.embedding(
    adata,
    basis="X_umap",
    color="cluster",
    size=20,
    legend_loc="on data",
    legend_fontoutline=3,
    ax=axs[0],
    show=False
)

# Right panel: UMAP base plot with starting-cell group highlighted in orange
sc.pl.embedding(
    adata,
    basis="X_umap",
    ax=axs[1],
    show=False,
    title="Start cells"
)
start_coords = adata.obsm["X_umap"][adata.obs["cluster"] == startgroup]
axs[1].scatter(
    start_coords[:, 0],
    start_coords[:, 1],
    s=15, color="orange", label=startgroup
)
axs[1].legend(frameon=False)

plt.tight_layout()
plt.savefig("ctumap.pdf")
plt.close()

# ── Step 7: Plot pseudotime and velocity streamlines ─────────────────────────────────
fig, axs = plt.subplots(ncols=2, nrows=1, figsize=(10, 4))

# Left panel: UMAP colored by pseudotime (red gradient)
sc.pl.embedding(
    adata,
    basis="X_umap",
    color="ptime",
    color_map="Reds",
    title="Pseudotime",
    ax=axs[0],
    show=False
)

# Right panel: cluster-colored UMAP with velocity streamlines overlaid
vf_ax = sc.pl.embedding(
    adata,
    basis="X_umap",
    color="cluster",
    palette=pal3,
    title="Velocity streamplot",
    legend_loc=None,
    alpha=0.4,
    size=30,
    ax=axs[1],
    show=False
)
vf_ax.streamplot(
    adata.uns["E_grid"][0],
    adata.uns["E_grid"][1],
    adata.uns["V_grid"][0],
    adata.uns["V_grid"][1],
    color="black",
    linewidth=1.5,
    density=1.8,
    arrowsize=1.2
)

plt.tight_layout()
plt.savefig("pseutimeumap.pdf")
plt.close()

