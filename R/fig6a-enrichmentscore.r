
# ── Load required libraries ───────────────────────────────────────────────────────────
library(Seurat)     # Single-cell data processing and module scoring
library(ggplot2)    # Plotting backend
library(ggpubr)     # Publication-ready plot utilities
library(ggsci)      # Scientific journal color palettes
library(tidyverse)  # Data wrangling (dplyr, tibble, etc.)
library(pheatmap)   # Clustered heatmap visualization

# ── Step 1: Load input data ───────────────────────────────────────────────────────────
# Load the Seurat object (neutrophil subset)
obj_neu <- readRDS(infile)

# Load functional gene set table (columns: gene, func)
df_neu <- read.table(
    "Neu_funcgeneset.txt",
    sep    = "\t",
    header = TRUE
)

# ── Step 2: Filter gene sets to genes present in the Seurat object ───────────────────
# Retain only genes that appear in the dataset's feature space
df_neu <- subset(df_neu, gene %in% rownames(obj_neu))

# Split gene vector by functional group to create a named list of gene sets
neu_geneset <- split(df_neu$gene, df_neu$func)

# ── Step 3: Add module scores for each functional gene set ────────────────────────────
# AddModuleScore appends one score column per gene set to the metadata
obj_neu <- AddModuleScore(
    object   = obj_neu,
    features = neu_geneset
)

# Rename the newly appended score columns to match gene-set names
n_sets     <- length(neu_geneset)
score_cols <- (ncol(obj_neu@meta.data) - n_sets + 1):ncol(obj_neu@meta.data)
colnames(obj_neu@meta.data)[score_cols] <- names(neu_geneset)

# ── Step 4: Compute average module score per subcluster ──────────────────────────────
# Extract barcode, subcluster label, and all module score columns
meta_sub <- obj_neu@meta.data[, c("barcode", "subcluster", names(neu_geneset))]

# Average each module score within each subcluster
avg_scores           <- aggregate(
    meta_sub[, names(neu_geneset)],
    by   = list(subcluster = meta_sub$subcluster),
    FUN  = mean
)
rownames(avg_scores) <- avg_scores$subcluster
avg_scores           <- avg_scores[, -1]          # Drop the grouping column

# Transpose so rows = gene sets, columns = subclusters (for heatmap orientation)
avg_mat <- t(avg_scores)

# ── Step 5: Generate and save the heatmap ────────────────────────────────────────────
# Row-scale scores to highlight relative enrichment across subclusters
p <- pheatmap(
    avg_mat,
    cluster_rows = FALSE,   # Preserve original gene-set row order
    cluster_cols = FALSE,   # Preserve original subcluster column order
    scale        = "row",   # Z-score scaling within each gene set (row)
    border_color = "white"  # White grid lines between cells
)

ggsave(
    filename = "neu_func_mean_heat.pdf",
    plot     = p,
    width    = 5,
    height   = 6
)

