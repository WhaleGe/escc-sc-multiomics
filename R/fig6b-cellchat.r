
# ── Parse command-line arguments ──────────────────────────────────────────────────────
args   <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]   # Output directory
inrds  <- args[2]   # Path to the input Seurat RDS file
group  <- args[3]   # Cell-grouping column in meta.data (e.g., celltype / subcluster)
pre    <- args[4]   # Filename prefix for saved outputs
ts1    <- args[5]   # First tissue label  (reference condition)
ts2    <- args[6]   # Second tissue label (comparison condition)

# ── Load required libraries ───────────────────────────────────────────────────────────
library(patchwork)       # Combine ggplot2 panels
library(Seurat)          # Single-cell data processing
library(ComplexHeatmap)  # Advanced heatmap visualization
library(CellChat)        # Cell-cell communication inference
setwd(outdir)

# ── Global options ────────────────────────────────────────────────────────────────────
options(
    stringsAsFactors       = FALSE,
    future.globals.maxSize = 5e10   # Increase memory limit for parallel processing
)

# Load the full Seurat object
objs <- readRDS(inrds)

# ── Step 1: Build per-tissue CellChat objects ─────────────────────────────────────────
# Iterate over both tissue labels; results are stored in a named list
object.list <- list()

for (ts in c(ts1, ts2)) {

    message("Processing tissue: ", ts)

    # Subset Seurat object to the current tissue
    obj <- subset(objs, Tissue == ts)

    # Library-size normalization
    obj <- NormalizeData(obj)

    # Ensure group labels are stored as plain character strings
    obj@meta.data[[group]] <- as.character(obj@meta.data[[group]])

    # ── 1a: Create CellChat object and assign human ligand-receptor database ──────────
    cellchat       <- createCellChat(object = obj, group.by = group)
    cellchat@DB    <- CellChatDB.human   # Use the built-in human CellChat database

    # Subset to expressed signaling molecules only
    cellchat <- subsetData(cellchat)

    # ── 1b: Identify over-expressed genes and interactions ────────────────────────────
    # Enable parallel processing with 10 workers
    future::plan("multiprocess", workers = 10)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # ── 1c: Compute communication probabilities and filter low-cell groups ────────────
    # population.size = FALSE: do not scale by cell-type abundance
    cellchat <- computeCommunProb(cellchat, population.size = FALSE)
    cellchat <- filterCommunication(cellchat, min.cells = 10)

    # ── 1d: Infer pathway-level communication and aggregate network ───────────────────
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    # ── 1e: Compute signaling centrality scores ───────────────────────────────────────
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

    # Store the processed CellChat object
    object.list[[ts]] <- cellchat
}

# ── Step 2: Harmonize cell-type labels and merge CellChat objects ─────────────────────
# Derive a unified set of cell-type labels covering both tissues
group.new <- unique(c(
    levels(object.list[[1]]@idents),
    levels(object.list[[2]]@idents)
))

# Lift each object so that all cell types are represented (pad with zeros if absent)
object.list[[1]] <- liftCellChat(object.list[[1]], group.new)
object.list[[2]] <- liftCellChat(object.list[[2]], group.new)

# Merge the two CellChat objects for comparative analysis
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
message("Merged CellChat object created successfully.")

# ── Step 3: Visualize increased signaling in ts2 ─────────────────────────────────────
# Bubble plot: shows LR pairs that are more active in ts2 (max.dataset = 2)
gg1 <- netVisual_bubble(
    cellchat,
    sources.use    = cl1,            # Source cell-type indices or names
    targets.use    = cl2,            # Target cell-type indices or names
    comparison     = c(1, 2),        # Compare dataset 1 (ts1) vs dataset 2 (ts2)
    max.dataset    = 2,              # Highlight interactions enriched in dataset 2
    title.name     = paste0("Increased signaling in ", ts2),
    angle.x        = 45,             # Rotate x-axis labels by 45 degrees
    remove.isolate = TRUE            # Remove LR pairs with no significant interactions
)

ggsave(
    filename = "LRcom_bubble_up.pdf",
    plot     = gg1,
    width    = 6,
    height   = 6
)

message("Bubble plot saved to: LRcom_bubble_up.pdf")