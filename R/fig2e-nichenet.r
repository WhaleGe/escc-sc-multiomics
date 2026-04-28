
# ── Load required libraries ───────────────────────────────────────────────────────────
library(readxl)       # Read Excel files
library(tidyverse)    # Data wrangling and visualization (includes dplyr, ggplot2, etc.)
library(Seurat)       # Single-cell / spatial data framework
library(ggplot2)      # Plotting
library(RColorBrewer) # Color palettes
library(nichenetr)    # NicheNet ligand-receptor analysis

# ── Step 0a: Load target gene set from pathway file ───────────────────────────────────
# Read pathway gene sets from Excel; extract unique individual gene symbols
df   <- read_xlsx("pathways.xlsx", sheet = 1)
gene <- unique(unlist(strsplit(df$geneID, split = "/")))

# ── Step 0b: Load and preprocess Seurat object ────────────────────────────────────────
# Subset to primary tumor (P_T) and metastatic tumor (M_T) spots/cells
merge      <- readRDS(inrds)
seurat_obj <- subset(merge, Tissue %in% c("P_T", "M_T"))

# Create a combined cell-type label: "<subclass>-<tissue>"
seurat_obj@meta.data$celltype_aggregate <- paste(
    seurat_obj@meta.data$subclass,
    seurat_obj@meta.data$Tissue,
    sep = "-"
)

# Set identity class to the combined cell-type label
celltype_id <- "celltype_aggregate"
seurat_obj  <- SetIdent(seurat_obj, value = seurat_obj[[celltype_id]])

# Normalize and scale expression data
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# ── Step 0c: Load NicheNet ligand-receptor network and target matrix ──────────────────
# Rows = target genes, Columns = ligands
ligand_target_matrix <- readRDS("00.database/ligand_target_matrix.rds")
lr_network           <- readRDS("00.database/lr_network.rds")

# Flag bona-fide interactions (exclude PPI-predicted entries)
lr_network <- lr_network %>%
    mutate(bonafide = !database %in% c("ppi_prediction", "ppi_prediction_go"))

# Standardize column names to "ligand" and "receptor"
lr_network <- lr_network %>%
    dplyr::rename(ligand = from, receptor = to) %>%
    distinct(ligand, receptor, bonafide)

organism <- "human"

# ── Step 1: Define the niches / microenvironments of interest ─────────────────────────
# Each niche specifies sender cell types and a receiver cell type
niches <- list(
    "MT_niche" = list(
        "sender" = c(
            "CD4T-M_T", "CD8T-M_T", "NK-M_T", "Bcells-M_T", "Plasma-M_T",
            "Macrophage-M_T", "Monocyte-M_T", "Mast-M_T", "DC-M_T",
            "Neutrophile-M_T", "Fibroblasts-M_T", "Pericytes-M_T", "Endothelial-M_T"
        ),
        "receiver" = c("m3-M_T")
    ),
    "PT_niche" = list(
        "sender" = c(
            "CD4T-P_T", "CD8T-P_T", "NK-P_T", "Bcells-P_T", "Plasma-P_T",
            "Macrophage-P_T", "Monocyte-P_T", "Mast-P_T", "DC-P_T",
            "Neutrophile-P_T", "Fibroblasts-P_T", "Pericytes-P_T", "Endothelial-P_T"
        ),
        "receiver" = c("m3-P_T")
    )
)

# ── Step 2: Calculate differential expression between niches ──────────────────────────
assay_oi <- "RNA"   # Assay to use (alternatives: "SCT", etc.)

# DE analysis for sender cell types (ligands only)
DE_sender <- calculate_niche_de(
    seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()),
    niches     = niches,
    type       = "sender",
    assay_oi   = assay_oi
)

# DE analysis for receiver cell types (receptors only)
DE_receiver <- calculate_niche_de(
    seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()),
    niches     = niches,
    type       = "receiver",
    assay_oi   = assay_oi
)

# Replace Inf / -Inf log2FC values with finite max / min to avoid downstream errors
DE_sender <- DE_sender %>%
    mutate(avg_log2FC = ifelse(
        avg_log2FC == Inf,  max(avg_log2FC[is.finite(avg_log2FC)]),
        ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)
    ))

DE_receiver <- DE_receiver %>%
    mutate(avg_log2FC = ifelse(
        avg_log2FC == Inf,  max(avg_log2FC[is.finite(avg_log2FC)]),
        ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)
    ))

# Process DE tables and combine sender-receiver pairs
expression_pct       <- 0.10
DE_sender_processed  <- process_niche_de(DE_sender,   niches, expression_pct, type = "sender")
DE_receiver_processed <- process_niche_de(DE_receiver, niches, expression_pct, type = "receiver")

specificity_score_LR_pairs <- "min_lfc"
DE_sender_receiver <- combine_sender_receiver_de(
    DE_sender_processed,
    DE_receiver_processed,
    lr_network,
    specificity_score = specificity_score_LR_pairs
)

# ── Step 3: Spatial DE (optional) ────────────────────────────────────────────────────
# Set both flags to FALSE if spatial region information is unavailable
include_spatial_info_sender   <- FALSE
include_spatial_info_receiver <- FALSE

# If no spatial info is available, create a neutral (mock) spatial tibble
if (!include_spatial_info_sender & !include_spatial_info_receiver) {
    spatial_info <- tibble(
        celltype_region_oi  = NA,
        celltype_other_region = NA
    ) %>% mutate(
        niche         = niches %>% names() %>% head(1),
        celltype_type = "sender"
    )
}

# Process sender spatial DE (or assign neutral spatial scores if unavailable)
if (include_spatial_info_sender) {
    sender_spatial_DE <- calculate_spatial_DE(
        seurat_obj   = seurat_obj %>% subset(features = lr_network$ligand %>% unique()),
        spatial_info = spatial_info %>% filter(celltype_type == "sender")
    )
    sender_spatial_DE_processed <- process_spatial_de(
        DE_table       = sender_spatial_DE,
        type           = "sender",
        lr_network     = lr_network,
        expression_pct = expression_pct,
        specificity_score = specificity_score_spatial
    )
    # Add neutral spatial scores for senders without spatial relevance
    sender_spatial_DE_others    <- get_non_spatial_de(niches, spatial_info, type = "sender", lr_network)
    sender_spatial_DE_processed <- sender_spatial_DE_processed %>%
        bind_rows(sender_spatial_DE_others) %>%
        mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
} else {
    # All senders receive a neutral spatial score
    sender_spatial_DE_processed <- get_non_spatial_de(niches, spatial_info, type = "sender", lr_network) %>%
        mutate(scaled_ligand_score_spatial = scale_quantile_adapted(ligand_score_spatial))
}

# Process receiver spatial DE (or assign neutral spatial scores if unavailable)
if (include_spatial_info_receiver) {
    receiver_spatial_DE <- calculate_spatial_DE(
        seurat_obj   = seurat_obj %>% subset(features = lr_network$receptor %>% unique()),
        spatial_info = spatial_info %>% filter(celltype_type == "receiver")
    )
    receiver_spatial_DE_processed <- process_spatial_de(
        DE_table       = receiver_spatial_DE,
        type           = "receiver",
        lr_network     = lr_network,
        expression_pct = expression_pct,
        specificity_score = specificity_score_spatial
    )
    # Add neutral spatial scores for receivers without spatial relevance
    receiver_spatial_DE_others    <- get_non_spatial_de(niches, spatial_info, type = "receiver", lr_network)
    receiver_spatial_DE_processed <- receiver_spatial_DE_processed %>%
        bind_rows(receiver_spatial_DE_others) %>%
        mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
} else {
    # All receivers receive a neutral spatial score
    receiver_spatial_DE_processed <- get_non_spatial_de(niches, spatial_info, type = "receiver", lr_network) %>%
        mutate(scaled_receptor_score_spatial = scale_quantile_adapted(receptor_score_spatial))
}

# ── Step 4: Calculate ligand activities and infer active ligand-target links ──────────
lfc_cutoff              <- 0.25
specificity_score_targets <- "min_lfc"

DE_receiver_targets <- calculate_niche_de_targets(
    seurat_obj     = seurat_obj,
    niches         = niches,
    lfc_cutoff     = lfc_cutoff,
    expression_pct = expression_pct,
    assay_oi       = assay_oi
)

DE_receiver_processed_targets <- process_receiver_target_de(
    DE_receiver_targets = DE_receiver_targets,
    niches              = niches,
    expression_pct      = expression_pct,
    specificity_score   = specificity_score_targets
)

# Define background gene set (all expressed target genes)
background <- DE_receiver_processed_targets %>% pull(target) %>% unique()

# Niche 1 (MT): use custom pathway gene set loaded at the top
geneset_niche1 <- gene

# Niche 2 (PT): use differentially expressed target genes above score threshold
geneset_niche2 <- DE_receiver_processed_targets %>%
    filter(
        receiver         == niches[[2]]$receiver,
        target_score     >= 0.75,
        target_significant == 1,
        target_present   == 1
    ) %>%
    pull(target) %>%
    unique()

top_n_target <- 250

niche_geneset_list <- list(
    "MT_niche" = list(
        "receiver"   = niches[[1]]$receiver,
        "geneset"    = geneset_niche1,
        "background" = background
    ),
    "PT_niche" = list(
        "receiver"   = niches[[2]]$receiver,
        "geneset"    = geneset_niche2,
        "background" = background
    )
)

ligand_activities_targets <- get_ligand_activities_targets(
    niche_geneset_list   = niche_geneset_list,
    ligand_target_matrix = ligand_target_matrix,
    top_n_target         = top_n_target
)

# ── Step 5: Compute scaled expression of ligands, receptors, and targets ──────────────
features_oi <- union(lr_network$ligand, lr_network$receptor) %>%
    union(ligand_activities_targets$target) %>%
    setdiff(NA)

# Build DotPlot to extract expression and fraction data
dotplot  <- suppressWarnings(
    Seurat::DotPlot(
        seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()),
        features = features_oi,
        assay    = assay_oi
    )
)

exprs_tbl <- dotplot$data %>%
    as_tibble() %>%
    rename(
        celltype           = id,
        gene               = features.plot,
        expression         = avg.exp,
        expression_scaled  = avg.exp.scaled,
        fraction           = pct.exp
    ) %>%
    mutate(fraction = fraction / 100) %>%
    select(celltype, gene, expression, expression_scaled, fraction) %>%
    distinct() %>%
    arrange(gene) %>%
    mutate(gene = as.character(gene))

# Split expression table by molecule type
exprs_tbl_ligand <- exprs_tbl %>%
    filter(gene %in% lr_network$ligand) %>%
    rename(
        sender                  = celltype,
        ligand                  = gene,
        ligand_expression       = expression,
        ligand_expression_scaled = expression_scaled,
        ligand_fraction         = fraction
    )

exprs_tbl_receptor <- exprs_tbl %>%
    filter(gene %in% lr_network$receptor) %>%
    rename(
        receiver                   = celltype,
        receptor                   = gene,
        receptor_expression        = expression,
        receptor_expression_scaled = expression_scaled,
        receptor_fraction          = fraction
    )

exprs_tbl_target <- exprs_tbl %>%
    filter(gene %in% ligand_activities_targets$target) %>%
    rename(
        receiver                  = celltype,
        target                    = gene,
        target_expression         = expression,
        target_expression_scaled  = expression_scaled,
        target_fraction           = fraction
    )

# Scale expression and fraction values within each table
exprs_tbl_ligand <- exprs_tbl_ligand %>%
    mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>%
    mutate(ligand_fraction_adapted = ligand_fraction) %>%
    mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct) %>%
    mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))

exprs_tbl_receptor <- exprs_tbl_receptor %>%
    mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled)) %>%
    mutate(receptor_fraction_adapted = receptor_fraction) %>%
    mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct) %>%
    mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))

# ── Step 6: Combine expression fractions with sender-receiver DE ──────────────────────
exprs_sender_receiver <- lr_network %>%
    inner_join(exprs_tbl_ligand,  by = "ligand") %>%
    inner_join(exprs_tbl_receptor, by = "receptor") %>%
    inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))

# Compute a composite ligand-scaled receptor expression-fraction score
ligand_scaled_receptor_expression_fraction_df <- exprs_sender_receiver %>%
    group_by(ligand, receiver) %>%
    mutate(
        rank_receptor_expression = dense_rank(receptor_expression),
        rank_receptor_fraction   = dense_rank(receptor_fraction)
    ) %>%
    mutate(
        ligand_scaled_receptor_expression_fraction = 0.5 * (
            (rank_receptor_fraction  / max(rank_receptor_fraction)) +
            (rank_receptor_expression / max(rank_receptor_expression))
        )
    ) %>%
    distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction, bonafide) %>%
    ungroup()

# ── Step 7: Prioritize ligand-receptor and ligand-target links ────────────────────────
# Each weight reflects the relative importance of each scoring component
prioritizing_weights <- c(
    "scaled_ligand_score"                          = 5,    # Ligand DE score (highest weight)
    "scaled_ligand_expression_scaled"              = 1,
    "ligand_fraction"                              = 1,
    "scaled_ligand_score_spatial"                  = 2,
    "scaled_receptor_score"                        = 0.5,
    "scaled_receptor_expression_scaled"            = 0.5,
    "receptor_fraction"                            = 1,
    "ligand_scaled_receptor_expression_fraction"   = 1,
    "scaled_receptor_score_spatial"                = 0,
    "scaled_activity"                              = 0,
    "scaled_activity_normalized"                   = 1,
    "bona_fide"                                    = 1
)

# Assemble all intermediate results into a single output list
output <- list(
    DE_sender_receiver                      = DE_sender_receiver,
    ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df,
    sender_spatial_DE_processed             = sender_spatial_DE_processed,
    receiver_spatial_DE_processed           = receiver_spatial_DE_processed,
    ligand_activities_targets               = ligand_activities_targets,
    DE_receiver_processed_targets           = DE_receiver_processed_targets,
    exprs_tbl_ligand                        = exprs_tbl_ligand,
    exprs_tbl_receptor                      = exprs_tbl_receptor,
    exprs_tbl_target                        = exprs_tbl_target
)

# Compute final prioritization tables (prioritization_score = overall priority score)
prioritization_tables <- get_prioritization_tables(output, prioritizing_weights)

# Quick inspection of top results for niche 1 receiver
prioritization_tables$prioritization_tbl_ligand_receptor %>%
    filter(receiver == niches[[1]]$receiver) %>% head(10)

prioritization_tables$prioritization_tbl_ligand_target %>%
    filter(receiver == niches[[1]]$receiver) %>% head(10)

# ── Step 8: Visualization ─────────────────────────────────────────────────────────────
receiver_oi <- "m3-M_T"

# -- 8a: Select top ligands and their best-matching receptor pairs --------------------
# Keep the single highest-scoring niche per ligand
top_ligand_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    select(niche, sender, receiver, ligand, receptor, prioritization_score) %>%
    group_by(ligand) %>%
    top_n(1, prioritization_score) %>%
    ungroup() %>%
    select(ligand, receptor, niche) %>%
    rename(top_niche = niche)   # Filter by ligand only

# Keep the single highest-scoring niche per ligand-receptor combination
top_ligand_receptor_niche_df <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    select(niche, sender, receiver, ligand, receptor, prioritization_score) %>%
    group_by(ligand, receptor) %>%
    top_n(1, prioritization_score) %>%
    ungroup() %>%
    select(ligand, receptor, niche) %>%
    rename(top_niche = niche)   # Filter by LR combination

# Select top 50 ligands per niche for the receiver of interest
ligand_prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    select(niche, sender, receiver, ligand, prioritization_score) %>%
    group_by(ligand, niche) %>%
    top_n(1, prioritization_score) %>%
    ungroup() %>%
    distinct() %>%
    inner_join(top_ligand_niche_df) %>%
    filter(niche == top_niche) %>%
    group_by(niche) %>%
    top_n(50, prioritization_score) %>%
    ungroup()

# Build the prioritized ligand-receptor table for the receiver of interest
filtered_ligands  <- ligand_prioritized_tbl_oi %>%
    filter(receiver == receiver_oi) %>%
    pull(ligand) %>%
    unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    filter(ligand %in% filtered_ligands) %>%
    select(niche, sender, receiver, ligand, receptor, ligand_receptor, prioritization_score) %>%
    distinct() %>%
    inner_join(top_ligand_receptor_niche_df) %>%
    group_by(ligand) %>%
    filter(receiver == receiver_oi) %>%
    top_n(2, prioritization_score) %>%
    ungroup()

# -- 8b: LFC heatmap for top LR pairs (minimum LFC compared to other niches) ----------
lfc_plot <- make_ligand_receptor_lfc_plot(
    receiver_oi, prioritized_tbl_oi,
    prioritization_tables$prioritization_tbl_ligand_receptor,
    plot_legend = FALSE, heights = NULL, widths = NULL
)
ggsave(filename = "m3-M_T_top100LR_lfcplot.pdf", plot = lfc_plot, width = 12, height = 20)
write.table(lfc_plot$data, "m3-M_T_top100LR_lfcplot.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -- 8c: Ligand expression, activity and target heatmap (top 50 ligands) --------------
exprs_activity_target_plot <- make_ligand_activity_target_exprs_plot(
    receiver_oi,
    prioritized_tbl_oi,
    prioritization_tables$prioritization_tbl_ligand_receptor,
    prioritization_tables$prioritization_tbl_ligand_target,
    output$exprs_tbl_ligand,
    output$exprs_tbl_target,
    lfc_cutoff,
    ligand_target_matrix,
    plot_legend = FALSE, heights = NULL, widths = NULL
)
ggsave("m3-M_T_top100LT.pdf",        plot = exprs_activity_target_plot$combined_plot, width = 20, height = 16)
ggsave("m3-M_T_top100LT_legend.pdf", plot = exprs_activity_target_plot$legends,       width = 10, height = 10)

# -- 8d: Ligand expression, activity and target heatmap (top 20 ligands) --------------
filtered_ligands <- ligand_prioritized_tbl_oi %>%
    filter(receiver == receiver_oi) %>%
    top_n(20, prioritization_score) %>%
    pull(ligand) %>%
    unique()

prioritized_tbl_oi <- prioritization_tables$prioritization_tbl_ligand_receptor %>%
    filter(ligand %in% filtered_ligands) %>%
    select(niche, sender, receiver, ligand, receptor, ligand_receptor, prioritization_score) %>%
    distinct() %>%
    inner_join(top_ligand_receptor_niche_df) %>%
    group_by(ligand) %>%
    filter(receiver == receiver_oi) %>%
    top_n(2, prioritization_score) %>%
    ungroup()

exprs_activity_target_plot <- make_ligand_activity_target_exprs_plot(
    receiver_oi,
    prioritized_tbl_oi,
    prioritization_tables$prioritization_tbl_ligand_receptor,
    prioritization_tables$prioritization_tbl_ligand_target,
    output$exprs_tbl_ligand,
    output$exprs_tbl_target,
    lfc_cutoff,
    ligand_target_matrix,
    plot_legend = FALSE, heights = NULL, widths = NULL
)
ggsave("m3-M_T_top20LT.pdf", plot = exprs_activity_target_plot$combined_plot, width = 20, height = 10)