
# ── Load required libraries ───────────────────────────────────────────────────────────
library(reshape2)   # dcast(): wide-format reshaping
library(stringr)    # str_extract(): regex-based string extraction
library(dplyr)      # Data manipulation (group_by, summarise, across, etc.)

# ── Parse command-line arguments ──────────────────────────────────────────────────────
# (Reserved for future use; no positional arguments are consumed at present)
args <- commandArgs(trailingOnly = TRUE)

# ── Step 1: Load and reshape inferCNV CNV-region predictions ─────────────────────────
# Read the annotated CNV-region file; retain cell_group_name, cytoband, and state
cluster_mut <- read.table(
    "EpiImmu_6P.pred_cnv_regions.annoband.txt",
    header = TRUE
)
cluster_mut <- cluster_mut[, c(1, 3, 7)]   # Keep: cell_group_name | cytoband | state

# Create a composite key for deduplication before aggregation
cluster_mut$id <- paste(cluster_mut$cell_group_name, cluster_mut$cytoband, sep = "_")

# Average numeric columns (state) and take the max of character columns (cell_group /
# cytoband) within each cell-group × cytoband combination
cluster_mut <- cluster_mut %>%
    group_by(id) %>%
    summarise(
        across(where(is.numeric),   mean),
        across(where(is.character), max)
    )
cluster_mut <- cluster_mut[, -1]   # Drop the composite key column

# Reshape to a wide matrix: rows = cytobands, columns = cell groups
mut_wide <- dcast(cluster_mut, cytoband ~ cell_group_name, value.var = "state")

# Fill missing CNV state values with 3 (diploid baseline)
mut_wide[is.na(mut_wide)] <- 3

# ── Step 2: Derive chromosomal arm labels and filter ambiguous entries ────────────────
# Extract arm label (e.g., "1p", "2q") from the cytoband string
mut_wide$arm <- str_extract(mut_wide$cytoband, "^\\d+(p|q)")

# Remove cytobands that span both p and q arms (centromeric / ambiguous regions)
mut_wide <- mut_wide[!grepl("p-q", mut_wide$arm), ]
mut_wide <- mut_wide[!grepl("q-p", mut_wide$arm), ]

# ── Step 3: Compute arm-level mean CNV state and centre on diploid (state = 3) ───────
# Average all cytobands belonging to the same arm
arm_mut <- mut_wide %>%
    group_by(arm) %>%
    summarise(across(where(is.numeric), mean))

# Z-score scale each arm row using diploid state (3) as the centre
# Result: positive values = amplification, negative values = deletion
scale_arm_mut <- t(scale(
    t(arm_mut[, -1]),
    center = rep(3, nrow(arm_mut))
))
rownames(scale_arm_mut) <- arm_mut$arm

# ── Step 4: Up-weight WGS-confirmed somatic mutation arms ────────────────────────────
# For each patient, multiply the scaled CNV values on mutated arms by 10 + 1 = 11.
# Arms without WGS-confirmed mutations retain their original scaled value (× 1).
WGS <- read.table("WGS.mutarm.txt", header = TRUE)

for (i in seq_len(nrow(WGS))) {

    pid      <- WGS$PID[i]
    mut_arms <- unlist(strsplit(WGS$MutArm[i], ","))   # Comma-separated arm list

    # Identify columns belonging to the current patient
    patient_cols <- grepl(pid, colnames(scale_arm_mut))

    # Build a per-arm weight vector: 11 for mutated arms, 1 for all others
    arm_weight <- (rownames(scale_arm_mut) %in% mut_arms) * 10 + 1

    # Apply the weight and replace NAs introduced by missing arms with 0
    scale_arm_mut[, patient_cols] <- scale_arm_mut[, patient_cols] * arm_weight
    scale_arm_mut[is.na(scale_arm_mut)] <- 0
}

# ── Step 5: Compute per-cell-group CNV score ─────────────────────────────────────────
# CNV score = mean of absolute scaled values across all arms (higher → more aneuploid)
cnv_score <- data.frame(
    score = apply(scale_arm_mut, 2, function(x) mean(abs(x)))
)

write.table(
    cnv_score,
    "EpiImmu_6P.cnvscore.txt",
    sep       = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote     = FALSE
)

# ── Step 6: Join CNV scores to cell groupings and classify tumour vs. normal ─────────
# Reload saved CNV scores (ensures reproducibility from this point onward)
cnv_score <- read.table(
    "EpiImmu_6P.cnvscore.txt",
    sep    = "\t",
    header = TRUE
)

# Load inferCNV cell-grouping table (maps individual cells to subcluster groups)
grouping <- read.table(
    "../result/17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings",
    sep    = "\t",
    header = TRUE
)

# Attach CNV score to each row using cell_group_name as the join key
grouping$cnvscore <- cnv_score[grouping$cell_group_name, "score"]

# Cell groups with no matching CNV score (e.g., immune / stromal groups) → score of 0
grouping[is.na(grouping$cnvscore), "cnvscore"] <- 0

# Classify each cell group as tumour ("tEpithelial") or normal ("nEpithelial")
# Threshold: CNV score ≥ 0.65 indicates substantial copy-number aberration
grouping$TNtype <- ifelse(
    grouping$cnvscore < 0.65,
    "nEpithelial",   # Normal epithelial
    "tEpithelial"    # Tumour epithelial
)

# Save the annotated grouping table
write.table(
    grouping,
    "EpiImmu_6P.cnvscore_TNtype.txt",
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
)

message("Done. Results written to EpiImmu_6P.cnvscore_TNtype.txt")