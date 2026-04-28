
# 1. Environment and Parameter Setup
# ------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
work_dir      <- args[1] # Output/Working directory
fragment_file <- args[2] # scATAC fragment file path
metadata_path <- args[3] # Metadata with cell-type labels
sample_id     <- args[4] # Current sample ID

library(rtracklayer)
library(Alleloscope)
setwd(work_dir)

# Load reference data
chr_size <- read.table("chrNameLength.txt", stringsAsFactors = FALSE, sep = '\t')
genome_size <- chr_size[1:22, ] # Focus on autosomes
bin_bed <- read.table("hg38arc.bin.500k.bed", header = FALSE, sep = "\t")

# 2. Bin-by-Cell Matrix Generation Function
# ------------------------------------------------------------------------------
generate_bin_matrix <- function(bin_bed, barcodes, fragment_path) {
    cat("Importing fragment file...\n")
    # Load fragments as GenomicRanges object
    frags <- import.bed(fragment_path, extraCols = c("type" = "character", "score" = "integer"))
    colnames(mcols(frags)) <- c("barcode", "dup_counts")
    
    # Filter for fragments within valid cells
    frags_in_cells <- frags[frags$barcode %in% barcodes]
    
    # Define query bins
    query_bins <- GRanges(paste0('chr', bin_bed[, 1]), 
                          IRanges(as.numeric(bin_bed[, 2]) + 1, as.numeric(bin_bed[, 3])))
    
    # Find overlaps between fragments and genomic bins
    overlaps <- findOverlaps(query_bins, frags_in_cells)
    overlaps_mat <- as.matrix(overlaps)
    
    cat("Generating contingency table...\n")
    # Build matrix: Rows = Bins, Columns = Cells
    bin_cell_mat <- table(overlaps_mat[, 1], match(frags_in_cells$barcode[overlaps_mat[, 2]], barcodes))
    return(bin_cell_mat)
}

# 3. CNV Processing
# ------------------------------------------------------------------------------
# Filter metadata for the specific sample
meta_data <- read.delim(metadata_path, header = TRUE, stringsAsFactors = FALSE)
sample_meta <- meta_data[meta_data$Sample == sample_id, ]
valid_barcodes <- rownames(sample_meta)

# Execute matrix generation
raw_counts <- generate_bin_matrix(bin_bed, valid_barcodes, fragment_file)

# Infer CNV and plot (Normal B/T cells used as reference)
cnv_plots <- plot_scATAC_cnv(
    raw_mat = as.matrix(raw_counts), 
    cell_type = as.matrix(sample_meta[, "Celltype", drop = FALSE]), 
    normal_lab = c("B cells", "T cells"), 
    size = genome_size, 
    plot_path = paste0(sample_id, "_cnv_profile.png")
)

# Export CNV Matrix
cnv_matrix <- as.data.frame(t(cnv_plots$plot_matrix))
write.table(cnv_matrix, paste0(sample_id, "_cnv_matrix.txt"), sep = "\t", quote = FALSE)