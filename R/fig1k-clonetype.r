# Load command line arguments
args       <- commandArgs(trailingOnly = TRUE)
inSTfile1  <- args[1]   # Path to the first Seurat object (e.g., Lymph Node)
inSTfile2  <- args[2]   # Path to the second Seurat object (e.g., Primary Tumor)
pre        <- args[3]   # Sample prefix for naming
day        <- args[4]   # Date string or identifier
outdir     <- args[5]   # Output directory

# ------------------------------------------------------------------------------
# 1. Environment Setup
# ------------------------------------------------------------------------------
# Load required libraries using pacman for efficiency
pacman::p_load(Seurat, ggplot2, SCEVAN, dplyr, pagoda2, patchwork)

# Create output directory if it doesn't exist and set as working directory
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
setwd(outdir)

# ------------------------------------------------------------------------------
# 2. Data Loading and Cell Selection
# ------------------------------------------------------------------------------
# Load Spatial Seurat objects
obj1 <- readRDS(inSTfile1)
obj2 <- readRDS(inSTfile2)

# Function to sample reference cells and select malignant cells
get_selected_cells <- function(seurat_obj) {
    # Identify target populations based on subtype metadata
    epi_cells <- colnames(seurat_obj)[seurat_obj$subtype == "Malignant"]
    fibro_cells <- colnames(seurat_obj)[seurat_obj$subtype == "Fibroblasts"]
    bcells <- colnames(seurat_obj)[seurat_obj$subtype == "Bcells"]
    
    # Randomly sample up to 500 cells for each reference type (Normal cells)
    set.seed(42) # Ensure reproducibility for sampling
    fibro_sampled <- sample(fibro_cells, size = min(500, length(fibro_cells)))
    bcells_sampled <- sample(bcells, size = min(500, length(bcells)))
    
    return(list(
        all_selected = c(epi_cells, fibro_sampled, bcells_sampled),
        refs = c(fibro_sampled, bcells_sampled)
    ))
}

# Process both objects for cell selection
selection1 <- get_selected_cells(obj1)
selection2 <- get_selected_cells(obj2)

# Subset and merge the objects
obj1.s <- subset(obj1, cells = selection1$all_selected)
obj2.s <- subset(obj2, cells = selection2$all_selected)
obj.combined <- merge(obj1.s, obj2.s)

# Extract raw counts and clean barcodes (remove Seurat's merge suffixes _1, _2)
count_matrix <- obj.combined@assays$Spatial@counts
colnames(count_matrix) <- gsub("_1$|_2$", "", colnames(count_matrix))

# Define the combined normal cell list for SCEVAN baseline
norm_names <- c(selection1$refs, selection2$refs)

# ------------------------------------------------------------------------------
# 3. SCEVAN Pipeline: CNA and Subclone Analysis
# ------------------------------------------------------------------------------
# Run the SCEVAN pipeline to infer CNA and identify subclones
# Note: par_cores is set to 20 for high-performance computing environments
results <- SCEVAN::pipelineCNA(
    count_matrix, 
    sample = pre, 
    norm_cell = norm_names, 
    FIXED_NORMAL_CELLS = TRUE, 
    par_cores = 20, 
    SUBCLONES = TRUE, 
    plotTree = TRUE
)

# Export the raw SCEVAN results
output_txt <- paste0(pre, "_SCEVAN_subclone_", day, ".txt")
write.table(results, output_txt, quote = FALSE, row.names = TRUE, col.names = TRUE, sep = '\t')

# ------------------------------------------------------------------------------
# 4. Result Integration and Final Export
# ------------------------------------------------------------------------------
# Create a new Seurat object to store the results and metadata
obj_final <- CreateSeuratObject(
    counts = count_matrix, 
    project = 'Spatial', 
    assay = 'Spatial', 
    min.cells = 0, 
    min.features = 0
)

# Calculate mitochondrial percentage
obj_final[["percent.mt"]] <- PercentageFeatureSet(obj_final, pattern = "^MT-")

# Clean up identity and synchronize metadata from the combined subset
obj_final$orig.ident <- gsub("_\\d+$", "", rownames(obj_final@meta.data))

# Transfer original metadata (excluding redundant standard columns)
meta_to_add <- obj.combined@meta.data %>% 
    select(-orig.ident, -nCount_Spatial, -nFeature_Spatial, -percent.mt)
obj_final <- AddMetaData(obj_final, metadata = meta_to_add)

# Add SCEVAN subclone and CNA results to metadata
obj_final <- AddMetaData(obj_final, metadata = results)

# Save the final processed object
output_rds <- paste0(pre, "_SCEVAN_subclone_", day, ".rds")
saveRDS(obj_final, output_rds)
