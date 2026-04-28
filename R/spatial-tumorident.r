
# ── Command-line arguments ────────────────────────────────────────────────────────────
# Usage: Rscript spatial-tumorident.r <input_rds_file>
# Example (for debugging, uncomment the line below):
# infile <- "spacet.rds"

args  <- commandArgs(trailingOnly = TRUE)
infile <- args[1]

# ── Load required libraries ───────────────────────────────────────────────────────────
library(SpaCET)   # Spatial cell-type deconvolution
library(Seurat)   # Single-cell / spatial data framework
library(dplyr)    # Data manipulation utilities

# ── Main workflow ─────────────────────────────────────────────────────────────────────
# Branch on whether the input file is already a SpaCET object (filename contains "spacet")
# or a raw Seurat object that needs to be converted first.

if (!grepl("spacet", basename(infile))) {

    # ── Branch 1: Input is a raw Seurat object ────────────────────────────────────────

    message("Input detected as a Seurat object. Creating SpaCET object...")

    # Load the Seurat object from disk
    obj.bulk <- readRDS(infile)

    # Extract the raw count matrix from the Spatial assay
    counts <- obj.bulk@assays$Spatial@counts

    # Remove named dimensions on column names to avoid downstream conflicts
    coln        <- colnames(counts)
    names(coln) <- NULL
    colnames(counts) <- coln

    # Create a SpaCET object using counts and spot spatial coordinates
    SpaCET_obj <- create.SpaCET.object(
        counts          = counts,
        spotCoordinates = obj.bulk@meta.data[, c("x", "y")],
        imagePath       = NA,          # No image path provided for Stereo-seq
        platform        = "Stereo-seq"
    )

    # Perform cell-type deconvolution for esophageal carcinoma (ESCA)
    SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType = "ESCA", coreNo = 8)

    # Save the SpaCET object; update filename prefix from "epi_" to "epi_spacet_"
    out_path <- gsub("epi_", "epi_spacet_", infile)
    message("Saving SpaCET object to: ", out_path)
    saveRDS(SpaCET_obj, out_path)

} else {

    # ── Branch 2: Input is already a SpaCET object ────────────────────────────────────

    message("Input detected as a SpaCET object. Running deconvolution directly...")

    # Load the existing SpaCET object from disk
    SpaCET_obj <- readRDS(infile)

    # Perform cell-type deconvolution for esophageal carcinoma (ESCA)
    SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType = "ESCA", coreNo = 8)

    # Save the result; update filename prefix from "epi_spacet_" to "epi_spacet_res_"
    out_path <- gsub("epi_spacet_", "epi_spacet_res_", infile)
    message("Saving deconvolution results to: ", out_path)
    saveRDS(SpaCET_obj, out_path)

}