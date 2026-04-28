
# Load required libraries
library(Seurat)
library(dplyr)
library(NMF)
library(reticulate)
library(cowplot)
library(ComplexHeatmap)
library(ggsci)
library(ggplot2)

# ------------------------------------------------------------------------------
# 1. Data Loading and Initial Filtering
# ------------------------------------------------------------------------------
# Load the predefined Seurat object
obj <- readRDS("tumor_dim15res0.5k20_Cluster.rds")

# Filter out mitochondrial, ribosomal, and non-coding RNA genes
gene <- rownames(obj)[grep("^MT-|^AC[0-9]|^AL[0-9]|^LNC|^LINC|^IGH|^IGK|^IGL|^RPS|^RPL", invert = TRUE, rownames(obj))]
obj <- subset(obj, features = gene)

# Split object by sample
samp.list <- SplitObject(obj, split.by = "Sample")

# Filter out samples containing fewer than 100 tumor cells
samp.list <- samp.list[lapply(samp.list, ncol) > 100]

# ------------------------------------------------------------------------------
# 2. Dimensionality Reduction (PCA) per Sample
# ------------------------------------------------------------------------------
pcaplot <- list()

for (i in 1:length(samp.list)) {
    # Extract valid genes for the current sample
    current_genes <- rownames(samp.list[[i]])[grep("^MT-|^AC[0-9]|^AL[0-9]|^LNC|^LINC|^IGH|^IGK|^IGL|^RPS|^RPL", invert = TRUE, rownames(samp.list[[i]]))]
    samp.list[[i]] <- subset(samp.list[[i]], features = current_genes)
    
    # Normalize, find variable features, and scale data
    samp.list[[i]] <- NormalizeData(samp.list[[i]]) %>% 
        FindVariableFeatures() %>% 
        ScaleData(do.center = TRUE)
    
    # Run PCA and generate Elbow plots
    samp.list[[i]] <- RunPCA(samp.list[[i]], npcs = 60)
    pcaplot[[i]] <- ElbowPlot(samp.list[[i]], ndims = 60)
}

# Export PCA Elbow plots
pdf("PCAPlot.pdf", 5, 2 * length(pcaplot))
plot_grid(plotlist = pcaplot, ncol = 1)
dev.off()

# ------------------------------------------------------------------------------
# 3. Clustering, UMAP, and Differentially Expressed Genes (DEGs)
# ------------------------------------------------------------------------------
# Define custom PCA dimensions and resolutions for each sample
npcs <- c(25, rep(30, 2), 20, 30, 30, 20, rep(30, 3))
res <- c(0.3, rep(0.4, 4), 0.5, rep(0.4, 4))

umapplot <- list()
heatdeg <- list()

for (i in 1:length(samp.list)) {
    # Find neighbors and perform clustering
    samp.list[[i]] <- FindNeighbors(samp.list[[i]], reduction = "pca", dims = 1:npcs[i], nn.eps = 0, k.param = 10)
    samp.list[[i]] <- FindClusters(samp.list[[i]], resolution = res[i], n.start = 10, algorithm = 1)
    
    # Run UMAP for visualization
    samp.list[[i]] <- RunUMAP(samp.list[[i]], dims = 1:npcs[i], verbose = FALSE, min.dist = 0.2)
    umapplot[[i]] <- DimPlot(samp.list[[i]], reduction = "umap", pt.size = 0.5, label = TRUE)
    
    # Identify DEGs for each cluster
    markers <- FindAllMarkers(samp.list[[i]], only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    filter_markers <- markers[markers$p_val < 0.05, ]
    top20 <- filter_markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
    
    # Export DEG lists
    write.table(filter_markers, paste0(names(samp.list)[i], ".deg.txt"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    
    # Generate heatmap for top 20 DEGs
    heatdeg[[i]] <- DoHeatmap(samp.list[[i]], features = top20$gene, size = 6) + 
        theme(text = element_text(size = 5))
}

# Export UMAP plots
pdf("Umap.pdf", 9, 4 * length(umapplot) / 2)
plot_grid(plotlist = umapplot, ncol = 2)
dev.off()

# Export Top 20 DEG heatmaps
pdf("Top20DEG.heat.pdf", 6, 6 * length(umapplot))
plot_grid(plotlist = heatdeg, ncol = 1)
dev.off()

save(samp.list, file = "samp.list.Rdata")

# ------------------------------------------------------------------------------
# 4. Preparation for Non-negative Matrix Factorization (NMF)
# ------------------------------------------------------------------------------
Input <- list()

for (i in 1:length(samp.list)) {
    samp <- names(samp.list)[i]
    counts <- as.matrix(samp.list[[i]]@assays$RNA@counts)
    
    # Select genes expressed in more than 10 cells
    gene <- rownames(counts)[apply(counts, 1, function(x) sum(x > 0) > 10)]
    
    # Scale data for NMF input and convert negative values to 0
    samp.list[[i]] <- ScaleData(samp.list[[i]], features = gene, scale.max = 10)
    mat <- as.matrix(samp.list[[i]]@assays$RNA@scale.data)
    mat[mat < 0] <- 0
    
    Input[[i]] <- mat
    write.table(mat, paste0(samp, ".nmf.input.txt"), row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
}

# ------------------------------------------------------------------------------
# 5. Run NMF Analysis
# ------------------------------------------------------------------------------
# Predefined NMF ranks per sample
nmfn <- c(6, 2, 4, 7, 10, 8, 6, 5, 5, 9)

GENES <- list()
CELLS <- list()

for (i in 1:length(nmfn)) {
    samp <- names(samp.list)[i]
    mat <- Input[[i]]
    
    # Execute NMF using the 'snmf/r' method
    res <- nmf(mat, rank = nmfn[i], method = "snmf/r")
    
    # Extract gene programs (basis)
    genes <- basis(res)
    program <- paste0(samp, ".p", 1:ncol(genes))
    colnames(genes) <- program
    GENES[[i]] <- genes
    names(GENES)[i] <- samp
    
    # Extract cell usage (coefficients)
    cells <- coef(res)
    rownames(cells) <- program
    CELLS[[i]] <- cells
    names(CELLS)[i] <- names(samp.list)[i]
}

save(GENES, file = "GENES.Rdata")
save(CELLS, file = "CELLS.Rdata")

# ------------------------------------------------------------------------------
# 6. Filter Programs and Evaluate Correlations
# ------------------------------------------------------------------------------
program_score <- c()
SCELLS <- list()

# Calculate standard deviation for each NMF program across cells
for (i in 1:length(CELLS)) {
    cells <- t(CELLS[[i]])
    scells <- cells
    SCELLS[[i]] <- scells
    
    if (i == 1) {
        program_score <- apply(scells, 2, sd)
    } else {
        program_score <- c(program_score, apply(scells, 2, sd))
    }
}

# Retain programs with SD >= 0.8
retainp <- names(program_score)[program_score >= 0.8]
olpGene <- Reduce(intersect, lapply(GENES, rownames))
sGENES <- lapply(GENES, function(x) x[rownames(x) %in% olpGene, ])
pmat <- do.call(cbind, sGENES)

# Extract Top 100 genes per NMF program
topn <- 100
igenelist <- list()

for (i in 1:length(sGENES)) {
    for (j in 1:ncol(sGENES[[i]])) {
        genes <- sort(sGENES[[i]][, j], decreasing = TRUE)[1:topn]
        genes <- data.frame(Load = genes, Gene = names(genes))
        colnames(genes)[1] <- colnames(sGENES[[i]])[j]
        
        if (j == 1) {
            metagenes <- genes
        } else {
            metagenes <- merge(metagenes, genes, by = "Gene", all = TRUE)
        }
    }
    metagenes[is.na(metagenes)] <- 0
    rownames(metagenes) <- metagenes$Gene
    igenelist[[i]] <- metagenes
}
save(igenelist, file = "TOP100.igenelist.Rdata")

# Calculate program correlation matrix based on union of top genes
pmat <- pmat[rownames(pmat) %in% unique(unlist(lapply(igenelist, rownames))), ]
cor <- cor(pmat)
save(cor, file = "TOP100.cor.Rdata")

# ------------------------------------------------------------------------------
# 7. Visualization of Program Correlation (Heatmap)
# ------------------------------------------------------------------------------
pcol <- pal_jco()(length(SCELLS))
names(pcol) <- names(GENES)

row_ha = rowAnnotation(PID = gsub("\\.\\S+$", "", rownames(cor)), col = list(PID = pcol))
top_ha = HeatmapAnnotation(PID = gsub("\\.\\S+$", "", rownames(cor)), col = list(PID = pcol))

# Cap correlation matrix values for visualization: max 0.8, min 0
test <- abs(cor)
test[test > 0.8] <- 0.8
test[test < 0] <- 0

p1 <- Heatmap(
    test, 
    name = "Cor", 
    left_annotation = row_ha, 
    top_annotation = top_ha, 
    show_row_names = TRUE,
    show_column_names = FALSE,  
    cluster_rows = TRUE, 
    cluster_columns = TRUE, 
    row_names_gp = gpar(fontsize = 8), 
    row_title_gp = gpar(fontsize = 10), 
    row_title_rot = 0, 
    col = colorRampPalette(rev(c("red4", "gold1", "white")))(10), 
    clustering_distance_rows = function(x) as.dist(1 - x), 
    clustering_distance_columns = function(x) as.dist(1 - x), 
    cell_fun = function(j, i, x, y, width, height, fill) {
        if (test[i, j] > 0.5 && test[i, j] < 0.8) {
            grid.text("*", x, y, gp = gpar(fontsize = 10))
        }
    }
)

pdf("cor.heat.pdf", 10, 8)
p1
dev.off()

# ------------------------------------------------------------------------------
# 8. Meta-program Definition and Top Gene Extraction
# ------------------------------------------------------------------------------
# Define specific meta-programs manually based on correlation hierarchy
m1 <- c("P28_T_0.p8", "P21_T_0.p4", "P26_T_0.p5")
m2 <- c("P21_T_0.p1", "P06_T_0.p2", "P22_T_0.p4")
m3 <- c("P27_T_0.p3", "P26_T_0.p4", "P08_T_0.p4", "P07_T_0.p2")
m4 <- c("P25_T_0.p6", "P06_T_0.p4", "P21_T_0.p5")
m5 <- c("P28_T_0.p2", "P25_T_0.p3", "P27_T_0.p2", "P06_T_0.p5", "P21_T_0.p3")
m6 <- c("P06_T_0.p1", "P21_T_0.p8", "P14_T_0.p6")
m7 <- c("P08_T_0.p3", "P28_T_0.p3", "P06_T_0.p6", "P27_T_0.p4", "P25_T_0.p1", "P22_T_0.p1", "P26_T_0.p3", "P21_T_0.p2", "P14_T_0.p4")
m8 <- c("P25_T_0.p4", "P21_T_0.p6", "P22_T_0.p6", "P26_T_0.p2", "P28_T_0.p9", "P08_T_0.p1", "P27_T_0.p5", "P06_T_0.p3", "P07_T_0.p1")

clusters <- list(m1, m2, m3, m4, m5, m6, m7, m8)
names(clusters) <- paste0("m", 1:8)
topn <- 100

genelist <- data.frame(matrix(nrow = topn, ncol = length(clusters)))
colnames(genelist) <- names(clusters)

# Approach 2: Aggregating all genes before extracting Top N (Active execution)
for (i in 1:length(clusters)) {
    for (j in 1:length(clusters[[i]])) {
        # Sort all genes for the respective program
        genes <- sort(pmat[, clusters[[i]][j]], decreasing = TRUE)
        genes <- data.frame(Load = genes, Gene = names(genes))
        colnames(genes)[1] <- clusters[[i]][j]
        
        # Merge gene loadings progressively
        if (j == 1) {
            metagenes <- genes
        } else {
            metagenes <- merge(metagenes, genes, by = "Gene", all = TRUE)
        }
    }
    
    # Fill missing values and set row names
    metagenes[is.na(metagenes)] <- 0
    rownames(metagenes) <- metagenes$Gene
    
    # Calculate cumulative score across aggregated programs and extract top N genes
    genelist[, i] <- names(sort(apply(metagenes[, 2:ncol(metagenes), drop = FALSE], 1, sum), decreasing = TRUE))[1:topn]
}

# Export final meta-program signatures
write.table(genelist, paste0("Top", topn, ".genelist.2.txt"), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)