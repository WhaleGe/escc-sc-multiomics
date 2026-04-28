
# Load command line arguments for flexible execution
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
pre    <- args[2]
day    <- args[3]

# Manual override for local/specific runs
outdir <- "./02.neu_score"
pre <- "Spatial_all_other"
day <- "260421"

# 1. Environment Setup
# ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(ggpubr)

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
setwd(outdir)

# 2. File Identification and Filtering
# ------------------------------------------------------------------------------
# List files for distance matrices and COSG scores
disfile <- list.files(pattern = "_m3_neu_dis.txt")
scorefile <- list.files(pattern = "_cosg.*txt")

# Exclude specific samples (e.g., Lymph Node or specific iterations)
exclude_patterns <- "LN_P|SstP03_T_0"
disfile <- disfile[!grepl(exclude_patterns, disfile)]
scorefile <- scorefile[!grepl(exclude_patterns, scorefile)]

# 3. Batch Processing of Neutrophil Subtypes
# ------------------------------------------------------------------------------
plist <- list()
target_types <- c("Neu", "Neu_CCL4", "Neu_DCN", "Neu_IFIT3", "Neu_S100A12", "Neu_SPP1")

for (k in target_types) {
    # Initialize an empty data frame for merging results
    combined_dat <- data.frame(sample = character(), dis = numeric(), 
                             score = numeric(), score_scale = numeric())
    
    for (i in 1:length(disfile)) {
        # Load distance and score datasets
        distance <- read.delim(disfile[i], sep = "\t", header = TRUE, check.names = FALSE)
        score_mat <- read.delim(scorefile[i], sep = "\t", header = TRUE, check.names = FALSE)
        
        # Min-Max normalization across all scores for the current sample
        score_mat <- score_mat %>% 
            mutate(across(everything(), ~ (. - min(.)) / (max(.) - min(.) + 1e-9)))

        # Define distance bins and sampling logic based on sample source
        is_sst <- grepl("Sst", disfile[i])
        bins <- if(is_sst) seq(10, 200, 10) else seq(50, 1000, 50)
        bin_width <- if(is_sst) 10 else 50
        
        sample_res <- data.frame(sample = character(), dis = numeric(), score = numeric())

        for (j in bins) {
            # Identify neutrophils within the specific distance bin from m3 cells
            target_neu <- apply(distance, 1, function(x) {
                colnames(distance)[which(j - bin_width <= x & x < j)]
            }) %>% unlist() %>% unique()
            
            # Extract scores if neutrophils are present in this bin
            if (length(target_neu) > 0) {
                tmp_df <- data.frame(
                    sample = gsub("_m3_neu_dis.txt", "", disfile[i]),
                    dis = j,
                    score = score_mat[target_neu, k],
                    row.names = NULL
                )
                sample_res <- rbind(sample_res, tmp_df)
            }
        }
        
        # Internal normalization for the specific sample's trend
        sample_res <- sample_res[complete.cases(sample_res), ]
        if (nrow(sample_res) > 0) {
            s_min <- min(sample_res$score)
            s_max <- max(sample_res$score)
            sample_res$score_scale <- (sample_res$score - s_min) / (s_max - s_min + 1e-9)
            
            # Adjust distance scale for Sst samples to match 50-1000 range
            if (is_sst) sample_res$dis <- 5 * sample_res$dis
            
            combined_dat <- rbind(combined_dat, sample_res)
        }
    }
    
    # 4. Statistical Analysis and Aggregation
    # ------------------------------------------------------------------------------
    combined_dat <- combined_dat[complete.cases(combined_dat), ]
    
    # Summarize mean, median, and SD per distance bin per sample
    dat_stat <- combined_dat %>% 
        group_by(sample, dis) %>% 
        summarise(
            count = n(),
            mean_score = mean(score, na.rm = TRUE),
            median_score = median(score, na.rm = TRUE),
            sd_score = sd(score, na.rm = TRUE),
            .groups = 'drop'
        ) %>% 
        group_by(sample) %>% 
        mutate(mean_scale = (mean_score - min(mean_score)) / (max(mean_score) - min(mean_score) + 1e-9)) %>% 
        as.data.frame()

    # Assign group (PT: Primary Tumor; MT: Metastatic Tumor)
    dat_stat$group <- ifelse(grepl("P06|P21|P22", dat_stat$sample), "PT", "MT")
    
    # Perform Wilcoxon test between MT and PT groups
    stat_test <- compare_means(mean_score ~ group, dat_stat, method = "wilcox.test")

    # 5. Visualization
    # ------------------------------------------------------------------------------
    plist[[k]] <- ggplot(data = dat_stat, aes(x = dis, y = mean_scale)) + 
        geom_point(aes(colour = group, shape = group, fill = group), size = 2) + 
        geom_smooth(aes(colour = group), span = 1, se = FALSE, method = "loess") +
        scale_color_manual(values = c("#8f7d78", "#f7b66f")) + 
        scale_x_continuous(
            limits = c(50, 1000), 
            breaks = seq(50, 1000, 50),
            labels = seq(50, 1000, 50) / 2 # Scale adjustment for physical distance
        ) +
        ylim(0, 1.1) + 
        theme_classic() + 
        labs(
            x = "Distance between neutrophils and m3 (um)", 
            y = paste0(k, " Mean Score"),
            title = paste0("Spatial Score Trend: ", k)
        ) + 
        annotate("text", x = 500, y = 1.05, label = paste0("MT vs PT, P.adj=", format.pval(stat_test$p.adj, digits = 3)), size = 4) +
        theme(
            axis.text.x = element_text(hjust = 1, vjust = 1, color = "black", angle = 45, size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.title = element_text(size = 14),
            plot.title = element_text(hjust = 0.5)
        )
}

# Export all plots to a single PDF
output_file <- paste0(pre, "_Neu_bym3dis1000_", day, ".pdf")
pdf(output_file, width = 7, height = 5)
for (p in plist) print(p)
dev.off()
