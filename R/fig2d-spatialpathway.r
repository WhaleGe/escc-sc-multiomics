
# Set working directory
# setwd("./04.malignant_pathway") # Adjusted based on environment

# Load required libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)

# ------------------------------------------------------------------------------
# 1. Data Loading and Preprocessing
# ------------------------------------------------------------------------------
# Load pathway statistics data
res.all <- read.table("Spatial_cb_m3-others_pathway_stat.txt", header = TRUE)

# Calculate -log10(p-value) and log2 fold change
# Added a small constant to p-value to avoid log(0)
res.all$logp <- -log10(res.all$p + 0.0001)
res.all$log2fc <- log(res.all$fc)

# Cap log2fc values between -1 and 1 for better visualization stability
res.all$log2fc <- ifelse(res.all$log2fc > 1, 1, res.all$log2fc)
res.all$log2fc <- ifelse(res.all$log2fc < (-1), -1, res.all$log2fc)

# Define Primary Tumor (PT) and Metastatic Tumor (MT) groups
# Grouping based on predefined sample names
res.all$group <- ifelse(res.all$sample %in% c("P06_T_0", "P21_T_0", "P22_T_0"), "PT", "MT")

# ------------------------------------------------------------------------------
# 2. Statistical Aggregation (Metastatic Group Focus)
# ------------------------------------------------------------------------------
# Subset data for Metastatic Tumor (MT) group
res.mt <- subset(res.all, group %in% c("MT"))

# Calculate comprehensive statistics per pathway
pathway_stats <- res.mt %>%
    group_by(pathway) %>%
    summarise(
        total_samples = n(),                                     # Total number of samples
        significant_count = sum(p < 0.05, na.rm = TRUE),        # Number of samples with p < 0.05
        significant_proportion = mean(p < 0.05, na.rm = TRUE),   # Proportion of significant samples
        significant_percentage = round(mean(p < 0.05, na.rm = TRUE) * 100, 2),
        mean_log2fc = round(mean(log2fc, na.rm = TRUE), 4),      # Mean log2FC
        median_log2fc = round(median(log2fc, na.rm = TRUE), 4),  # Median log2FC
        sd_log2fc = round(sd(log2fc, na.rm = TRUE), 4),          # Standard deviation of log2FC
        .groups = 'drop'
    ) %>%
    # Sort pathways by their mean log2FC for visualization ordering
    arrange(mean_log2fc)

# Reorder factor levels for the Y-axis based on the calculated mean log2fc
pathway_stats$pathway <- factor(pathway_stats$pathway, levels = pathway_stats$pathway)

# ------------------------------------------------------------------------------
# 3. Visualization: Bar Plot
# ------------------------------------------------------------------------------
p1 <- ggplot(pathway_stats, aes(x = mean_log2fc, y = pathway)) + 
    # Create bar plot where fill represents the proportion of significant samples
    geom_bar(aes(fill = significant_proportion), width = 0.8, stat = "identity") +
    # Use a blue-to-red gradient for significance levels
    scale_fill_gradientn(colors = c("blue", "red"), name = "Sig. Proportion") +
    # Fixed x-axis limits to match the log2fc capping
    xlim(-1, 1) +
    labs(
        x = "Mean log2 Fold Change",
        y = "Biological Pathway",
        title = "Pathway Activity in MT Samples"
    ) +
    theme_classic() + 
    theme(
        text = element_text(family = "sans", color = "grey30", size = 12),
        axis.text.y = element_text(size = 10),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(), 
        axis.ticks = element_line(linewidth = 0.6, color = "gray30"),
        axis.ticks.length = unit(1.5, units = "mm"),
        plot.margin = unit(x = c(0.2, 0.5, 0.2, 0.1), units = "inches"),
        legend.position = "right"
    )

# Save the plot as PDF
pdf("Spatial_cb_m3-others_pathway_stat_bar.pdf", 12, 5)
print(p1)
dev.off()


