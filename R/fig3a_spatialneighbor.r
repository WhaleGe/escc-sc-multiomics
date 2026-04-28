
# Load command line arguments for flexible execution
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
pre    <- args[2]
day    <- args[3]

# Define output and naming parameters
outdir <- "./03.m3_neighbor"
pre <- "Spatial_cb"
day <- "250603"

# 1. Load Environment and Packages
# ------------------------------------------------------------------------------
# pacman::p_load ensures all packages are installed and loaded
pacman::p_load(Seurat, cowplot, ggplot2, dplyr, ggpubr, ggsci, fmsb, tibble)

# Create output directory if it doesn't exist and set as working directory
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
setwd(outdir)

# Define color palette for different samples and groups
pals <- c("#9e8e8a", "#f8c082", "#f7b66f", "#f8c082", "#9e8e8a", "#9e8e8a", "#8f7d78", "#9e8e8a", "#f8c082")

# 2. Data Acquisition and Merging
# ------------------------------------------------------------------------------
# List files for neighborhood proportions and p-values
file_prop <- list.files(pattern = "dd.prop.scale.txt")
file_pval <- list.files(pattern = "dd.p.txt")

# Initialize data frames with the first sample's data (targeting the 'm3' program)
# Using check.names=FALSE to preserve original sample naming formats
dat_s <- read.delim(file_prop[1], sep = "\t", header = TRUE, check.names = FALSE)
dat_s <- dat_s[, "m3", drop = FALSE]

dat_p <- read.delim(file_pval[1], sep = "\t", header = TRUE, check.names = FALSE)
dat_p <- dat_p[, "m3", drop = FALSE]

# Iteratively merge 'm3' neighborhood data from all detected files
for (i in 2:length(file_prop)) {
    # Process proportion files
    dat_s_temp <- read.delim(file_prop[i], sep = "\t", header = TRUE, check.names = FALSE)
    dat_s_temp <- dat_s_temp[, "m3", drop = FALSE]
    dat_s <- cbind.data.frame(dat_s, dat_s_temp)
    
    # Process p-value files
    dat_p_temp <- read.delim(file_pval[i], sep = "\t", header = TRUE, check.names = FALSE)
    dat_p_temp <- dat_p_temp[, "m3", drop = FALSE]
    dat_p <- cbind.data.frame(dat_p, dat_p_temp)
}

# Assign sample IDs as column names (extracted from filenames)
sample_ids <- c("P04_T_0", "P06_T_0", "P21_T_0", "P22_T_0", "P26_T_0", "P27_T_0", "P28_T_0")
colnames(dat_s) <- sample_ids
colnames(dat_p) <- sample_ids

# 3. Statistical Grouping and Median Calculation
# ------------------------------------------------------------------------------
# Transpose data and assign Primary Tumor (PT) vs Metastatic Tumor (MT) status
dat_st <- t(dat_s) %>% 
    as.data.frame() %>% 
    mutate(group = c("MT", "PT", "PT", "PT", "MT", "MT", "MT"))

# Calculate median neighborhood scores for each group (PT vs MT)
dat_st_m <- aggregate(dat_st[, 1:20], by = list(dat_st$group), FUN = "median") %>% 
    column_to_rownames("Group.1")

# 4. Radar Chart Preparation
# ------------------------------------------------------------------------------
# Define the maximum and minimum boundaries for the radar chart axes
var_max <- rep(2.5, 20)
var_min <- rep(-2.5, 20)

# Combine boundaries, individual samples, and group medians into a single data frame
# Row order: 1.Max, 2.Min, 3-9.Samples, 10-11.Group Medians
dat_plot <- rbind(var_max, var_min, as.data.frame(t(dat_s)), dat_st_m)

# 5. Visualization: Radar Chart
# ------------------------------------------------------------------------------
pdf(paste0(pre, "_m3_neighborhood_radar1_", day, ".pdf"), width = 8, height = 8)

radarchart(
    dat_plot,
    axistype = 1, 
    seg = 10,                          # Number of segments for the axis
    pty = c(rep(32, 7), 16, 16),       # Point symbol types
    pcol = pals,                       # Line/Point colors
    plty = 7,                          # Line type (dash/dotted)
    plwd = c(rep(0.75, 7), 3, 3),      # Line width (thicker lines for group medians)
    cglcol = "grey",                   # Background grid color
    cglty = 1,                         # Grid line type
    axislabcol = "grey",               # Axis label color
    cglwd = 1.5,                       # Grid line width
    caxislabels = seq(-2.5, 2.5, 0.5), # Axis labels from min to max
    vlcex = 0.8,                       # Text size for variables
    title = "m3 neighborhood"
)

# Add legend for sample identification
legend(
    x = -1.42, y = -0.8, 
    legend = rownames(dat_plot[c(3:9), ]), 
    bty = "n", pch = 20, col = pals, 
    text.col = pals, cex = 1, pt.cex = 1
)
dev.off()

