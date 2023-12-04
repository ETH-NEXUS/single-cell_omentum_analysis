###################################################
## File name: cohort_analysis_omentum.R 
## Adapted from: analyse_scRNA_cohort_gexOnly.R
## Author: Anne Bertolini
## Date created: September 2020
## R Version: 3.5.1
###################################################

## GENERAL:
## This script combines GEX results of different samples into one object.
## Those results can then be visualised together.
## No integration or normalisation step is included.

library(SingleCellExperiment)
library(ggplot2)
library(scMerge)
library(uwot)
library(rhdf5)
library(stringr)
library(patchwork)
library(optparse)
library(limma)
library(cowplot)
library(reshape2)
library(ggbeeswarm)
library(igraph)

cat("\n\n\nPrint sessionInfo():\n\n")
print(sessionInfo())

# parse command line arguments
option_list <- list(
  make_option("--rdsFile", type = "character", help = "RDS file with merged single cell data."),
  make_option("--outName", type = "character", help = "Sample name that will be prefix of all plot names."),
  make_option("--outdir", type = "character", help = "Path to the outdir, in which files will be written.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""

###################################
### General plotting parameters ###
###################################

# aspect ratio of UMAP plots
UMAP_aspect_ratio <- "1"

###################################
###   Read in and format data   ###
###################################

combined = readRDS(opt$rdsFile)
colData_df <- as.data.frame(colData(combined))

cat("# plot UMAP clustering\n")

df_labels <- data.frame(matrix(ncol = 3, nrow = 0))
myNames <- c("center_of_mass_x", "center_of_mass_y", "clusterID")
colnames(df_labels) <- myNames

clusterIDs = unique(colData(combined)$umap_cohort_cl)
for (clust in clusterIDs){
  print(clust)
  subsetClust = combined[, colData(combined)$umap_cohort_cl == clust]
  xMed = median(reducedDim(subsetClust)$GEX_umap1)
  yMed = median(reducedDim(subsetClust)$GEX_umap2)
  newRow <- c(xMed,yMed,clust)
  df_labels <- rbind(df_labels,newRow)
}

df_labels
colnames(df_labels) <- myNames
dim(df_labels)

ggplot(colData_df, aes(x = reducedDim(combined)$GEX_umap1, y = reducedDim(combined)$GEX_umap2, color = factor(colData(combined)$umap_cohort_cl))) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20, title="Cluster IDs")) +
  geom_label(data=df_labels, inherit.aes = FALSE, aes(x = center_of_mass_x, y = center_of_mass_y, label = clusterID)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Umap_clustering_labels.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")
