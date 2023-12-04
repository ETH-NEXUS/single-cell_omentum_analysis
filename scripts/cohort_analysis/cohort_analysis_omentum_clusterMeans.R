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

# have matrix with normcounts
# remove all genes that have 0 counts in all cells
normcounts_all.zero.removed <- normcounts(combined)
mask_all_zero <- apply(normcounts_all.zero.removed, 1, sum) > 0
normcounts_all.zero.removed <- normcounts_all.zero.removed[mask_all_zero, ]
stopifnot(length(normcounts_all.zero.removed[, 1]) == sum(apply(normcounts_all.zero.removed, 1, sum) > 0))

# get mean expression per umap cluster

all_umap_clusters <- unique(colData(combined)$umap_cohort_cl)
print("umap clusters:")
print(all_umap_clusters)

all_cluster_means <- as.data.frame(rownames(combined))

names(all_cluster_means) <- "gene"
all_cluster_means$gene <- as.character(all_cluster_means$gene)

# Use relative counts normalisation before calculating mean expression per cluster
# get mean expression for each umap cluster
for (cluster_ID in all_umap_clusters) {
  print(paste0("Calculating mean expression for cluster ", cluster_ID))
  col_header <- paste0("clust_", cluster_ID)
  print(col_header)
  subset_cluster <- combined[, combined$umap_cohort_cl == cluster_ID]
  #subset_cluster <- subset(combined, combined$umap_cohort_cl == cluster_ID)
  print(paste0("Number of cells in cluster: ", dim(subset_cluster)[2]))
  stopifnot(all_cluster_means$gene == rownames(subset_cluster))
  all_cluster_means[col_header] <- apply(assay(subset_cluster, "normcounts"), 1, mean)
  all_cluster_means[col_header] <- round(all_cluster_means[col_header], digits = 4)
  all_cluster_means[col_header][all_cluster_means[col_header] < 0] <- 0
}
# write mean expression table
filename_mean_expr <- paste0(opt$outdir, opt$outName, ".mean_expression_per_umap_cluster.tsv")
print(filename_mean_expr)
write.table(all_cluster_means, file = filename_mean_expr, quote = FALSE, sep = "\t", row.names = FALSE)


