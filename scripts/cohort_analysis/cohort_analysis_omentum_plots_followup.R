###################################################
## File name: cohort_analysis_omentum_plots_followUp.R 
## R Version: 3.5.1
###################################################

## GENERAL:
## Load cohort SCE object and perform follow-up analysis + plotting

library(Seurat)
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
  make_option("--metadata_table", type = "character", help = "Table with columns RDS_files (paths to input RDS files that contain SCE objects with GEX and ADT results), h5_files (paths to h5 files that contain GEX most variable genes), and sample_name."),
  make_option("--colour_config", type = "character", help = "Colour config for visualising the cell types. Sometimes union of different cell type colour config files is required."),
  make_option("--selected_genes", type = "character", help = "Table with selected genes that will be shown in gene expression plots."),
  make_option("--sampleName", type = "character", help = "Sample name that will be prefix of all plot names."),
  make_option("--toggle_label", type = "logical", action = "store_true", default = TRUE,  help = "Set gene plot title labels to include user-defined gene aliases. (On by default.)"),
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

# create sub directories to group plots
dir_gene_expr <- opt$outdir %&% "gene_expression/"
dir.create(dir_gene_expr, showWarnings = FALSE)
path_gene_expr <- opt$outdir %&% "gene_expression/" %&% opt$sampleName

dir_qc <- opt$outdir %&% "qc_plots/"
dir.create(dir_qc, showWarnings = FALSE)
path_qc_plots <- opt$outdir %&% "qc_plots/" %&% opt$sampleName

###################################
###   Read in and format data   ###
###################################

# read in metadata table
cat("\nRead in metadata table:\n")

metadata <- read.table("gfb_omentum_2020/git/gfb_omentum_2020/omentum_required/metadata_omentum_cohort.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
str(metadata)

# read in colour config
cat("\nRead in Color-Config\n")
config <- read.csv("gfb_omentum_2020/git/gfb_omentum_2020/omentum_required/2021-04_reanalysis/colour_config_omentum_2021-05.txt", sep = "\t", stringsAsFactors = FALSE)
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(config$colour == unique(config$colour))
ct.color <- c(config$colour, "grey50", "black")
names(ct.color) <- c(config$cell_type, "uncertain", "unknown")
print(ct.color)

cohort_sce = readRDS("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/cohort/OmentumCohort.merged_cohort_nn30_mD03.RDS")
str(cohort_sce)
###############################################
combined = readRDS(rdsFile)

#metadata(combined) = c(metadata(combined), list(umap_modularity=modularity(km)))

colData(combined)$site  = metadata_table$site
for sample in mySamples{
  mySite = dict(sample)
  combined$site[sce$sample_name == sample] = mySite
}

colData_df <- as.data.frame(colData(combined))
umap_gex_df <- as.data.frame(reducedDim(combined, "umap_gex"))
# add umap coordinates to the meta data table of cells
colData_df <- plyr::join(colData_df, umap_gex_df, by = "barcodes")
######################################

# get mean expression per umap cluster

all_umap_clusters <- levels(combined[[]]$umap_cohort_cl)
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
  subset_cluster <- subset(x = combined, subset = umap_cohort_cl == cluster_ID)
  print(paste0("Number of cells in cluster: ", dim(subset_cluster)[2]))
  stopifnot(all_cluster_means$gene == rownames(subset_cluster))
  all_cluster_means[col_header] <- apply(assay(subset_cluster, "normcounts"), 1, mean)
  all_cluster_means[col_header] <- round(all_cluster_means[col_header], digits = 4)
}
# write mean expression table
filename_mean_expr <- paste0(opt$outdir, opt$sampleName, ".mean_expression_per_umap_cluster.tsv")
print(filename_mean_expr)
write.table(all_cluster_means, file = filename_mean_expr, quote = FALSE, sep = "\t", row.names = FALSE)


