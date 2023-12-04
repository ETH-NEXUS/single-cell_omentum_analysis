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
path_qc_plots <- opt$outdir %&% "qc_plots/" %&% opt$sampleName

###################################
###   Read in and format data   ###
###################################

# read in metadata table
cat("\nRead in metadata table:\n")
metadata <- read.table(opt$metadata_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
print(metadata)

# read in data from all samples of the cohort
cat("\nRead in one SCE object per sample:\n")
all_sce_files <- metadata$RDS_files
cell_desc <- character()

# make a list of all gex SCE objects
all_sce_objects <- lapply(seq_along(all_sce_files), function(x) {
  cat("\n\nRead in sample ", x, "\n")
  cat(all_sce_files[x], "\n")
  my_sce <- readRDS(all_sce_files[x])
  my_sce@metadata$sample_name <- metadata$sample_name[x]
  my_sce@metadata$file_name <- basename(all_sce_files[x])
  cat(my_sce@metadata$sample_name, "\n")
  colData(my_sce)$sampleID <- metadata$Sample_ID[x]
  colData(my_sce)$file_name <- basename(all_sce_files[x])
  colData(my_sce)$Date_of_Exp <- metadata$Date_of_Exp[x]
  colData(my_sce)$Cohort_Type_Old <- metadata$Cohort_Type_Old[x]
  colData(my_sce)$Cohort_Type <- metadata$Cohort_Type[x]
  colData(my_sce)$Cohort_Type_Highlevel <- metadata$Cohort_Type_Highlevel[x]
  colData(my_sce)$Exp <- metadata$Exp[x]
  colData(my_sce)$Patient <- metadata$Patient[x]
  colData(my_sce)$Site <- metadata$Site[x]
  colData(my_sce)$Tissue_Region_old <- metadata$Tissue_Region_old[x]
  colData(my_sce)$Tissue_Region <- metadata$Tissue_Region[x]
  colData(my_sce)$Diagnosis <- metadata$Diagnosis[x]
  colData(my_sce)$Primary_Relapse <- metadata$Primary_Relapse[x]
  colData(my_sce)$FIGO <- metadata$FIGO[x]
  colData(my_sce)$Grade <- metadata$Grade[x]
  colData(my_sce)$Subtype <- metadata$Subtype[x]
  colData(my_sce)$Age <- metadata$Age[x]
  colData(my_sce)$RBC_Lysis <- metadata$RBC_Lysis[x]
  colData(my_sce)$DCR <- metadata$DCR[x]
  colData(my_sce)$Viability <- metadata$Viability[x]
  colData(my_sce)$Clump_level <- metadata$Clump_level[x]
  # give all cell barcodes a sample specific prefix
  colnames(my_sce) <- paste0(x, "-", colnames(my_sce))
  colData(my_sce)$barcodes <- paste0(x, "-", colData(my_sce)$barcodes)
  stopifnot(colnames(my_sce) == colData(my_sce)$barcodes)
  stopifnot(rownames(colData(my_sce)) == colData(my_sce)$barcodes)
  my_sce
})

# read in colour config
cat("\nRead in Color-Config\n")
cat("\n", opt$colour_config, "\n")
config <- read.csv(opt$colour_config, sep = "\t", stringsAsFactors = FALSE)
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(config$colour == unique(config$colour))
ct.color <- c(config$colour, "grey50", "black")
names(ct.color) <- c(config$cell_type, "uncertain", "unknown")
print(ct.color)

# get list of all cell description columns (colData) in all samples to check if
# GEX cell description data is identical between samples
print("Checking colData")
# get character vector with all colData headers of all samples
all_colData_headers <- character()
for (sample in seq_along(all_sce_objects)) {
  all_colData_headers <- c(all_colData_headers, names(colData(all_sce_objects[[sample]])))
}

# get frequency of each colData header
table_headers <- table(all_colData_headers)
# identify columns that are not represented in all samples
header_not_all_samples <- table_headers != length(all_sce_objects)
headers_remove <- table_headers[header_not_all_samples]
cat("\n\n\nWARNING: removing colData columns that are not found in all samples:\n\n")
print(headers_remove)
names_headers_remove <- names(headers_remove)

# remove those columns from all samples
for (sample in seq_along(all_sce_objects)) {
  for (header in names_headers_remove) {
    print(header)
    print(str(header))
    colData(all_sce_objects[[sample]])[header] <- NULL
  }
}


# check if now all colData columns are in all samples
headers_after_filter <- character()
for (sample in seq_along(all_sce_objects)) {
  print(sample)
  headers_after_filter <- c(headers_after_filter, names(colData(all_sce_objects[[sample]])))
}
# get frequency of each colData header
table_headers_after <- table(headers_after_filter)
# identify columns that are not represented in all samples
header_not_all_after <- table_headers_after != length(all_sce_objects)
cat("\n\nCheck if the cell description columns in all samples are identical.\n")
stopifnot(sum(header_not_all_after) == 0)
unique_cell_desc <- unique(headers_after_filter)

###   combine all SCE objects   ###
combined <- scMerge::sce_cbind(all_sce_objects,
                      exprs = c("counts", "normcounts", "pearson_resid"),
                      method = "union",
                      cut_off_batch = 0.0,
                      cut_off_overall = 0.0,
                      colData_names = unique_cell_desc)
# find duplicated barcodes
barcodes_no_prefix <- gsub(pattern = "[^-]+-(.*)", "\\1", colnames(combined))
indices_duplicated <- which(duplicated(barcodes_no_prefix))
barcodes_duplicated <- barcodes_no_prefix[indices_duplicated]
indices_ALL_duplicated <- which(barcodes_no_prefix %in% barcodes_duplicated)
mask_ALL_duplicated <- barcodes_no_prefix %in% barcodes_duplicated
cat("\n\nDuplicated barcodes:\n")
print(table(mask_ALL_duplicated))
cat("\n")

# somehow R turns hyphens in column headers into dots (when using as.data.frame())
# Therefore all hyphens are replaced by underscores.
names_without_hyphen <- gsub("-", "_", names(colData(combined)))
colData_df <- as.data.frame(colData(combined))
names(colData_df) <- names_without_hyphen
# add column with duplicated barcodes
colData_df$duplicated_barcode <- mask_ALL_duplicated

# get union of samples' most variable genes in the GEX pipeline
# only those will be used for calculation of UMAP

var_genes_files <- metadata$h5_files

var_genes_list <- lapply(seq_along(var_genes_files), function(x) {
  cat("\n\nRead in variable genes sample ", x, "\n")
  cat(var_genes_files[x], "\n")
  var_genes_temp <- rhdf5::h5read(var_genes_files[x], "gene_attrs/gene_names")
  var_genes_temp
})
var_genes <- unique(unlist(var_genes_list))
# have matrix with GEX pearson residuals of all variable genes
residuals_variableGenes <- assay(combined, "pearson_resid")[rownames(combined) %in% var_genes, ]
# rows must contain observations (cells)
residuals_variableGenes <- t(residuals_variableGenes)

###########################
###   Calculate UMAPS   ###
###########################

set.seed(3792)
cat("\nGet GEX based UMAP embedding.\n")
# use union of all highly variable genes to calculate UMAP
umap_gex <- uwot::umap(residuals_variableGenes,
                       n_neighbors = 30, pca = 50, spread = 1, min_dist = 0.3, ret_nn = T)
cat("\n\n umap output:\n")
print(str(umap_gex))
# reformat UMAP coordinates
umap_gex_df <- as.data.frame(umap_gex$embedding)
names(umap_gex_df) <- c("GEX_umap1", "GEX_umap2")
umap_gex_df$barcodes <- rownames(residuals_variableGenes)
# make sure all cell dimensions are the same
stopifnot(combined$barcodes == umap_gex_df$barcodes)
# add UMAP coordinates to SCE object
reducedDim(combined, "umap_gex") <- umap_gex_df
# add umap coordinates to the meta data table of cells
colData_df <- plyr::join(colData_df, umap_gex_df, by = "barcodes")

#######################################
###        Composition plots        ###
#######################################

# match colours to cell type levels
id.final.ct <- match(levels(as.factor(colData_df$celltype_final)), names(ct.color))
id.major.ct <- match(levels(as.factor(colData_df$celltype_major)), names(ct.color))
celltype_list <- list(id.final.ct, id.major.ct)
names(celltype_list) <- c("celltype_final", "celltype_major")

# plot cell related data to characterize the cohort
cat("\n\nGenerate plots that characterise cohort.\n\n")
umap1 <- "GEX_umap1"
umap2 <- "GEX_umap2"

# plot sample IDs
cat("# plot Site\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Site)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Site.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Region\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Tissue_Region)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Tissue_Region.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# Exp\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = as.character(Exp))) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10, title="Experiment")) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_Experiment.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

print("Finished plotting.")
