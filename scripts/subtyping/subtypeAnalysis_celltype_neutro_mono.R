###################################################
## File name: cohort_analysis_omentum.R 
## Adapted from: analyse_scRNA_cohort_gexOnly.R
## Author: Anne Bertolini
## Date created: September 2020
## R Version: 3.5.1
###################################################

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
  make_option("--rdsFile", type = "character", help = "RDS file with merged single cell data."),
  make_option("--selected_genes", type = "character", help = "Table with selected genes that will be shown in gene expression plots."),
  make_option("--colour_config", type = "character", help = "Colour config for visualising the cell types. Sometimes union of different cell type colour config files is required."),
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
dir_gene_expr <- opt$outdir  %&% "/gene_expression/"
dir.create(dir_gene_expr, showWarnings = FALSE)
path_gene_expr <- opt$outdir %&% "/gene_expression/" %&% opt$sampleName

###################################
###   Read in and format data   ###
###################################

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

#read in full RDS
subset_combined = readRDS(opt$rdsFile)
colData_subset <- as.data.frame(colData(subset_combined))

# match colours to cell type levels
id.final.ct <- match(levels(as.factor(colData_subset$celltype_final)), names(ct.color))
id.major.ct <- match(levels(as.factor(colData_subset$celltype_major)), names(ct.color))
celltype_list <- list(id.final.ct, id.major.ct)
names(celltype_list) <- c("celltype_final", "celltype_major")

#######################################
###        Composition plots        ###
#######################################

# plot cell related data to characterize the cohort
  
cat("\n\nGenerate plots that characterise cohort.\n\n")
umap1 <- "umap1"
umap2 <- "umap2"

# plot first cell types
cat("# plot first cell type\n")
ggplot(colData_subset, aes(x = reducedDim(subset_combined)$umap1, y = reducedDim(subset_combined)$umap2, color = celltype_major)) +
scale_color_manual(name = "Cell type",
                     values = ct.color[id.major.ct],
                     drop = F) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20)) +
  theme(aspect.ratio = UMAP_aspect_ratio,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.key.width = unit(0.5, "line"),
        legend.key.height = unit(0.5, "line"),
        legend.spacing.y = unit(0.5, "line"),
        legend.spacing.x = unit(0.5, "line")
        )
filename <- paste0(opt$outdir, opt$sampleName, "_celltypes_major.png")
ggsave(filename = filename, width = 35, height = 20, dpi = 300, units = "cm")

# plot final cell types
cat("# plot final cell type\n")
ggplot(colData_subset, aes(x = reducedDim(subset_combined)$umap1, y = reducedDim(subset_combined)$umap2, color = celltype_final)) +
  scale_color_manual(name = "Cell type",
                     values = ct.color[id.final.ct],
                     drop = F) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20)) +
  theme(aspect.ratio = UMAP_aspect_ratio,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.key.width = unit(0.5, "line"),
        legend.key.height = unit(0.5, "line"),
        legend.spacing.y = unit(0.5, "line"),
        legend.spacing.x = unit(0.5, "line")
        )
filename <- paste0(opt$outdir, opt$sampleName, "_celltypes_final.png")
ggsave(filename = filename, width = 35, height = 20, dpi = 300, units = "cm")
