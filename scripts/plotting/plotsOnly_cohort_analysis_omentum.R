###################################################
## File name: plotOnly_cohort_analysis_omentum.R 
## Adapted from: cohort_analysis_omentum.R
## Author: Anne Bertolini
## Date created: September 2020
## R Version: 3.5.1
###################################################

## GENERAL:
## This script reads in an existing SCE object and creates plots

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
  make_option("--rdsfile", type = "character", help = "RDS file that contains the sce object of the combined cohort samples."),
  make_option("--colour_config", type = "character", help = "Colour config for visualising the cell types. Sometimes union of different cell type colour config files is required."),
  make_option("--sampleName", type = "character", help = "Sample name that will be prefix of all plot names."),
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

combined = readRDS(opt$rdsfile)
colData_df <- as.data.frame(colData(combined))

#get embedding coordinates

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

# subset cohort type and plot cell types

subset_ <- combined[, combined$umap_cohort_cl == cluster_ID]


# plot sample IDs
cat("# plot sample IDs\n")
#ggplot(colData_df, aes(x = GEX_umap1, y = GEX_umap2)) +
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = sampleID)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 4, shape = 15), nrow = 18)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_samples.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot cohort type\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Cohort_Type)) +
#MiPr: was sollen die Doppelklammern? Besser: siehe vorherigen plot
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Cohort_Type.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Patient\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Patient)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Patient.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Type\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Tissue_Type)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Tissue_Type.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Region\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Tissue_Region)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Tissue_Region.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Diagnosis\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Diagnosis)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Diagnosis.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot UMAP clustering\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = factor(umap_cohort_cl))) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Umap_clustering.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

# NOTE: Not for all samples available
cat("# plot Primary_Relapse\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Primary_Relapse)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Primary_Relapse.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

# NOTE: Not for all samples available
cat("# FIGO\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = FIGO)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_FIGO.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

# plot file names (sanity check)
cat("# plot file names\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = file_name)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_fileNames.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# RBC_Lysis\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = RBC_Lysis)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_RBC_Lysis.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# Date_of_Exp\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Date_of_Exp)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_Date_of_Exp.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# Exp\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = Exp)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_Exp.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# DCR\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = DCR)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_DCR.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

# plot cell cycle phase
cat("# plot cell cycle phase\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = cycle_phase)) +
  geom_point(size = 0.2) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_cellCycle.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

# plot first cell types
cat("# plot first cell type\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = celltype_major)) +
scale_color_manual(name = "Cell type",
                     values = ct.color[id.major.ct],
                     drop = F) +
  geom_point(size = 0.2) +
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

# plot second cell types
cat("# plot second cell type\n")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = celltype_final)) +
  scale_color_manual(name = "Cell type",
                     values = ct.color[id.final.ct],
                     drop = F) +
  geom_point(size = 0.2) +
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

# set several continuous values to be plotted
cat("\n# set several continuous values to be plotted\n")
values_to_plot <- c("fractionMT", "log_umi", "n_gene")
print(values_to_plot)
# plot various continuous values
for (continuousValue in values_to_plot) {
  cat("\nplot", continuousValue, "\n")
  ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]], color = colData_df[[continuousValue]])) +
    geom_point(size = 0.2) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_light(base_size = 15) +
    theme(aspect.ratio = UMAP_aspect_ratio) +
    scale_colour_viridis_c(name = continuousValue)
  filename <- paste0(path_qc_plots, "_", continuousValue, ".png")
  ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")
}

# plot duplicated barcodes
# colours for plotting duplicates
palet <- ifelse(colData_df$duplicated_barcode, "red", "grey")
ggplot(colData_df, aes(x = colData_df[[umap1]], y = colData_df[[umap2]])) +
  geom_point(size = 0.2, color = "grey") +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  theme(aspect.ratio = UMAP_aspect_ratio) +
  geom_point(data = colData_df[which(colData_df$duplicated_barcode), ],
             aes(x = colData_df[which(colData_df$duplicated_barcode), umap1], y = colData_df[which(colData_df$duplicated_barcode), umap2]),
             size = 2, colour = "red") +
  geom_text(data = colData_df[which(colData_df$duplicated_barcode), ],
            aes(x = colData_df[which(colData_df$duplicated_barcode), umap1],
                y = colData_df[which(colData_df$duplicated_barcode), umap2],
                label = barcodes),
            size = 2) +
  theme(legend.position = "none") +
  ggtitle("Duplicated barcodes")
filename <- paste0(opt$outdir, opt$sampleName, "_duplicatedBarcodes.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

#######################################
###            GSVA plots           ###
#######################################
cat("\n\nGSVA plots\n")
print(Sys.time())


mask_GSVA_cols <- grepl("GSVA_", names(colData_df))
gsva_counts <- subset(colData_df, select = names(colData_df)[mask_GSVA_cols])
gsva_counts$barcodes <- colData_df$barcodes
# remove the "GSVA_" tag
colnames(gsva_counts) <- gsub(pattern = "[^_]+_(.*)", "\\1", colnames(gsva_counts))
gsva_melted <- reshape2::melt(gsva_counts, id.vars = "barcodes")

# read in pathway list
pathway.list.all <- read.table(opt$pathway_dictionary, head = T, sep = "\t", stringsAsFactors = F)
# keep only pathways that are in matrix
pathway.list <- pathway.list.all[pathway.list.all$Pathway %in% colnames(gsva_counts), ]
# have list per group
pathway.list <- tapply(pathway.list$Pathway, pathway.list$Group, c)
pathway.list.all <- tapply(pathway.list.all$Pathway, pathway.list.all$Group, c)
# sort gene names alphabetically for better overview in plot
pathway.list <- lapply(pathway.list, sort)
pathway.list.all <- lapply(pathway.list.all, sort)

# Trim outlier values and have one scale for all
gsva_melted$value_trimmed <- gsva_melted$value
limit_min <- quantile(gsva_melted$value, prob = 0.01)
gsva_melted$value_trimmed[gsva_melted$value <= limit_min] <- limit_min
print("Number of data points set to limit_min:")
print(table(gsva_melted$value_trimmed == limit_min))
limit_max <- quantile(gsva_melted$value, prob = 0.99)
gsva_melted$value_trimmed[gsva_melted$value >= limit_max] <- limit_max
print("Number of data points set to limit_max:")
print(table(gsva_melted$value_trimmed == limit_max))
print("summary(gsva_melted$value_trimmed) after trimmming the outliers:")
print(summary(gsva_melted$value_trimmed))
# get breaks
gsva_melted$value_trimmed <- round(gsva_melted$value_trimmed, digits = 2)
break_low <- min(gsva_melted$value_trimmed)
break_high <- max(gsva_melted$value_trimmed)
label_low <- paste("<=", break_low)
print(label_low)
label_high <- paste(">=", break_high)
print(label_high)
print("summary(gsva_melted$value_trimmed) after rounding the values to two digits:")
print(summary(gsva_melted$value_trimmed))

# get chunks of gene sets
all_genesets <- colnames(gsva_counts)
chunks_genesets <- split(all_genesets, ceiling(seq_along(all_genesets) / chunksize_genesets))
# get umap coordinates of the integrated cells
umap_coordinates <- subset(colData_df, select = c("GEX_umap1", "GEX_umap2"))
umap_coordinates$barcodes <- colData_df$barcodes

# set plotting theme for gsva plots
theme_set(theme_grey())
# generate plot per chunk of gene sets
for (chunk in seq_along(pathway.list.all)) {
  cat("\n\nChunk: ", chunk, "\n")
  cat("Plot chunk of gene sets:\n")
  print(pathway.list.all[[chunk]])
  sub_gsva_melted <- gsva_melted[gsva_melted$variable %in% pathway.list.all[[chunk]], ]
  gsva_data_plotting <- plyr::join(sub_gsva_melted, umap_coordinates, type = "inner")
  cat("Check if the cell barcodes in both melted gsva objects are identical.\n")
  stopifnot(gsva_data_plotting$barcodes == sub_gsva_melted$barcodes)
  # plot
  p_genesets <- ggplot(gsva_data_plotting, aes(x = GEX_umap1, y = GEX_umap2)) +
    geom_point(aes(color = value_trimmed), size = 1) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme(panel.background = element_rect(fill = "grey80"),
          strip.text.x = element_text(size = 10.5),
          aspect.ratio = UMAP_aspect_ratio) +
    scale_color_distiller(name = "", palette = "RdBu",
                          breaks = c(break_low, 0, break_high),
                          labels = c(label_low, "0.00", label_high)) +
    facet_wrap(~variable, ncol = columns_genesets)
  p_genesets
  filename_gsva_chunk <- paste0(path_pathways, ".integrated_gsva_umap__", names(pathway.list.all)[chunk], ".png")
  ggsave(filename = filename_gsva_chunk,
         width = 40, height = 7.5 * ceiling(length(pathway.list.all[[chunk]]) / columns_genesets),
         dpi = 300, units = "cm")
}


#######################################
###      Gene expression plots      ###
#######################################

# have matrix with normcounts
# remove all genes that have 0 counts in all cells
normcounts_all.zero.removed <- normcounts(combined)
mask_all_zero <- apply(normcounts_all.zero.removed, 1, sum) > 0
normcounts_all.zero.removed <- normcounts_all.zero.removed[mask_all_zero, ]
stopifnot(length(normcounts_all.zero.removed[, 1]) == sum(apply(normcounts_all.zero.removed, 1, sum) > 0))

# plot expression of marker genes in UMAP
# read in gene list
gene.list.all <- read.table(opt$selected_genes, head = T, sep = "\t", stringsAsFactors = F)
# read protein/alias names of genes for alternate plot titles
if (opt$toggle_label) {
  # use regex to parse gene aliases in "notes" col of list file w/ format: ";ALIAS:alias_name;"
  gene.list.protein <- tapply(gene.list.all$notes, gene.list.all$SYMBOL, function(x) regmatches(x, regexec(";ALIAS:(.*?);", x))[[1]][2])
  #function for converting 'gene name' -> 'gene name/protein name'
  convert_gene_name <- function(gene_list) {
    gene_list <- unique(gene_list)
    out_list <- c()
    for (gene_symbol in gene_list) {
      gene_alias <- gene.list.protein[[gene_symbol]][1]
      if (!is.na(gene_alias) && gene_alias != gene_symbol) {
        final_name <- paste(gene_symbol, "/", gene_alias, sep = "")
      } else {
        final_name <- gene_symbol
      }
      out_list <- c(out_list, final_name)
    }
    return(out_list)
  }
}

# keep only genes that are in matrix
gene.list <- gene.list.all[gene.list.all$SYMBOL %in% rownames(normcounts_all.zero.removed), ]
# have list per group
gene.list <- tapply(gene.list$SYMBOL, gene.list$group, c)
gene.list.all <- tapply(gene.list.all$SYMBOL, gene.list.all$group, c)
# sort gene names alphabetically for better overview in plot
gene.list <- lapply(gene.list, sort)
gene.list.all <- lapply(gene.list.all, sort)
#cat("\nstr(gene.list) for expression plotting:\n")
#print(str(gene.list))
# get indices of genes in matrix
list_groups <- lapply(seq_along(gene.list), function(x) {
  indices <- match(gene.list[[x]], rownames(normcounts_all.zero.removed))
  indices
})
names(list_groups) <- names(gene.list)
# sort by position in matrix, do sanity check with old list
list_groups_sorted <- lapply(list_groups, sort)
list_groups_orig <- limma::ids2indices(gene.list, rownames(normcounts_all.zero.removed), remove.empty = F)
print("str(list_groups_sorted):")
print(str(list_groups_sorted))
#stopifnot(all.equal(list_groups_sorted, list_groups_orig))

#####   get cell type coloured UMAP and legend as a reference   #####
fontsize <- theme(axis.text = element_text(size = 3), axis.title = element_text(size = 3))
theme_set(theme_bw(4) + fontsize)
# plot UMAP (based on highly variable genes), colours = final cell type, as reference for expression plots
p_final_ref <- ggplot(colData_df, aes(x = GEX_umap1, y = GEX_umap2, color = celltype_final)) +
  scale_color_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  geom_point(size = 0.0001) +
  theme(aspect.ratio = 1) +
  #coord_fixed(ratio = 1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "none")

p_legend_ref <- ggplot(colData_df, aes(x = GEX_umap1, y = GEX_umap2, color = celltype_final)) +
  scale_color_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  geom_point(size = 1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 0.8, shape = 15), nrow = 17)) +
  theme(legend.key.width = unit(0.05, "line"),
        legend.key.height = unit(0.001, "line"),
        #legend.text = element_text(size = 3.5),
        legend.spacing.y = unit(0.001, "line"),
        legend.spacing.x = unit(0.25, "line"))
legend_ref <- cowplot::get_legend(p_legend_ref)

sample_ref <- ggplot(colData_df, aes(x = GEX_umap1, y = GEX_umap2, color = Tissue_Type)) +
  #scale_color_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F)+
  geom_point(size = 0.0001) +
  theme(aspect.ratio = 1) +
  #coord_fixed(ratio = 1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "none")

p_sample_legend <- ggplot(colData_df, aes(x = GEX_umap1, y = GEX_umap2, color = Tissue_Type)) +
  geom_point(size = 1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 0.8, shape = 15), nrow = 17)) +
  theme(legend.key.width = unit(0.05, "line"),
        legend.key.height = unit(0.001, "line"),
        #legend.text = element_text(size = 3.5),
        legend.spacing.y = unit(0.001, "line"),
        legend.spacing.x = unit(0.25, "line"))
sample_legend_ref <- cowplot::get_legend(p_sample_legend)

ref_empty <- ggplot() + theme_void()


#####  loop over groups of marker genes, plot expression on UMAP with one output file per group   #####
for (group in seq(length(gene.list.all))) {
  group_name <- names(gene.list.all[group])
  cat("Plot group:\n")
  cat(group_name, "\n")
  # get all genes that have zero counts in all cells
  mask.not.expressed <- !(gene.list.all[[group]] %in% rownames(normcounts_all.zero.removed[list_groups[[group_name]], , drop = F]))
  genes.not.expressed <- sort(gene.list.all[[group]][mask.not.expressed])
  if (opt$toggle_label) {
    genes.not.expressed <- convert_gene_name(genes.not.expressed)
  }
  # split number of genes per line based on # of characters of total line
  if (nchar(paste(genes.not.expressed[1:5], collapse = "")) > 32) {
    split_by <- 3
  } else {
    split_by <- 5
  }
  if (length(genes.not.expressed) > 0) {
    genes.not.expressed.split <- split(genes.not.expressed, ceiling(seq_along(genes.not.expressed) / split_by))
  } else {
    genes.not.expressed.split <- NULL
  }
  genes.not.expressed.split.collapsed <- lapply(seq_along(genes.not.expressed.split), function(x) {
    paste(genes.not.expressed.split[[x]], collapse = ", ")
  })
  # generate plot that consists of text and lists all genes that were not found in this sample
  if (length(genes.not.expressed.split.collapsed) > 0) {
    text1 <- paste(genes.not.expressed.split.collapsed, collapse = ",\n")
  } else {
    text1 <- "-"
  }
  text2 <- paste0("For the following genes \n no expression was detected in this sample:\n", text1)
  plot_missed <- ggplot() +
    annotate("text", x = 1, y = 1, label = text2, size = 1.2) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          aspect.ratio = 1)
  # check if there are genes expressed and a submatrix can be generated
  if (length(list_groups[[group_name]]) > 0) {
    # get matrix with all cells and all genes of the respective group
    matrix_group <- normcounts_all.zero.removed[list_groups[[group_name]], , drop = F]
  } else {
    matrix_group <- matrix(0, nrow = 1, ncol = ncol(normcounts(combined)))  # NOTE: was my_sce
  }
  # get empty list that will be filled with expression plots per gene
  plot_gene <- list()
  # check first if any of the groups' genes are expressed in this sample
  if (length(list_groups[[group_name]]) > 0) {
    # loop over the expressed genes of the current group
    for (gene in seq(length(list_groups[[group_name]]))) {
      # empty list that will be filled with 1) expression plot and 2) violin plot of one gene
      merged_plots <- list()
      # combine cell_attributes with expression of current gene
      cell_attributes_gene <- colData_df
      cell_attributes_gene$normcounts <- matrix_group[gene, ]
      # set all values below 0 to 0 for the plotting
      cell_attributes_gene$normcounts[cell_attributes_gene$normcounts < 0] <- 0
      #print("table(cell_attributes_gene$normcounts < 0):")
      #print(table(cell_attributes_gene$normcounts < 0))
      # maximum count found for current gene
      max_count <- max(cell_attributes_gene$normcounts)
      # upper limit of gene expression colour scale is either maximum count or 3
      upper_limit <- max(3, max_count)
      # match gene symbol w/ conventional/protein name for plot labels
      gene_symbol <- rownames(matrix_group)[gene]
      if (opt$toggle_label) {
        legend <- convert_gene_name(gene_symbol)
      } else {
        legend <- gene_symbol
      }
      # plot gene expression
      merged_plots[["expr"]] <- ggplot(cell_attributes_gene, aes(x = GEX_umap1, y = GEX_umap2)) +
        geom_point(aes(color = normcounts), size = rel(0.001)) +
        xlab("UMAP 1") + ylab("UMAP 2") +
        #        scale_color_distiller(name="", palette = "Spectral", direction = -1,
        #        scale_color_viridis(option = "plasma", name = "",
        scale_color_gradientn(name = "", colours = c("slateblue4", "royalblue1",
                                                     "aquamarine3", "khaki", 383,
                                                     "sienna1", "orangered4"),
                              limits = c(0, max(3, upper_limit)),
                              breaks = c(floor(upper_limit / 3), round(2 * (upper_limit / 3)), upper_limit)) +
        coord_fixed(ratio = 1) +
        ggtitle(legend) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 4),
              plot.margin = margin(1, 1, 1, 1, "pt"),
              plot.background = element_rect(fill = "white",colour = "black", size = 0.2),
              #legend.text = element_text(size = 8),
              legend.key.size = unit(0.3, "line"),
              legend.text = element_text(size = rel(0.7)),
              legend.margin = margin(0.2, 0.2, 0.2, 0.2),
              axis.text = element_text(size = rel(0.6)),
              aspect.ratio = 1)
      merged_plots[["expr"]]
      # plot violin plot
      merged_plots[["violin"]] <- ggplot(cell_attributes_gene, aes(x = factor(umap_cohort_cl), y = normcounts)) +
        xlab("Umap cluster ID") +
        scale_fill_brewer(name = "", palette = "Set2") + scale_color_brewer(name = "", palette = "Set2") +
        geom_violin(scale = "count", alpha = 0.2, size = 0.2) +
        ggbeeswarm::geom_quasirandom(alpha = 0.2, groupOnX = T, size = 0.05) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 3),
              plot.margin = margin(3, 5, 6, 1),
              plot.background = element_rect(
                fill = "white",
                colour = "white",
                size = 1),
              axis.text = element_text(size = 3),
              axis.title = element_text(size = 3.5))
      # combine expression plot and violin plot of current gene to one plot
      plot_gene[[gene]] <- cowplot::plot_grid(plotlist = merged_plots, ncol = 1, rel_heights = c(2.1, .9))
      # add combined plot to list that aggregates plots of all genes of the current group
      plot_gene[[gene]]
    }
  }
  names(plot_gene) <- rownames(matrix_group)
  plot_gene <- c(list(plot_missed, p_final_ref, legend_ref, sample_ref, sample_legend_ref, ref_empty), plot_gene)
  nr.cols <- min(6, 2 * (length(list_groups[[group_name]]) + 6))
  nr.rows <- ceiling((length(list_groups[[group_name]]) + 6) / nr.cols)
  plots_group <- plot_grid(plotlist = plot_gene, ncol = nr.cols)
  # save png for each group
  ggsave(path_gene_expr %&% "." %&% group_name %&% "_gene_expression.png",
         plots_group, dpi = 600, width = 20, height = 4 * nr.rows, units = "cm")
  cat("Gene expression plot of group -", group_name, "- plotted.\n\n")
}

print("Finished plotting.")
# get mean expression per umap cluster

#all_umap_clusters <- levels(combined[[]]$umap_cohort_cl)
#print("umap clusters:")
#print(all_umap_clusters)
#all_cluster_means <- as.data.frame(rownames(combined))
#names(all_cluster_means) <- "gene"
#all_cluster_means$gene <- as.character(all_cluster_means$gene)
# Use relative counts normalisation before calculating mean expression per cluster
# get mean expression for each seurat cluster
#for (cluster_ID in all_umap_clusters) {
#  print(paste0("Calculating mean expression for cluster ", cluster_ID))
#  col_header <- paste0("clust_", cluster_ID)
#  print(col_header)
#  subset_cluster <- subset(x = combined, subset = umap_cohort_cl == cluster_ID)
#  print(paste0("Number of cells in cluster: ", dim(subset_cluster)[2]))
#  stopifnot(all_cluster_means$gene == rownames(subset_cluster))
#  all_cluster_means[col_header] <- apply(X = GetAssayData(subset_cluster, assay = "SCT", slot = "data"), 1, mean)
#  all_cluster_means[col_header] <- round(all_cluster_means[col_header], digits = 4)
#}
# write mean expression table
#filename_mean_expr <- paste0(opt$outdir, opt$sampleName, ".mean_expression_per_umap_cluster.tsv")
#print(filename_mean_expr)
#write.table(all_cluster_means, file = filename_mean_expr, quote = FALSE, sep = "\t", row.names = FALSE)
