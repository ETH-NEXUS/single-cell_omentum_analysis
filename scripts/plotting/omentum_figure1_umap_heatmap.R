###################################################
## File name: omentum_figure1_umap_heatmap.R 
## Date created: August 2021
## R Version: 3.5.1
###################################################

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
  make_option("--colour_config", type = "character", help = "Colour config for visualising the cell types. Sometimes union of different cell type colour config files is required."),
  make_option("--selected_genes", type = "character", help = "Table with selected genes that will be shown in gene expression plots."),
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
#path_gene_expr <- opt$outdir %&% opt$sampleName

###################################
###   Read in and format data   ###
###################################

combined = readRDS(opt$rdsFile)

# read in color config
cat("\nRead in Color-Config\n")
cat("\n", opt$colour_config, "\n")
config <- read.csv(opt$colour_config, sep = "\t", stringsAsFactors = FALSE)
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(config$colour == unique(config$colour))
ct.color <- c(config$colour, "grey50", "black")
names(ct.color) <- c(config$cell_type, "uncertain", "unknown")
print(ct.color)

###################################
###   #relabel cell type, also based on cluster ID
###################################

unique(colData(combined)$celltype_final)

colData(combined)$celltype_figure1 <- as.character(colData(combined)$celltype_final)
combined$celltype_figure1[combined$celltype_figure1 == "Fibroblasts"] = "Mesenchymal.progenitor"
combined$celltype_figure1[combined$umap_cohort_cl == "10"] = "Plasmacytoid.dendritic.cells"
combined$celltype_figure1[combined$umap_cohort_cl == "6"] = "Lymphatic.endothelial.cells"
# update in Novemeber 2021, also change cluster 15 myeloid cells to Macrophages and change Bcell subtypes
combined$celltype_figure1[combined$umap_cohort_cl == "15" & combined$celltype_figure1 == "Myeloid.cells"] = "Macrophages"
head(combined$celltype_figure1[combined$umap_cohort_cl == "15" & combined$celltype_figure1 == "Myeloid.cells"])
length(combined$celltype_figure1[combined$umap_cohort_cl == "15" & combined$celltype_figure1 == "Myeloid.cells"])
combined$celltype_figure1[combined$celltype_figure1 == "B.cells.naive"] = "B.cells"
combined$celltype_figure1[combined$celltype_figure1 == "B.cells.memory"] = "B.cells"

unique(colData(combined)$celltype_figure1)

# match colors to cell type levels
id.final.ct <- match(levels(as.factor(combined$celltype_figure1)), names(ct.color))
id.major.ct <- match(levels(as.factor(combined$celltype_major)), names(ct.color))
celltype_list <- list(id.final.ct, id.major.ct)
names(celltype_list) <- c("celltype_final", "celltype_major")

colData_df <- as.data.frame(colData(combined))

# plot final cell types
cat("# plot final cell type labelling \n")
ggplot(colData_df, aes(x = reducedDim(combined)$GEX_umap1, y = reducedDim(combined)$GEX_umap2, color = celltype_figure1)) +
  scale_color_manual(name = "Cell type",
                     values = ct.color[id.final.ct],
                     drop = F) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20)) +
  theme(aspect.ratio = 1,
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.key.width = unit(0.5, "line"),
        legend.key.height = unit(0.5, "line"),
        legend.spacing.y = unit(0.5, "line"),
        legend.spacing.x = unit(0.5, "line")
  )
filename <- paste0(opt$outdir, opt$outName, "_celltypes_curated_labels.png")
ggsave(filename = filename, width = 35, height = 20, dpi = 300, units = "cm")


###################################
###  cell type marker UMAPS
###################################

# have matrix with normcounts
# remove all genes that have 0 counts in all cells
normcounts_all.zero.removed <- normcounts(combined)
mask_all_zero <- apply(normcounts_all.zero.removed, 1, sum) > 0
normcounts_all.zero.removed <- normcounts_all.zero.removed[mask_all_zero, ]
stopifnot(length(normcounts_all.zero.removed[, 1]) == sum(apply(normcounts_all.zero.removed, 1, sum) > 0))

# plot expression of marker genes in UMAP
# read in gene list
gene.list.all <- read.table(opt$selected_genes, head = T, sep = "\t", stringsAsFactors = F)

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

##### reference UMAPs: celltype manually curated and clustering   #####
fontsize <- theme(axis.text = element_text(size = 3), axis.title = element_text(size = 3))
theme_set(theme_bw(4) + fontsize)
# plot UMAP (based on highly variable genes), colours = final cell type, as reference for expression plots
p_final_ref <- ggplot(colData_df, aes(x = reducedDim(combined)$GEX_umap1, y = reducedDim(combined)$GEX_umap2, color = celltype_figure1)) +
  scale_color_manual(name = "Cell type", values = ct.color[id.final.ct], drop = F) +
  geom_point(size = 0.00001) +
  theme(aspect.ratio = 1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "none")

p_legend_ref <- ggplot(colData_df, aes(x = reducedDim(combined)$GEX_umap1, y = reducedDim(combined)$GEX_umap2, color = celltype_figure1)) +
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

sample_ref <- ggplot(colData_df, aes(x = reducedDim(combined)$GEX_umap1, y = reducedDim(combined)$GEX_umap2, color = factor(umap_cohort_cl))) +
  geom_point(size = 0.00001) +
  theme(aspect.ratio = 1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme(legend.position = "none")

p_sample_legend <- ggplot(colData_df, aes(x = reducedDim(combined)$GEX_umap1, y = reducedDim(combined)$GEX_umap2, color = factor(umap_cohort_cl))) +
  geom_point(size = 1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 2, shape = 15), nrow = 10, title="Cluster IDs")) +
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
    annotate("text", x = 1, y = 1, label = text2, size = 1.5) +
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
      legend <- gene_symbol
      # plot gene expression
      merged_plots[["expr"]] <- ggplot(cell_attributes_gene, aes(x = reducedDim(combined)$GEX_umap1, y = reducedDim(combined)$GEX_umap2)) +
        geom_point(aes(color = normcounts), size = rel(0.0001)) +
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
  filename_group <- paste0(opt$outdir, opt$outName, ".", group_name , "_gene_expression.png")
  ggsave(filename = filename_group, plots_group, dpi = 600, width = 20, height = 4 * nr.rows, units = "cm")
  cat("Gene expression plot of group -", group_name, "- plotted.\n\n")
}
