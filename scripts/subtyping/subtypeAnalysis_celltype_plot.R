###################################################
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

path_qc_plots <- opt$outdir %&% "/qc_plots/" %&% opt$sampleName
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

#read in full RDS
subtype_sce = readRDS(opt$rdsFile)
colData_subset <- as.data.frame(colData(subtype_sce))

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
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = celltype_major)) +
scale_color_manual(name = "Cell type",
                     values = ct.color[id.major.ct],
                     drop = F) +
  geom_point(size = 0.3) +
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
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = celltype_final)) +
  scale_color_manual(name = "Cell type",
                     values = ct.color[id.final.ct],
                     drop = F) +
  geom_point(size = 0.3) +
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


### 
### UMAP VISUALIZATIONS
###

# plot sample IDs
cat("# plot sample IDs\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = sampleID)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 15), nrow = 18)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_samples.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot cohort type\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Cohort_Type)) +
  #MiPr: was sollen die Doppelklammern? Besser: siehe vorherigen plot
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Cohort_Type.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot cohort type old\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Cohort_Type_Old)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Cohort_Type_Old.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot cohort type highlevel\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Cohort_Type_Highlevel)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Cohort_Type_Highlevel.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Patient\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Patient)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 17)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Patient.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Type\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Tissue_Type)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Tissue_Type.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Region Old\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Tissue_Region_old)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Tissue_Region_old.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Region\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Tissue_Region)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Tissue_Region.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Diagnosis\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Diagnosis)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Diagnosis.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot UMAP clustering\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = factor(umap_subtype_cl))) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20, title="Cluster IDs")) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Umap_clustering.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

clusterIDs = unique(colData(subtype_sce)$umap_subtype_cl)
clust_length = length(clusterIDs)

df_labels <- data.frame(matrix(ncol = clust_length, nrow = 0))
myNames <- c("center_of_mass_x", "center_of_mass_y", "clusterID")
colnames(df_labels) <- myNames

for (clust in clusterIDs){
  print(clust)
  subsetClust = subtype_sce[, colData(subtype_sce)$umap_subtype_cl == clust]
  xMed = median(reducedDim(subsetClust)$umap1)
  yMed = median(reducedDim(subsetClust)$umap2)
  newRow <- c(xMed,yMed,clust)
  df_labels <- rbind(df_labels,newRow)
}

df_labels
colnames(df_labels) <- myNames
dim(df_labels)

ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = factor(umap_subtype_cl))) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20, title="Cluster IDs")) +
  geom_label(data=df_labels, inherit.aes = FALSE, aes(x = center_of_mass_x, y = center_of_mass_y, label = clusterID)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$sampleName, "_Umap_clustering_labels.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# RBC_Lysis\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = RBC_Lysis)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_RBC_Lysis.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# Date_of_Exp\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Date_of_Exp)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 18)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_Date_of_Exp.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# Exp\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = as.character(Exp))) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10, title="Experiment")) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_Exp.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# DCR\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = DCR)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_DCR.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

# plot cell cycle phase
cat("# plot cell cycle phase\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = cycle_phase)) +
  geom_point(size = 0.3) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_cellCycle.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")


# set several continuous values to be plotted
cat("\n# set several continuous values to be plotted\n")
values_to_plot <- c("fractionMT", "log_umi", "n_gene")
print(values_to_plot)
# plot various continuous values
for (continuousValue in values_to_plot) {
  cat("\nplot", continuousValue, "\n")
  ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = colData_subset[[continuousValue]])) +
    geom_point(size = 0.3) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_light(base_size = 15) +
    theme(aspect.ratio = UMAP_aspect_ratio) +
    scale_colour_viridis_c(name = continuousValue)
  filename <- paste0(path_qc_plots, "_", continuousValue, ".png")
  ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")
}
