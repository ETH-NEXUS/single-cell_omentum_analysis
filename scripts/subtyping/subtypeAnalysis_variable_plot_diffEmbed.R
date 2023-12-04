###################################################
## File name: subtypeAnalysis_subset_and_normalize.R
## Date created: July 2021
## R Version: 4
###################################################

## GENERAL:
## for the single-cell omentum project; subset the cohort into 4 different celltype subsets based on cell type assignment and clustering

library(SingleCellExperiment)
library(ggplot2)
library(scMerge)
library(uwot)
library(patchwork)
library(optparse)
library(reshape2)
library(igraph)
library(ggrepel)
library(scran)
#library(dbscan)
library(dynamicTreeCut)
library(RColorBrewer)
library(pheatmap)

cat("\n\n\nPrint sessionInfo():\n\n")
print(sessionInfo())

# parse command line arguments
option_list <- list(
  make_option("--rdsFile", type = "character", help = "RDS file with merged single cell data."),
  make_option("--outName", type = "character", help = "Sample name that will be prefix of all plot names."),
  make_option("--file_genes_lily", type = "character", help = "File with genes selected by Lily."),
  make_option("--outdir", type = "character", help = "Path to the outdir, in which files will be written.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""
UMAP_aspect_ratio <- "1"

# create sub directories to group plots
dir_gene_expr <- opt$outdir  %&% "/gene_expression/"
dir.create(dir_gene_expr, showWarnings = FALSE)
path_gene_expr <- opt$outdir %&% "/gene_expression/" %&% opt$outName

dir_qc <- opt$outdir %&% "/qc_plots/"
dir.create(dir_qc, showWarnings = FALSE)
path_qc_plots <- opt$outdir %&% "/qc_plots/" %&% opt$outName

###################################
###   Read in and format data   ###
###################################

#mesothelial
#rdsFile = "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/subtype_analysis/subset_normalize/Omentum.mesothelial.subtype.RDS"
#outdir = "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/subtype_analysis/mesothelial/"
#outName = "Mesothelial.cells"
#subtype_sce = readRDS(rdsFile)

#test
#rdsFile_test = "gfb_omentum_2020/analysis/tests/testCohort/TestCohort.mesothelial.subtype.RDS"
#outDir_test = "gfb_omentum_2020/analysis/tests/testCohort/"
#file_genes_lily = "gfb_omentum_2020/data/data_from_lilys_analysis/genes_MMT_Lily_2021-06-30.txt"
#subtype_sce = readRDS(rdsFile_test)

subtype_sce = readRDS(opt$rdsFile)
genes_lily <- read.table(opt$file_genes_lily, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# four different cell type subsets are interesting:
# Fibroblasts; Mesothelial.cells ; T.NK.cells; HGSOC

####### first compute the highly variable genes again

# TODO: does it need the log? Not possible with the normcount assay!
meanvartrend = modelGeneVar(normcounts(subtype_sce))
hvg = getTopHVGs(meanvartrend, fdr.threshold = 0.05)
length(hvg)
head(hvg)

gene_desc_subtype = meanvartrend
gene_desc_subtype$symbol = rownames(meanvartrend)
gene_desc_subtype$hvg = (gene_desc_subtype$symbol %in% hvg)

#show hvgs
ggplot(as.data.frame(gene_desc_subtype), aes(x=mean, y=total)) + 
  geom_point(aes(color=hvg)) + geom_smooth() +
  geom_line(aes(x=mean, y=tech), color="darkgreen", size=2) 
filename <- paste0(opt$outdir, opt$outName, "_hvg_mean_var.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

max_hvg_length = 500 # TEST WITH 500 INSTEAD OF 1000
# pick max number of hvg's
dtest = gene_desc_subtype[hvg, ]
dtest = dtest[order(dtest$bio, decreasing = T), ]    
dtest = dtest[seq_len(min(nrow(dtest), max_hvg_length)), ]
gene_desc_subtype[setdiff(rownames(gene_desc_subtype), rownames(dtest)), ]$hvg = FALSE
hvg = rownames(gene_desc_subtype)[which(gene_desc_subtype$hvg)]

# CALCULATE UMAP EMBEDDING
# CHANGED TO NORMCOUNT ASSAY! ALSO MIN_DIST = 0.01 and neighbors = 50
residuals_variableGenes <- assay(subtype_sce, "normcounts")[rownames(subtype_sce) %in% hvg, ]
residuals_variableGenes <- t(residuals_variableGenes)

set.seed(3792)
cat("\nGet UMAP embedding.\n")
# use union of all highly variable genes to calculate UMAP
umap_gex_subtype <- uwot::umap(residuals_variableGenes,
                       n_neighbors = 50, pca = 50, spread = 1, min_dist = 0.01, ret_nn = T)

umap_gex_subtype_df <- as.data.frame(umap_gex_subtype$embedding)
names(umap_gex_subtype_df) <- c("umap1", "umap2")
umap_gex_subtype_df$barcodes <- rownames(residuals_variableGenes)
# make sure all cell dimensions are the same
stopifnot(subtype_sce$barcodes == umap_gex_subtype_df$barcodes)
# add UMAP coordinates to SCE object
reducedDim(subtype_sce, "subtype_umap") <- umap_gex_subtype_df

# PERFORM CLUSTERING

metadata(subtype_sce) = c(metadata(subtype_sce), list(subtype_umap_dist=umap_gex_subtype$nn$euclidean))
## nearest neighbors
uu.nn = umap_gex_subtype$nn$euclidean$idx
uu.nn[uu.nn==0] = 1
## weights as 1-distance
wgt = 1 - umap_gex_subtype$nn$euclidean$dist/min(max(umap_gex_subtype$nn$euclidean$dist),1e5)
wgt[wgt < 0] = 0

wgt = melt(wgt, varnames = c("i", "var"), value.name = "x")
adj = melt(uu.nn, varnames = c("i", "var"), value.name="j", id.var=1)
adj = merge(wgt, adj, sort=F)
adj = adj[, c("i", "j", "x")]
adj = sparseMatrix(i=adj$i, j=adj$j, x=adj$x)
# convert to weighted graph
g = graph.adjacency(adj, mode="max", weighted = T, diag = F)
km = igraph::cluster_louvain(g)
ph.membership = km$membership
names(ph.membership) = km$names
colData(subtype_sce)$umap_subtype_cl = km$membership
#colData_df$umap_cohort_cl = km$membership
metadata(subtype_sce) = c(metadata(subtype_sce), list(umap_modularity_subtype=modularity(km)))

unique(colData(subtype_sce)$umap_subtype_cl)

# save sce object in RDS file
#saveRDS(subtype_sce, outdir %&% outName %&% ".subtype_sce.RDS")
saveRDS(subtype_sce, opt$outdir %&% opt$outName %&% ".subtype_sce.RDS")

### 
### CLUSTER MEAN EXPRESSION
###

normcounts_all.zero.removed <- normcounts(subtype_sce)
mask_all_zero <- apply(normcounts_all.zero.removed, 1, sum) > 0
normcounts_all.zero.removed <- normcounts_all.zero.removed[mask_all_zero, ]
stopifnot(length(normcounts_all.zero.removed[, 1]) == sum(apply(normcounts_all.zero.removed, 1, sum) > 0))

all_umap_clusters <- unique(colData(subtype_sce)$umap_subtype_cl)
print("umap clusters:")
print(all_umap_clusters)

all_cluster_means <- as.data.frame(rownames(subtype_sce))

names(all_cluster_means) <- "gene"
all_cluster_means$gene <- as.character(all_cluster_means$gene)

# Use relative counts normalisation before calculating mean expression per cluster
# get mean expression for each umap cluster
for (cluster_ID in all_umap_clusters) {
  print(paste0("Calculating mean expression for cluster ", cluster_ID))
  col_header <- paste0("clust_", cluster_ID)
  print(col_header)
  subset_cluster <- subtype_sce[, subtype_sce$umap_subtype_cl == cluster_ID]
  print(paste0("Number of cells in cluster: ", dim(subset_cluster)[2]))
  stopifnot(all_cluster_means$gene == rownames(subset_cluster))
  
  all_cluster_means[col_header] <- apply(assay(subset_cluster, "normcounts"), 1, mean)
  all_cluster_means[col_header] <- round(all_cluster_means[col_header], digits = 4)
  all_cluster_means[col_header][all_cluster_means[col_header] < 0] <- 0
}
# write mean expression table
filename_mean_expr <- paste0(opt$outdir, opt$outName, ".mean_expression_per_subtype_umap_cluster.tsv")
write.table(all_cluster_means, file = filename_mean_expr, quote = FALSE, sep = "\t", row.names = FALSE)

### 
### UMAP VISUALIZATIONS
###

colData_subset <- as.data.frame(colData(subtype_sce))

# plot sample IDs
cat("# plot sample IDs\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = sampleID)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 3, shape = 15), nrow = 18)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_samples.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot cohort type\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Cohort_Type)) +
  #MiPr: was sollen die Doppelklammern? Besser: siehe vorherigen plot
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Cohort_Type.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot cohort type old\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Cohort_Type_Old)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Cohort_Type_Old.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot cohort type highlevel\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Cohort_Type_Highlevel)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Cohort_Type_Highlevel.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Patient\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Patient)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 17)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Patient.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Type\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Tissue_Type)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Tissue_Type.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Region Old\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Tissue_Region_old)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Tissue_Region_old.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Tissue_Region\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Tissue_Region)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Tissue_Region.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot Diagnosis\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Diagnosis)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Diagnosis.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# plot UMAP clustering\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = factor(umap_subtype_cl))) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20, title="Cluster IDs")) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Umap_clustering.png")
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
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 20, title="Cluster IDs")) +
  geom_label(data=df_labels, inherit.aes = FALSE, aes(x = center_of_mass_x, y = center_of_mass_y, label = clusterID)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(opt$outdir, opt$outName, "_Umap_clustering_labels.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# RBC_Lysis\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = RBC_Lysis)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_RBC_Lysis.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# Date_of_Exp\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = Date_of_Exp)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 18)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_Date_of_Exp.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# Exp\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = as.character(Exp))) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10, title="Experiment")) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_Exp.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

cat("# DCR\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = DCR)) +
  geom_point(size = 0.1) +
  xlab("UMAP 1") + ylab("UMAP 2") +
  theme_light(base_size = 15) +
  guides(colour = guide_legend(override.aes = list(size = 5, shape = 15), nrow = 10)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
filename <- paste0(path_qc_plots, "_DCR.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

# plot cell cycle phase
cat("# plot cell cycle phase\n")
ggplot(colData_subset, aes(x = reducedDim(subtype_sce)$umap1, y = reducedDim(subtype_sce)$umap2, color = cycle_phase)) +
  geom_point(size = 0.1) +
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
    geom_point(size = 0.1) +
    xlab("UMAP 1") + ylab("UMAP 2") +
    theme_light(base_size = 15) +
    theme(aspect.ratio = UMAP_aspect_ratio) +
    scale_colour_viridis_c(name = continuousValue)
  filename <- paste0(path_qc_plots, "_", continuousValue, ".png")
  ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")
}


### 
### HEATMAP CLUSTER MEAN EXPRESSION
### GENES FROM LILY
###


head(all_cluster_means)
length(genes_lily$gene)
subset_lilys_genes = subset(all_cluster_means, gene %in% genes_lily$gene)
length(subset_lilys_genes$gene)
rownames(subset_lilys_genes) = subset_lilys_genes$gene

dropCol = c("gene")
subset_lilys_genes_df =  subset_lilys_genes[,!(names(subset_lilys_genes) %in% dropCol)]
hm_lily = pheatmap(subset_lilys_genes_df, scale="none", cluster_rows = T, clustering_method = "ward.D2",
         show_rownames = T, show_colnames = T, fontsize=8, treeheight_col = 50)

filename <- paste0(opt$outdir, opt$outName, "_heatmap_lilys_genes.png")
ggsave(filename = filename, hm_lily, width = 20, height = 30, dpi = 300, units = "cm")

subset_lily_transform = subset_lilys_genes_df + 1
subset_lilys_genes_log = log10(subset_lily_transform)

hm_lily_log =  pheatmap(subset_lilys_genes_log, scale="none", cluster_rows = T, clustering_method = "ward.D2",
         show_rownames = T, show_colnames = T, fontsize=8, treeheight_col = 50)

filename <- paste0(opt$outdir, opt$outName, "_heatmap_lilys_genes_log.png")
ggsave(filename = filename, hm_lily_log, width = 20, height = 30, dpi = 300, units = "cm")


### 
### HEATMAP CLUSTER MEAN EXPRESSION
### TOP 50 and TOP 100 HVG
###

dtest50 = dtest[seq_len(min(nrow(dtest), 50)), ]
subset_hvg50 = subset(all_cluster_means, gene %in% rownames(dtest50))

rownames(subset_hvg50) = subset_hvg50$gene

dropCol = c("gene")
subset_hvg50_df =  subset_hvg50[,!(names(subset_hvg50) %in% dropCol)]
subset_hvg50_transform = subset_hvg50_df + 1
subset_hvg50_log = log10(subset_hvg50_transform)

hm_hvg50 = pheatmap(subset_hvg50_df, scale="none", cluster_rows = T, clustering_method = "ward.D2",
         show_rownames = T, show_colnames = T, fontsize=8, treeheight_col = 50)
filename <- paste0(opt$outdir, opt$outName, "_heatmap_top50_hvg.png")
ggsave(filename = filename, hm_hvg50, width = 20, height = 40, dpi = 300, units = "cm")

hm_hvg50_log = pheatmap(subset_hvg50_log, scale="none", cluster_rows = T, clustering_method = "ward.D2",
         show_rownames = T, show_colnames = T, fontsize=8, treeheight_col = 50)
filename <- paste0(opt$outdir, opt$outName, "_heatmap_top50_hvg_log.png")
ggsave(filename = filename, hm_hvg50_log, width = 20, height = 40, dpi = 300, units = "cm")

dtest100 = dtest[seq_len(min(nrow(dtest), 100)), ]
subset_hvg100 = subset(all_cluster_means, gene %in% rownames(dtest100))

rownames(subset_hvg100) = subset_hvg100$gene

dropCol = c("gene")
subset_hvg100_df =  subset_hvg100[,!(names(subset_hvg100) %in% dropCol)]
subset_hvg100_transform = subset_hvg100_df + 1
subset_hvg100_log = log10(subset_hvg100_transform)

hm_hvg100 = pheatmap(subset_hvg100_df, scale="none", cluster_rows = T, clustering_method = "ward.D2",
         show_rownames = T, show_colnames = T, fontsize=8, treeheight_col = 50)
filename <- paste0(opt$outdir, opt$outName, "_heatmap_top100_hvg.png")
ggsave(filename = filename, hm_hvg100, width = 20, height = 60, dpi = 300, units = "cm")

hm_hvg100_log = pheatmap(subset_hvg100_log, scale="none", cluster_rows = T, clustering_method = "ward.D2",
         show_rownames = T, show_colnames = T, fontsize=8, treeheight_col = 50)
filename <- paste0(opt$outdir, opt$outName, "_heatmap_top100_hvg_log.png")
ggsave(filename = filename, hm_hvg100_log, width = 20, height = 60, dpi = 300, units = "cm")


