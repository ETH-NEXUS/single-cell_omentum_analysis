####
# Script to plot a heatmap with cluster to gene average counts, tcells only!
###

library(ggplot2)
library(dynamicTreeCut)
library(RColorBrewer)
library(pheatmap)

# 8 main cohort tcell clusters: 3, 11, 12, 17, 19, 20, 22, 23
# read in cluster means
cluster_avg <- as.data.frame(read.csv("gfb_omentum_2020/git/gfb_omentum_2020/omentum_required/2021-08-celltype-subtypes/OmentumCohort.mean_expression_per_umap_cluster_labeled.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))

# read in selected genes
marker_genes <- as.data.frame(read.table("gfb_omentum_2020/git/gfb_omentum_2020/omentum_required/2021-08-celltype-subtypes/tcell_marker.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE))

head(cluster_avg)
head(marker_genes)

length(marker_genes$Gene)
length(unique(marker_genes$Gene))

#check which genes are missing
subset(marker_genes, !(Gene %in% cluster_avg$gene))

subset_genes = subset(cluster_avg, gene %in% marker_genes$Gene)
length(subset_genes$gene) #134
rownames(subset_genes) = subset_genes$gene

names(subset_genes)

keepCol = c("T.NK_3","T.NK_11","T.NK_12","T.NK_17","T.NK_19","T.NK_20","T.NK_22","T.NK_23")
subset_genes_df =  subset_genes[,(names(subset_genes) %in% keepCol)]
head(subset_genes_df)

#normal
hm =  pheatmap(t(subset_genes_df), scale="none", cluster_rows = T, clustering_method = "ward.D2",
                   show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50) 

filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/figure2_tcells/", "Tcells_MainCohort_Clusters", "_heatmap_figure2_tcell_marker.png")
ggsave(filename = filename, hm, width = 60, height = 20, dpi = 300, units = "cm")

#log
subset_genes_transform = subset_genes_df + 1
subset_genes_log = log10(subset_genes_transform)

hm_log =  pheatmap(t(subset_genes_log), scale="none", cluster_rows = T, clustering_method = "ward.D2",
                   show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50) 

filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/figure2_tcells/", "Tcells_MainCohort_Clusters", "_heatmap_figure2_tcell_marker_log.png")
ggsave(filename = filename, hm_log, width = 60, height = 20, dpi = 300, units = "cm")

####
####check all-zero genes
####
subset_genes_df_noZero = subset(subset_genes_df, rowSums(subset_genes_df) > 0)
length(rownames(subset_genes_df_noZero)) #129

######
###second heatmap on the subtype-derived clusters (tcell subtype subsetted)
######

cluster_avg_subset <- as.data.frame(read.csv("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/subtypes_hvg500/tcells/T_NK_cells.mean_expression_per_subtype_umap_cluster.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
head(cluster_avg_subset)

#check which genes are missing
subset(marker_genes, !(Gene %in% cluster_avg_subset$gene))

subset_genes_tcellsub = subset(cluster_avg_subset, gene %in% marker_genes$Gene)
length(subset_genes_tcellsub$gene)
rownames(subset_genes_tcellsub) = subset_genes_tcellsub$gene

names(subset_genes_tcellsub)

dropCol = c("gene")
subset_genes_tcellsub_df =  subset_genes_tcellsub[,!(names(subset_genes_tcellsub) %in% dropCol)]

head(subset_genes_tcellsub_df)

hm_tcellsub =  pheatmap(t(subset_genes_tcellsub_df), scale="none", cluster_rows = T, clustering_method = "ward.D2",
                            show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)
filename_subset <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/figure2_tcells/", "Tcells_SubCohort_Clusters", "_heatmap_figure2_tcell_marker.png")
ggsave(filename = filename_subset, hm_tcellsub, width = 60, height = 20, dpi = 300, units = "cm")

subset_genes_tcellsub_transform = subset_genes_tcellsub_df + 1
subset_genes_tcellsub_log = log10(subset_genes_tcellsub_transform)

hm_tcellsub_log =  pheatmap(t(subset_genes_tcellsub_log), scale="none", cluster_rows = T, clustering_method = "ward.D2",
                   show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)
filename_subset <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/figure2_tcells/", "Tcells_SubCohort_Clusters", "_heatmap_figure2_tcell_marker_log.png")
ggsave(filename = filename_subset, hm_tcellsub_log, width = 60, height = 20, dpi = 300, units = "cm")

####
####check all-zero genes
####
subset_genes_tcellsub_df_noZero = subset(subset_genes_tcellsub_df, rowSums(subset_genes_tcellsub_df) > 0)
length(rownames(subset_genes_tcellsub_df_noZero)) #118

