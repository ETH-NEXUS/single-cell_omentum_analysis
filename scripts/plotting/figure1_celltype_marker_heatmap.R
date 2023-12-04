####
# Script to plot a heatmap with cluster to gene average counts
###

library(ggplot2)
library(dynamicTreeCut)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(viridis)
library(som)

# read in cluster means

cluster_avg <- as.data.frame(read.csv("gfb_omentum_2020/git/gfb_omentum_2020/omentum_required/2021-08-celltype-subtypes/OmentumCohort.mean_expression_per_umap_cluster_labeled.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))

#### 
# write out ranked gene list per cluster
####

out_geneexp = "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/ranked_gene_expression/"
head(cluster_avg)
colnames(cluster_avg)
for (clusterID in colnames(cluster_avg)){
  if (clusterID != "gene"){
    print(clusterID)
    rank_clust = cluster_avg
    ranked_df = rank_clust[order(rank_clust[[clusterID]], decreasing = T), ]
    newTable = data.frame("rank" = seq_along(1:length(ranked_df$gene)), "gene" = ranked_df$gene, "mean_exp" = ranked_df[[clusterID]])
    filename_rank <- paste0(out_geneexp, clusterID, ".rankedGenes.tsv")
    print(filename_rank)
    write.table(newTable, file = filename_rank, quote = FALSE, sep = "\t", row.names = FALSE)
  }
}

rank_clust = cluster_avg
head(rank_clust)
colnames(rank_clust)
dtest = rank_clust[order(rank_clust$Mesenchymal.progenitor_2, decreasing = T), ]    
head(dtest)
newTable = data.frame("rank" = seq_along(1:length(dtest$gene)), "gene" = dtest$gene, "mean_exp" = dtest$Mesenchymal.progenitor_2)
head(newTable)
filename_rank <- paste0(out_geneexp, "Mesenchymal.progenitor_2", ".rankedGenes.tsv")
print(filename_rank)
write.table(newTable, file = filename_rank, quote = FALSE, sep = "\t", row.names = FALSE)

dtest = dtest[seq_len(min(nrow(dtest), max_hvg_length)), ]

# read in selected genes

marker_genes<- as.data.frame(read.table("gfb_omentum_2020/git/gfb_omentum_2020/omentum_required/2021-08-celltype-subtypes/selected_genes_celltype_marker.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE))

head(cluster_avg)
head(marker_genes)
anno_names = colnames(cluster_avg)[2:length(colnames(cluster_avg))]
anno_names
anno_names = str_replace(anno_names, "(r_H)", "r-H")
anno_names = str_replace(anno_names, "([_1234567890])", "")
anno_names = str_replace(anno_names, "([1234567890])", "")
anno_names = str_replace(anno_names, "([1234567890])", "")
anno_names
# colour annotation
annot_colors = c("green4","grey50","magenta","plum1","mediumpurple1", "chocolate1","red2","gold","cyan2","black","darkcyan","lightsalmon","blue2", "darkviolet")

ann_col = data.frame(celltype = sort(unique(anno_names)))
ann_col_2 = data.frame(celltype = factor(anno_names))
ann_col_2
rownames(ann_col_2) = colnames(cluster_avg)[2:length(colnames(cluster_avg))]
levels(ann_col_2$celltype)

anno_cols = list(celltype = annot_colors[seq_along(sort(levels(ann_col_2$celltype)))])

names(anno_cols[[1]]) = sort(levels(ann_col_2$celltype))

length(marker_genes$SYMBOL) # contains THY1 twice!
sort(marker_genes$SYMBOL)
#sort(subset_genes$gene)
subset_genes = subset(cluster_avg, gene %in% marker_genes$SYMBOL)
length(subset_genes$gene)
rownames(subset_genes) = subset_genes$gene

dropCol = c("gene")
subset_genes_df =  subset_genes[,!(names(subset_genes) %in% dropCol)]

subset_genes_df_zscore = normalize(subset_genes_df, byrow=TRUE)
colnames(subset_genes_df_zscore) = colnames(subset_genes_df)

#first heatmap (zscaled, no log, default colours) is the one preferred by Ulli and Francis
hm = pheatmap(subset_genes_df_zscore, scale="none", cluster_rows = T, clustering_method = "ward.D2", #color = viridis(n=8, option = "turbo"), 
              annotation_col = ann_col_2, annotation_colors = anno_cols,
              show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)
filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/figure1_heatmap/", "OmentumCohort", "_heatmap_figure1_ct_marker_2021-11_zScaled_annorow.png")
ggsave(filename = filename, hm, width = 25, height = 30, dpi = 300, units = "cm")

hm = pheatmap(subset_genes_df_zscore, scale="none", cluster_rows = T, clustering_method = "ward.D2", color = viridis(n=4, option = "viridis"), 
              annotation_col = ann_col_2, annotation_colors = anno_cols,
              show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)
filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/figure1_heatmap/", "OmentumCohort", "_heatmap_figure1_ct_marker_2021-11_zScaled_annorow_viridis.png")
ggsave(filename = filename, hm, width = 25, height = 30, dpi = 300, units = "cm")

hm = pheatmap(subset_genes_df_zscore, scale="none", cluster_rows = T, clustering_method = "ward.D2", color = viridis(n=5, option = "viridis"), 
              annotation_col = ann_col_2, annotation_colors = anno_cols,
              show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)
filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/figure1_heatmap/", "OmentumCohort", "_heatmap_figure1_ct_marker_2021-11_zScaled_annorow_viridis5.png")
ggsave(filename = filename, hm, width = 25, height = 30, dpi = 300, units = "cm")

hm = pheatmap(subset_genes_df_zscore, scale="none", cluster_rows = T, clustering_method = "ward.D2", color = viridis(n=6, option = "viridis"), 
              annotation_col = ann_col_2, annotation_colors = anno_cols,
              show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)
filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/figure1_heatmap/", "OmentumCohort", "_heatmap_figure1_ct_marker_2021-11_zScaled_annorow_viridis6.png")
ggsave(filename = filename, hm, width = 25, height = 30, dpi = 300, units = "cm")

hm = pheatmap(subset_genes_df_zscore, scale="none", cluster_rows = T, clustering_method = "ward.D2", color = viridis(n=7, option = "viridis"), 
              annotation_col = ann_col_2, annotation_colors = anno_cols,
              show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)
filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/figure1_heatmap/", "OmentumCohort", "_heatmap_figure1_ct_marker_2021-11_zScaled_annorow_viridis7.png")
ggsave(filename = filename, hm, width = 25, height = 30, dpi = 300, units = "cm")



hm = pheatmap(subset_genes_df, scale="row", cluster_rows = T, clustering_method = "ward.D2", #color = viridis(n=8, option = "turbo"), 
              annotation_col = ann_col_2, annotation_colors = anno_cols,
                   show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)
filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/figure1_heatmap/", "OmentumCohort", "_heatmap_figure1_ct_marker_2021-11_rowScaled_annorow.png")
ggsave(filename = filename, hm, width = 25, height = 30, dpi = 300, units = "cm")

subset_genes_transform = subset_genes_df + 1
subset_genes_log = log10(subset_genes_transform)

hm_log =  pheatmap(subset_genes_log, scale = "row", cluster_rows = T, clustering_method = "ward.D2",
                   annotation_col = ann_col_2, annotation_colors = anno_cols,
                        show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)

filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/figure1_heatmap/", "OmentumCohort", "_heatmap_figure1_ct_marker_log_2021-11_rowScaled_annorow.png")
ggsave(filename = filename, hm_log, width = 25, height = 30, dpi = 300, units = "cm")

# have a second heatmap without VIM, because it is dominating across clusters
subset_noVim <- subset_genes[ !(subset_genes$gene %in% c("VIM")), ]
length(subset_noVim$gene)

subset_noVim_df =  subset_noVim[,!(names(subset_noVim) %in% dropCol)]

subset_noVim_transform = subset_noVim_df + 1
subset_noVim_log = log10(subset_noVim_transform)

hm_noVIM_log =  pheatmap(subset_noVim_log, scale="row", cluster_rows = T, clustering_method = "ward.D2",
                         annotation_col = ann_col_2, annotation_colors = anno_cols,
                   show_rownames = T, show_colnames = T, fontsize=13, treeheight_col = 50)

filename <- paste0("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-11_followUp/figure1_heatmap/", "OmentumCohort", "_heatmap_figure1_ct_marker_noVIM_log_2021-11_rowscaled_annorow.png")
ggsave(filename = filename, hm_noVIM_log, width = 25, height = 30, dpi = 300, units = "cm")
