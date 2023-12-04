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
library(rhdf5)
library(stringr)
library(patchwork)
library(optparse)
library(limma)
library(cowplot)
library(reshape2)
library(ggbeeswarm)
library(igraph)
library(sctransform)

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
###   Read in and format data   ###
###################################

combined = readRDS(opt$rdsFile)

# first update according to latest clinical metadata

colData(combined)$Tissue_Region_old <- colData(combined)$Tissue_Region
combined$Tissue_Region[combined$Tissue_Region == "Omentum"] = "benign"
combined$Tissue_Region[combined$Tissue_Region == "Distal"] = "tumor_distant"
combined$Tissue_Region[combined$Tissue_Region == "Next"] = "tumor_peritumoral"
combined$Tissue_Region[combined$Tissue_Region == "Tumor"] = "tumor_metastasis"

unique(combined$Tissue_Region)
unique(combined$Tissue_Region_old)

# four different cell type subsets are interesting:
# Fibroblasts; Mesothelial.cells ; T.NK.cells; HGSOC

# define the cluster IDs that should be taken into account for the different cell types

take_meso = c(1, 4, 13, 18, 24)
take_tumor = c(14)
take_fibro = c(2, 8)
take_tcell = c(3, 11, 12, 17, 19, 20, 22, 23)

#first take the desired clusters
meso_sub1 <- combined[, combined$umap_cohort_cl %in% take_meso]
tumor_sub1 <- combined[, combined$umap_cohort_cl %in% take_tumor]
fibro_sub1 <- combined[, combined$umap_cohort_cl %in% take_fibro]
tcell_sub1 <- combined[, combined$umap_cohort_cl %in% take_tcell]

# then take the desired cell type
meso_sub <- meso_sub1[, colData(meso_sub1)$celltype_final == "Mesothelial.cells"]
tumor_sub <- tumor_sub1[, colData(tumor_sub1)$celltype_final == "HGSOC"]
fibro_sub <- fibro_sub1[, colData(fibro_sub1)$celltype_final == "Fibroblasts"]
tcell_sub <- tcell_sub1[, colData(tcell_sub1)$celltype_final == "T.NK.cells"]

#check; must be "TRUE" only
unique(meso_sub$celltype_final == "Mesothelial.cells")
unique(tumor_sub$celltype_final == "HGSOC")
unique(fibro_sub$celltype_final == "Fibroblasts")
unique(tcell_sub$celltype_final == "T.NK.cells")

# take the celltype subset and perform an additional sctransform normalization, on the counts assay

#Fibroblasts
print("Fibroblasts")
genecounts_fibro = rowSums(counts(fibro_sub))
cellcounts_fibro = colSums(counts(fibro_sub))
checkcounts_fibro = apply(counts(fibro_sub), 2, function(x) length(which(x>0)))

fibro_sub_2 = fibro_sub[!is.na(genecounts_fibro) & genecounts_fibro > 100, cellcounts_fibro > 0 & checkcounts_fibro > 0]  # some genes removed

vst_out_fibro = sctransform::vst(counts(fibro_sub_2), cell_attr = colData(fibro_sub_2), method="nb_fast", 
                                latent_var = c('log_umi'), latent_var_nonreg = c("g2m_score", "s_score"), 
                                return_gene_attr = T, return_cell_attr = T)

print("Start performing sctransform::smooth_via_pca: ")
y_smooth_fibro = sctransform::smooth_via_pca(vst_out_fibro$y, do_plot = TRUE)

print("Start performing sctransform::correct: ")
dat_cor_fibro = sctransform::correct(vst_out_fibro, data = y_smooth_fibro,
                                    do_round = TRUE, do_pos = FALSE)

# Plot mean variance plot
mvp_fibro <- ggplot(vst_out_fibro$gene_attr, aes(log10(gmean), log10(residual_variance))) +
  geom_point(alpha=0.3, shape=16) +
  #geom_density_2d(size = 0.3)
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".fibroblasts.mean_vs_variance_plot.png", mvp_fibro, dpi = 300)
mvp2_fibro <- ggplot(vst_out_fibro$gene_attr, aes(log10(gmean), residual_variance)) +
  geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3) +
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".fibroblasts.mean_vs_variance_plot_2.png", mvp2_fibro, dpi = 300)

## make SingleCellExperiment object and save to disk
vst_out_fibro$gene_attr$gene_ids = rownames(vst_out_fibro$gene_attr)
rowData(fibro_sub_2)$gene_ids = rownames(fibro_sub_2)

gene_desc = merge(vst_out_fibro$gene_attr, rowData(fibro_sub_2), sort=F, by = "gene_ids")
idx = match(rownames(dat_cor_fibro), rownames(fibro_sub_2))
idy = match(colnames(dat_cor_fibro), colnames(fibro_sub_2))
sce_fibro = SingleCellExperiment(assays = list(counts=fibro_sub_2[idx, idy],
                                               normcounts=dat_cor_fibro,
                                               pearson_resid=vst_out_fibro$y),
                                 colData = colData(fibro_sub_2),
                                 rowData = gene_desc)

saveRDS(sce_fibro, opt$outdir %&% opt$outName %&% ".fibroblasts.subtype.RDS")


#T.NK.cells
print("T.NK.cells")
genecounts_tcell = rowSums(counts(tcell_sub))
cellcounts_tcell = colSums(counts(tcell_sub))
checkcounts_tcell = apply(counts(tcell_sub), 2, function(x) length(which(x>0)))

tcell_sub_2 = tcell_sub[!is.na(genecounts_tcell) & genecounts_tcell > 100, cellcounts_tcell > 0 & checkcounts_tcell > 0]  # some genes removed

vst_out_tcell = sctransform::vst(counts(tcell_sub_2), cell_attr = colData(tcell_sub_2), method="nb_fast", 
                                latent_var = c('log_umi'), latent_var_nonreg = c("g2m_score", "s_score"), 
                                return_gene_attr = T, return_cell_attr = T)

print("Start performing sctransform::smooth_via_pca: ")
y_smooth_tcell = sctransform::smooth_via_pca(vst_out_tcell$y, do_plot = TRUE)

print("Start performing sctransform::correct: ")
dat_cor_tcell = sctransform::correct(vst_out_tcell, data = y_smooth_tcell,
                                    do_round = TRUE, do_pos = FALSE)

# Plot mean variance plot
mvp_tcell <- ggplot(vst_out_tcell$gene_attr, aes(log10(gmean), log10(residual_variance))) +
  geom_point(alpha=0.3, shape=16) +
  #geom_density_2d(size = 0.3)
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".tcells.mean_vs_variance_plot.png", mvp_tcell, dpi = 300)
mvp2_tcell <- ggplot(vst_out_tcell$gene_attr, aes(log10(gmean), residual_variance)) +
  geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3) +
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".tcells.mean_vs_variance_plot_2.png", mvp2_tcell, dpi = 300)

## make SingleCellExperiment object and save to disk
vst_out_tcell$gene_attr$gene_ids = rownames(vst_out_tcell$gene_attr)
rowData(tcell_sub_2)$gene_ids = rownames(tcell_sub_2)

gene_desc = merge(vst_out_tcell$gene_attr, rowData(tcell_sub_2), sort=F, by = "gene_ids")
idx = match(rownames(dat_cor_tcell), rownames(tcell_sub_2))
idy = match(colnames(dat_cor_tcell), colnames(tcell_sub_2))
sce_tcell = SingleCellExperiment(assays = list(counts=tcell_sub_2[idx, idy],
                                               normcounts=dat_cor_tcell,
                                               pearson_resid=vst_out_tcell$y),
                                 colData = colData(tcell_sub_2),
                                 rowData = gene_desc)

saveRDS(sce_tcell, opt$outdir %&% opt$outName %&% ".t_nk_cells.subtype.RDS")
