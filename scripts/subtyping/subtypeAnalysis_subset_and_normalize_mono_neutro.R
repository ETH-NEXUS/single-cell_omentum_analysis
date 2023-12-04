###################################################
## File name: subtypeAnalysis_subset_and_normalize.R
## Date created: July 2021
## R Version: 4
###################################################

## GENERAL:
##Monocytes and Neutrophils

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

# Monocytes and neutrophil cells - here only subset by cluster level!

# define the cluster IDs that should be taken into account for the different cell types

take_neutro = c(5)
take_mono = c(15)

#take the desired clusters
neutro_sub <- combined[, combined$umap_cohort_cl %in% take_neutro]
mono_sub <- combined[, combined$umap_cohort_cl %in% take_mono]

# take the celltype subset and perform an additional sctransform normalization, on the counts assay

#Neutrophils
print("Neutrophils")
genecounts_neutro = rowSums(counts(neutro_sub))
cellcounts_neutro = colSums(counts(neutro_sub))
checkcounts_neutro = apply(counts(neutro_sub), 2, function(x) length(which(x>0)))

neutro_sub_2 = neutro_sub[!is.na(genecounts_neutro) & genecounts_neutro > 40, cellcounts_neutro > 0 & checkcounts_neutro > 0]  # some genes removed

vst_out_neutro = sctransform::vst(counts(neutro_sub_2), cell_attr = colData(neutro_sub_2), method="nb_fast", 
                                latent_var = c('log_umi'), latent_var_nonreg = c("g2m_score", "s_score"), 
                                return_gene_attr = T, return_cell_attr = T)

print("Start performing sctransform::smooth_via_pca: ")
y_smooth_neutro = sctransform::smooth_via_pca(vst_out_neutro$y, do_plot = TRUE)

print("Start performing sctransform::correct: ")
dat_cor_neutro = sctransform::correct(vst_out_neutro, data = y_smooth_neutro,
                                    do_round = TRUE, do_pos = FALSE)

# Plot mean variance plot
mvp_neutro <- ggplot(vst_out_neutro$gene_attr, aes(log10(gmean), log10(residual_variance))) +
  geom_point(alpha=0.3, shape=16) +
  #geom_density_2d(size = 0.3)
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".neutrophils.mean_vs_variance_plot.png", mvp_neutro, dpi = 300)
mvp2_neutro <- ggplot(vst_out_neutro$gene_attr, aes(log10(gmean), residual_variance)) +
  geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3) +
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".neutrophils.mean_vs_variance_plot_2.png", mvp2_neutro, dpi = 300)

## make SingleCellExperiment object and save to disk
vst_out_neutro$gene_attr$gene_ids = rownames(vst_out_neutro$gene_attr)
rowData(neutro_sub_2)$gene_ids = rownames(neutro_sub_2)

gene_desc = merge(vst_out_neutro$gene_attr, rowData(neutro_sub_2), sort=F, by = "gene_ids")
idx = match(rownames(dat_cor_neutro), rownames(neutro_sub_2))
idy = match(colnames(dat_cor_neutro), colnames(neutro_sub_2))
sce_neutro = SingleCellExperiment(assays = list(counts=neutro_sub_2[idx, idy],
                                               normcounts=dat_cor_neutro,
                                               pearson_resid=vst_out_neutro$y),
                                 colData = colData(neutro_sub_2),
                                 rowData = gene_desc)

saveRDS(sce_neutro, opt$outdir %&% opt$outName %&% ".neutrophils.subtype.RDS")


#Monocytes
print("Monocytes")
genecounts_mono = rowSums(counts(mono_sub))
cellcounts_mono = colSums(counts(mono_sub))
checkcounts_mono = apply(counts(mono_sub), 2, function(x) length(which(x>0)))

mono_sub_2 = mono_sub[!is.na(genecounts_mono) & genecounts_mono > 40, cellcounts_mono > 0 & checkcounts_mono > 0]  # some genes removed

vst_out_mono = sctransform::vst(counts(mono_sub_2), cell_attr = colData(mono_sub_2), method="nb_fast", 
                                latent_var = c('log_umi'), latent_var_nonreg = c("g2m_score", "s_score"), 
                                return_gene_attr = T, return_cell_attr = T)

print("Start performing sctransform::smooth_via_pca: ")
y_smooth_mono = sctransform::smooth_via_pca(vst_out_mono$y, do_plot = TRUE)

print("Start performing sctransform::correct: ")
dat_cor_mono = sctransform::correct(vst_out_mono, data = y_smooth_mono,
                                    do_round = TRUE, do_pos = FALSE)

# Plot mean variance plot
mvp_mono <- ggplot(vst_out_mono$gene_attr, aes(log10(gmean), log10(residual_variance))) +
  geom_point(alpha=0.3, shape=16) +
  #geom_density_2d(size = 0.3)
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".monocytes.mean_vs_variance_plot.png", mvp_mono, dpi = 300)
mvp2_mono <- ggplot(vst_out_mono$gene_attr, aes(log10(gmean), residual_variance)) +
  geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3) +
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".monocytes.mean_vs_variance_plot_2.png", mvp2_mono, dpi = 300)

## make SingleCellExperiment object and save to disk
vst_out_mono$gene_attr$gene_ids = rownames(vst_out_mono$gene_attr)
rowData(mono_sub_2)$gene_ids = rownames(mono_sub_2)

gene_desc = merge(vst_out_mono$gene_attr, rowData(mono_sub_2), sort=F, by = "gene_ids")
idx = match(rownames(dat_cor_mono), rownames(mono_sub_2))
idy = match(colnames(dat_cor_mono), colnames(mono_sub_2))
sce_mono = SingleCellExperiment(assays = list(counts=mono_sub_2[idx, idy],
                                               normcounts=dat_cor_mono,
                                               pearson_resid=vst_out_mono$y),
                                 colData = colData(mono_sub_2),
                                 rowData = gene_desc)

saveRDS(sce_mono, opt$outdir %&% opt$outName %&% ".monocytes.subtype.RDS")
