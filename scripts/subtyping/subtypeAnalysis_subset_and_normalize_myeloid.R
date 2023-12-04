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

take_mast = c(25)
take_neutro= c(5)
take_mono = c(15)
take_dendr = c(10)

all_subset <- combined[, combined$umap_cohort_cl %in% c(take_mast, take_neutro, take_mono, take_dendr)]
all_subset$subtype <- ""
all_subset$subtype[all_subset$umap_cohort_cl %in% take_mast] <- "Mast.cells"
all_subset$subtype[all_subset$umap_cohort_cl %in% take_neutro] <- "Neutrophils"
all_subset$subtype[all_subset$umap_cohort_cl %in% take_mono] <- "Monocytes.macrophages"
all_subset$subtype[all_subset$umap_cohort_cl %in% take_dendr] <- "Plasmacytoid.dendritic"

print("working on the combined subset Myeloid")
genecounts_sub = rowSums(counts(all_subset))
cellcounts_sub = colSums(counts(all_subset))
checkcounts_sub = apply(counts(all_subset), 2, function(x) length(which(x>0)))

all_subset_2 = all_subset[!is.na(genecounts_sub) & genecounts_sub > 100, cellcounts_sub > 0 & checkcounts_sub > 0]  # some genes removed

vst_out_sub = sctransform::vst(counts(all_subset_2), cell_attr = colData(all_subset_2), method="nb_fast", 
                                latent_var = c('log_umi'), latent_var_nonreg = c("g2m_score", "s_score"), 
                                return_gene_attr = T, return_cell_attr = T)

print("Start performing sctransform::smooth_via_pca: ")
y_smooth_sub = sctransform::smooth_via_pca(vst_out_sub$y, do_plot = TRUE)

print("Start performing sctransform::correct: ")
dat_cor_sub = sctransform::correct(vst_out_sub, data = y_smooth_sub,
                                    do_round = TRUE, do_pos = FALSE)

# Plot mean variance plot
mvp_sub <- ggplot(vst_out_sub$gene_attr, aes(log10(gmean), log10(residual_variance))) +
  geom_point(alpha=0.3, shape=16) +
  #geom_density_2d(size = 0.3)
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".subset_myeloid.mean_vs_variance_plot.png", mvp_sub, dpi = 300)
mvp2_sub <- ggplot(vst_out_sub$gene_attr, aes(log10(gmean), residual_variance)) +
  geom_point(alpha=0.3, shape=16) +
  geom_density_2d(size = 0.3) +
  geom_smooth()
ggsave(opt$outdir %&% opt$outName %&% ".subset_myeloid.mean_vs_variance_plot_2.png", mvp2_sub, dpi = 300)

## make SingleCellExperiment object and save to disk
vst_out_sub$gene_attr$gene_ids = rownames(vst_out_sub$gene_attr)
rowData(all_subset_2)$gene_ids = rownames(all_subset_2)

gene_desc = merge(vst_out_sub$gene_attr, rowData(all_subset_2), sort=F, by = "gene_ids")
idx = match(rownames(dat_cor_sub), rownames(all_subset_2))
idy = match(colnames(dat_cor_sub), colnames(all_subset_2))
sce_sub = SingleCellExperiment(assays = list(counts=counts(all_subset_2[idx, idy]),
                                               normcounts=dat_cor_sub,
                                               pearson_resid=vst_out_sub$y),
                                 colData = colData(all_subset_2),
                                 rowData = gene_desc)

saveRDS(sce_sub, opt$outdir %&% opt$outName %&% ".subset_myeloid.subtype.RDS")

