###################################################
## Date created: November 2021
## R Version: 4
###################################################

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
library(plyr)
library(dplyr)

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
max_hvg_length = 1000
UMAP_aspect_ratio <- "1"

###################################
###   Read in and format data   ###
###################################

subtype_sce = readRDS(opt$rdsFile)

colData(subtype_sce)$Tissue_Region_old <- colData(subtype_sce)$Tissue_Region
subtype_sce$Tissue_Region[subtype_sce$Tissue_Region == "Omentum"] = "benign"
subtype_sce$Tissue_Region[subtype_sce$Tissue_Region == "Distal"] = "tumor_distant"
subtype_sce$Tissue_Region[subtype_sce$Tissue_Region == "Next"] = "tumor_peritumoral"
subtype_sce$Tissue_Region[subtype_sce$Tissue_Region == "Tumor"] = "tumor_metastasis"

all_umap_clusters <- unique(colData(subtype_sce)$umap_subtype_cl)
print("umap clusters:")
print(all_umap_clusters)

all_clust_tissue <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_tissue = c("clustID", "Tissue_Region", "cells_Tissue_Region", "frac_Tissue_Region")
colnames(all_clust_tissue) <- myNames_tissue

all_clust_cohort <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cohort = c("clustID", "Cohort_Highlevel", "cells_Cohort_Highlevel", "frac_Cohort_Highlevel")
colnames(all_clust_cohort) <- myNames_cohort

all_clust_site <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_site = c("clustID", "Site", "cells_Site", "frac_Site")
colnames(all_clust_site) <- myNames_site

all_clust_patient <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_patient = c("clustID", "Patient", "cells_Patient", "frac_Patient")
colnames(all_clust_patient) <- myNames_patient

all_clust_celltype <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_celltype = c("clustID", "Celltype", "cells_Celltype", "frac_Celltype")
colnames(all_clust_celltype) <- myNames_celltype

all_clust_celltype_mj <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_celltype_mj = c("clustID", "Celltype_Major", "cells_Celltype_Major", "frac_Celltype_Major")
colnames(all_clust_celltype_mj) <- myNames_celltype_mj

for (cluster_ID in all_umap_clusters) {
  print(cluster_ID)
  subset_cluster <- subtype_sce[, subtype_sce$umap_subtype_cl == cluster_ID]
  print(paste0("Number of cells in cluster: ", dim(subset_cluster)[2]))
  total_cells = dim(subset_cluster)[2]
  
  print("Celltype major")
  subset_cluster_celltype_mj <- unique(colData(subset_cluster)$celltype_major)
  print(subset_cluster_celltype_mj)
  for (celltype_mj in subset_cluster_celltype_mj){
    subset_celltype_mj <- subset_cluster[, colData(subset_cluster)$celltype_major == celltype_mj]
    celltype_cells_mj = dim(subset_celltype_mj)[2]
    print(celltype_cells_mj)
    frac_celltype_mj = (celltype_cells_mj/total_cells) * 100.0
    print(frac_celltype_mj)
    temp = data.frame(cluster_ID,celltype_mj,celltype_cells_mj,frac_celltype_mj)
    colnames(temp) = myNames_celltype_mj
    all_clust_celltype_mj = rbind(all_clust_celltype_mj, temp)
  }
  
  print("Celltype final")
  subset_cluster_celltype <- unique(colData(subset_cluster)$celltype_final)
  print(subset_cluster_celltype)
  for (celltype in subset_cluster_celltype){
    subset_celltype <- subset_cluster[, colData(subset_cluster)$celltype_final == celltype]
    celltype_cells = dim(subset_celltype)[2]
    print(celltype_cells)
    frac_celltype = (celltype_cells/total_cells) * 100.0
    print(frac_celltype)
    temp = data.frame(cluster_ID,celltype,celltype_cells,frac_celltype)
    colnames(temp) = myNames_celltype
    all_clust_celltype = rbind(all_clust_celltype, temp)
  }
  
  print("Region")
  subset_cluster_tissue_regions <- unique(colData(subset_cluster)$Tissue_Region)
  print(subset_cluster_tissue_regions)
  for (tissue in subset_cluster_tissue_regions){
    subset_tissue <- subset_cluster[, colData(subset_cluster)$Tissue_Region == tissue]
    tissue_cells = dim(subset_tissue)[2]
    print(tissue_cells)
    frac_tissue = (tissue_cells/total_cells) * 100.0
    print(frac_tissue)
    temp = data.frame(cluster_ID,tissue,tissue_cells,frac_tissue)
    colnames(temp) = myNames_tissue
    all_clust_tissue = rbind(all_clust_tissue, temp)
  }
  
  print("Cohort")
  subset_cluster_cohortType <- unique(colData(subset_cluster)$Cohort_Type_Highlevel)
  print(subset_cluster_cohortType)
  for (cohort in subset_cluster_cohortType){
    subset_cohort <- subset_cluster[, colData(subset_cluster)$Cohort_Type_Highlevel == cohort]
    cohort_cells = dim(subset_cohort)[2]
    print(cohort_cells)
    frac_cohort = (cohort_cells/total_cells) * 100.0
    print(frac_cohort)
    temp = data.frame(cluster_ID,cohort,cohort_cells,frac_cohort)
    colnames(temp) = myNames_cohort
    all_clust_cohort = rbind(all_clust_cohort, temp)
  }
  
  print("Site")
  subset_cluster_site <- unique(colData(subset_cluster)$Tissue_Type)
  print(subset_cluster_site)
  for (site in subset_cluster_site){
    subset_site <- subset_cluster[, colData(subset_cluster)$Tissue_Type == site]
    site_cells = dim(subset_site)[2]
    print(site_cells)
    frac_site = (site_cells/total_cells) * 100.0
    print(frac_site)
    temp = data.frame(cluster_ID,site,site_cells,frac_site)
    colnames(temp) = myNames_site
    all_clust_site = rbind(all_clust_site, temp)
  }
  
  print("Patient")
  subset_cluster_patient <- unique(colData(subset_cluster)$Patient)
  print(subset_cluster_patient)
  for (patient in subset_cluster_patient){
    print(patient)
    subset_patient <- subset_cluster[, colData(subset_cluster)$Patient == patient]
    patient_cells = dim(subset_patient)[2]
    print(patient_cells)
    frac_patient = (patient_cells/total_cells) * 100.0
    print(frac_patient)
    temp = data.frame(cluster_ID,patient,patient_cells,frac_patient)
    colnames(temp) = myNames_patient
    all_clust_patient = rbind(all_clust_patient, temp)
  }
}
print("Tissue df")
all_clust_tissue
print("cohort df")
all_clust_cohort
print("Site df")
all_clust_site
print("Patient df")
all_clust_patient
print("Celltype major df")
all_clust_celltype_mj
print("Celltype final df")
all_clust_celltype

filename_celltype <- paste0(opt$outdir, opt$outName, ".celltype_final_comp_cluster.tsv")
write.table(all_clust_celltype, file = filename_celltype, quote = FALSE, sep = "\t", row.names = FALSE)

filename_celltype_major <- paste0(opt$outdir, opt$outName, ".celltype_major_comp_cluster.tsv")
write.table(all_clust_celltype_mj, file = filename_celltype_major, quote = FALSE, sep = "\t", row.names = FALSE)

filename_tissue <- paste0(opt$outdir, opt$outName, ".tissue_comp_cluster.tsv")
write.table(all_clust_tissue, file = filename_tissue, quote = FALSE, sep = "\t", row.names = FALSE)

filename_cohort <- paste0(opt$outdir, opt$outName, ".cohort_comp_cluster.tsv")
write.table(all_clust_cohort, file = filename_cohort, quote = FALSE, sep = "\t", row.names = FALSE)

filename_site <- paste0(opt$outdir, opt$outName, ".site_comp_cluster.tsv")
write.table(all_clust_site, file = filename_site, quote = FALSE, sep = "\t", row.names = FALSE)

filename_patient <- paste0(opt$outdir, opt$outName, ".patient_comp_cluster.tsv")
write.table(all_clust_patient, file = filename_patient, quote = FALSE, sep = "\t", row.names = FALSE)


ggplot(all_clust_celltype, aes(fill=Celltype, y=frac_Celltype, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Cluster ID", y = "Fraction") +
  theme(axis.text.x = element_text(size = 13),
	axis.title.x = element_text(size = 18),
	axis.title.y = element_text(size = 18),
	axis.text.y= element_text(size = 13),
	legend.text = element_text(size = 13),
	legend.title = element_text(size = 18))

filename <- paste0(opt$outdir, opt$outName, "_clusters_celltype_final.waterfall.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")

ggplot(all_clust_celltype_mj, aes(fill=Celltype_Major, y=frac_Celltype_Major, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Cluster ID", y = "Fraction") +
  theme(axis.text.x = element_text(size = 13),
	axis.title.x = element_text(size = 18),
	axis.title.y = element_text(size = 18),
	axis.text.y= element_text(size = 13),
	legend.text = element_text(size = 13),
	legend.title = element_text(size = 18))

filename <- paste0(opt$outdir, opt$outName, "_clusters_celltype_major.waterfall.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")

all_clust_tissue_ord <- arrange(all_clust_tissue[all_clust_tissue$Tissue_Region == "benign",], desc(frac_Tissue_Region))
clustTissue_ordering = union(c(all_clust_tissue_ord$clustID),c(unique(all_clust_tissue$clustID)))

ggplot(all_clust_tissue, aes(fill=Tissue_Region, y=frac_Tissue_Region, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_x_discrete(limits = factor(clustTissue_ordering)) +
  labs(x ="Cluster ID", y = "Fraction") + 
  theme(axis.text.x = element_text(size = 13),
	axis.title.x = element_text(size = 18),
	axis.title.y = element_text(size = 18),
	axis.text.y= element_text(size = 13),
	legend.text = element_text(size = 13),
	legend.title = element_text(size = 18))

filename <- paste0(opt$outdir, opt$outName, "_clusters_tissue.waterfall.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")

all_clust_cohort_ord <- arrange(all_clust_cohort[all_clust_cohort$Cohort_Highlevel == "Benign_Omentum",], desc(frac_Cohort_Highlevel))
clustCohort_ordering = union(c(all_clust_cohort_ord$clustID),c(unique(all_clust_cohort$clustID)))

ggplot(all_clust_cohort, aes(fill=Cohort_Highlevel, y=frac_Cohort_Highlevel, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(limits = factor(clustCohort_ordering)) +
  labs(x ="Cluster ID", y = "Fraction") + 
  theme(axis.text.x = element_text(size = 13),
	axis.title.x = element_text(size = 18),
	axis.title.y = element_text(size = 18),
	axis.text.y= element_text(size = 13),
	legend.text = element_text(size = 13),
	legend.title = element_text(size = 18))

filename <- paste0(opt$outdir, opt$outName, "_clusters_cohort.waterfall.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")

all_clust_site_ord <- arrange(all_clust_site[all_clust_site$Site == "Majus",], desc(frac_Site))
clustSite_ordering = union(c(all_clust_site_ord$clustID),c(unique(all_clust_site$clustID)))

ggplot(all_clust_site, aes(fill=Site, y=frac_Site, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") +
  scale_x_discrete(limits = factor(clustSite_ordering)) +
  labs(x ="Cluster ID", y = "Fraction") + 
  theme(axis.text.x = element_text(size = 13),
	axis.title.x = element_text(size = 18),
	axis.title.y = element_text(size = 18),
	axis.text.y= element_text(size = 13),
	legend.text = element_text(size = 13),
	legend.title = element_text(size = 18))

filename <- paste0(opt$outdir, opt$outName, "_clusters_site.waterfall.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")

patient_colours = c(brewer.pal(12, "Paired"),brewer.pal(8, "Accent"))
length(unique(all_clust_patient$Patient))
patient_colours[seq_along(unique(all_clust_patient$Patient))]

ggplot(all_clust_patient, aes(fill=Patient, y=frac_Patient, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = patient_colours[seq_along(unique(all_clust_patient$Patient))]) +
  labs(x ="Cluster ID", y = "Fraction") + 
  theme(axis.text.x = element_text(size = 13),
	axis.title.x = element_text(size = 18),
	axis.title.y = element_text(size = 18),
	axis.text.y= element_text(size = 13),
	legend.text = element_text(size = 13),
	legend.title = element_text(size = 18))

filename <- paste0(opt$outdir, opt$outName, "_clusters_patient.waterfall.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")

