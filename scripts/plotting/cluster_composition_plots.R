###################################################
## Date created: July 2021
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

all_umap_clusters <- unique(colData(subtype_sce)$umap_cohort_cl)
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
  subset_cluster <- subtype_sce[, subtype_sce$umap_cohort_cl == cluster_ID]
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

ggplot(all_clust_celltype, aes(fill=Celltype, y=frac_Celltype, x=clustID)) + 
  geom_bar(position="stack", stat="identity")
filename <- paste0(opt$outdir, opt$outName, "_clusters_celltype_final.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

ggplot(all_clust_celltype_mj, aes(fill=Celltype_Major, y=frac_Celltype_Major, x=clustID)) + 
  geom_bar(position="stack", stat="identity")
filename <- paste0(opt$outdir, opt$outName, "_clusters_celltype_major.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

ggplot(all_clust_tissue, aes(fill=Tissue_Region, y=frac_Tissue_Region, x=clustID)) + 
  geom_bar(position="stack", stat="identity")
filename <- paste0(opt$outdir, opt$outName, "_clusters_tissue.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

ggplot(all_clust_cohort, aes(fill=Cohort_Highlevel, y=frac_Cohort_Highlevel, x=clustID)) + 
  geom_bar(position="stack", stat="identity")
filename <- paste0(opt$outdir, opt$outName, "_clusters_cohort.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

ggplot(all_clust_site, aes(fill=Site, y=frac_Site, x=clustID)) + 
  geom_bar(position="stack", stat="identity")
filename <- paste0(opt$outdir, opt$outName, "_clusters_site.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

ggplot(all_clust_patient, aes(fill=Patient, y=frac_Patient, x=clustID)) + 
  geom_bar(position="stack", stat="identity")
filename <- paste0(opt$outdir, opt$outName, "_clusters_patient.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 300, units = "cm")

#########
# provide celltype composition per tissue_type, tissue_region, cohort_highlevel, patient, sample
#########

all_patients <- unique(colData(subtype_sce)$Patient)
print("patients:")
print(all_patients)

# celltype final
all_cts_patient <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_patient = c("Patient", "Celltype", "cells_Celltype", "frac_Celltype")
colnames(all_cts_patient) <- myNames_cts_patient
# celltype major
all_cts_mj_patient <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_mj_patient = c("Patient", "Celltype_Major", "cells_Celltype_Major", "frac_Celltype_Major")
colnames(all_cts_mj_patient) <- myNames_cts_mj_patient

for (patient in all_patients) {
  print(patient)
  subset_patient <- subtype_sce[, colData(subtype_sce)$Patient == patient]
  print(paste0("Number of cells in patient subset: ", dim(subset_patient)[2]))
  total_cells = dim(subset_patient)[2]
  
  print("Celltype major")
  subset_patient_ct_mj <- unique(colData(subset_patient)$celltype_major)
  print(subset_patient_ct_mj)
  for (celltype_mj in subset_patient_ct_mj){
    subset_celltype_mj <- subset_patient[, colData(subset_patient)$celltype_major == celltype_mj]
    celltype_cells_mj = dim(subset_celltype_mj)[2]
    print(celltype_cells_mj)
    frac_celltype_mj = (celltype_cells_mj/total_cells) * 100.0
    print(frac_celltype_mj)
    temp = data.frame(patient,celltype_mj,celltype_cells_mj,frac_celltype_mj)
    colnames(temp) = myNames_cts_mj_patient
    all_cts_mj_patient = rbind(all_cts_mj_patient, temp)
  }
  print("Celltype final")
  subset_patient_ct <- unique(colData(subset_patient)$celltype_final)
  print(subset_patient_ct)
  for (celltype in subset_patient_ct){
    subset_celltype <- subset_patient[, colData(subset_patient)$celltype_final == celltype]
    celltype_cells = dim(subset_celltype)[2]
    print(celltype_cells)
    frac_celltype = (celltype_cells/total_cells) * 100.0
    print(frac_celltype)
    temp = data.frame(patient,celltype,celltype_cells,frac_celltype)
    colnames(temp) = myNames_cts_patient
    all_cts_patient = rbind(all_cts_patient, temp)
  }
}

all_tissue_types <- unique(colData(subtype_sce)$Tissue_Type)
print("Tissue_Type:")
print(all_tissue_types)

# celltype final
all_cts_tissue_type <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_tissue_type = c("Tissue_Type", "Celltype", "cells_Celltype", "frac_Celltype")
colnames(all_cts_tissue_type) <- myNames_cts_tissue_type
# celltype major
all_cts_mj_tissue_type <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_mj_tissue_type = c("Tissue_Type", "Celltype_Major", "cells_Celltype_Major", "frac_Celltype_Major")
colnames(all_cts_mj_tissue_type) <- myNames_cts_mj_tissue_type

for (tissue_type in all_tissue_types) {
  print(tissue_type)
  subset_tissue_type <- subtype_sce[, colData(subtype_sce)$Tissue_Type == tissue_type]
  print(paste0("Number of cells in tissue_type subset: ", dim(subset_tissue_type)[2]))
  total_cells = dim(subset_tissue_type)[2]
  
  print("Celltype major")
  subset_tissue_type_ct_mj <- unique(colData(subset_tissue_type)$celltype_major)
  print(subset_tissue_type_ct_mj)
  for (celltype_mj in subset_tissue_type_ct_mj){
    subset_celltype_mj <- subset_tissue_type[, colData(subset_tissue_type)$celltype_major == celltype_mj]
    celltype_cells_mj = dim(subset_celltype_mj)[2]
    print(celltype_cells_mj)
    frac_celltype_mj = (celltype_cells_mj/total_cells) * 100.0
    print(frac_celltype_mj)
    temp = data.frame(tissue_type,celltype_mj,celltype_cells_mj,frac_celltype_mj)
    colnames(temp) = myNames_cts_mj_tissue_type
    all_cts_mj_tissue_type = rbind(all_cts_mj_tissue_type, temp)
  }
  print("Celltype final")
  subset_tissue_type_ct <- unique(colData(subset_tissue_type)$celltype_final)
  print(subset_tissue_type_ct)
  for (celltype in subset_tissue_type_ct){
    subset_celltype <- subset_tissue_type[, colData(subset_tissue_type)$celltype_final == celltype]
    celltype_cells = dim(subset_celltype)[2]
    print(celltype_cells)
    frac_celltype = (celltype_cells/total_cells) * 100.0
    print(frac_celltype)
    temp = data.frame(tissue_type,celltype,celltype_cells,frac_celltype)
    colnames(temp) = myNames_cts_tissue_type
    all_cts_tissue_type = rbind(all_cts_tissue_type, temp)
  }
}

all_tissue_regions <- unique(colData(subtype_sce)$Tissue_Region)
print("Tissue_Region:")
print(all_tissue_regions)

# celltype final
all_cts_tissue_region <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_tissue_region = c("Tissue_Region", "Celltype", "cells_Celltype", "frac_Celltype")
colnames(all_cts_tissue_region) <- myNames_cts_tissue_region
# celltype major
all_cts_mj_tissue_region <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_mj_tissue_region = c("Tissue_Region", "Celltype_Major", "cells_Celltype_Major", "frac_Celltype_Major")
colnames(all_cts_mj_tissue_region) <- myNames_cts_mj_tissue_region

for (tissue_region in all_tissue_regions) {
  print(tissue_region)
  subset_tissue_region <- subtype_sce[, colData(subtype_sce)$Tissue_Region == tissue_region]
  print(paste0("Number of cells in tissue_region subset: ", dim(subset_tissue_region)[2]))
  total_cells = dim(subset_tissue_region)[2]
  
  print("Celltype major")
  subset_tissue_region_ct_mj <- unique(colData(subset_tissue_region)$celltype_major)
  print(subset_tissue_region_ct_mj)
  for (celltype_mj in subset_tissue_region_ct_mj){
    subset_celltype_mj <- subset_tissue_region[, colData(subset_tissue_region)$celltype_major == celltype_mj]
    celltype_cells_mj = dim(subset_celltype_mj)[2]
    print(celltype_cells_mj)
    frac_celltype_mj = (celltype_cells_mj/total_cells) * 100.0
    print(frac_celltype_mj)
    temp = data.frame(tissue_region,celltype_mj,celltype_cells_mj,frac_celltype_mj)
    colnames(temp) = myNames_cts_mj_tissue_region
    all_cts_mj_tissue_region = rbind(all_cts_mj_tissue_region, temp)
  }
  print("Celltype final")
  subset_tissue_region_ct <- unique(colData(subset_tissue_region)$celltype_final)
  print(subset_tissue_region_ct)
  for (celltype in subset_tissue_region_ct){
    subset_celltype <- subset_tissue_region[, colData(subset_tissue_region)$celltype_final == celltype]
    celltype_cells = dim(subset_celltype)[2]
    print(celltype_cells)
    frac_celltype = (celltype_cells/total_cells) * 100.0
    print(frac_celltype)
    temp = data.frame(tissue_region,celltype,celltype_cells,frac_celltype)
    colnames(temp) = myNames_cts_tissue_region
    all_cts_tissue_region = rbind(all_cts_tissue_region, temp)
  }
}

all_cohorts <- unique(colData(subtype_sce)$Cohort_Type_Highlevel)
print("Cohort_Type_Highlevel:")
print(all_cohorts)

# celltype final
all_cts_cohort <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_cohort = c("Cohort", "Celltype", "cells_Celltype", "frac_Celltype")
colnames(all_cts_cohort) <- myNames_cts_cohort
# celltype major
all_cts_mj_cohort <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_mj_cohort = c("Cohort", "Celltype_Major", "cells_Celltype_Major", "frac_Celltype_Major")
colnames(all_cts_mj_cohort) <- myNames_cts_mj_cohort

for (cohort in all_cohorts) {
  print(cohort)
  subset_cohort <- subtype_sce[, colData(subtype_sce)$Cohort_Type_Highlevel == cohort]
  print(paste0("Number of cells in cohort subset: ", dim(subset_cohort)[2]))
  total_cells = dim(subset_cohort)[2]
  
  print("Celltype major")
  subset_cohort_ct_mj <- unique(colData(subset_cohort)$celltype_major)
  print(subset_cohort_ct_mj)
  for (celltype_mj in subset_cohort_ct_mj){
    subset_celltype_mj <- subset_cohort[, colData(subset_cohort)$celltype_major == celltype_mj]
    celltype_cells_mj = dim(subset_celltype_mj)[2]
    print(celltype_cells_mj)
    frac_celltype_mj = (celltype_cells_mj/total_cells) * 100.0
    print(frac_celltype_mj)
    temp = data.frame(cohort,celltype_mj,celltype_cells_mj,frac_celltype_mj)
    colnames(temp) = myNames_cts_mj_cohort
    all_cts_mj_cohort = rbind(all_cts_mj_cohort, temp)
  }
  print("Celltype final")
  subset_cohort_ct <- unique(colData(subset_cohort)$celltype_final)
  print(subset_cohort_ct)
  for (celltype in subset_cohort_ct){
    subset_celltype <- subset_cohort[, colData(subset_cohort)$celltype_final == celltype]
    celltype_cells = dim(subset_celltype)[2]
    print(celltype_cells)
    frac_celltype = (celltype_cells/total_cells) * 100.0
    print(frac_celltype)
    temp = data.frame(cohort,celltype,celltype_cells,frac_celltype)
    colnames(temp) = myNames_cts_cohort
    all_cts_cohort = rbind(all_cts_cohort, temp)
  }
}

all_samples <- unique(colData(subtype_sce)$sampleID)
print("sampleID:")
print(all_samples)

# celltype final
all_cts_sample <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_sample = c("Sample", "Celltype", "cells_Celltype", "frac_Celltype")
colnames(all_cts_sample) <- myNames_cts_sample
# celltype major
all_cts_mj_sample <- data.frame(matrix(ncol = 4, nrow = 0))
myNames_cts_mj_sample = c("Sample", "Celltype_Major", "cells_Celltype_Major", "frac_Celltype_Major")
colnames(all_cts_mj_sample) <- myNames_cts_mj_sample

for (sample in all_samples) {
  print(sample)
  subset_sample <- subtype_sce[, colData(subtype_sce)$sampleID == sample]
  print(paste0("Number of cells in sample subset: ", dim(subset_sample)[2]))
  total_cells = dim(subset_sample)[2]
  
  print("Celltype major")
  subset_sample_ct_mj <- unique(colData(subset_sample)$celltype_major)
  print(subset_sample_ct_mj)
  for (celltype_mj in subset_sample_ct_mj){
    subset_celltype_mj <- subset_sample[, colData(subset_sample)$celltype_major == celltype_mj]
    celltype_cells_mj = dim(subset_celltype_mj)[2]
    print(celltype_cells_mj)
    frac_celltype_mj = (celltype_cells_mj/total_cells) * 100.0
    print(frac_celltype_mj)
    temp = data.frame(sample,celltype_mj,celltype_cells_mj,frac_celltype_mj)
    colnames(temp) = myNames_cts_mj_sample
    all_cts_mj_sample = rbind(all_cts_mj_sample, temp)
  }
  print("Celltype final")
  subset_sample_ct <- unique(colData(subset_sample)$celltype_final)
  print(subset_sample_ct)
  for (celltype in subset_sample_ct){
    subset_celltype <- subset_sample[, colData(subset_sample)$celltype_final == celltype]
    celltype_cells = dim(subset_celltype)[2]
    print(celltype_cells)
    frac_celltype = (celltype_cells/total_cells) * 100.0
    print(frac_celltype)
    temp = data.frame(sample,celltype,celltype_cells,frac_celltype)
    colnames(temp) = myNames_cts_sample
    all_cts_sample = rbind(all_cts_sample, temp)
  }
}

# write data frames to files

# patient
# celltype final
filename_celltype <- paste0(opt$outdir, opt$outName, ".celltype_final_comp_patient.tsv")
write.table(all_cts_patient, file = filename_celltype, quote = FALSE, sep = "\t", row.names = FALSE)
#celltype major
filename_celltype_major <- paste0(opt$outdir, opt$outName, ".celltype_major_comp_patient.tsv")
write.table(all_cts_mj_patient, file = filename_celltype_major, quote = FALSE, sep = "\t", row.names = FALSE)

# tissue_type
# celltype final
filename_celltype <- paste0(opt$outdir, opt$outName, ".celltype_final_comp_tissue_type.tsv")
write.table(all_cts_tissue_type, file = filename_celltype, quote = FALSE, sep = "\t", row.names = FALSE)
#celltype major
filename_celltype_major <- paste0(opt$outdir, opt$outName, ".celltype_major_comp_tissue_type.tsv")
write.table(all_cts_mj_tissue_type, file = filename_celltype_major, quote = FALSE, sep = "\t", row.names = FALSE)

# tissue_region
# celltype final
filename_celltype <- paste0(opt$outdir, opt$outName, ".celltype_final_comp_tissue_region.tsv")
write.table(all_cts_tissue_region, file = filename_celltype, quote = FALSE, sep = "\t", row.names = FALSE)
#celltype major
filename_celltype_major <- paste0(opt$outdir, opt$outName, ".celltype_major_comp_tissue_region.tsv")
write.table(all_cts_mj_tissue_region, file = filename_celltype_major, quote = FALSE, sep = "\t", row.names = FALSE)

# cohort
# celltype final
filename_celltype <- paste0(opt$outdir, opt$outName, ".celltype_final_comp_cohort.tsv")
write.table(all_cts_cohort, file = filename_celltype, quote = FALSE, sep = "\t", row.names = FALSE)
#celltype major
filename_celltype_major <- paste0(opt$outdir, opt$outName, ".celltype_major_comp_cohort.tsv")
write.table(all_cts_mj_cohort, file = filename_celltype_major, quote = FALSE, sep = "\t", row.names = FALSE)

# sample
# celltype final
filename_celltype <- paste0(opt$outdir, opt$outName, ".celltype_final_comp_sample.tsv")
write.table(all_cts_sample, file = filename_celltype, quote = FALSE, sep = "\t", row.names = FALSE)
#celltype major
filename_celltype_major <- paste0(opt$outdir, opt$outName, ".celltype_major_comp_sample.tsv")
write.table(all_cts_mj_sample, file = filename_celltype_major, quote = FALSE, sep = "\t", row.names = FALSE)


