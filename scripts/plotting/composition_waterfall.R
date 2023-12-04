###################################################
## Date created: July 2021
## R Version: 4
###################################################
library(ggplot2)
library(scMerge)
library(uwot)
library(patchwork)
library(optparse)
library(reshape2)
library(igraph)
library(ggrepel)
library(dynamicTreeCut)
library(RColorBrewer)
library(pheatmap)
library(plyr)
library(dplyr)

cat("\n\n\nPrint sessionInfo():\n\n")
print(sessionInfo())

path_clusterComp = "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/cluster_composition/"
path_celltypeComp = "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/"

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")

###################################
###   Read in and format data   ###
###################################

# for clusters show: celltype, celltype major, tissue_region, tissue_type, cohort, patient
# show celltypes (final and major) also for patient, sample, tissue_type, tissue_region, and cohort

clustCohort <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/cluster_composition/OmentumCohort.cohort_comp_cluster.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
clustRegion <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/cluster_composition/OmentumCohort.tissue_comp_cluster.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
clustSite <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/cluster_composition/OmentumCohort.site_comp_cluster.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
clustPatient <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/cluster_composition/OmentumCohort.patient_comp_cluster.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
clustCT_final <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/cluster_composition/OmentumCohort.celltype_final_comp_cluster.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE))
clustCT_major <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/cluster_composition/OmentumCohort.celltype_major_comp_cluster.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE))

patientCT_final <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_final_comp_patient.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
patientCT_major <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_major_comp_patient.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))

sampleCT_final <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_final_comp_sample_short.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
sampleCT_major <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_major_comp_sample_short.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))

cohortCT_final <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_final_comp_cohort.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
cohortCT_major <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_major_comp_cohort.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))

regionCT_final <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_final_comp_tissue_region.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
regionCT_major <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_major_comp_tissue_region.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))

siteCT_final <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_final_comp_tissue_type.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))
siteCT_major <- as.data.frame(read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_composition_plots/celltype_composition/OmentumCohort.celltype_major_comp_tissue_type.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE))

#########
## TISSUE REGION
#########

#unsorted:
ggplot(clustRegion, aes(fill=Tissue_Region, y=frac_Tissue_Region, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") + 
  xlab("cluster ID") + 
  ylab("fraction")

#subset(clustRegion, clustRegion$Tissue_Region == "benign")$frac_Tissue_Region
clustRegion_ordered <- arrange(clustRegion[clustRegion$Tissue_Region == "benign",], desc(frac_Tissue_Region))

ggplot(clustRegion, aes(fill=Tissue_Region, y=frac_Tissue_Region, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_x_discrete(limits = factor(clustRegion_ordered$clustID)) +
  labs(x ="Cluster ID", y = "Fraction") + 
  theme(axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_clusterComp, "OmentumCohort_clusters_tissue.png")
ggsave(filename = filename, width = 30, height = 15, dpi = 600, units = "cm")

#########
## COHORT
#########

#unsorted:
ggplot(clustCohort, aes(fill=Cohort_Highlevel, y=frac_Cohort_Highlevel, x=clustID)) + 
  geom_bar(position="stack", stat="identity")

clustCohort_ordered <- arrange(clustCohort[clustCohort$Cohort_Highlevel == "Benign_Omentum",], desc(frac_Cohort_Highlevel))

ggplot(clustCohort, aes(fill=Cohort_Highlevel, y=frac_Cohort_Highlevel, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_x_discrete(limits = factor(clustCohort_ordered$clustID)) +
  labs(x ="Cluster ID", y = "Fraction") + 
  theme(axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_clusterComp, "OmentumCohort_clusters_cohort.png")
ggsave(filename = filename, width = 30, height = 15, dpi = 600, units = "cm")

#########
## SITE
#########

#unsorted:
ggplot(clustSite, aes(fill=Site, y=frac_Site, x=clustID)) + 
  geom_bar(position="stack", stat="identity")

#clustSite_ordered <- arrange(clustSite[clustSite$Site == "Gastric",], desc(frac_Site))
clustSite_ordered <- arrange(clustSite[clustSite$Site == "Majus",], desc(frac_Site))

ggplot(clustSite, aes(fill=Site, y=frac_Site, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") + 
  scale_x_discrete(limits = factor(clustSite_ordered$clustID)) +
  labs(x ="Cluster ID", y = "Fraction") + 
  theme(axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_clusterComp, "OmentumCohort_clusters_site.png")
ggsave(filename = filename, width = 30, height = 15, dpi = 600, units = "cm")

#########
## PATIENT
#########

#unsorted, here sorting is probably not necessary:
patient_colours = c(brewer.pal(12, "Paired"),brewer.pal(8, "Accent"))
#Patient = patient_colours[seq_along(levels(ann_col$Patient))])
length(unique(clustPatient$Patient))
patient_colours[seq_along(unique(clustPatient$Patient))]

ggplot(clustPatient, aes(fill=Patient, y=frac_Patient, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Cluster ID", y = "Fraction") +
  #scale_fill_manual(values = c("#1B9E77", "#D95F02" ,"#7570B3" , "#66A61E"))
  scale_fill_manual(values = patient_colours[seq_along(unique(clustPatient$Patient))]) +
  theme(axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_clusterComp, "OmentumCohort_clusters_patient.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")

#########
## CELLTYPE FINAL
#########

config <- read.csv("gfb_omentum_2020/git/gfb_omentum_2020/omentum_required/2021-04_reanalysis/colour_config_omentum_2021-05.txt", sep = "\t", stringsAsFactors = FALSE)
config$cell_type <- gsub(pattern = "([^_]+)_.*", "\\1", config$cell_type)
stopifnot(length(config$colour) == length(config$cell_type))
stopifnot(config$colour == unique(config$colour))
ct.color <- c(config$colour, "grey50", "black")
names(ct.color) <- c(config$cell_type, "uncertain", "unknown")
print(ct.color)

#unsorted, here sorting is probably not necessary:
ggplot(clustCT_final, aes(fill=Celltype, y=frac_Celltype, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Cluster ID", y = "Fraction") +
  scale_fill_manual(values = ct.color[seq_along(unique(clustCT_final$Celltype))]) +
  theme(axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_clusterComp, "OmentumCohort_clusters_celltype_final.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")


#########
## CELLTYPE MAJOR
#########

subset_config_major = subset(config, config$cell_type %in% clustCT_major$Celltype_Major)
ct.color_major <- c(subset_config_major$colour, "grey50", "black")
names(ct.color_major) <- c(subset_config_major$cell_type, "uncertain", "unknown")
print(ct.color_major)

#unsorted, here sorting is probably not necessary:
ggplot(clustCT_major, aes(fill=Celltype_Major, y=frac_Celltype_Major, x=factor(clustID))) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Cluster ID", y = "Fraction") +
  scale_fill_manual(values = ct.color_major[seq_along(unique(clustCT_major$Celltype_Major))]) +
  theme(axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_clusterComp, "OmentumCohort_clusters_celltype_major.png")
ggsave(filename = filename, width = 30, height = 20, dpi = 600, units = "cm")


#################### 
### metadata celltype composition
####################

#########
## COHORT
#########

head(cohortCT_final)
#unsorted:
ggplot(cohortCT_final, aes(fill=Celltype, y=frac_Celltype, x=Cohort)) + 
  geom_bar(position="stack", stat="identity")

#ordering 
unique(cohortCT_final$Cohort)
cohort_order = c("Benign_Omentum", "Omentum_without_Metastasis", "Omentum_with_Metastasis")
factor(cohort_order)

ggplot(cohortCT_final, aes(fill=Celltype, y=frac_Celltype, x=Cohort)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Cohort", y = "Fraction") +
  scale_fill_manual(values = ct.color[seq_along(unique(cohortCT_final$Celltype))]) +
  scale_x_discrete(limits = cohort_order) +
  theme(axis.text.x = element_text(angle = 90,size = 13), #(angle = 75,size = 13,vjust = 0.5)
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_cohort_celltype_final.png")
ggsave(filename = filename, width = 15, height = 25, dpi = 600, units = "cm")


head(cohortCT_major)
#ordering 
unique(cohortCT_major$Cohort)
cohort_order = c("Benign_Omentum", "Omentum_without_Metastasis", "Omentum_with_Metastasis")
factor(cohort_order)

ggplot(cohortCT_major, aes(fill=Celltype_Major, y=frac_Celltype_Major, x=Cohort)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Cohort", y = "Fraction") +
  scale_fill_manual(values = ct.color_major[seq_along(unique(cohortCT_major$Celltype_Major))]) +
  scale_x_discrete(limits = cohort_order) +
  theme(axis.text.x = element_text(angle = 90,size = 13), #(angle = 75,size = 13,vjust = 0.5)
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_cohort_celltype_major.png")
ggsave(filename = filename, width = 15, height = 25, dpi = 600, units = "cm")


#########
## TISSUE REGION
#########

head(regionCT_final)
#unsorted:
ggplot(regionCT_final, aes(fill=Celltype, y=frac_Celltype, x=Tissue_Region)) + 
  geom_bar(position="stack", stat="identity")

#ordering 
unique(regionCT_final$Tissue_Region)
tissue_order = c("benign", "tumor_distant", "tumor_peritumoral","tumor_metastasis")
factor(tissue_order)

ggplot(regionCT_final, aes(fill=Celltype, y=frac_Celltype, x=Tissue_Region)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Tissue Region", y = "Fraction") +
  scale_fill_manual(values = ct.color[seq_along(unique(regionCT_final$Celltype))]) +
  scale_x_discrete(limits = tissue_order) +
  theme(axis.text.x = element_text(angle = 90,size = 13), #(angle = 75,size = 13,vjust = 0.5)
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_tissue_region_celltype_final.png")
ggsave(filename = filename, width = 18, height = 25, dpi = 600, units = "cm")


head(regionCT_major)
ggplot(regionCT_major, aes(fill=Celltype_Major, y=frac_Celltype_Major, x=Tissue_Region)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Tissue Region", y = "Fraction") +
  scale_fill_manual(values = ct.color_major[seq_along(unique(regionCT_major$Celltype_Major))]) +
  scale_x_discrete(limits = tissue_order) +
  theme(axis.text.x = element_text(angle = 90,size = 13), #(angle = 75,size = 13,vjust = 0.5)
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_tissue_region_celltype_major.png")
ggsave(filename = filename, width = 18, height = 25, dpi = 600, units = "cm")


#########
## SITE
#########

head(siteCT_final)
#unsorted:
ggplot(siteCT_final, aes(fill=Celltype, y=frac_Celltype, x=Tissue_Type)) + 
  geom_bar(position="stack", stat="identity")

#ordering 
unique(regionCT_final$Tissue_Region)
site_order = c("Majus", "Minus", "Gastric")
factor(site_order)

ggplot(siteCT_final, aes(fill=Celltype, y=frac_Celltype, x=Tissue_Type)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Site", y = "Fraction") +
  scale_fill_manual(values = ct.color[seq_along(unique(siteCT_final$Celltype))]) +
  scale_x_discrete(limits = site_order) +
  theme(axis.text.x = element_text(size = 13), #(angle = 75,size = 13,vjust = 0.5)
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_site_celltype_final.png")
ggsave(filename = filename, width = 18, height = 25, dpi = 600, units = "cm")


head(siteCT_major)
ggplot(siteCT_major, aes(fill=Celltype_Major, y=frac_Celltype_Major, x=Tissue_Type)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Site", y = "Fraction") +
  scale_fill_manual(values = ct.color_major[seq_along(unique(siteCT_major$Celltype_Major))]) +
  scale_x_discrete(limits = site_order) +
  theme(axis.text.x = element_text(size = 13), #(angle = 75,size = 13,vjust = 0.5)
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_site_celltype_major.png")
ggsave(filename = filename, width = 15, height = 25, dpi = 600, units = "cm")

#########
## Patient
#########

head(patientCT_final)
#unsorted:
ggplot(patientCT_final, aes(fill=Celltype, y=frac_Celltype, x=Patient)) + 
  geom_bar(position="stack", stat="identity")

unique(patientCT_final$Patient)
# order by cohort highlevel
patient_order = c("B345","B349","B427","B480","B340","B419","B479","B468","B436","B438","B453","B454","B486","B497","B500")
length(patient_order)

ggplot(patientCT_final, aes(fill=Celltype, y=frac_Celltype, x=Patient)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Patient", y = "Fraction") +
  scale_fill_manual(values = ct.color[seq_along(unique(patientCT_final$Celltype))]) +
  scale_x_discrete(limits = patient_order) +
  theme(axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_patient_celltype_final.png")
ggsave(filename = filename, width = 35, height = 25, dpi = 600, units = "cm")


head(patientCT_major)
ggplot(patientCT_major, aes(fill=Celltype_Major, y=frac_Celltype_Major, x=Patient)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Patient", y = "Fraction") +
  scale_fill_manual(values = ct.color_major[seq_along(unique(patientCT_major$Celltype_Major))]) +
  scale_x_discrete(limits = patient_order) +
  theme(axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_patient_celltype_major.png")
ggsave(filename = filename, width = 30, height = 25, dpi = 600, units = "cm")

#########
## Samples
#########

head(sampleCT_final)
#unsorted:
ggplot(sampleCT_final, aes(fill=Celltype, y=frac_Celltype, x=Sample)) + 
  geom_bar(position="stack", stat="identity")

length(unique(sampleCT_final$Sample))
# order by cohort highlevel and by tissue region (majus, Minus, Gastric)
sample_order = c("B345_Majus","B345_Minus","B345_Gastric","B349_Majus","B349_Minus","B349_Gastric","B427_Majus","B427_Minus","B427_Gastric","B480_Majus","B340_Majus","B340_Minus","B419_Majus","B419_Minus","B419_Gastric","B479_Majus","B468_Majus_Tumor_BL","B436_Tumor","B436_Next_Majus","B436_Distal_Gastric","B438_Tumor_Next_Majus","B453_Tumor","B453_Next_Majus","B453_Distal_Minus","B454_Tumor_BL","B454_Next_Majus","B454_Distal_Gastric","B486_Tumor","B486_Next_Majus","B486_Distal_Gastric","B497_Tumor","B497_Next_Majus","B497_Distal_Minus","B500_Tumor","B500_Next_Majus","B500_Distal_Minus")
length(sample_order)

ggplot(sampleCT_final, aes(fill=Celltype, y=frac_Celltype, x=Sample)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Sample", y = "Fraction") +
  scale_fill_manual(values = ct.color[seq_along(unique(sampleCT_final$Celltype))]) +
  scale_x_discrete(limits = sample_order) +
  theme(axis.text.x = element_text(angle = 90, size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_sample_celltype_final.png")
ggsave(filename = filename, width = 50, height = 30, dpi = 600, units = "cm")


head(sampleCT_major)
ggplot(sampleCT_major, aes(fill=Celltype_Major, y=frac_Celltype_Major, x=Sample)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x ="Sample", y = "Fraction") +
  scale_fill_manual(values = ct.color_major[seq_along(unique(sampleCT_major$Celltype_Major))]) +
  scale_x_discrete(limits = sample_order) +
  theme(axis.text.x = element_text(angle = 90,size = 13),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.y= element_text(size = 13),
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 18))

filename <- paste0(path_celltypeComp, "OmentumCohort_sample_celltype_major.png")
ggsave(filename = filename, width = 50, height = 30, dpi = 600, units = "cm")

