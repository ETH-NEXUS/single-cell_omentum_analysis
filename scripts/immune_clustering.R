'%&%' = function(a,b) paste(a,b,sep="")
library(pheatmap)
library(dynamicTreeCut)
library(RColorBrewer)

path = "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/single_samples/"
outdir = "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-07_celltype_composition/"

fn_clin = "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/cohort/metadata_omentum_cohort_2021-05.txt"


clin = read.table(fn_clin, head=T, sep="\t")
table(clin$Cohort_Type)
table(clin$Cohort_Type_Highlevel)
table(clin$Site)
table(clin$Tissue_Region)
table(clin$Patient)
table(clin$Diagnosis)

str(clin)
fn = list.files(path %&% "composition/composition_table", 
                pattern = "atypical.*graph_cell.*\\.txt$", full.names = T)
samples = gsub(".*composition_table/(.*).genes_.*", "\\1", fn)
samples
# test
td = read.table(fn[1], header = T)
idx = which(td$Cluster == "Sum")
td[idx, seq_len(ncol(td)-3)]

# load all data
dd_l = vector("list", length(fn))
names(dd_l) = samples
for(ii in seq_along(fn)){
  td = read.table(fn[ii], header = T)
  idx = which(td$Cluster == "Sum")
  td = td[idx, seq(2, ncol(td)-3)]
  td = apply(td, 2, as.numeric)
  td = data.frame(celltype = names(td), counts = as.numeric(td))
  td$sample = samples[ii]
  dd_l[[ii]] = td
}
dd = do.call(rbind, dd_l)
str(dd)
dd = dd[which(dd$sample %in% clin$sample_name), ]

dd_allTypes = dd

dd = dd[which(!(dd$celltype %in% c("uncertain", "unknown"))), ]

dd


# inspect
all_celltypes = unique(dd$celltype)

###### ONLY FOR MAJOR

dd$ct_major = dd$celltype
idx = grep("^B\\.", dd$celltype)
dd$ct_major[idx] = "B.cells"
idx = grep("^T(\\.|c)", dd$celltype)
dd$ct_major[idx] = "T.cells"
idx = grep("Mesothe", dd$celltype)
dd$ct_major[idx] = "Mesothelial.cells"
idx = grep("ndothel", dd$celltype)
dd$ct_major[idx] = "Endothelial.cells"
idx = grep("ibrobl", dd$celltype)
dd$ct_major[idx] = "Fibroblasts"
idx = grep("NK", dd$celltype)
dd$ct_major[idx] = "NK.cells"
idx = grep("yeloi", dd$celltype)
dd$ct_major[idx] = "Myeloid.cells"
idx = grep("HGSOC", dd$celltype)
dd$ct_major[idx] = "Tumor.cells"
table(dd$ct_major) # number of sample with this celltype present

mct = reshape2::acast(dd, ct_major ~ sample, value.var = "counts", fun.aggregate = sum)

################ END ONLY FOR MAJOR
mct = reshape2::acast(dd, celltype ~ sample, value.var = "counts", fun.aggregate = sum)
mct[is.na(mct)] = 0

total = colSums(mct)
mct_norm = sweep(mct, 2, total, "/")
mct_vst = asin(sqrt(mct_norm))

hc = hclust(dist(t(mct_norm)), method = "ward.D2")
dtc = cutreeHybrid(hc, as.matrix(dist(t(mct_vst))))
names(dtc$labels) = colnames(mct_vst)

ann_col = data.frame(cluster = factor("c" %&% (dtc$labels)))
levels(ann_col$cluster)

rownames(ann_col) = colnames(mct_vst)
idx = match(colnames(mct_vst), clin$sample_name)

table(clin$Cohort_Type)
table(clin$Cohort_Type_Highlevel)
table(clin$Site)
table(clin$Tissue_Region)
table(clin$Patient)
table(clin$Diagnosis)

ann_col$Cohort_Type = factor(clin$Cohort_Type[idx])  #1
ann_col$Cohort_Type_Highlevel = factor(clin$Cohort_Type_Highlevel[idx])  #2
ann_col$Site = factor(clin$Site[idx]) #3
ann_col$Tissue_Region = factor(clin$Tissue_Region[idx]) #4
ann_col$Patient = factor(clin$Patient[idx]) #5
ann_col$Diagnosis = factor(clin$Diagnosis[idx]) #6


patient_colours = c(brewer.pal(12, "Paired"),brewer.pal(8, "Accent"))
annot_colors = list(#cluster = c("orangered", "skyblue"),
                    Cohort_Type = brewer.pal(9, "Set1")[seq_along(levels(ann_col$Cohort_Type))],
                    Cohort_Type_Highlevel = brewer.pal(9, "Set1")[seq_along(levels(ann_col$Cohort_Type_Highlevel))],
                    Site = brewer.pal(8, "Dark2")[seq_along(levels(ann_col$Site))],
                    Tissue_Region = brewer.pal(9, "Pastel1")[seq_along(levels(ann_col$Tissue_Region))],
                    Diagnosis = brewer.pal(12, "Paired")[seq_along(levels(ann_col$Diagnosis))],
                    Patient = patient_colours[seq_along(levels(ann_col$Patient))])

#names(annot_colors[[1]]) = levels(ann_col$cluster)
names(annot_colors[[1]]) = levels(ann_col$Cohort_Type)
names(annot_colors[[2]]) = levels(ann_col$Cohort_Type_Highlevel)
names(annot_colors[[3]]) = levels(ann_col$Site)
names(annot_colors[[4]]) = levels(ann_col$Tissue_Region)
names(annot_colors[[5]]) = levels(ann_col$Diagnosis)
names(annot_colors[[6]]) = levels(ann_col$Patient)

drop <- c("cluster")
ann_col_sub = ann_col[,!(names(ann_col) %in% drop)]

hm = pheatmap(mct_norm, scale="none", cluster_rows = F, clustering_method = "ward.D2",
              annotation_col = ann_col_sub, annotation_colors = annot_colors, 
              show_rownames = T, show_colnames = T, fontsize=8, treeheight_col = 70)
ggplot2::ggsave(outdir %&% "celltype_composition_omentum_cohort_hm_scale_none.png", hm$gtable, width = 36, height = 22, units = "cm")


############# NOT REQUIRED
mct_norm_hm = mct_norm
mct_norm_hm[mct_norm_hm > 0.95] 

head(mct_norm_hm)
head(mct_norm)

hm_cutted = pheatmap(mct_norm_hm, scale="none", cluster_rows = F, clustering_method = "ward.D2",
              annotation_col = ann_col_sub, annotation_colors = annot_colors, 
              show_rownames = T, show_colnames = T, fontsize=8, treeheight_col = 70)
ggplot2::ggsave(outdir %&% "celltype_composition_omentum_cohort_hm_cut_scale_none.png", hm_cutted$gtable, width = 36, height = 22, units = "cm")


# write.table(mct_norm, path %&% "allmel_cellfractions.txt", sep="\t", col.names=NA)
