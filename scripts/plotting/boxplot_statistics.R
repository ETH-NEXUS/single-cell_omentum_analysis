library(zeallot)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(dplyr)
library(MASS)
library(multcomp)
library(broom)
library(ggpubr)
theme_set(theme_bw(base_size=10))
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}

output_path="./"

# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
theme_set(theme_bw(base_size=10))

tmp=readRDS("OmentumCohort.merged_cohort_nn30_mD03.RDS")
tmp_l=unique(tmp$sampleID)

tmp_list = list()
aaa=list.files("single_samples-full_object")
aaa=grep("_complete.RDS", aaa, v=T)
ttt = gsub("_complete.RDS", "",aaa)
j=1
for(i in 1:length(aaa)){
 if(!ttt[i]%in%tmp_l){next()}
 tmp_list[[j]]<-readRDS(paste0("single_samples-full_object/",aaa[i]))
 print(paste(i, "of", length(aaa)))
 names(tmp_list)[j]=ttt[i]
 j=j+1
}

#Cluster-based celltypes
colData(tmp)$celltype <- "PLACEHOLDER"
tmp$celltype[tmp$umap_cohort_cl == "14"] = "HGSOC"#
tmp$celltype[tmp$umap_cohort_cl %in% c("1","4","13","18")] = "Mesothelial"#
tmp$celltype[tmp$umap_cohort_cl %in% c("2","8")] = "Mesenchymal.progenitor"#
tmp$celltype[tmp$umap_cohort_cl == "24"] = "Transdifferentiating"#
tmp$celltype[tmp$umap_cohort_cl == "9"] = "Endothelial"#
tmp$celltype[tmp$umap_cohort_cl == "6"] = "Lymphatic.endothelial"#
tmp$celltype[tmp$umap_cohort_cl == "15"] = "Monocytes.macrophages"
tmp$celltype[tmp$umap_cohort_cl == "5"] = "Neutrophils"#
tmp$celltype[tmp$umap_cohort_cl == "10"] = "Plasmacytoid.dendritic"#
tmp$celltype[tmp$umap_cohort_cl == "25"] = "Mast.cells"#
tmp$celltype[tmp$umap_cohort_cl %in% c("3","11","12","17","19","20","22","23")] = "T.NK"#
tmp$celltype[tmp$umap_cohort_cl == "16"] = "B.cells"#
tmp$celltype[tmp$umap_cohort_cl == "21"] = "Plasma"#
tmp$celltype[which(tmp$celltype == "PLACEHOLDER")] = "Noise"#

myct=table(tmp$celltype)

f_make_annot = function(test_vars){
  annot = as.data.frame(d_clin[, test_vars])
  rownames(annot) = d_clin$tupro_id
  for(mm in seq_len(ncol(annot))){
    annot[, mm] = factor(annot[, mm])
  }
  annot_colors = vector("list", length(test_vars))
  names(annot_colors) = test_vars
  for(nn in seq_along(annot_colors)){
    # nn = 1
    this_pal =  brewer.pal(palette_set$pal_length[nn], 
                           palette_set$pal_name[nn])[seq_along(levels(annot[, nn]))]
    annot_colors[[test_vars[nn]]] = this_pal
    names(annot_colors[[nn]]) = levels(annot[, nn])
  }
  return(list(annot, annot_colors))
}
palette_set = data.frame(pal_name = c("Paired", "Dark2", "Set1", "Accent", 
                                      "Set3", "Dark2", "Paired", "Accent"),
                         pal_length = c(12, 8, 8, 8, 12, 8, 12, 8))

##get the single samples
d_clin <- read.delim("metadata_omentum_cohort.txt")
d_clin <- d_clin[which(d_clin$Sample_ID%in%tmp_l),]
colnames(d_clin)[which(colnames(d_clin)=="Site")]<-"Tissue_Type"
colnames(d_clin)[which(colnames(d_clin)=="Sample_ID")]<-"sample_name_short"

d_clin$tupro_id=d_clin$sample_name_short
d_clin$Cohort<-d_clin$Cohort_Type_Highlevel
d_clin$Cohort[which(d_clin$Cohort%in%c("Benign_Omentum","Omentum_without_Metastasis"))]="non_metastatic_omentum"
d_clin$Cohort[which(d_clin$Cohort=="Omentum_with_Metastasis")]="metastatic_omentum"

test_vars=c("Tissue_Region")
c(annot, annot_colors) %<-% f_make_annot(test_vars)
colnames(annot)[1]<-"Tissue_Region"
# Correct the colors with the new colors defined in the color table
annot_colors$Tissue_Region = c("benign"="#1E88E5","tumor_distant"="#004D40","tumor_peritumoral"="#FFC107", "tumor_metastasis"="#D81B60")
annot_colors$Cohort_Type_Highlevel = c("Benign_Omentum"="#648FFF","Omentum_with_Metastasis"="#FE3000","Omentum_without_Metastasis"="#FFB000")
annot_colors$Cohort = c("metastatic_omentum"="#4D0092","non_metastatic_omentum"="#1AFF1A")
# Calculate the celltype numbers for each sample. cann the table mct
ct=unique(tmp$celltype)
mct=NULL
for(i in 1:length(tmp_list)){
aaa=table(tmp_list[[i]]$celltype)
bbb=rep(0, length(ct))
names(bbb)=as.character(ct)
for(j in 1:length(aaa)){
ccc = which(names(bbb)==names(aaa)[j])
bbb[ccc]=aaa[j]
}
mct=rbind(mct, bbb)
}

rownames(mct)=names(tmp_list)


annot$tupro_id = rownames(annot)

id_match = match(rownames(mct), annot$tupro_id)
annot = annot[id_match[!is.na(id_match)], ]
mct = mct[as.character(annot$tupro_id), ]
annot$total = rowSums(mct)
df = melt(mct, varnames = c("tupro_id", "celltype"), value.name = "cell_nr")
df = left_join(df, annot)
df$prop = df$cell_nr / df$total
df$logtotal = log(df$total)

############################
## Negative binomial models
## Global test: for each variable, is there an association 
## between any level of the variable and any cluster?

for(ii in seq_along(test_vars)){
  rfit = data.frame(test_var = test_vars[ii], pval=1)
  # groupvar = test_vars[1]
  groupvar = test_vars[ii]
  td = df[!is.na(df[, groupvar]), ]
  min_pats = 3*length(unique(td$celltype))
  tab = table(td[, groupvar], useNA="ifany")
  tab = tab[tab >=min_pats]
  if(length(tab) >= 2){
    td = td[td[, groupvar] %in% names(tab), ]
    formel = "cell_nr ~ tupro_id + logtotal + celltype*" %&% groupvar
    tfit = MASS::glm.nb(as.formula(formel), data=td, control=glm.control(maxit=100))
    formel = "cell_nr ~ tupro_id + logtotal + celltype + " %&% groupvar
    tfit0 = MASS::glm.nb(as.formula(formel), data=td, control=glm.control(maxit=100))
    rfit$pval = anova(tfit, tfit0)$"Pr(Chi)"[2]
    td[,groupvar] <- factor(td[,groupvar], levels=c("benign","tumor_distant","tumor_peritumoral","tumor_metastasis"))
    ggplot(td, aes(x=td[, groupvar], y=(cell_nr)/total, fill = td[, groupvar])) +
      geom_boxplot() + geom_jitter(width=0.2, show.legend = F) +
      xlab("") + ylab("Proportion") + scale_y_log10() + 
      theme(axis.text.x = element_blank()) + 
      scale_fill_manual(values = annot_colors[[groupvar]], name="") +
      facet_wrap(~celltype, nrow=1)
    ggsave(output_path %&% "NB__" %&% groupvar %&% "__celltype_boxplot_log_" %&% test_vars[ii] %&% ".pdf", 
           width = 50, height = 25, units = "cm", dpi = 300)
	rfit$padj = p.adjust(rfit$pval, "BH")
	fname = output_path %&% "NB__celltype_condi__glm_global_LRT_" %&% test_vars[ii] %&% ".txt"
	write.table(rfit[order(rfit$padj), ], fname, sep="\t", row.names = F)
	fname = output_path %&% "NB__" %&% groupvar %&% "__celltype_boxplot_" %&% test_vars[ii] %&% "_table.tsv"
	write.table(td, fname, row.names=F, sep="\t")
  }
}

## Global test: for each cluster, is there an association 
## with any level of the variable?
allfits <- grep("NB__celltype_condi__glm_global_LRT_", list.files(), v=T)
for(ii in seq_len(length(allfits))){
rfit = read.delim(allfits[ii])
rfit = rfit[rfit$padj < 0.05, ]
if(nrow(rfit)==0){next()}
celltypes = unique(df$celltype)
dres = data.frame(NULL)
  this_var = rfit$test_var#[ii]
  td = df[!is.na(df[, this_var]), ]
  min_pats = 3*length(unique(td$celltype))
  tab = table(td[, this_var], useNA="ifany")
  tab = tab[tab >=min_pats]
  var_levels = names(tab)
  td = td[td[, this_var] %in% var_levels, ]
  td[, this_var] = factor(td[, this_var])
  if(length(tab) >= 2){
    td$condi = td[, this_var]
    glm_lrt = function(X) {
      tfit = try(glm.nb(cell_nr ~ logtotal + condi, data=X, control=glm.control(maxit=100)), silent = T)
      tfit0 = try(glm.nb(cell_nr ~ logtotal, data=X, init.theta = tfit$theta, control=glm.control(maxit=100)), silent = T)
      if(class(tfit)[1] == "negbin" & class(tfit0)[1] == "negbin"){
        return(anova(tfit, tfit0)$"Pr(Chi)"[2])
      } else {
        return(1)
      }
    }
    these_cts = unique(td$celltype)
    tmp = vector("list", length(these_cts))
    names(tmp) = these_cts
    for(jj in seq_along(these_cts)){
      ttd = td[td$celltype == these_cts[jj], c("cell_nr", "logtotal", "condi")]
      tmp[[jj]] = glm_lrt(ttd)
    }
    clu_ktest_pval = melt(tmp)
    names(clu_ktest_pval) = c("pval", "celltype")
    clu_ktest_pval$variable = this_var
    clu_ktest_pval$comparison = "NegBin_glm"
    tres = clu_ktest_pval[, c("celltype", "variable", "comparison", "pval")]
    dres = rbind(dres, tres)
dres$pval[dres$pval == 0] = min(dres$pval[dres$pval > 0])
mres = acast(dres, celltype~variable, value.var = "pval")
dres$padj = p.adjust(dres$pval, "BH")
fname = output_path %&% "NB__celltype_condi__NegBin_LRT_" %&% this_var %&% ".txt"
write.table(dres[order(dres$padj), ], fname, sep="\t", row.names = F)
  }
}

## all pairwise comparisons
alldres <- grep("NB__celltype_condi__NegBin_LRT_", list.files(), v=T)
for(ii in 1:length(alldres)){
dout = data.frame(NULL)
dres=read.delim(alldres[ii])
dres = dres[dres$padj < 0.05, ]
if(nrow(dres)==0){next()}
for(jj in seq_len(nrow(dres))){
  this_ct = dres$celltype[jj]
  this_var = dres$variable[jj]
  td = df[!is.na(df[, this_var]), ]
  td = td[td$celltype == this_ct, ]
  min_pats = 3
  tab = table(td[, this_var], useNA="ifany")
  tab = tab[tab >= min_pats]
  var_levels = names(tab)
  td = td[td[, this_var] %in% var_levels, ]
  td[, this_var] = factor(td[, this_var])
  if(length(tab) >= 2){
    td$condi = factor(td[, this_var])
    td$condi = factor(gsub("-", "__", td$condi))
    tapply(td$prop, td$condi, mean)
    tapply(td$prop, td$condi, sd)
    td$logtotal = log(td$total)
    tfit = glm.nb(cell_nr ~ logtotal + condi, data=td)
    summary(tfit)$coef
    tmp = glht(tfit, linfct = mcp(condi = "Tukey"))
    tmp = summary(tmp, test = adjusted("none"))
    tmp = as.data.frame(tidy(tmp))
    tmp$term = this_var
    tmp$comparison = "pairwise_NB_glm"
    tmp$celltype = this_ct
    dout = rbind(dout, tmp[, -3])
}}
dout$p.value[dout$p.value == 0] = min(dout$p.value[dout$p.value > 0])
dout$padj = p.adjust(dout$p.value, "BH")
fname = output_path %&% "NB__celltype_condi_NegBin_pairwise_" %&% this_var %&% ".txt"
write.table(dout[order(dout$padj), ], fname, sep="\t", row.names = F)
  
}

doutfiles <- grep("pairwise", list.files(), v=T)
groups1 = list(c("benign", "benign","benign","tumor_distant","tumor_distant","tumor_peritumoral"), c("Benign_Omentum","Benign_Omentum","Omentum_without_Metastasis"), c("non_metastatic_omentum"))
groups2 = list(c("tumor_distant","tumor_peritumoral","tumor_metastasis","tumor_peritumoral","tumor_metastasis","tumor_metastasis"), c("Omentum_without_Metastasis","Omentum_with_Metastasis","Omentum_with_Metastasis"), c("metastatic_omentum"))
repnum=list(c(6),c(6),c(6),c(6))
ypos=list(c(rep(c(1,1.6,2.2, 2.8),6)), c(1), c(rep(c(1.2),36)))

for(ii in seq_along(test_vars)){
  mydout <- grep(test_vars[ii] %&% ".txt", doutfiles, v=T)
  if(length(mydout)==0){next()}
  mydout <- read.delim(mydout)
  mydout <- mydout[which(mydout$padj<0.05),]
  mydoutgroup1 <- unlist(lapply(strsplit(mydout$contrast," - "), function(x)x[2]))
  mydoutgroup2 <- unlist(lapply(strsplit(mydout$contrast," - "), function(x)x[1]))
  mydoutuid <- paste0(mydout$celltype,mydoutgroup1,mydoutgroup2)
  mydoutuid2 <- paste0(mydout$celltype,mydoutgroup2,mydoutgroup1)
  rfit = data.frame(test_var = test_vars[ii], pval=1)
  groupvar = test_vars[ii]
  td = df[!is.na(df[, groupvar]), ]
  min_pats = 3*length(unique(td$celltype))
  tab = table(td[, groupvar], useNA="ifany")
  tab = tab[tab >=min_pats]
  if(length(tab) >= 2){
    td = td[td[, groupvar] %in% names(tab), ]
    formel = "cell_nr ~ tupro_id + logtotal + celltype*" %&% groupvar
    tfit = MASS::glm.nb(as.formula(formel), data=td, control=glm.control(maxit=100))
    formel = "cell_nr ~ tupro_id + logtotal + celltype + " %&% groupvar
    tfit0 = MASS::glm.nb(as.formula(formel), data=td, control=glm.control(maxit=100))
    rfit$pval = anova(tfit, tfit0)$"Pr(Chi)"[2]
  annotation_df <- data.frame(
    celltype = td$celltype, y="Proportion",
    group1 = groups1[[ii]], group2=groups2[[ii]], p=1, p.adj=1, p.format=1, p.adj.signif="", y.position=NA)
  anndfuid <- paste0(annotation_df$celltype, annotation_df$group1, annotation_df$group2)
  for(kk in 1:length(mydoutuid)){
    tmp=which(anndfuid %in% mydoutuid[kk])
    if(length(tmp)==0){
      tmp=which(anndfuid %in% mydoutuid2[kk])
    }
    annotation_df$p.adj.signif[tmp] <- format(mydout$padj[kk],scientific=T,digits=2)
  }
  annoloc_4pos=which(annotation_df$p.adj.signif!="")
  annotation_df$y.position[annoloc_4pos] <- ypos[[ii]]
   
    ggplot(td, aes(x=td[, groupvar], y=(cell_nr+1)/total, fill = td[, groupvar])) +
      geom_boxplot() + geom_jitter(width=0.2, show.legend = F) +
      xlab("") + ylab("Proportion") + scale_y_log10() +
      theme(axis.text.x = element_blank()) +
      scale_fill_manual(values = annot_colors[[groupvar]], name="") +
      facet_wrap(~celltype, nrow=1) + stat_pvalue_manual(annotation_df, label="p.adj.signif")
    ggsave(output_path %&% "NB__" %&% groupvar %&% "__celltype_boxplot_" %&% test_vars[ii] %&% "significance.pdf",
           width = 30, height = 25, units = "cm", dpi = 300)
  }
}

