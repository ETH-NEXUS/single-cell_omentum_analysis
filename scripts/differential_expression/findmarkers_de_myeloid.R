###################################################
## File name: findmarkers_de_myeloid.R
## Author: Matteo Carrara
###################################################

## GENERAL:
## This script performs a full DE analysis of myeloid cells using findMarkers
library(Seurat)
library(ggplot2)
library(ggrepel)


tmp=readRDS("omentum_myeloid_.subtype_sce.RDS")
seurat_integrated=as.Seurat(tmp, counts="counts", data="normcounts")
Idents(seurat_integrated)=seurat_integrated$umap_subtype_cl

tmp=unique(Idents(seurat_integrated))

markers_gex_pairwise=list()

aaa=unique(seurat_integrated$umap_subtype_cl)

for(k in 1:length(aaa)){
  message(paste("Working on cluster:", aaa[k]))
  for(j in 1:length(aaa)){
    if(k==j){next()}
    message(paste0("contrast ", aaa[k], " vs ", aaa[j]))
    markers_gex_pairwise=FindMarkers(object=seurat_integrated, ident.1=aaa[k], ident.2=aaa[j],  assay="originalexp", logfc.threshold=0, min.cells.group=1,min.cells.feature=1,min.pct=0,only.pos=FALSE)
    tmpname = paste0("cl_", aaa[k], "_vs_cl_", aaa[j], ".RDS")
    saveRDS(markers_gex_pairwise, tmpname)
}
}

myfiles=list.files()
myfiles=grep("_vs_", myfiles,v=T)

for(i in 1:length(myfiles)){
  mypvalthr=0.01
  myfcthr=1.5
  mygenes=readRDS(paste0(myfiles[i]))
  mygenes$status="not"
  mygenes$status[(mygenes$avg_log2FC > myfcthr | mygenes$avg_log2FC < -myfcthr) & mygenes$p_val < mypvalthr] = "signif&strong"
  mygenes$status[(mygenes$avg_log2FC > -myfcthr & mygenes$avg_log2FC < myfcthr) & mygenes$p_val < mypvalthr] = "signif"
  mygenes$status[(mygenes$avg_log2FC > myfcthr | mygenes$avg_log2FC < -myfcthr) & mygenes$p_val > mypvalthr] = "strong"
  mygenes$label=NA
  mygenes$label[mygenes$status=="signif&strong"]=rownames(mygenes)[mygenes$status=="signif&strong"]
  p = ggplot(data=mygenes, aes(x=avg_log2FC, y=-log10(p_val), col=status, label=label)) +geom_point() +theme_minimal() +
        geom_hline(yintercept=-log10(mypvalthr)) +
        geom_vline(xintercept=c(myfcthr, -myfcthr)) +
        scale_color_manual(values=c("black","blue", "red","darkgreen"))
  p
  ggsave(paste0("volcano/",myfiles[i], "_volcano.png"))
}

myfiles_gex = myfiles
mygenes_group_gex=NULL
for(i in 1:length(myfiles_gex)){
  mygenes=readRDS(paste0(myfiles_gex[i]))
  mygenes_group_gex=c(mygenes_group_gex,rownames(mygenes)[1:10])
}
mygenes_group_gex=unique(mygenes_group_gex)
myexpression_gex =AverageExpression(seurat_integrated, assays="originalexp", features=mygenes_group_gex,slot="counts")
myexpression_gex=myexpression_gex$originalexp
myexpression_gex_n =AverageExpression(seurat_integrated, assays="originalexp", features=mygenes_group_gex,slot="data")
myexpression_gex_n=myexpression_gex_n$originalexp
library(RColorBrewer)
hmcol=rev(colorRampPalette(brewer.pal(8,"RdBu"))(127))
this_annot=data.frame(umap_subtype_cl=rep(NA, ncol(myexpression_gex)))
rownames(this_annot)= colnames(myexpression_gex)
this_annot$umap_subtype_cl= unlist(lapply(strsplit(rownames(this_annot), " "), function(x)x[1]))


annot_colors=list(
  umap_subtype_cl=c(
"1"= "#f8766d",
"2"= "#e38900",
"3"= "#c49a00",
"4"= "#99a802",
"5"= "#53b400",
"6"= "#00bc57",
"7"= "#02c094",
"8"= "#00bfc4",
"9"= "#01b6eb",
"10"="#03a5ff",
"11"="#a58aff",
"12"="#df70f8",
"13"="#fb61d7",
"14"="#ff66a8"
  )
)


library(pheatmap)
library(RColorBrewer)
library(viridis)
library(ggplot2)

phm <- pheatmap(t(myexpression_gex), color = hmcol, scale = "column", fontsize = 12,
                             show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                             annotation_row = this_annot, annotation_colors = annot_colors,
                             main = "Top differentially expressed GEX average raw expression",
                             fontsize_col=8)
              ggsave("Myeloid_DE_heatmap_gex_top10_noselection.pdf", phm$gtable,
                     width = 70, height = 20, units = "cm")

phm <- pheatmap(t(myexpression_gex_n), color = hmcol, scale = "column", fontsize = 12,
                             show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                             annotation_row = this_annot, annotation_colors = annot_colors,
                             main = "Top differentially expressed GEX average normalized expression",
                             fontsize_col=8)
              ggsave("Myeloid_DE_heatmap_gex_top10_noselection_normcounts.pdf", phm$gtable,
                     width = 70, height = 20, units = "cm")

save.image(file="session.RData")
#selection of genes that are unique to the cluster
#we select significant genes that appear DE in at least 10%, 20%, 50%, 75%, 100% of the pairwise contrasts for that cluster
#For each cluster, load in a list, separately, all contrasts and filter for sig and strong
#For each filtered gene, check in how many filtered DE tables it appears
#if it passes the threshold, add it to the list
#save the final list
#Create an hm with the top genes of each final list

myclusters=tmp

find_gene_presence <- function(x, mytmpgene){
  tmpnum = length(which(rownames(x)==mytmpgene))
  return(tmpnum)
}

gene=list()
myrank=NULL
for(i in 1:length(myclusters)){
 mycontrasts=list.files("pairwise")
 mycontrasts=grep(paste0("cl_", myclusters[i], "_vs_"), mycontrasts, v=T)

 contrast_list = NULL
 gene_presence = NULL
 for(j in 1:length(mycontrasts)){
  contrast_list[[j]]=readRDS(paste0("pairwise/",mycontrasts[j]))
 }

 sig = lapply(contrast_list, function(x){x[x$p_val_adj <1e-39,]})
 sig_str = lapply(sig, function(x){x[abs(x$avg_log2FC) > 2 ,]})
 gene_names = unique(unlist(lapply(sig_str, rownames)))
 for(j in 1:length(gene_names)){
  mytmpgene=gene_names[j]
  tmpgenenum=sum(unlist(lapply(sig_str, find_gene_presence, mytmpgene)))
  tmpgenefrac=tmpgenenum/length(myclusters)
  gene_presence = rbind(gene_presence, data.frame(name=mytmpgene,num_contrasts=tmpgenenum, frac_contrasts=tmpgenefrac))
 }
 write.csv(gene_presence, paste0("gene_presence_clust_", myclusters[i], ".csv"))

 gene[[i]] = gene_presence[which(gene_presence$frac_contrasts>0.859),]


 message(nrow(gene_presence[which(gene_presence$frac_contrasts>0.85),]))


}

mygenes_group_gex=unique(unlist(lapply(gene,function(x)x[,1])))
length(mygenes_group_gex)
#251

myexpression_gex =AverageExpression(seurat_integrated, assays="originalexp", features=mygenes_group_gex,slot="counts")
myexpression_gex=myexpression_gex$originalexp
myexpression_gex_n =AverageExpression(seurat_integrated, assays="originalexp", features=mygenes_group_gex,slot="data")
myexpression_gex_n=myexpression_gex_n$originalexp
library(RColorBrewer)
hmcol=rev(colorRampPalette(brewer.pal(8,"RdBu"))(127))

this_annot=data.frame(umap_subtype_cl=rep(NA, ncol(myexpression_gex)))
rownames(this_annot)= colnames(myexpression_gex)
this_annot$umap_subtype_cl= unlist(lapply(strsplit(rownames(this_annot), " "), function(x)x[1]))

annot_colors=list(
  umap_subtype_cl=c(
  "1"= "#f8766d",
"2"= "#e38900",
"3"= "#c49a00",
"4"= "#99a802",
"5"= "#53b400",
"6"= "#00bc57",
"7"= "#02c094",
"8"= "#00bfc4",
"9"= "#01b6eb",
"10"="#03a5ff",
"11"="#a58aff",
"12"="#df70f8",
"13"="#fb61d7",
"14"="#ff66a8"
)
)


library(pheatmap)
library(RColorBrewer)
library(viridis)
library(dendsort)
library(ggplot2)

phm <- pheatmap(t(myexpression_gex), color = hmcol, scale = "column", fontsize = 12,
                             show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                             annotation_row = this_annot, annotation_colors = annot_colors,
                             main = "Top differentially expressed GEX average raw expression",
                             fontsize_col = 4)
              ggsave("DE_heatmap_gex_085_select.pdf", phm$gtable,
                     width = 60, height = 20, units = "cm")

phm <- pheatmap(t(myexpression_gex_n), color = hmcol, scale = "column", fontsize = 12,
                             show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                             annotation_row = this_annot, annotation_colors = annot_colors,
                             main = "Top differentially expressed GEX average normalized expression",
                             fontsize_col = 4)
              ggsave("DE_heatmap_gex_085_select_normcounts.pdf", phm$gtable,
                     width = 60, height = 20, units = "cm")

####
phm <- pheatmap(t(myexpression_gex), color = hmcol, scale = "row", fontsize = 12,
                             show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                             annotation_row = this_annot, annotation_colors = annot_colors,
                             main = "Top differentially expressed GEX average raw expression",
                             fontsize_col = 4)
              ggsave("DE_heatmap_gex_085_select_rowscale.pdf", phm$gtable,
                     width = 60, height = 20, units = "cm")

phm <- pheatmap(t(myexpression_gex_n), color = hmcol, scale = "row", fontsize = 12,
                             show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                             annotation_row = this_annot, annotation_colors = annot_colors,
                             main = "Top differentially expressed GEX average normalized expression",
                             fontsize_col = 4)
              ggsave("DE_heatmap_gex_085_select_normcounts_rowscale.pdf", phm$gtable,
                     width = 60, height = 20, units = "cm")

phm <- pheatmap(t(myexpression_gex), color = hmcol, scale = "none", fontsize = 12,
                             show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                             annotation_row = this_annot, annotation_colors = annot_colors,
                             main = "Top differentially expressed GEX average raw expression",
                             fontsize_col = 4)
              ggsave("DE_heatmap_gex_085_select_nonescale.pdf", phm$gtable,
                     width = 60, height = 20, units = "cm")  

phm <- pheatmap(t(myexpression_gex_n), color = hmcol, scale = "none", fontsize = 12,
                             show_colnames = TRUE, show_rownames = FALSE, annotation_names_col = TRUE,
                             annotation_row = this_annot, annotation_colors = annot_colors,
                             main = "Top differentially expressed GEX average normalized expression",
                             fontsize_col = 4)
              ggsave("DE_heatmap_gex_085_select_normcounts_nonescale.pdf", phm$gtable,
                     width = 60, height = 20, units = "cm")








