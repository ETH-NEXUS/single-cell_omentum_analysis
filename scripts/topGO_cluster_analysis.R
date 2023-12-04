###
##Analysis of clusters, topGO
###

library(tidyverse)
library(org.Hs.eg.db)
library(topGO) #Depends: BioGenerics, graph, Biobase, GO.db, AnnotationDbi, SparseM
library(Rgraphviz)
library(UpSetR)
library(ggplot2)

path <- "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/GO_term_analysis/" #path for output files
ontol <- "BP" #ontology of interest (BP: biological process; MF: molecular function; CC: cellular component)

#Prepare: gene name list and gene/GO mapping
x <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

geneGO2ID <- annFUN.org(whichOnto = ontol, mapping = "org.Hs.eg.db")
geneID2GO <- inverseList(geneGO2ID)
geneNames <- names(geneID2GO)

# define functions

getCluster_avg_geneNum <- function(subtype_name,cluster_avgs){
  print(subtype_name)
  print(colnames(cluster_avgs))
  for (clust in colnames(cluster_avgs)){   
    if (clust != 'gene'){
      print(clust)
      # extract genes of interest based on expression above average
      avg_clust = mean(cluster_avgs[[clust]])
      print(avg_clust)
      clust_sort = cluster_avgs[order(cluster_avgs[[clust]], decreasing = T), ]
      print(length(subset(clust_sort, (clust_sort[[clust]]) > avg_clust)$gene))
    }
  }
}

analyze_topGO <- function(subtype_name,cluster_avgs,geneNames,geneID2GO,ontol,xx,path){
  print(subtype_name)
  print(colnames(cluster_avgs))
  for (clust in colnames(cluster_avgs)){   
   if (clust != 'gene'){
    print(clust)
    
    # extract genes of interest based on expression above average
    avg_clust = mean(cluster_avgs[[clust]])
    print(avg_clust)
    clust_sort = cluster_avgs[order(cluster_avgs[[clust]], decreasing = T), ]
    print(length(subset(clust_sort, (clust_sort[[clust]]) > avg_clust)$gene))
    clust_genes = subset(clust_sort, (clust_sort[[clust]]) > avg_clust)$gene
    
    #Create a named vector: geneList. 
    #If gene in the list: 1 else 0
    clust_geneL <- clust_genes[clust_genes %in% names(xx)]
    clust_geneL_ID <- rep(0,length(clust_geneL))
    for (i in 1:length(clust_geneL)){
      clust_geneL_ID[i] <- unlist(xx[which(clust_geneL[i] == names(xx))])
    }
    clust_geneList <- factor(as.integer(geneNames %in% clust_geneL_ID))
    names(clust_geneList) <- geneNames
    
    #TopGO analysis
    clust_GOdata <- new("topGOdata",
                         ontology = ontol,
                         allGenes = clust_geneList,
                         annot = annFUN.gene2GO,
                         gene2GO = geneID2GO)
    clust_resultFis <- runTest(clust_GOdata,
                                algorithm = "elim", #"classic" for classical enrichment analysis
                                statistic = "fisher")
    clust_resultAdj <- clust_resultFis
    clust_resultAdj@score <- p.adjust(clust_resultFis@score,method = "BH")
    
    #GO term table
    clust_allRes <- GenTable(clust_GOdata,
                              elim = clust_resultFis,
                              p_adj = clust_resultAdj,
                              topNodes = 20)
    write_tsv(clust_allRes,paste(path,"topgo_geneTable_", subtype_name, "_", clust, ".tsv",sep=""))
    
    #Sig nodes plot
    png(paste(path,"topgo_SigNodes_", subtype_name, "_", clust, ".png",sep=""),
        width = 25, height = 25, units = "cm",res=300)
    showSigOfNodes(clust_GOdata, score(clust_resultFis), firstSigNodes = 5, useInfo = 'all')
    dev.off()
  }
}
}

# read in cluster means, neutrophils
neutro_avg <- as.data.frame(read.csv("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/subtypes_hvg500/neutrophils/Neutrophils.mean_expression_per_subtype_umap_cluster.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
head(neutro_avg)
rownames(neutro_avg) = neutro_avg$gene
colnames(neutro_avg)

analyze_topGO("neutrophils",neutro_avg,geneNames,geneID2GO,ontol,xx,path)
getCluster_avg_geneNum("neutrophils",neutro_avg)

#tumor
tumor_avg <- as.data.frame(read.csv("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/subtypes_hvg500/HGSOC/HGSOC.mean_expression_per_subtype_umap_cluster.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
head(tumor_avg)
rownames(tumor_avg) = tumor_avg$gene
colnames(tumor_avg)

analyze_topGO("HGSOC",tumor_avg,geneNames,geneID2GO,ontol,xx,path)
getCluster_avg_geneNum("HGSOC",tumor_avg)

#fibroblasts
fibro_avg <- as.data.frame(read.csv("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/subtypes_hvg500/fibroblasts/Fibroblasts.mean_expression_per_subtype_umap_cluster.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
head(fibro_avg)
rownames(fibro_avg) = fibro_avg$gene
colnames(fibro_avg)

analyze_topGO("fibroblasts",fibro_avg,geneNames,geneID2GO,ontol,xx,path)
getCluster_avg_geneNum("fibroblasts",fibro_avg)

#Mesothelial
mesothelial_avg <- as.data.frame(read.csv("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/subtypes_hvg500/mesothelial/Mesothelial.cells.mean_expression_per_subtype_umap_cluster.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
head(mesothelial_avg)
rownames(mesothelial_avg) = mesothelial_avg$gene
colnames(mesothelial_avg)

analyze_topGO("Mesothelial",mesothelial_avg,geneNames,geneID2GO,ontol,xx,path)
getCluster_avg_geneNum("Mesothelial",mesothelial_avg)

#Monocytes
monocytes_avg <- as.data.frame(read.csv("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/subtypes_hvg500/monocytes/Monocytes.mean_expression_per_subtype_umap_cluster.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
head(monocytes_avg)
rownames(monocytes_avg) = monocytes_avg$gene
colnames(monocytes_avg)

analyze_topGO("Monocytes",monocytes_avg,geneNames,geneID2GO,ontol,xx,path)
getCluster_avg_geneNum("Monocytes",monocytes_avg)

#T_NK_cells
tcells_avg <- as.data.frame(read.csv("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/subtypes_hvg500/tcells/T_NK_cells.mean_expression_per_subtype_umap_cluster.tsv", sep = "\t", stringsAsFactors = FALSE, header = TRUE))
head(tcells_avg)
rownames(tcells_avg) = tcells_avg$gene
colnames(tcells_avg)

analyze_topGO("T_NK_cells",tcells_avg,geneNames,geneID2GO,ontol,xx,path)
getCluster_avg_geneNum("T_NK_cells",tcells_avg)


#########
# UpSetR plot
# NOTE: function was not possible because plots were not saved
#########

#neutrophils
input_neutro <- read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/GO_term_analysis/neutrophils/clusters_to_sets.txt", sep="\t", header = TRUE)
head(input_neutro)
rownames(input_neutro) <- input_neutro$Sets

  png(paste(path,"upsetr_neutrophils", ".png",sep=""), width = 70, height = 50, res = 300, units="cm")
  upset(data=input_neutro, nsets=100, nintersects=NA, order.by=c("freq"), 
        #intersections = list("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7"),
        matrix.color = 3,        # color of dots and lines
        shade.color = 3,        # shading in table
        shade.alpha = 0.25,        # Transparency of shading in matrix, 0 is no shading
        keep.order = FALSE,
        main.bar.color = "grey20",
        line.size = 0.5,        # lines connecting dots
        point.size = 5,
        mainbar.y.label = "Number of GO terms",
        #mainbar.y.max = "100" #?
        sets.x.label = "Input GO terms per cluster",        # label of sets barplots
        #sets.bar.color = "#56B4E9"
        set_size.show = FALSE,
        #att.color = 5,    #?
        #group.by = "degree", #?
        matrix.dot.alpha = 1,        # transparency of dots not in set, 1 is strongest
        sets.bar.color = "grey38",
        #set_size.scale_max = 1,
        mb.ratio = c(0.7, 0.3),        # ratio between main bar plot and matrix plot
        text.scale = c(5,3,2,2,3,3)
        # Can be a universal scale, or a vector containing individual scales in the following format: 
        # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  )
  dev.off()


#fibroblasts
input_fibro <- read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/GO_term_analysis/fibroblasts/clusters_to_sets.txt", sep="\t", header = TRUE)
head(input_fibro)
rownames(input_fibro) <- input_fibro$Sets

  png(paste(path,"upsetr_fibroblasts", ".png",sep=""), width = 70, height = 50, res = 300, units="cm")
  upset(data=input_fibro, nsets=100, nintersects=NA, order.by=c("freq"), 
        #intersections = list("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7"),
        matrix.color = 3,        # color of dots and lines
        shade.color = 3,        # shading in table
        shade.alpha = 0.25,        # Transparency of shading in matrix, 0 is no shading
        keep.order = FALSE,
        main.bar.color = "grey20",
        line.size = 0.5,        # lines connecting dots
        point.size = 5,
        mainbar.y.label = "Number of GO terms",
        #mainbar.y.max = "100" #?
        sets.x.label = "Input GO terms per cluster",        # label of sets barplots
        #sets.bar.color = "#56B4E9"
        set_size.show = FALSE,
        #att.color = 5,    #?
        #group.by = "degree", #?
        matrix.dot.alpha = 1,        # transparency of dots not in set, 1 is strongest
        sets.bar.color = "grey38",
        #set_size.scale_max = 1,
        mb.ratio = c(0.7, 0.3),        # ratio between main bar plot and matrix plot
        text.scale = c(5,3,2,2,3,3)
        # Can be a universal scale, or a vector containing individual scales in the following format: 
        # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  )
  dev.off()


# HGSOC
input_tumor <- read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/GO_term_analysis/HGSOC/clusters_to_sets.txt", sep="\t", header = TRUE)
head(input_tumor)
rownames(input_tumor) <- input_tumor$Sets

  png(paste(path,"upsetr_HGSOC", ".png",sep=""), width = 70, height = 50, res = 300, units="cm")
  upset(data=input_tumor, nsets=100, nintersects=NA, order.by=c("freq"), 
        #intersections = list("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7"),
        matrix.color = 3,        # color of dots and lines
        shade.color = 3,        # shading in table
        shade.alpha = 0.25,        # Transparency of shading in matrix, 0 is no shading
        keep.order = FALSE,
        main.bar.color = "grey20",
        line.size = 0.5,        # lines connecting dots
        point.size = 5,
        mainbar.y.label = "Number of GO terms",
        #mainbar.y.max = "100" #?
        sets.x.label = "Input GO terms per cluster",        # label of sets barplots
        #sets.bar.color = "#56B4E9"
        set_size.show = FALSE,
        #att.color = 5,    #?
        #group.by = "degree", #?
        matrix.dot.alpha = 1,        # transparency of dots not in set, 1 is strongest
        sets.bar.color = "grey38",
        #set_size.scale_max = 1,
        mb.ratio = c(0.7, 0.3),        # ratio between main bar plot and matrix plot
        text.scale = c(5,3,2,2,3,3)
        # Can be a universal scale, or a vector containing individual scales in the following format: 
        # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  )
  dev.off()


# Mesothelial
input_meso <- read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/GO_term_analysis/mesothelial/clusters_to_sets.txt", sep="\t", header = TRUE)
head(input_meso)
rownames(input_meso) <- input_meso$Sets

  png(paste(path,"upsetr_mesothelial", ".png",sep=""), width = 70, height = 50, res = 300, units="cm")
  upset(data=input_meso, nsets=100, nintersects=NA, order.by=c("freq"), 
        #intersections = list("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7"),
        matrix.color = 3,        # color of dots and lines
        shade.color = 3,        # shading in table
        shade.alpha = 0.25,        # Transparency of shading in matrix, 0 is no shading
        keep.order = FALSE,
        main.bar.color = "grey20",
        line.size = 0.5,        # lines connecting dots
        point.size = 5,
        mainbar.y.label = "Number of GO terms",
        #mainbar.y.max = "100" #?
        sets.x.label = "Input GO terms per cluster",        # label of sets barplots
        #sets.bar.color = "#56B4E9"
        set_size.show = FALSE,
        #att.color = 5,    #?
        #group.by = "degree", #?
        matrix.dot.alpha = 1,        # transparency of dots not in set, 1 is strongest
        sets.bar.color = "grey38",
        #set_size.scale_max = 1,
        mb.ratio = c(0.7, 0.3),        # ratio between main bar plot and matrix plot
        text.scale = c(5,3,2,2,3,3)
        # Can be a universal scale, or a vector containing individual scales in the following format: 
        # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  )
  dev.off()


# Monocytes
input_mono <- read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/GO_term_analysis/monocytes/clusters_to_sets.txt", sep="\t", header = TRUE)
head(input_mono)
rownames(input_mono) <- input_mono$Sets

  png(paste(path,"upsetr_monocytes", ".png",sep=""), width = 70, height = 50, res = 300, units="cm")
  upset(data=input_mono, nsets=100, nintersects=NA, order.by=c("freq"), 
        #intersections = list("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7"),
        matrix.color = 3,        # color of dots and lines
        shade.color = 3,        # shading in table
        shade.alpha = 0.25,        # Transparency of shading in matrix, 0 is no shading
        keep.order = FALSE,
        main.bar.color = "grey20",
        line.size = 0.5,        # lines connecting dots
        point.size = 5,
        mainbar.y.label = "Number of GO terms",
        #mainbar.y.max = "100" #?
        sets.x.label = "Input GO terms per cluster",        # label of sets barplots
        #sets.bar.color = "#56B4E9"
        set_size.show = FALSE,
        #att.color = 5,    #?
        #group.by = "degree", #?
        matrix.dot.alpha = 1,        # transparency of dots not in set, 1 is strongest
        sets.bar.color = "grey38",
        #set_size.scale_max = 1,
        mb.ratio = c(0.7, 0.3),        # ratio between main bar plot and matrix plot
        text.scale = c(5,3,2,2,3,3)
        # Can be a universal scale, or a vector containing individual scales in the following format: 
        # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  )
  dev.off()


# T_NK_cells
input_tcell <- read.table("gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/GO_term_analysis/t_nk_cells/clusters_to_sets.txt", sep="\t", header = TRUE)
head(input_tcell)
rownames(input_tcell) <- input_tcell$Sets

  png(paste(path,"upsetr_tcells", ".png",sep=""), width = 90, height = 70, res = 300, units="cm")
  upset(data=input_tcell, nsets=100, nintersects=NA, order.by=c("freq"), 
        #intersections = list("X0", "X1", "X2", "X3", "X4", "X5", "X6", "X7"),
        matrix.color = 3,        # color of dots and lines
        shade.color = 3,        # shading in table
        shade.alpha = 0.25,        # Transparency of shading in matrix, 0 is no shading
        keep.order = FALSE,
        main.bar.color = "grey20",
        line.size = 0.6,        # lines connecting dots
        point.size = 6,
        mainbar.y.label = "Number of GO terms",
        #mainbar.y.max = "100" #?
        sets.x.label = "Input GO terms per cluster",        # label of sets barplots
        #sets.bar.color = "#56B4E9"
        set_size.show = FALSE,
        #att.color = 5,    #?
        #group.by = "degree", #?
        matrix.dot.alpha = 1,        # transparency of dots not in set, 1 is strongest
        sets.bar.color = "grey38",
        #set_size.scale_max = 1,
        mb.ratio = c(0.7, 0.3),        # ratio between main bar plot and matrix plot
        text.scale = c(8,5,3,3,5,5)
        # Can be a universal scale, or a vector containing individual scales in the following format: 
        # c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  )
  dev.off()






#############
# test
#############

c3_sort = neutro_avg[order(neutro_avg$clust_3, decreasing = T), ]    
head(c3_sort)

median(neutro_avg$clust_3)
min(neutro_avg$clust_3)
max(neutro_avg$clust_3)
avg_c3 = mean(neutro_avg$clust_3)
length(subset(c3_sort, (c3_sort$clust_3) > avg_c3)$gene) # 367 greater than average
neutro_genes = subset(c3_sort, (c3_sort$clust_3) > avg_c3)$gene

#Create a named vector: geneList. 
#If gene in the list: 1 else 0
neutro_geneL <- neutro_genes[neutro_genes %in% names(xx)]
neutro_geneL_ID <- rep(0,length(neutro_geneL))
for (i in 1:length(neutro_geneL)){
  neutro_geneL_ID[i] <- unlist(xx[which(neutro_geneL[i] == names(xx))])
}
neutro_geneList <- factor(as.integer(geneNames %in% neutro_geneL_ID))
names(neutro_geneList) <- geneNames

#TopGO analysis
neutro_GOdata <- new("topGOdata",
              ontology = ontol,
              allGenes = neutro_geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)
neutro_resultFis <- runTest(neutro_GOdata,
                     algorithm = "elim", #"classic" for classical enrichment analysis
                     statistic = "fisher")
neutro_resultAdj <- neutro_resultFis
neutro_resultAdj@score <- p.adjust(neutro_resultFis@score,method = "BH")

#GO term table
neutro_allRes <- GenTable(neutro_GOdata,
                   elim = neutro_resultFis,
                   p_adj = neutro_resultAdj,
                   topNodes = 20)
write_tsv(neutro_allRes,paste(path,"topgo_geneTable_neutrophils_clust3.tsv",sep=""))

#Sig nodes plot
png(paste(path,"topgo_SigNodes_neutrophils_clust3.png",sep=""),
    width = 25, height = 25, units = "cm",res=300)
showSigOfNodes(neutro_GOdata, score(neutro_resultFis), firstSigNodes = 5, useInfo = 'all')
dev.off()

