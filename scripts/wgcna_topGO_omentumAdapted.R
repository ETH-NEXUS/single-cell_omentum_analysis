#####
# GO term analysis, adapted from Lars' script wgcna_toGO.R
#####

library(tidyverse)
library(org.Hs.eg.db)
library(topGO) #Depends: BioGenerics, graph, Biobase, GO.db, AnnotationDbi, SparseM
library(Rgraphviz)

#Input data
#Input data geneL: vector with gene names
geneL <- c("ADAR","CEBPB","CASP1","MS4A1","AP2M1",#for demonstration purposes (can be removed)
           "CETN2","TPP1","CASP4","CANX","DAGLA",
           "FMNL1","MRPL49","CAMLG","KYAT1","CCT")

path <- "gfb_omentum_2020/analysis/2021-04_cohort_reanalysis/2021-08_fig1_2_subtypes/GO_term_analysis/test/" #path for output files
ontol <- "BP" #ontology of interest (BP: biological process; MF: molecular function; CC: cellular component)

#Prepare: gene name list and gene/GO mapping
x <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

geneGO2ID <- annFUN.org(whichOnto = ontol, mapping = "org.Hs.eg.db")
geneID2GO <- inverseList(geneGO2ID)
geneNames <- names(geneID2GO)

#Create a named vector: geneList. 
#If gene in the list: 1 else 0
geneL <- geneL[geneL %in% names(xx)]
geneL_ID <- rep(0,length(geneL))
for (i in 1:length(geneL)){
  geneL_ID[i] <- unlist(xx[which(geneL[i] == names(xx))])
}
geneList <- factor(as.integer(geneNames %in% geneL_ID))
names(geneList) <- geneNames

#TopGO analysis
GOdata <- new("topGOdata",
              ontology = ontol,
              allGenes = geneList,
              annot = annFUN.gene2GO,
              gene2GO = geneID2GO)
resultFis <- runTest(GOdata,
                     algorithm = "elim", #"classic" for classical enrichment analysis
                     statistic = "fisher")
resultAdj <- resultFis
resultAdj@score <- p.adjust(resultFis@score,method = "BH")

#GO term table
allRes <- GenTable(GOdata,
                   elim = resultFis,
                   p_adj = resultAdj,
                   topNodes = 20)
write_tsv(allRes,paste(path,"topgo_geneTable.tsv",sep=""))

#Sig nodes plot
png(paste(path,"topgo_SigNodes.png",sep=""),
    width = 25, height = 25, units = "cm",res=300)
showSigOfNodes(GOdata, score(resultFis), firstSigNodes = 5, useInfo = 'all')
dev.off()
