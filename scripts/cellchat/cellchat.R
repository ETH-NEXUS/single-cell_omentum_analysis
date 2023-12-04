###################################################
## File name: cellchat.R
## Author: Matteo Carrara
###################################################

## GENERAL:
## This script performs a full cell-cell-communication analysis using cellchat
library(SingleCellExperiment)
library(optparse)
library(NMF)
library(ggalluvial)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(ggplot2)
library(yaml)

cat("\n\n\nPrint sessionInfo():\n\n")
print(sessionInfo())

option_list <- list(
  make_option("--rdsFile", type = "character", help = "RDS file with merged single cell data."),
  make_option("--subset_var", type = "character", default="none", help = "the variable to subset on, if any."),
  make_option("--subset_values", type="character", default="none", help="the values of the variable to subset on, separated by comma"),
  make_option("--group_by", type="character", help="the variable to group by. Usually the celltype"),
  make_option("--color_yaml", type="character", help="the yaml file with the necessary color configuration"),
  make_option("--pathways", type="character", default="none", help="pathways of interest to plot, separated by comma"),
  make_option("--n_patterns_out", type="numeric", default=3,help="Number of patterns to use for global OUTGOING communication patterns. Defaults to 3, but can be adapted based on the signaling_score_outgoing plot produced in the first run"),
  make_option("--n_patterns_in", type="numeric", default=3,help="Number of patterns to use for global INCOMING communication patterns. Defaults to 3, but can be adapted based on the signaling_score_incoming plot produced in the first run"),
  make_option("--run_cluster_infer", type="character", default="no", help="Run the cluster number infer step. This is skipped by default as it should be ran beforehand to set n_patterns_in and n_patterns_out properly"),
  make_option("--min_cells", type="numeric", default=10, help="Mininum number of cells in a group to consider the group for the analysis"),
  make_option("--base_outdir", type="character", default="output", help="Name of the base output directory where all results are stored")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if((opt$run_cluster_infer) == "yes"){
  message("run_cluster_infer is set. The script will plot the results under output/latent_patterns/qc and quit. Please review the plots, decide the values to assign to n_patterns_out and n_patterns_in accordingly and re-run the procedure without run_cluster_infer")
}

"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""

message("Reading the main object. This may take some time...")
combined <- readRDS(opt$rdsFile)
str(combined)

# Rename the tissue regions to be more user-friendly
combined$Tissue_Region[combined$Tissue_Region=="Omentum"] <- "benign"
combined$Tissue_Region[combined$Tissue_Region=="Tumor"] <- "tumor_metastasis"
combined$Tissue_Region[combined$Tissue_Region=="Next"] <- "tumor_peritumoral"
combined$Tissue_Region[combined$Tissue_Region=="Distal"] <- "tumor_distant"
combined$Tissue_Region <- factor(combined$Tissue_Region, levels=c("benign", "tumor_distant", "tumor_peritumoral","tumor_metastasis"))

# Assign the cluster-based celltypes
combined$celltype="PLACEHOLDER"
combined$celltype[combined$umap_cohort_cl == "14"] = "HGSOC"#
combined$celltype[combined$umap_cohort_cl %in% c("1","4","13","18")] = "Mesothelial"#
combined$celltype[combined$umap_cohort_cl %in% c("2","8")] = "Mesenchymal.progenitor"#
combined$celltype[combined$umap_cohort_cl == "24"] = "Transdifferentiating"#
combined$celltype[combined$umap_cohort_cl == "9"] = "Endothelial"#
combined$celltype[combined$umap_cohort_cl == "6"] = "Lymphatic.endothelial"#
combined$celltype[combined$umap_cohort_cl == "15"] = "Monocytes.macrophages"
combined$celltype[combined$umap_cohort_cl == "5"] = "Neutrophils"#
combined$celltype[combined$umap_cohort_cl == "10"] = "Plasmacytoid.dendritic"#
combined$celltype[combined$umap_cohort_cl == "25"] = "Mast.cells"#
combined$celltype[combined$umap_cohort_cl %in% c("3","11","12","17","19","20","22","23")] = "T.NK"#
combined$celltype[combined$umap_cohort_cl == "16"] = "B.cells"#
combined$celltype[combined$umap_cohort_cl == "21"] = "Plasma"#
combined$celltype[which(combined$celltype == "PLACEHOLDER")] = "Noise"#

#We are not interested in the noise here
combined=combined[,which(combined$celltype!="Noise")]

message("Squashing negative normalized counts to 0 for use in cellchat")
normcounts(combined)[which(normcounts(combined)<0)] <- 0

message("Reading the color config")
myallcolors <- read_yaml(opt$color_yaml)
mycelltypecols <- unlist(myallcolors[[which(names(myallcolors)==opt$group_by)]])

all.data.input <- normcounts(combined)
all.meta <- as.data.frame(colData(combined))
rownames(all.meta) <- combined$barcodes

message("creating output dir")
dir.create(opt$base_outdir)

if(opt$subset_var != "none"){
  meta_col_subset <- which(colnames(all.meta)==opt$subset_var)
  if(opt$subset_values=="none"){
    stop("Error: provided a subset_var but no subset_values")
  }
  meta_values <- unlist((strsplit(opt$subset_values,",")))
  values_forname <- paste(meta_values, collapse="-")
}else{
  meta_col_subset <- "all"
  values_forname <- "all"

}

if(opt$pathways!="none"){
  message("will analyse manual pathways")
pathways <- unlist((strsplit(opt$pathways,",")))
}else{
  message("will analyse only top pathways")
}

dir.create(opt$base_outdir%&%"/"%&%opt$subset_var%&% "__" %&% values_forname)
outdir=opt$base_outdir%&%"/"%&%opt$subset_var%&% "__" %&% values_forname %&% "/"
message(outdir)

message("subsetting the object following user-provided options")
cell.use <- rownames(all.meta)[all.meta[,meta_col_subset] %in% meta_values]
data.input <- all.data.input[, cell.use]
meta <- all.meta[cell.use, ]
 
message("Creating the cellchat object")
cellchat <- createCellChat(object = data.input, meta = meta, group.by = opt$group_by)#"Cohort_Type_Highlevel")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = opt$group_by)#Cohort_Type_Highlevel")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))

#Redorder the colors vector to match the order of the Idents, because the plot functions cannot match using named vectors and rely on position in the vector instead
mycelltypecols <- mycelltypecols[order(factor(names(mycelltypecols), levels=levels(cellchat@idents)))]

message("loading the cellchat database part relevant for the analysis")
CellChatDB <- CellChatDB.human
pdf(outdir%&%"database_categories.pdf")
showDatabaseCategory(CellChatDB)
dev.off()
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
message("Starting main analysis. This may take some time..")
future::plan("multisession", workers = 5)
options(future.globals.maxSize= 891289600)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

##To check if the default method works we can test genes of interest:
#computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.25)
#and if it doesn't work, use trim-0.1 or trim 0.05
# if that works, use type = "truncatedMean" and trim = 0.1 in computeCommunProb

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)

#check cell-cell comm of interest:
message("saving the table of significant CCC events")
df.net <- subsetCommunication(cellchat) #returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
write.table(df.net, paste0(outdir,"/ligand-receptor_interactions.tsv"),sep="\t")
# can subset for spefic cell groups with sources.use and targets.use. Can subset by genes using signaling
auto_pathways <- unique(df.net$pathway_name)
df.net.types <- unique(df.net[,which(colnames(df.net)%in%c("pathway_name", "annotation"))])
write.table(df.net.types, paste0(outdir,"/pathway_interactions_types.tsv"), sep="\t", row.names=F)
write.table(table(df.net.types$annotation)/sum(table(df.net.types$annotation))*100, paste0(outdir,"/pathway_interactions_types_proportions.tsv"))
df.net.types <- data.frame("type"=names(table(df.net.types$annotation)),"values"=table(df.net.types$annotation))
df.net.types$type <- factor(df.net.types$type, levels=c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"))
df.net.types <- df.net.types %>% 
  arrange(desc(type)) %>%
  mutate(prop = values.Freq / sum(df.net.types$values.Freq) *100) 
df.net.types$prop=as.character(round(df.net.types$prop,digits=1)) %&% "%"

message("plotting pathways categories")
pdf(outdir%&%"pathway_interactions_categories.pdf")
ggplot(df.net.types, aes(x = "", y = values.Freq, fill = type)) +
  geom_col(color = "black") +
  geom_text(aes(label = prop),
            position = position_stack(vjust = 0.5)) +
  coord_polar(theta = "y") + theme_void()
dev.off()

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat) #USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use

dir.create(outdir %&% "/latent_patterns/")
dir.create(outdir %&% "/latent_patterns/qc")
if(opt$run_cluster_infe=="yes"){
  message("Running procedure to infer the number of incoming clusters. This will take a long time and the procedure will stop to allow review of output and a proper definition of n_patterns_out and n_patterns_in")
  pdf(outdir %&% "/latent_patterns/qc/" %&% "signaling_score_outgoing.pdf")
  print(selectK(cellchat, pattern = "outgoing"))
  dev.off()
  pdf(outdir %&% "/latent_patterns/qc/" %&% "signaling_score_incoming.pdf")
  print(selectK(cellchat, pattern = "incoming"))
  dev.off()
}else{

message("plotting the aggregate network")
groupSize <- as.numeric(table(cellchat@idents))
pdf(paste0(outdir,"aggregate_network.pdf"), width=10, height=10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions", color.use=mycelltypecols)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength", color.use=mycelltypecols)
dev.off()

message("plotting the network by group")
mat <- cellchat@net$weight
pdf(paste0(outdir, "/by_group_network.pdf"), width=20, height=20)
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i], color.use=mycelltypecols)
}
dev.off()

message("writing to disk the list of all unique significant pathways")
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# Add the custom pathways, if any
if(opt$pathways!="none"){
    pathways.show.all <- c(pathways.show.all, pathways)
}
write(pathways.show.all, outdir%&%"all_significant_pathways.txt")

# check the order of cell identity to set suitable vertex.receiver
myct=levels(cellchat@idents)
vertex.receiver = seq(1,length(myct))

message("plotting the cell-level information")
cellchat <- rankNetPairwise(cellchat)
dir.create(paste0(outdir,"/chord_sent"))
dir.create(paste0(outdir,"/chord_received"))
dir.create(paste0(outdir,"/bubble"))
for(i in 1:length(myct)){
  #all sending from X
  if(length(which(df.net$source==myct[i]))!=0){
    pdf(paste0(outdir,"chord_sent/chord_diagram_sentfrom_",myct[i],".pdf"), width=10, height=10)
    netVisual_chord_gene(cellchat, sources.use = i, targets.use = 1:length(myct), lab.cex = 0.5,legend.pos.y = 0,legend.pos.x = 0, color.use=mycelltypecols)
    dev.off()
  }
  if(length(which(df.net$target==myct[i]))!=0){
    #all received by X
    pdf(paste0(outdir,"/chord_received/chord_diagram_receivedby_",myct[i],".pdf"), width=10, height=10)
    netVisual_chord_gene(cellchat, sources.use = 1:length(myct), targets.use = i, lab.cex = 0.5,legend.pos.y = 0, legend.pos.x = 0,color.use=mycelltypecols)
    dev.off()
  }
  if(length(which(df.net$source==myct[i]))!=0){
    # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
    pdf(paste0(outdir,"/bubble/bubble_plot_",myct[i],".pdf"), width=20, height=20)
    print(netVisual_bubble(cellchat, sources.use = i, targets.use = 1:length(myct), remove.isolate = FALSE, font.size=15 ))
    dev.off()
    # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
  }
}

message("plotting centrality and signaling role")
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#visualize the dominant senders (sources) and receivers (targets) in a 2D space using scatter plot.
pdf(outdir%&%"/signaling_role.pdf")
print(netAnalysis_signalingRole_scatter(cellchat, color.use=mycelltypecols))
dev.off()
#which signals contributing most to outgoing or incoming signaling of certain cell groups.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf(outdir%&%"/signaling_role_hm_outgoing.pdf", width=10, height=10)
print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", color.use=mycelltypecols, width=20, height=20))
dev.off()
pdf(outdir%&%"/signaling_role_hm_incoming.pdf", width=10, height=10)
print(netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", color.use=mycelltypecols, width=20, height=20))
dev.off()

# CUSTOMIZATION
message("Pathway deep-dive")
dir.create(paste0(outdir,"/pathways"))
for(i in 1:length(pathways.show.all)){
  message(pathways.show.all[i] %&% " " %&% i %&% " of " %&% length(pathways.show.all))
  dir.create(outdir%&%"/pathways/"%&%pathways.show.all[i])
  pathdir <- outdir%&%"/pathways/"%&%pathways.show.all[i]

  cwd <- getwd()
  setwd(pathdir)
  #  Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "circle", out.format="pdf", height=15, color.use=mycelltypecols)
  setwd(cwd)

  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathdir, "/", pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 7, height = 2, units = 'in', dpi = 300)
  # show all the significant interactions (L-R pairs) associated with certain signaling pathways
  pdf(pathdir%&%"/chord_"%&%pathways.show.all[i]%&%".pdf", width=15, height=15)
  netVisual_chord_gene(cellchat, sources.use = 1:length(myct), targets.use = 1:length(myct), signaling = pathways.show.all[i], color.use=mycelltypecols,
    legend.pos.x=0, legend.pos.y=0)
  dev.off()
  #Plot gene expression of all genes in a specific pathway
  pdf(pathdir%&%"/expression_"%&%pathways.show.all[i]%&%".pdf", width=15, height=15)
  print(plotGeneExpression(cellchat, signaling = pathways.show.all[i], color.use=mycelltypecols))
  dev.off()
  # Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
  pdf(pathdir%&%"/centrality_hm_"%&%pathways.show.all[i]%&%".pdf")
  netAnalysis_signalingRole_network(cellchat, signaling = pathways.show.all[i], width = 8, height = 2.5, font.size = 10, color.use=mycelltypecols)
  dev.off()
  # Signaling role analysis on the cell-cell communication networks of interest
  #visualize the dominant senders (sources) and receivers (targets) in a 2D space using scatter plot.
  pdf(pathdir%&%"/signaling_role_"%&%pathways.show.all[i]%&%".pdf")
  print(netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show.all[i], color.use=mycelltypecols))
  dev.off()
  
}

#Global communication patterns, i.e. how multiple cells and pathways coordinate together
#OUTGOING
message("Plotting latent patterns")
newplotdir <- outdir %&% "/latent_patterns/"
dir.create(newplotdir)
#Look at when the score starts to drop, this is the best number of patterns
nPatterns = opt$n_patterns_out
#The commands are sensitive to colors in mycelltypecols relative to filtered-out or absent celltypes. We need to filter pre-emptively
mycelltypecols_bck  <- mycelltypecols
tmp <- table(cellchat@meta$celltype)
tmp <- names(tmp)[which(tmp>opt$min_cells)]
mycelltypecols <- as.vector(mycelltypecols[which(names(mycelltypecols)%in%tmp)])
names(mycelltypecols) <- tmp
#plot
pdf(newplotdir%&%"/communication_pattern_outgoing_hm.pdf", width=15, height=15)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns, color.use=mycelltypecols, height=20)
dev.off()
pdf(newplotdir%&%"/communication_patterns_outgoing_alluvial.pdf", width=15, height=15)
print(netAnalysis_river(cellchat, pattern = "outgoing", color.use=mycelltypecols))
dev.off()
pdf(newplotdir%&%"/communication_patterns_outgoing_dotplot.pdf", width=15, height=15)
print(netAnalysis_dot(cellchat, pattern = "outgoing", color.use=mycelltypecols, font.size=15))
dev.off()
#INCOMING
#Look at when the score starts to drop, this is the best number of patterns
nPatterns = opt$n_patterns_out
pdf(newplotdir%&%"/communication_pattern_incoming_hm.pdf", width=15, height=15)
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns, color.use=mycelltypecols, height=20)
dev.off()
pdf(newplotdir%&%"/communication_patterns_incoming_alluvial.pdf", width=15, height=15)
print(netAnalysis_river(cellchat, pattern = "incoming", color.use=mycelltypecols))
dev.off()
pdf(newplotdir%&%"/communication_patterns_incoming_dotplot.pdf", width=15, height=15)
print(netAnalysis_dot(cellchat, pattern = "incoming", color.use=mycelltypecols, font.size=15))
dev.off()

#Manifold classification
#pathway similarity grouping
message("Plotting pathway grouping")
newplotdir <- outdir %&% "pathway_grouping"
dir.create(newplotdir)
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf(newplotdir%&%"/functional_grouping.pdf", width=15, height=15)
print(netVisual_embedding(cellchat, type = "functional", label.size = 3.5))
dev.off()
pdf(newplotdir%&%"/functional_grouping_zoom.pdf", width=15, height=15)
print(netVisual_embeddingZoomIn(cellchat, type = "functional", nCol = 2))
dev.off()
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
pdf(newplotdir%&%"/structural_grouping.pdf", width=15, height=15)
print(netVisual_embedding(cellchat, type = "structural", label.size = 3.5))
dev.off()
pdf(newplotdir%&%"/structural_grouping_zoom.pdf", width=15, height=15)
print(netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2))
dev.off()

message("writing the cellchat object to disk")
saveRDS(cellchat, outdir%&%"/cellchat.RDS")
}#end of else, allowing to skip the entire procedure if run_cluster_infer is set

message("End of analysis")
