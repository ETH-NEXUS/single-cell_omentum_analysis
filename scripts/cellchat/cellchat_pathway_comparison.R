###################################################
## File name: cellchat_pathway_comparison.R
## Author: Matteo Carrara
###################################################

## GENERAL:
## This script calculates pairwise comparisons of pathways obtained from a cellchat analysis

library(SingleCellExperiment)
library(optparse)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(ggplot2)
library(yaml)

cat("\n\n\nPrint sessionInfo():\n\n")
print(sessionInfo())

option_list <- list(
  make_option("--rdsFiles", type = "character", help = "The two RDS files with the cellchat data, separated by comma. Only two are accepted because this is a pairwise comparison"),
  make_option("--color_yaml", type="character", help="the yaml file with the necessary color configuration"),
  make_option("--group_by", type="character", help="the variable the cellchat object was grouped by. Usually the celltype"),
  make_option("--comparison_names", type="character", help="Names to assign to each rdsFile/element of comparison. Needs to be in the same order as the RDS files"),
  make_option("--subset_var", type="character", help="Name of the metadata variable used to subset. It must correspond to the column name containing the variable on which to compare"),
  make_option("--pos_dataset", type="character", help="Name of one of the comparisons in comparison_names to use as positive. Positive fold changes will refer to positive regulation in the selected pathway")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""

message("Reading the color config")
myallcolors <- read_yaml(opt$color_yaml)
mycelltypecols <- unlist(myallcolors[[which(names(myallcolors)==opt$group_by)]])

message("creating output dir")
dir.create("./output_comparison")

cellchatfiles <- unlist((strsplit(opt$rdsFiles,",")))
if(length(cellchatfiles)>2){stop("Please provide exactly 2 cellchat files")}
cellchatnames <- unlist((strsplit(opt$comparison_names,",")))
if(length(cellchatnames)>2){stop("Please provide exactly 2 names for the cellchat files")}

cellchatlist <- list()
for(i in 1:length(cellchatfiles)){
    cellchatlist[[i]]<-readRDS(cellchatfiles[i])
    names(cellchatlist)[i]<-cellchatnames[i]
}

# Define the cell labels to lift up
group.new = levels(cellchatlist[[1]]@idents)
cellchatlist[[2]] <- liftCellChat(cellchatlist[[2]], group.new)
#>The CellChat object will be lifted up using the cell labels B.cells, Endothelial, HGSOC, Lymphatic.endothelial, Mast.cells, Mesenchymal.progenitor, Mesothelial, Monocytes.macrophages, Neutrophils, Plasma, Plasmacytoid.dendritic, T.NK, Transdifferentiating
#> Update slots object@net, object@netP, object@idents in a single dataset...
object.list <- list(a = cellchatlist[[1]], b = cellchatlist[[2]])
names(object.list)<-cellchatnames
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
#> Warning in mergeCellChat(object.list, add.names = names(object.list),
#> cell.prefix = TRUE): Prefix cell names!
#> The cell barcodes in merged 'meta' is  rep1_AAACCTGCACCAACCG rep1_AAACGGGAGCCGATTT rep1_AAACGGGAGTATCGAA rep1_AAACGGGCATCTCCCA rep1_AAAGATGCACTTGGAT rep1_AAAGATGCAGTTCATG
#> Warning in mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE): The cell barcodes in merged 'meta' is different from those in the used data matrix.
#>               We now simply assign the colnames in the data matrix to the rownames of merged 'mata'!
#> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 'idents', 'var.features' , 'DB', and 'LR'.

myct=levels(cellchat@idents$joint)

outdir <- "./output_comparison/" %&% cellchatnames[1] %&% "-" %&% cellchatnames[2] %&% "/"
message("saving in ", outdir)
dir.create(outdir)

message("comparing probabilities...")
dir.create(outdir %&% "prob_compare")
newout <- outdir %&% "prob_compare/"

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg<-gg1 + gg2
pdf(newout %&% "cell-level_heatmap.pdf", width=15, height=15)
print(gg)
dev.off()

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg<-gg1 + gg2
pdf(newout %&% "pathway_information_flow.pdf", width=15, height=15)
print(gg)
dev.off()

dir.create(newout %&% "bubble")
for(i in 1:length(myct)){
  pdf(newout %&% "bubble/" %&% myct[i] %&% "_LR_comm_prob_comparison_bubble_out.pdf",  width=15, height=15)
  gg1 <- try({netVisual_bubble(cellchat, sources.use = c(myct[i]), targets.use = c(1:length(myct)),  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", cellchatnames[2]), angle.x = 45, remove.isolate = T)},silent=TRUE)
  if(is(gg1)=="try-error"){gg1 <- NULL}
  gg2 <- try({netVisual_bubble(cellchat, sources.use = c(myct[i]), targets.use = c(1:length(myct)),  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in ", cellchatnames[2]), angle.x = 45, remove.isolate = T)},silent=TRUE)
  if(is(gg2)=="try-error"){gg2 <- NULL}
  gg<-gg1 + gg2
  print(gg)
  dev.off()
  pdf(newout %&% "bubble/" %&% myct[i] %&% "_LR_comm_prob_comparison_bubble_in.pdf",  width=15, height=15)
  try({gg1 <- netVisual_bubble(cellchat, sources.use = c(1:length(myct)), targets.use = myct[i],  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", cellchatnames[2]), angle.x = 45, remove.isolate = T)},silent=TRUE)
  if(is(gg1)=="try-error"){gg1 <- NULL}
  try({gg2 <- netVisual_bubble(cellchat, sources.use = c(1:length(myct)), targets.use = myct[i],  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in ", cellchatnames[2]), angle.x = 45, remove.isolate = T)},silent=TRUE)
  if(is(gg2)=="try-error"){gg2 <- NULL}
  gg<-gg1 + gg2
  print(gg)
  dev.off()
}

message("comparing differential gene expression")
dir.create(outdir %&% "de_compare")
newout<-outdir %&% "de_compare/"

message("setting the reference dataset for positive fold changes as: ", cellchatnames[2], ". A positive fold change will mean upregulation in this dataset")
pos.dataset<-opt$pos_dataset
features.name<-paste(cellchatnames[2], cellchatnames[1], sep="_")
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = opt$subset_var, pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in the condition
net.up <- subsetCommunication(cellchat, net = net, datasets = pos.dataset,ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in the other condition, i.e.,downregulated in the main condition
net.down <- subsetCommunication(cellchat, net = net, datasets = cellchatnames[1],ligand.logFC = -0.1, receptor.logFC = -0.1)

#Since the signaling genes in the net.up and net.down might be complex with multi-subunits, we can do further deconvolution to obtain the individual signaling genes.
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

myct=levels(cellchat@idents$joint)

dir.create(newout %&% "bubble")
pairLR.use.up = net.up[, "interaction_name", drop = F]
pairLR.use.down = net.down[, "interaction_name", drop = F]

for(i in 1:length(myct)){
  try({gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = myct[i], targets.use = c(1:length(myct)), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", pos.dataset))},silent=TRUE)
  if(is(gg1)=="try-error"){gg1 <- NULL}
  try({gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = myct[i], targets.use = c(1:length(myct)), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", pos.dataset))},silent=TRUE)
  if(is(gg2)=="try-error"){gg2 <- NULL}
  gg<-gg1 + gg2
  pdf(newout %&% "bubble/bubble_differential_" %&% myct[i] %&% "_outgoing.pdf",  width=15, height=15)
  print(gg)
  dev.off()

  try({gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:length(myct)), targets.use = myct[1], comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", pos.dataset))},silent=TRUE)
  if(is(gg1)=="try-error"){gg1 <- NULL}
  try({gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:length(myct)), targets.use = myct[1], comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", pos.dataset))},silent=TRUE)
  if(is(gg2)=="try-error"){gg2 <- NULL}
  gg<-gg1 + gg2
  pdf(newout %&% "bubble/bubble_differential_" %&% myct[i] %&% "_incoming.pdf",  width=15, height=15)
  print(gg)
  dev.off()
}
dir.create(newout %&% "chords")
# Chord diagram
for(i in 1:length(myct)){
  message("chords for ", myct[i])
  pdf(newout %&% "chords/chords_differential_" %&% myct[i] %&% "_outgoing.pdf",  width=15, height=15)
  par(mfrow = c(1,2), xpd=TRUE)
  try({netVisual_chord_gene(object.list[[2]], sources.use = myct[i], targets.use = c(1:length(myct)), slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", pos.dataset),color.use=mycelltypecols)},silent=TRUE)
  try({netVisual_chord_gene(object.list[[1]], sources.use = myct[i], targets.use = c(1:length(myct)), slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", pos.dataset),color.use=mycelltypecols)},silent=TRUE)
  dev.off()
  pdf(newout %&% "chords/chords_differential_" %&% myct[i] %&% "_incoming.pdf",  width=15, height=15)
  par(mfrow = c(1,2), xpd=TRUE)
  try({netVisual_chord_gene(object.list[[2]], sources.use = c(1:length(myct)), targets.use = myct[i], slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", pos.dataset),color.use=mycelltypecols)},silent=TRUE)
  try({netVisual_chord_gene(object.list[[1]], sources.use = c(1:length(myct)), targets.use = myct[i], slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", pos.dataset),color.use=mycelltypecols)},silent=TRUE)
  dev.off()
}

#remove NAs for correct visualisation
net.up.bck<-net.up
net.down.bck<-net.down
net.down$receptor.logFC[which(is.na(net.down$receptor.logFC))]<-0
net.down$receptor.pvalues[which(is.na(net.down$receptor.pvalues))]<-1
net.down$receptor.pct.1[which(is.na(net.down$receptor.pct.1))]<-0
net.down$receptor.pct.2[which(is.na(net.down$receptor.pct.2))]<-0
net.up$receptor.logFC[which(is.na(net.up$receptor.logFC))]<-0
net.up$receptor.pvalues[which(is.na(net.up$receptor.pvalues))]<-1
net.up$receptor.pct.1[which(is.na(net.up$receptor.pct.1))]<-0
net.up$receptor.pct.2[which(is.na(net.up$receptor.pct.2))]<-0
# visualize the enriched ligands in the first condition
pdf(newout %&% "/wordcloud_enriched_down_" %&% pos.dataset %&% ".pdf")
print(computeEnrichmentScore(net.down, species = 'human',color.use=mycelltypecols))
dev.off()
# visualize the enriched ligands in the second condition
pdf(newout %&% "/wordcloud_enriched_up_" %&% pos.dataset %&% ".pdf")
print(computeEnrichmentScore(net.up, species = 'human',color.use=mycelltypecols))
dev.off()

#We can plot the gene expression distribution of signaling genes related to L-R pairs or signaling pathway using a Seurat wrapper function plotGeneExpression.
#cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("NL", "LS")) # set factor level
#plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)

message("done")

