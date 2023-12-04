###################################################
## File name: cell_level_pw_comparison.R
## Author: Matteo Carrara
###################################################

## GENERAL:
## This script generate pathway comparison results at cell level using cellchat

library(SingleCellExperiment)
library(optparse)
library(CellChat)
library(yaml)
library(ggbreak)
options(stringsAsFactors = FALSE)


option_list <- list(
make_option("--rdsFiles", type = "character", help = "All RDS files with the cellchat object saved during the comaprison analysis, separated by comma."),
  make_option("--color_yaml", type="character", help="the yaml file with the necessary color configuration"),
  make_option("--comparison_names", type="character", help="Names to assign to each rdsFile/element of comparison. Needs to be in the same order as the RDS files"),
  make_option("--top", type="numeric", help="Number of top genes to select for each comparison")
  make_option("--padj", type="numeric", help="p_adjusted threshold for significance")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

#opt=list()
#opt$group_by="Tissue_Region"
#opt$rdsFiles=c("../benign-tumor_distant/cellchat_obj.RDS","../benign-tumor_metastasis/cellchat_obj.RDS","../benign-tumor_peritumoral/cellchat_obj.RDS","../tumor_distant-tumor_metastasis/cellchat_obj.RDS","../tumor_distant-tumor_peritumoral/cellchat_obj.RDS","../tumor_peritumoral-tumor_metastasis/cellchat_obj.RDS")
#opt$comparison_names=c("benign-distant", "benign-metastasis", "benign-peritumoral", "distant-metastasis", "distant-peritumoral", "peritumoral-metastasis")
#opt$color_yaml="../../../color_config_full.yaml"
#opt$top=1000000

"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""

message("Reading the color config")
myallcolors <- read_yaml(opt$color_yaml)
mycelltypecols <- unlist(myallcolors[[which(names(myallcolors)==opt$group_by)]])

cellchatfiles <- unlist((strsplit(opt$rdsFiles,",")))
cellchatnames <- unlist((strsplit(opt$comparison_names,",")))
if(length(cellchatnames)!=length(cellchatfiles)){stop("Please provide the same number of cellchat files and names")}

cellchatlist <- list()
gglist <- list()
datalist <- NULL
for(i in 1:length(cellchatfiles)){
    cellchatlist[[i]]<-readRDS(cellchatfiles[i])
    names(cellchatlist)[i]<-cellchatnames[i]
    gglist[[i]]<-rankNet(cellchatlist[[i]], mode = "comparison", stacked = F, do.stat = TRUE)
    names(gglist)[i]<-cellchatnames[i]

    tmpdata <- gglist[[i]]$data
    tmpdata$comparison <- cellchatnames[i]
    datalist <- rbind(datalist, tmpdata)
}

final <- NULL
for(k in 1:length(cellchatlist)){
    myct=levels(cellchatlist[[k]]@idents$joint)
    for(i in 1:length(myct)){

        incr_out <- try({netVisual_bubble(cellchatlist[[k]], sources.use = c(myct[i]), targets.use = c(1:length(myct)),  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", cellchatnames[2]), angle.x = 45, remove.isolate = T)},silent=TRUE)
        if(is(incr_out)=="try-error"){incr_out <- NULL}else{incr_out <- incr_out$data}
        decr_out <- try({netVisual_bubble(cellchatlist[[k]], sources.use = c(myct[i]), targets.use = c(1:length(myct)),  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in ", cellchatnames[2]), angle.x = 45, remove.isolate = T)},silent=TRUE)
        if(is(decr_out)=="try-error"){decr_out <- NULL}else{decr_out <- decr_out$data}
        incr_in <- try({netVisual_bubble(cellchatlist[[k]], sources.use = c(1:length(myct)), targets.use = myct[i],  comparison = c(1, 2), max.dataset = 2, title.name = paste0("Increased signaling in ", cellchatnames[2]), angle.x = 45, remove.isolate = T)},silent=TRUE)
        if(is(incr_in)=="try-error"){incr_in <- NULL}else{incr_in <- incr_in$data}
        decr_in <- try({netVisual_bubble(cellchatlist[[k]], sources.use = c(1:length(myct)), targets.use = myct[i],  comparison = c(1, 2), max.dataset = 1, title.name = paste0("Decreased signaling in ", cellchatnames[2]), angle.x = 45, remove.isolate = T)},silent=TRUE)
        if(is(decr_in)=="try-error"){decr_in <- NULL}else{decr_in <- decr_in$data}
        tmp <- rbind(incr_out, incr_in, decr_out, decr_in)
        tmp <- tmp[which(tmp$source==myct[i]),]
        tmp <- tmp[which(tmp$target==myct[i]),]
        if(!is.null(tmp)){
            tmp$comparison <- cellchatnames[k]
        }
        final <- rbind(final, tmp)
    }

}
final <- final[which(!duplicated(final)),]
write.table(final, "pw_comparison_from_bubbleplot.tsv",sep="\t")
