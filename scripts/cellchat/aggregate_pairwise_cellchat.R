###################################################
## File name: aggregate_pairwise_cellchat.R 
## Author: Matteo Carrara
###################################################

## GENERAL:
## This script gathers a set of pairwise cellchat cell-cell-communication comparisons in RDS files
## and generates aggregated results

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

#multiple comparison correction
datalist$name <- as.character(datalist$name)
datalist$group <- as.character(datalist$group)
datalist$padj <- p.adjust(datalist$pvalues, method="bonferroni")
#datalist_nonzero <- datalist[which(datalist$contribution!=0),]
#datalist_sig <- datalist_nonzero[which(datalist_nonzero$padj<opt$padj),]
tmpsig <- datalist[which(datalist$padj<opt$padj),]

# take the top N for each pariwise. Define top as those with lower padj. If there are too many at 0, prioritize those with larger contribution
mytop <- NULL
for(i in 1:length(cellchatnames)){
  tmp <- datalist[which(datalist$comparison==cellchatnames[i]),]
  tmp_sig <- tmp[order(tmp$padj),]
  tmp_sig <- tmp_sig[which(tmp_sig$padj<opt$padj),]
  tmp_sig <- tmp_sig[which(tmp_sig$contribution!=0),]

  if(length(which(tmp_sig$padj==0))>opt$top){
    tmp_zero <- tmp_sig[which(tmp_sig$padj==0),]
    tmp_zero <- tmp_zero[order(tmp_zero$contribution, decreasing=TRUE),]
  }else{
    tmp_zero <- tmp
  }

  if(opt$top>=1000000){
    topnum=nrow(tmp_zero)
  }else{
    topnum=opt$top
  }
  mytop_part <- tmp_zero[1:topnum,]
  if(nrow(mytop_part)!=topnum){stop("wrong number of top pw")}
  #now find the other part if not already present
  for(j in 1:nrow(mytop_part)){
    if(length(which(rownames(mytop_part)%in%rownames(mytop_part)[j]))==1){
      myline <- tmp[which(tmp$name==mytop_part$name[j]),]
      
      if(nrow(myline)>2){stop("found unexpected many matches")}
      mytop <- rbind(mytop, myline)
      message("i=",i, "j=",j," dim=",nrow(mytop))
    }
  }
  
}

mycols <- c("benign"="#1E88E5","tumor_distant"="#004D40","tumor_peritumoral"="#FFC107","tumor_metastasis"="#D81B60")
mycomp <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00","#CC79A7")
names(mycomp) <- unique(datalist$comparison)

myf <- datalist[which(datalist$name%in%tmpsig$name),]
myf <- myf[order(myf$name),]


mycollist<-list()
mysigname <- sort(unique(tmpsig$name))
for(i in 1:length(mysigname)){
  tmp <- tmpsig[which(tmpsig$name==mysigname[i]),]
  tmp <- sort(unique(tmp$comparison))
  if(length(tmp)>6){stop()}
  for(j in 1:length(tmp)){
    if(i==1){mycollist[[j]]=mycomp[tmp[j]]}else{
    mycollist[[j]]=c(mycollist[[j]],mycomp[tmp[j]])}
  }
  if(j<6){
    for(k in (j+1):6){
      if(i==1){mycollist[[k]]="#FFFFFF"}else{
      mycollist[[k]]=c(mycollist[[k]],"#FFFFFF")
    }
    }

  }
}

ggplot(data=myf, aes(x=name, y=contribution.scaled, fill=group)) +
  geom_bar(stat="identity", position=position_dodge(0.7), width=0.65)+
  scale_fill_manual(values=mycols) + #facet_grid(~comparison, scales="free_x") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text = element_text(size = 18), plot.margin=margin(1,1,6,1, "cm"), axis.title.x = element_blank(), legend.position=c(0.5,-0.8), legend.direction="horizontal") +
  annotate("point", x=c(1:length(mysigname)), y=-0.6, col=mycollist[[1]], size=4) +
  annotate("point", x=c(1:length(mysigname)), y=-0.66, col=mycollist[[2]], size=4) +
  annotate("point", x=c(1:length(mysigname)), y=-0.72, col=mycollist[[3]], size=4) +
  annotate("point", x=c(1:length(mysigname)), y=-0.78, col=mycollist[[4]], size=4) +
    annotate("point", x=c(1:length(mysigname)), y=-0.84, col=mycollist[[5]], size=4) +
  annotate("point", x=c(1:length(mysigname)), y=-0.9, col=mycollist[[6]], size=4) +

  scale_y_continuous(expand=c(0,0))+
    coord_cartesian(ylim = c(0, 1.55), clip="off") +
    geom_point(data = . %>% 
               slice(1:6) %>% 
               # optional: list colours in the desired label order
               mutate(col = forcats::fct_inorder(names(mycomp))),
             aes(colour = col), 
             alpha = 0) +
             scale_color_manual(name = "Comparison",
                     values = mycomp,
                     guide = guide_legend(override.aes = list(alpha = 1, shape = 15, size = 5)))
    

  ggsave("pariwise_ccc_aggregation.pdf", width=25, height=8)

write.table(myf, "ccc_aggregation_table.tsv", sep="\t")
