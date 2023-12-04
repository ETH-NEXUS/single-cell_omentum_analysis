
################################################################################
## cell type classification - GSVA
## NOTE: this is an independent piece of analysis. 
##       all output files, including plots, will already be handled in this script
################################################################################

lby = c("optparse", "scran", "reshape2", "uwot", "GSVA", "aroma.light",
        "ggplot2", "pheatmap","RColorBrewer", "limma", "viridis")
resp = lapply(lby, require, character.only=T, warn.conflicts=F, quietly=T)
if(!all(unlist(resp))) stop("Could not load one or more packages")
rm(resp, lby)

# options(error=function()traceback(2))
options(stringsAsFactors = FALSE)
fontsize = theme(axis.text=element_text(size=9), axis.title=element_text(size=11))
theme_set(theme_bw(12) + fontsize)
col.pal = colorRampPalette( rev(brewer.pal(11, "RdBu")) )(255)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# convenience function for string concatenation
'%&%' = function(a,b) paste(a,b,sep="")

# command line arguments are parsed
option_list = list(
  make_option("--SCE", type = "character", help = "Path to sce object file with input data (sce_celltypes_noatypical.RDS)."),
  make_option("--outputDirec", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--geneset", type = "character", help = "Hallmark geneset library gmt file."),
  make_option("--sampleName", type = "character", help = "Sample identifier. Attached to each output name.")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

path = opt$outputDirec %&% opt$sampleName

print(opt$SCE)

path = opt$outputDirec %&% opt$sampleName


################################################################################
## main code starts here
################################################################################
## load input data
sce_data = readRDS(opt$SCE)

################################################################################
#
# signature expression painting on tSNE plots (phenograph)
#
# load ref.gene.list
# make sure all columns are captured when reading in the gene set lists in gmt format
tmp = readLines(opt$geneset)
tmp = lapply(tmp, function(x) strsplit(x, "\\\t")[[1]])
names(tmp) = sapply(tmp, function(x) x[1])
# remove gene set names from gene list
gset = sapply(tmp, function(x) x[-1])
# remove gene set description from gene list
gset = sapply(gset, function(x) x[-1])
names(gset) = gsub("HALLMARK_", "", names(gset))

# select subset gene expression matrix
idxs = ids2indices(gset, rowData(sce_data)$SYMBOL)
all.genes = unique(as.numeric(unlist(idxs)))
             
t.m = assay(sce_data, "pearson_resid")[all.genes,]
             
# estimate geneset-sample matrix from gene-sample matrix
rgsa = gsva(t.m, gset, method="gsva")

# save sce object for gsva-scores instead of gene counts for later use 
sce_gsva = SingleCellExperiment(assays=list(gsva = rgsa), 
                                colData = colData(sce_data),
                                metadata = gset)
saveRDS(sce_gsva, path %&% ".sce_gsva.RDS")
# optional: save geneset-cell score matrix.
# write.table(rgsa, path %&% ".geneset_scores.tsv", sep="\t", col.names = NA)

