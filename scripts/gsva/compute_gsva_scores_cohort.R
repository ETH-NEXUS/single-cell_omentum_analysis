################################################################################
## Run GSVA
#################################################################################
# Run GSVA on a merged SCE object
# Script does not generate any plots.
# Project: scOmentum 
# May 2021

lby = c("optparse", "scran", "reshape2", "uwot", "GSVA", "aroma.light",
        "ggplot2", "pheatmap","RColorBrewer", "limma", "viridis")
resp = lapply(lby, require, character.only=T, warn.conflicts=F, quietly=T)
if(!all(unlist(resp))) stop("Could not load one or more packages")
rm(resp, lby)

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
  make_option("--outDir", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--geneset", type = "character", help = "Geneset library gmt file."),
  make_option("--outName", type = "character", help = "Sample identifier. Attached to each output name."),
  make_option("--genesetName", type = "character", help = "Gene set identifier. Attached to each output name.")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

print(opt$SCE)

#path = opt$outputDirec %&% opt$sampleName


################################################################################
## main code starts here
################################################################################

## load input data
sce_data = readRDS(opt$SCE)
print(sce_data)
cat("\n\n\nGenes:\n")
print(head(rownames(sce_data)))

################################################################################

# make sure all columns are captured when reading in the gene set lists in gmt format
tmp = readLines(opt$geneset)
tmp = lapply(tmp, function(x) strsplit(x, "\\\t")[[1]])
names(tmp) = sapply(tmp, function(x) x[1])
# remove gene set names from gene list
gset = sapply(tmp, function(x) x[-1])
# remove gene set description from gene list
gset = sapply(gset, function(x) x[-1])
names(gset)
#names(gset) = gsub("HALLMARK_", "", names(gset))

# select subset gene expression matrix
idxs = ids2indices(gset, rownames(sce_data))
all.genes = unique(as.numeric(unlist(idxs)))
t.m = assay(sce_data, "normcounts")[all.genes,]

# estimate geneset-sample matrix from gene-sample matrix
# GSVA
rgsa = gsva(t.m, gset, method="gsva")
# reformat results (for plotting)
t.plot = melt(rgsa)
names(t.plot) = c("gene.set", "barcodes", "value")
print("str(t.plot):")
print(str(t.plot))
print("summary(t.plot):")
print(summary(t.plot))

# Trim outlier values
t.plot$value_limited <- t.plot$value
perc_1 <- quantile(t.plot$value, prob = 0.01)
perc_99 <- quantile(t.plot$value, prob = 0.99)
t.plot$value_limited[t.plot$value_limited <= perc_1] <- perc_1
print("Number of data points set to perc_1:")
print(table(t.plot$value_limited == perc_1))
t.plot$value_limited[t.plot$value_limited >= perc_99] <- perc_99
print("Number of data points set to perc_99:")
print(table(t.plot$value_limited == perc_99))
print("str(t.plot$value_limited):")
print(str(t.plot$value_limited))
print("summary(t.plot$value_limited) after trimmming the outliers:")
print(summary(t.plot$value_limited))
t.plot$value_limited <- round(t.plot$value_limited, digits = 2)
break_low <- min(t.plot$value_limited)
print(break_low)
break_high <- max(t.plot$value_limited)
print(break_high)
label_low <- paste("<=", break_low)
print(label_low)
label_high <- paste(">=", break_high)
print(label_high)
print("summary(t.plot$value_limited) after rounding the values to two digits:")
print(summary(t.plot$value_limited))

# write results into file
filename <- paste0(opt$outDir, opt$outName, "_", opt$genesetName , "_GSVA.tsv")
write.table(t.plot, file = filename, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
