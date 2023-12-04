suppressPackageStartupMessages({
    library(SingleCellExperiment)
  library(dittoSeq)
})
# parse command line arguments
 option_list <- list(
   make_option("--cohort_object", type = "character", help = "Integrated Seurat Object"),
   make_option("--outdir", type = "character", help = "Full path to output directory."),
   make_option("--outName", type = "character", help = "Prefix name of output files."),
   make_option("--genes", type="character", help = "genes to plot in the DotPlot.")
 )
 opt_parser <- OptionParser(option_list = option_list)
 opt <- parse_args(opt_parser)


# convenience function for string concatenation
"%&%" <- function(a, b) paste(a, b, sep = "")
empty_string <- ""


###################################
###   Set plotting Parameters   ###
###################################

# If one of the values plot_aspect_ratio, chunksize, or number of columns is changed in the following section
# the other values should be checked as well.
# Also, check physical size of the generated output file in the respective plotting section of the script.

# IMPORTANT: when generating UMAP plots with the function RunUMAP() from Seurat the parameters "spread", "min.dist" and "n.neighbors"
# have a great impact on the resulting plot. They should be checked and optimised for each project.

# for all UMAP plots the same aspect ratio is used to have a comparable visualisation
plot_aspect_ratio <- 0.7

### The selected genes and gene sets are plotted in chunks so that the file size is not too large:

# chunksize_genes determines the number of genes shown in one expression plot
chunksize_genes <- 9
# columns_genes determines the number of columns the gene expression plots are arranged in.
## E.g. "columns_genes <- 1" arranges all plots in the output file below one another in a single column.
columns_genes <- 3

# chunk size of gene sets in GSVA plots
chunksize_genesets <- 12
# number of columns in GSVA plots
columns_genesets <- 4

# override the theme change of cowplot
theme_set(theme_grey())
print(Sys.time())


########################
###   Read in data   ###
########################


seurat_integrated <- readRDS(opt$cohort_object)
print(seurat_integrated)

set.seed(2)

epitopes <- read.delim(opt$genes, header=FALSE)
epitopes <- epitopes[,1]#sort(epitopes[,1])

epitopes[(!epitopes %in% rownames(seurat_integrated))]
epitopes <- epitopes[(epitopes %in% rownames(seurat_integrated))]

seurat_integrated$umap_subtype_cl=factor(as.character(seurat_integrated$umap_subtype_cl), levels=as.character(c(5,2,6,7,12,1,3,11,13,8,14,9, 4, 10)))

mydotplot <- dittoDotPlot(seurat_integrated, vars = epitopes, group.by = "umap_subtype_cl") +
  ylab("Cluster") + xlab("Gene") + scale_color_gradient2(low="#1f66ac", mid="grey", high="#af050e")
mydotplot$labels$colour="relative expression"
ggsave(filename= paste0(opt$outdir, opt$outName, "_DotPlot_subtype_cl.pdf"),
  width = 40, height = 20, dpi = 600, units = "cm", plot=mydotplot)

mydotplot <- dittoDotPlot(seurat_integrated, vars = epitopes, group.by = "umap_subtype_cl", scale=FALSE) +
  ylab("Cluster") + xlab("Gene") + scale_color_gradient2(low="#1f66ac", mid="grey", high="#af050e")
mydotplot$labels$colour="relative expression"
ggsave(filename= paste0(opt$outdir, opt$outName, "_DotPlot_subtype_cl_unscaled.pdf"),
  width = 40, height = 20, dpi = 600, units = "cm", plot=mydotplot)

