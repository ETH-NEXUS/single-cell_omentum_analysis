#####################
## Plot GSVA results
## Use R 4
#####################

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(optparse)

# give out session Info
cat("\n\n\nPrint sessionInfo:\n\n")
print(sessionInfo())
cat("\n\n\n\n")

# command line arguments are parsed
option_list = list(
  make_option("--rdsFile", type = "character", help = "Path to sce object file with input data (sce_celltypes_noatypical.RDS)."),
  make_option("--outdir", type = "character", help = "Path to the directory where output files will be written."),
  make_option("--gsva_results", type = "character", help = ""),
  make_option("--sample_name", type = "character", help = "Sample identifier. Attached to each output name.")
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


###   General parameters   ###
UMAP_aspect_ratio <- "1"


###   Read in UMAP results   ###
#umap_table <- read.csv(opt$umap_results, sep = "\t")

combined = readRDS(opt$rdsFile)
colData_df <- as.data.frame(colData(combined))
umap_gex_df <- as.data.frame(reducedDim(combined, "umap_gex"))
# add umap coordinates to the meta data table of cells
colData_df <- plyr::join(colData_df, umap_gex_df, by = "barcodes")

print(head(colData_df))


###   Read in GSVA results   ###
gsva_table <- read.csv(opt$gsva_results, sep = "\t")

#define column and rows
number_genesets = length(unique(sort(gsva_table$gene.set)))
print(number_genesets)
print("Before:")
print(unique(sort(gsva_table$gene.set)))

gsva_table$gene.set <- gsub("_", " ", gsva_table$gene.set)
print("After:")
print(unique(sort(gsva_table$gene.set)))

###   Join Results
md <- dplyr::full_join(colData_df, gsva_table, by = "barcodes")
stopifnot(unique(sort(md$barcodes)) == unique(sort(gsva_table$barcodes)))
stopifnot(unique(sort(md$barcodes)) == unique(sort(colData_df$barcodes)))

### get min and max for the legend
break_low <- min(md$value_limited)
print(break_low)
break_high <- max(md$value_limited)
print(break_high)
label_low <- paste("<=", break_low)
print(label_low)
label_high <- paste(">=", break_high)
print(label_high)

cat("\n\nsummary(md)$value_limited:\n")
print(summary(md$value_limited))

###   Plot GSVA   ###


#nr.cols <- min(6, number_genesets)
nr.cols = 6
nr.rows <- floor((number_genesets + 6) / nr.cols)
print(nr.rows)

pp <- ggplot(md, aes(x = GEX_umap1, y = GEX_umap2)) +
  geom_point(aes(color = value_limited), size = 0.1) +
  xlab("umap-1") + ylab("umap-2") +
  theme(panel.background = element_rect(fill = "grey80")) +
  scale_color_distiller(name = "", palette = "RdBu",
                        breaks = c(break_low, 0, break_high),
                        labels = c(label_low, "0.00", label_high)) +
  facet_wrap(~gene.set, ncol = 6, labeller = label_wrap_gen(width = 25, multi_line = TRUE)) +
  theme(aspect.ratio = UMAP_aspect_ratio)
# save plot
ggsave(paste0(opt$outdir, opt$sample_name, "_gsva_umap.png"), pp,
       dpi = 300, width = 50, height = (7.5 * nr.rows), units = "cm")
