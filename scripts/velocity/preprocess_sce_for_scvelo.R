################################################################################
## preprocess_sce_for_scvelo.R
## Author: Matteo Carrara (carrara@nexus.ethz.ch)
## Creation date: August 2022
## Description: processed, normalized counts and metadata from the pipeline come in an RDS file containing a SingleCellExperiment object. This is not compatible with scVelo and is exceedingly complex to reliably load in python as-is. This script convers the SingleCellExperiment object in a scvelo-friendly anndata object and saves it in an h5ad file, which is easier to load in python
################################################################################

library(zellkonverter)

parser <- argparse::ArgumentParser()
parser$add_argument("--scefile", type = "character", help = "Path to the RDS file containing the SingleCellExperiment object with the single cell analysis results")
parser$add_argument("--outile", type = "character", help = "Path and filename of the output h5ad file"
parser$add_argument("--distances", type = "character", help = "Name of the reducedDims element that store the distances of interest")

print("Reading the RDS file for the conversion of the SCE object into an AnnData object")
sce <- readRDS(args$scefile)
data_class <- class(sce)

if(data_class != "SingleCellExperiment"){
    stop("FATAL: Tried to conver in anndata an object that does not appear to be of class SingleCellExperiment")
}

print("The content matches the SCE specifications. Converting in AnnData and saving to disk")
sce_anndata <- SCE2AnnData(sce, verbose=T)

#Warning: anndata needs to be loaded after convertion to avoid conflicts
library(anndata)
write_h5ad(sce_anndata, args$outfile)

