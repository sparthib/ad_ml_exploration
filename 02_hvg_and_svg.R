library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(BiocParallel)
library(doParallel)
library(here)
devtools::install_github("cellgeni/sceasy")
library(sceasy)
library(reticulate)

sceasy::convertFormat(spe_postqc, from="sce", to="anndata",
                      outFile='filename.h5ad')
