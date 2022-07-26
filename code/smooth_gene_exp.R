library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(doParallel)
library(here)
library(sessioninfo)
library(dplyr)
#BiocManager::install("BayesSpace")
library(BayesSpace)



spe_harmony <- readRDS(here::here("input_data",
                                  paste0("spe_harmony_wholegenome.rds")))



spe <-readRDS(here::here("input_data",
                       paste0("spe_wholegenome_postqc.rds")))

spe<- scuttle::logNormCounts(spe)


colData(spe)$col <- colData(spe_harmony)$col
colData(spe)$row <- colData(spe_harmony)$row

rm(spe_harmony)

neighbors_list <- BayesSpace:::.find_neighbors(spe, "Visium")


for(i in seq_len(length(neighbors_list))){

    neighbors_list[[i]] <-  neighbors_list[[i]] + 1

}


for (i in seq_len(length(neighbors_list))){


}
