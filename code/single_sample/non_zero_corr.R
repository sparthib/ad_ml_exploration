#library(sgejobs)
# sgejobs::job_single('non_zero_corr',
#                     create_shell = TRUE,
#                     create_logdir = TRUE, task_num = 10,
#                     queue = "bluejay",
#                     cores = 5L,
#                     command = "Rscript non_zero_corr.R",
#                     memory = "20G")
#



##This script computes all combos of correlation values for each of the 10 samples.
library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(doParallel)
library(BiocParallel)
library(here)
library(sessioninfo)
library(dplyr)

#### cores ####


#numCores <- makeCluster(detectCores(), type='PSOCK') # grabs max available
#switch to 4 for now
# numCores <- 4
#
# options('mc.cores' = numCores)
# registerDoParallel(numCores)

#### load data ####


spe_postqc <-
    readRDS(here::here("input_data",
                       paste0("spe_wholegenome_postqc.rds")))

spe_postqc <- scuttle::logNormCounts(spe_postqc)
spe_counts <- logcounts(spe_postqc)

sample_ids <- c(
    "S1_A1_Br3874" ,
    "S1_B1_Br3854",
    "S1_C1_Br3873" ,
    "S1_D1_Br3880" ,
    "S2_A1_Br3874" ,
    "S2_B1_Br3854",
    "S2_C1_Br3873",
    "S2_D1_Br3880",
    "S3_A1_Br3874" ,
    "S3_D1_Br3880"
)
#### take only first sample ####

# for(i in sample_ids){
#     dir.create(here("plots",
#                     "02_sample_subset",
#                     i))
# }

s = as.numeric(Sys.getenv("SGE_TASK_ID"))
print(s)
ix <- colData(spe_postqc)$sample_id_short == sample_ids[s]
spe_sub <- spe_postqc[, ix]


PAbeta <- as.matrix(colData(spe_sub)$PAbeta)
names(PAbeta) <- rownames(colData(spe_sub))

mat2 = PAbeta[PAbeta !=0]
i = 1
res_list <- list()


# for (i in seq_len(nrow(logcounts(spe_sub)))){
#     #mat1 : for a given gene, all spots with non-zero expression
#     mat1 <- as.matrix(logcounts(spe_sub))[i,][as.matrix(logcounts(spe_sub))[i,] != 0]
#     if(length(mat1)!=0){
#         #if there remain some spots with non-zero expression
#         #samecols <- intersect(names(mat1),names(mat2))
#         res_list[[length(res_list)+1]] <- merge(mat1,mat2, by = 0)
#         gene_names <- c(gene_names, rownames(logcounts(spe_sub))[i])
#
#
#     }
# }



fun <- function(i) {
    mat1 <- as.matrix(logcounts(spe_sub))[i,][as.matrix(logcounts(spe_sub))[i,] != 0]
    merge(mat1,mat2, by = 0)


}

t0 = Sys.time()
res_list <- bplapply(1:nrow(logcounts(spe_sub)), fun, BPPARAM = MulticoreParam(5))

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
names(res_list) <- rownames(logcounts(spe_sub))



output_dir <- here("corr_outputs",
                   sample_ids[s])
saveRDS(res_list, paste0(output_dir, "/non_zero_PAbeta.RDS"))

#
# PAbeta <- as.matrix(colData(spe_sub)$PAbeta)
# PpTau <- as.matrix(colData(spe_sub)$PpTau)
# NAbeta <- as.matrix(colData(spe_sub)$NAbeta)
# NpTau <- as.matrix(colData(spe_sub)$NpTau)
#
# n_genes <- nrow(logcounts(spe_sub))
# spe_counts <- logcounts(spe_sub)
#
#
