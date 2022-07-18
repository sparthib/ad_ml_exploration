# library(sgejobs)
# sgejobs::job_single(
#             "calculate_correlation",
#             create_shell = TRUE,
#             queue = "bluejay",
#             cores = 15L,
#             memory = "10G",
#             command = "Rscript 01_calculate_correlation.R",
#             create_logdir = TRUE
#         )
library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(BiocParallel)
library(doParallel)
library(here)
library(sessioninfo)
#BiocManager::install("gpuMagic")

#numCores <- makeCluster(detectCores(), type='PSOCK') # grabs max available
#switch to 4 for now
numCores <- 4

options('mc.cores' = numCores)
registerDoParallel(numCores) #specifies that you want 8
#workers working every time you do a foreach loop

#EXAMPLE
#res <- foreach(i = 1:5 ) %dopar% {
   # i*i
#}

# i specifies the sequence of your iteration
# %dopar% means run in parallel
# assign everything to res

#test time

# system.time(
#     res <- foreach(i = 1:100000) %dopar%
#         {
#             i*i
#
#         }
# )

#for dopar
# user  system elapsed
# 20.316   2.427  24.209

#for do
# user  system elapsed
# 0.101   0.002   0.105

##got warning using the do instead of dopar

# total user time was 0.005 seconds

n_genes = 27853
# snowFORK <- SnowParam(workers = numCores, type = "FORK")

spe_postqc <-
    readRDS(
        here::here(
            "input_data",
            paste0("spe_wholegenome_postqc.rds")
        )
    )

spe_postqc <- scuttle::logNormCounts(spe_postqc)


PAbeta <- as.matrix(colData(spe_postqc)$PAbeta)
PpTau <- as.matrix(colData(spe_postqc)$PpTau)



#beta = pvec(1:10, corr_for_counts, spe = spe_postqc, mc.cores = numCores)
#
# beta = apply(X = as.matrix(counts(spe_postqc)[1:2000,]),
#             MARGIN = 1,
#             FUN = corr_for_counts)


spe_counts <- logcounts(spe_postqc)
rm(spe_postqc)


t0 = Sys.time()
res <- foreach(i = 1:27853,
               .combine = rbind,
               .multicombine = TRUE,
               .inorder = FALSE,
               .packages = c('data.table', 'doParallel')) %dopar% {
                   cor(as.matrix(spe_counts[i,]), PpTau, method = 'pearson')
               }

t1 = Sys.time()
##


# > beta = pvec(1:n_genes, corr_for_counts, spe = spe_postqc, mc.cores = numCores)


difftime(t1, t0, unit="secs")
# Time difference of 3314.799 secs

print(res[1:10])

PpTau_res <- cbind(rowData(spe_postqc)$gene_id,
                    res)
saveRDS(PpTau_res,
        here("corr_outputs",
             "pearson",
             "pearson_PpTau.RDS"))







