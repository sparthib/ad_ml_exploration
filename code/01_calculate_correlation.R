# library(sgejobs)
# sgejobs::job_single(
#             "calculate_correlation",
#             create_shell = TRUE,
#             queue = "bluejay",
#             cores = 15L,
#             memory = "20G",
#             command = "Rscript 01_calculate_correlation.R",
#             create_logdir = TRUE
#         )
library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(BiocParallel)
library(doParallel)
library(here)
library(sessionInfo)
#BiocManager::install("gpuMagic")

numCores <- makeCluster(detectCores(), type='PSOCK') # grabs max available

options('mc.cores' = numCores)
registerDoParallel(numCores)


n_genes = 27853
# snowFORK <- SnowParam(workers = numCores, type = "FORK")

spe_postqc <-
    readRDS(
        here::here(
            "input_data",
            paste0("spe_wholegenome_postqc.rds")
        )
    )

pAbeta <- as.matrix(colData(spe_postqc)$PAbeta)
corr_for_counts <- function(x){
    res <- cor(
        x,
        pAbeta, method = "spearman")
    res

    }

t0 = Sys.time()
#beta = pvec(1:10, corr_for_counts, spe = spe_postqc, mc.cores = numCores)
#
# beta = apply(X = as.matrix(counts(spe_postqc)[1:2000,]),
#             MARGIN = 1,
#             FUN = corr_for_counts)

spe_counts <- counts(spe_postqc)

res <- foreach(i = seq_len(n_genes),
               .combine = rbind,
               .multicombine = TRUE,
               .inorder = FALSE,
               .packages = c('data.table', 'doParallel')) %dopar% {
                   cor(as.matrix(spe_counts[i,]), pAbeta, method = 'pearson')
               }


# > beta = pvec(1:n_genes, corr_for_counts, spe = spe_postqc, mc.cores = numCores)

t1 = Sys.time()
difftime(t1, t0, unit="secs")

print(res[1:10])



# corrs <-apply(as.matrix(counts(spe_S1_A1_Br3874)[i])
#       , 1, function(x) cor(x, as.numeric(as.matrix(colData(spe_S1_A1_Br3874)$NAbeta))))
