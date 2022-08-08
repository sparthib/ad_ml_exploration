library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(BiocParallel)
library(doParallel)
library(here)
library(scran)
library(nnSVG)



spe_postqc <-
    readRDS(
        here::here(
            "input_data",
            paste0("spe_wholegenome_postqc.rds")
        )
    )

spe_postqc <- logNormCounts(spe_postqc)


####sample id names ####
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
print("sample id names listed")


#### get task id and subset spe ####
s = as.numeric(Sys.getenv("SGE_TASK_ID"))
ix <- colData(spe_postqc)$sample_id_short == sample_ids[s]
spe_sub <- spe_postqc[, ix]

print(sample_ids[s])


##### fit mean-variance relationship ####
dec <- modelGeneVar(spe_postqc)

# visualize mean-variance relationship
# fit <- metadata(dec)
# plot(fit$mean, fit$var,
#      xlab = "mean of log-expression", ylab = "variance of log-expression")
# curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)


saveRDS(top_hvgs, 'corr_outputs/top_hvgs.rds')

####

# run nnSVG once per sample-part within LC annotated regions
# and store lists of top SVGs
