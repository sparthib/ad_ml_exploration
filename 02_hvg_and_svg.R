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

# fit mean-variance relationship
dec <- modelGeneVar(spe_postqc)

# visualize mean-variance relationship
fit <- metadata(dec)
plot(fit$mean, fit$var,
     xlab = "mean of log-expression", ylab = "variance of log-expression")
curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)

saveRDS()

####








