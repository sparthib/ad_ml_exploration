# library(sgejobs)
# sgejobs::job_single('hvgs',
#                     create_shell = TRUE,
#                     create_logdir = FALSE, task_num = 10,
#                     queue = "bluejay",
#                     command = "Rscript hvgs.R",
#                     memory = "20G")



library(SpatialExperiment)
library(doParallel)
library(here)
library(dplyr)
library(scuttle)
library(scran)
library(sessioninfo)
library(nnSVG)
here()

#### load spe object and calc log counts ####
spe_postqc <-
    readRDS(here::here("input_data",
                       "spe_wholegenome_postqc.rds"))

print(class(spe_postqc))

spe_postqc <- scuttle::logNormCounts(spe_postqc)
print("logcounts calculated ")

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

# fit mean-variance relationship
dec <- scran::modelGeneVar(spe_sub)

# visualize mean-variance relationship
# fit <- metadata(dec)
# plot(fit$mean, fit$var,
#      xlab = "mean of log-expression", ylab = "variance of log-expression")
# curve(fit$trend(x), col = "dodgerblue", add = TRUE, lwd = 2)

top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)

output_dir <- here('corr_outputs',
                   sample_ids[s])
saveRDS(top_hvgs, paste0(output_dir, "/top_hvgs.RDS"))


