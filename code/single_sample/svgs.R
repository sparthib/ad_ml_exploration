# library(sgejobs)
# sgejobs::job_single('svgs',
#                     create_shell = TRUE,
#                     create_logdir = FALSE, task_num = 10,
#                     cores = 5L,
#                     queue = "bluejay",
#                     command = "Rscript svgs.R",
#                     memory = "20G")
#


library(SpatialExperiment)
library(doParallel)
library(here)
library(dplyr)
library(scuttle)
library(sessioninfo)
library(nnSVG)
library(BiocParallel)
library(doParallel)
here()

#### load spe object and calc log counts ####
spe_postqc <-
    readRDS(here::here("input_data",
                       "spe_wholegenome_postqc.rds"))

print(class(spe_postqc))

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


# run nnSVG filtering for mitochondrial gene and low-expressed genes
spe_sub <- filter_genes(
    spe_sub,
    filter_genes_ncounts = 2,
    filter_genes_pcspots = 1
)

# re-calculate logcounts after filtering
spe_sub <- logNormCounts(spe_sub)


# run nnSVG
set.seed(123)
spe_sub <- nnSVG(spe_sub, n_threads = 5)

# store results
res_list <- rowData(spe_sub)


# ------------
# save results
# ------------

# save nnSVG results
output_dir <- here('corr_outputs',
                   sample_ids[s])
saveRDS(res_list, paste0(output_dir, "/top_svgs.RDS"))
