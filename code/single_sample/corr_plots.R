# library(sgejobs)
# sgejobs::job_single('corr_plots',
#                     create_shell = TRUE,
#                     create_logdir = FALSE, task_num = 10,
#                     queue = "bluejay",
#                     cores = 4L,
#                     command = "Rscript corr_plots.R",
#                     memory = "20G")



library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(doParallel)
library(here)
library(sessioninfo)
library(dplyr)
library(scuttle)

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

#### load spe object and calc log counts ####
spe_postqc <-
    readRDS(here::here("input_data",
                       paste0("spe_wholegenome_postqc.rds")))
print("spe object loaded")

spe_postqc <- scuttle::logNormCounts(spe_postqc)
spe_counts <- logcounts(spe_postqc)

print("log counts for spe calculated")

#### get task id and subset spe ####
s = as.numeric(Sys.getenv("SGE_TASK_ID"))
print(s)
ix <- colData(spe_postqc)$sample_id_short == sample_ids[s]
spe_sub <- spe_postqc[, ix]

print(sample_ids[s])


#load corr values
output_dir <- here("corr_outputs",
                   sample_ids[s])

res <- readRDS(paste0(output_dir, "/sorted_corrs.RDS"))

names <- c ("PAbeta_pearson",
            "PpTau_pearson",
            "NAbeta_pearson",
            "NpTau_pearson",
            "PAbeta_spearman",
            "PpTau_spearman",
            "NAbeta_spearman",
            "NpTau_spearman")

names(res) <- names

dir.create(here("plots","02_sample_subset", sample_ids[s]))

y_labs <- c("PAbeta", "PpTau", "NAbeta", "NpTau",
            "PAbeta", "PpTau", "NAbeta", "NpTau")

for(i in seq_len(8)){
    n = nrow(res[[i]])
    gene_ids <- res[[i]]$gene_id[c(1:100, (n-99):n)] #get top and bottom 100 genes
    ix <- which(rowData(spe_sub)$gene_id %in% gene_ids)
    spe <- spe_sub[ix,]

    output_dir <- here("plots","02_sample_subset", sample_ids[s])
    pdf(file = paste0(output_dir, "/", names[i], ".pdf"))

    for (j in seq_len(200)){
        par(mfrow=c(1,1))

        if(i == 1 || i == 5){
            plot(as.vector(logcounts(spe)[j,]), as.vector(colData(spe)$PAbeta),
                 xlab = paste0(rowData(spe)$gene_name[j], " ",
                               rowData(spe)$gene_id[j]," ", "expression"),
                 ylab = paste0(y_labs[i]))
        }

        #dev.off()

        if(i == 2 || i == 6){
            plot(as.vector(logcounts(spe)[j,]), as.vector(colData(spe)$PpTau),
                 xlab = paste0(rowData(spe)$gene_name[j], " ",
                               rowData(spe)$gene_id[j]," ", "expression"),
                 ylab = paste0(y_labs[i]))

        }
        #dev.off()


        if(i == 3 || i == 7){
            plot(as.vector(logcounts(spe)[j,]), as.vector(colData(spe)$NAbeta),
                 xlab = paste0(rowData(spe)$gene_name[j], " ",
                               rowData(spe)$gene_id[j]," ", "expression"),
                 ylab = paste0(y_labs[i]))
        }

        #dev.off()


        if(i == 4 || i == 8){
            plot(as.vector(logcounts(spe)[j,]), as.vector(colData(spe)$NpTau),
                 xlab = paste0(rowData(spe)$gene_name[j], " ",
                               rowData(spe)$gene_id[j]," ", "expression"),
                 ylab = paste0(y_labs[i]))

        }
        #dev.off()

    }

    dev.off()

}







