library(sgejobs)
sgejobs::job_single('all_samples',
                    create_shell = TRUE,
                    create_logdir = TRUE, task_num = 10,
                    queue = "bluejay",
                    cores = 4L,
                    memory = "20G")

library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(doParallel)
library(here)
library(sessioninfo)
library(dplyr)

#### cores ####


#numCores <- makeCluster(detectCores(), type='PSOCK') # grabs max available
#switch to 4 for now
numCores <- 4

options('mc.cores' = numCores)
registerDoParallel(numCores)

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
PpTau <- as.matrix(colData(spe_sub)$PpTau)
NAbeta <- as.matrix(colData(spe_sub)$NAbeta)
NpTau <- as.matrix(colData(spe_sub)$NpTau)

n_genes <- nrow(logcounts(spe_sub))
spe_counts <- logcounts(spe_sub)


#### Pearson Correlations ####

#PpTau
t0 = Sys.time()
res <- foreach(
    i = 1:n_genes,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel')
) %dopar% {
    cor(as.matrix(spe_counts[i, ]), PpTau, method = 'pearson')
}

t1 = Sys.time()
difftime(t1, t0, unit = "secs")

PpTau_pearson_res <- cbind(rowData(spe_sub)$gene_id,
                           res)

#PAbeta
t0 = Sys.time()
res <- foreach(
    i = 1:n_genes,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel')
) %dopar% {
    cor(as.matrix(spe_counts[i, ]), PAbeta, method = 'pearson')
}

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
#Time difference of 100.7417 secs
PAbeta_pearson_res <- cbind(rowData(spe_sub)$gene_id,
                            res)

#NpTau
t0 = Sys.time()
res <- foreach(
    i = 1:n_genes,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel')
) %dopar% {
    cor(as.matrix(spe_counts[i, ]), NpTau, method = 'pearson')
}

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
#Time difference of 103.4176 secs

NpTau_pearson_res <- cbind(rowData(spe_sub)$gene_id,
                           res)

#NAbeta
t0 = Sys.time()
res <- foreach(
    i = 1:n_genes,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel')
) %dopar% {
    cor(as.matrix(spe_counts[i, ]), NAbeta, method = 'pearson')
}

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
#Time difference of 105.4904 secs
NAbeta_pearson_res <- cbind(rowData(spe_sub)$gene_id,
                            res)



#### Spearman Correlations ####

#PpTau
t0 = Sys.time()
res <- foreach(
    i = 1:n_genes,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel')
) %dopar% {
    cor(as.matrix(spe_counts[i, ]), PpTau, method = 'spearman')
}

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
#Time difference of 112.114 secs

PpTau_spearman_res <- cbind(rowData(spe_sub)$gene_id,
                            res)

#PAbeta
t0 = Sys.time()
res <- foreach(
    i = 1:n_genes,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel')
) %dopar% {
    cor(as.matrix(spe_counts[i, ]), PAbeta, method = 'spearman')
}

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
#Time difference of 108.9839 secs
PAbeta_spearman_res <- cbind(rowData(spe_sub)$gene_id,
                             res)

#NpTau
t0 = Sys.time()
res <- foreach(
    i = 1:n_genes,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel')
) %dopar% {
    cor(as.matrix(spe_counts[i, ]), NpTau, method = 'spearman')
}

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
#Time difference of 109.7761 secs
NpTau_spearman_res <- cbind(rowData(spe_sub)$gene_id,
                            res)

#NAbeta
t0 = Sys.time()
res <- foreach(
    i = 1:n_genes,
    .combine = rbind,
    .multicombine = TRUE,
    .inorder = FALSE,
    .packages = c('data.table', 'doParallel')
) %dopar% {
    cor(as.matrix(spe_counts[i, ]), NAbeta, method = 'spearman')
}

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
# Time difference of 108.3586 secs
NAbeta_spearman_res <- cbind(rowData(spe_sub)$gene_id,
                             res)


#### clean and sort ####

clean_table <- function(df) {
    df <- as.data.frame(df)
    colnames(df) <- c("gene_id", "correlation")
    df$correlation <- as.numeric(df$correlation)
    df <- na.omit(df)
    df <- df |> dplyr::arrange(desc(correlation))
}

x <- list(
    PAbeta_pearson_res,
    PpTau_pearson_res,
    NAbeta_pearson_res ,
    NpTau_pearson_res,
    PAbeta_spearman_res,
    PpTau_spearman_res,
    NAbeta_spearman_res,
    NpTau_spearman_res
)

names <- c ("PAbeta_pearson",
            "PpTau_pearson",
            "NAbeta_pearson",
            "NpTau_pearson",
            "PAbeta_spearman",
            "PpTau_spearman",
            "NAbeta_spearman",
            "NpTau_spearman")

y_labs <- c("PAbeta", "PpTau", "NAbeta", "NpTau",
            "PAbeta", "PpTau", "NAbeta", "NpTau")
res <- lapply(x, clean_table)

dir.create(here("corr_outputs",
                sample_ids[s]))
output_dir <- here("corr_outputs",
                   sample_ids[s])
names(res) <- names
saveRDS(res, paste0(output_dir, "/sorted_corrs.RDS"))


### Get logcounts for top 100 genes.
for(i in seq_len(8)){
    n = nrow(res[[i]])
    gene_ids <- res[[i]]$gene_id[c(1:100, (n-99):n)] #get top and bottom 100 genes
    ix <- which(rowData(spe_sub)$gene_id %in% gene_ids)
    spe <- spe_sub[ix,]

    output_dir <- here("plots",
                       "02_sample_subset",
                       "sample_1")

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







