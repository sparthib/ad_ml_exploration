library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(doParallel)
library(BiocParallel)
library(here)
library(sessioninfo)
library(dplyr)

here()

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

s = 3
#s = as.numeric(Sys.getenv("SGE_TASK_ID"))


res_list <- readRDS(here('corr_outputs',
                         sample_ids[s],
                         'non_zero_PAbeta.RDS'))

spe_postqc <-
    readRDS(here::here("input_data",
                       paste0("spe_wholegenome_postqc.rds")))


print(s)
ix <- colData(spe_postqc)$sample_id_short == sample_ids[s]
spe_sub <- spe_postqc[, ix]



PAbeta <- as.matrix(colData(spe_sub)$PAbeta)
names(PAbeta) <- rownames(colData(spe_sub))

mat2 = PAbeta[PAbeta !=0]


non_zero_genes <- c()
non_zero_corrs <- c()
for( i in seq_len(length(res_list))){
    if(nrow(res_list[[i]]) >=10) {
        print(i)
        non_zero_genes <- c(non_zero_genes,names(res_list)[i] )
        non_zero_corrs <- c( non_zero_corrs,
                             cor(res_list[[i]][,'x'],res_list[[i]][,'y'], method = 'pearson'))
    }
}

non_zero_pabeta_df <- as.data.frame(cbind(non_zero_genes, non_zero_corrs))
colnames(non_zero_pabeta_df) <- c("gene_id", "corr_value")

non_zero_pabeta_df <- na.omit(non_zero_pabeta_df)
non_zero_pabeta_df$corr_value <- as.numeric(non_zero_pabeta_df$corr_value)
non_zero_pabeta_df <- non_zero_pabeta_df |> arrange(corr_value, desc = TRUE)




plot(res_list$ENSG00000269246[,'x'], res_list$ENSG00000269246[,'y'],
     )

output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "non_zero_scatter_PAbeta.pdf"))
for( i in seq_len(length(res_list))){

    if(nrow(res_list[[i]]) >=10) {
        par(mfrow=c(1,1))
        non_zero_genes <- c(non_zero_genes,names(res_list)[i] )
        plot(res_list[[i]][,'x'],res_list[[i]][,'y'],
             xlab=names(res_list)[i], ylab = "PAbeta")
    }
}
dev.off()



#bin them by 0.1 corr values (Absolute)
#within bin, rank by most spots.
#plot ones with most spots
# filter ones with largest range of y-axis


