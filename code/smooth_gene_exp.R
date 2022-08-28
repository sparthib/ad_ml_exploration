library(SpatialExperiment)
library(spatialLIBD)
library(doParallel)
library(here)
library(dplyr)
library(scuttle)
library(sessioninfo)
library(nnSVG)
library(BiocParallel)

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


spe_harmony <- readRDS(here::here("input_data",
                                  paste0("spe_harmony_wholegenome.rds")))



spe <-readRDS(here::here("input_data",
                       paste0("spe_wholegenome_postqc.rds")))


colData(spe)$col <- colData(spe_harmony)$col
colData(spe)$row <- colData(spe_harmony)$row
rm(spe_harmony)

ix <- colData(spe)$sample_id_short == sample_ids[s]
spe<- spe[, ix]
#dim(spe)
# [1] 27853  2665
spe <- nnSVG::filter_genes(spe,
                    filter_mito = FALSE)
spe<- scuttle::logNormCounts(spe)
#dim(spe)
# [1] 1990 2665

neighbors_list <- BayesSpace:::.find_neighbors(spe, "Visium")
# Neighbors were identified for 2664 out of 2665 spots.
# > neighbors_list[[1]]
# [1]  281  343  751 1705 2594
for(i in seq_len(length(neighbors_list))){

    neighbors_list[[i]]  <-  neighbors_list[[i]] + 1

}

# > neighbors_list[[1]]
# [1]  282  344  752 1706 2595

neighbor_vals <- list()

for (i in seq_len(length(neighbors_list))){
    if(length(neighbors_list[[i]]) >= 1){
        neighbor_vals[[i]] <- as.matrix(logcounts(spe)[, c(neighbors_list[[i]], i)])
    }
    else{ neighbor_vals[[i]] <- as.matrix(logcounts(spe)[,i])}
    neighbor_vals[[i]] <- 2 ** neighbor_vals[[i]]
    neighbor_vals[[i]] <- rowMeans(neighbor_vals[[i]])
    neighbor_vals[[i]] <- log(neighbor_vals[[i]], base = 2 )

}

smoothened_log_counts <- matrix(unlist(neighbor_vals), ncol = 2665, nrow = 1990)
rm(neighbor_vals)




rownames(smoothened_log_counts) <- rownames(spe)
colnames(smoothened_log_counts) <- colnames(spe)



saveRDS(smoothened_log_counts, here('input_data',
                                    'sample3_smoothened_log_counts.RDS'))


non_zero_genes_smoothened <-rowSums(smoothened_log_counts != 0)

non_zero_genes_raw <-rowSums(logcounts(spe)!= 0)
hist(non_zero_genes_raw, breaks=26)



####PAbeta spearman#######

PAbeta <- as.matrix(colData(spe)$PAbeta)
names(PAbeta) <- rownames(colData(spe))

mat2 = PAbeta[PAbeta !=0]


##for each gene, make a list of all spots that have the gene expressed and PAbeta

res_list <- list()
for( s in seq_len(nrow(smoothened_log_counts))){


    mat_nonzero <- which(smoothened_log_counts[s,] != 0, arr.ind = T)

   res_list[[s]]<- merge( smoothened_log_counts[s,][mat_nonzero], mat2, by = 0)


}

names(res_list) <- rownames(logcounts(spe))



output_dir <- here("corr_outputs",
                   sample_ids[3])
saveRDS(res_list, paste0(output_dir, "/smooth_non_zero_PAbeta.RDS"))


# > nrow(res_list[["ENSG00000159251"]])
# [1] 50
# > nrow(res_list[["ENSG00000175899"]])
# [1] 598
# > nrow(res_list[["ENSG00000133392"]])
# [1] 135
# > nrow(res_list[["ENSG00000107796"]])
# [1] 200
# > nrow(res_list[["ENSG00000159251"]])
# [1] 50
# > nrow(res_list[["ENSG00000198467"]])
# [1] 190

#Compare with non-smoothened version
# > nrow(raw_res_list[["ENSG00000159251"]])
# > nrow(raw_res_list[["ENSG00000175899"]])
# [1] 14
# [1] 166
# > nrow(raw_res_list[["ENSG00000133392"]])
# > [1] 39
# nrow(raw_res_list[["ENSG00000107796"]])
# [1] 55
# > nrow(raw_res_list[["ENSG00000159251"]])
# > [1] 14
# nrow(raw_res_list[["ENSG00000198467"]])
# [1] 47


raw_res_list <- readRDS(paste0(output_dir, "/non_zero_PAbeta.RDS"))


non_zero_genes <- c()
non_zero_corrs <- c()

nrows <- c()

for(i in seq_len(length(raw_res_list))) {
    nrows[i] <- nrow(raw_res_list[[i]])
}

nrows_df <- as.data.frame(cbind(names(raw_res_list), nrows))
colnames(nrows_df) <- c("gene_id", "nrows")
nrows_df$nrows <- as.numeric(nrows_df$nrows)

for( i in seq_len(length(raw_res_list))){
        non_zero_genes <- c(non_zero_genes,names(raw_res_list)[i] )
        non_zero_corrs <- c( non_zero_corrs,
                             cor(raw_res_list[[i]][,'x'],raw_res_list[[i]][,'y'], method = 'pearson'))
    }



######Corr df

non_zero_pabeta_df <- as.data.frame(cbind(non_zero_genes, non_zero_corrs))
colnames(non_zero_pabeta_df) <- c("gene_id", "corr_value")

non_zero_pabeta_df <- na.omit(non_zero_pabeta_df)
non_zero_pabeta_df$corr_value <- as.numeric(non_zero_pabeta_df$corr_value)

non_zero_pabeta_df$abs_corr <- abs(non_zero_pabeta_df$corr_value)
non_zero_pabeta_df <- non_zero_pabeta_df |> arrange(corr_value, desc = TRUE)

non_zero_pabeta_df <- non_zero_pabeta_df |> mutate(bin = case_when(
    abs_corr >= 0 & abs_corr < 0.1 ~ '1',
    abs_corr >= 0.1 & abs_corr < 0.2 ~ '2',
    abs_corr >= 0.2 & abs_corr < 0.3 ~ '3',
    abs_corr >= 0.3 & abs_corr < 0.4 ~ '4',
    abs_corr >= 0.4 & abs_corr < 0.5 ~ '5',
    abs_corr >= 0.5 & abs_corr < 0.6 ~ '6',
    abs_corr >= 0.6 & abs_corr < 0.7 ~ '7',
    abs_corr >= 0.7 & abs_corr < 0.8 ~ '8',
    abs_corr >= 0.8 & abs_corr < 0.9 ~ '9',
    abs_corr >= 0.9 & abs_corr <= 1.0  ~ '10',

))


non_zero_pabeta_df <- merge(non_zero_pabeta_df, nrows_df, by = "gene_id")

non_zero_groups <- non_zero_pabeta_df  |>
    group_by(bin)  |> top_n(n = 50, nrows) |> arrange(by = abs_corr, desc = TRUE)

# gene_id         corr_value abs_corr bin   nrows
# <chr>                <dbl>    <dbl> <chr> <dbl>
#     1 ENSG00000150938      0.416    0.416 5       470
# 2 ENSG00000172403      0.418    0.418 5       238
# 3 ENSG00000265972      0.421    0.421 5       508
# 4 ENSG00000188783      0.428    0.428 5       288
# 5 ENSG00000122786      0.444    0.444 5       438
# 6 ENSG00000118523      0.456    0.456 5       272
# 7 ENSG00000101335      0.458    0.458 5       445
# 8 ENSG00000149591      0.465    0.465 5       327
# 9 ENSG00000076706      0.482    0.482 5       261
# 10 ENSG00000162614      0.490    0.490 5        88
# 11 ENSG00000177469      0.492    0.492 5       207
# 12 ENSG00000077157      0.499    0.499 5       384
# arrange(dep_delay, .by_group = TRUE)

# gene_id         corr_value abs_corr bin   nrows
# <chr>                <dbl>    <dbl> <chr> <dbl>
#     1 ENSG00000090097      0.400    0.400 5        56
# 2 ENSG00000230590      0.401    0.401 5        48
# 3 ENSG00000150347      0.402    0.402 5        77
# 4 ENSG00000138036      0.405    0.405 5        48
# 5 ENSG00000150938      0.408    0.408 5       107
# 6 ENSG00000102024      0.409    0.409 5        80
# 7 ENSG00000187479      0.409    0.409 5        74
# 8 ENSG00000100029      0.410    0.410 5        41
# 9 ENSG00000159884      0.411    0.411 5        53
# 10 ENSG00000170776      0.412    0.412 5        50
# # … with 41 more rows

output_dir <- here("plots","02_sample_subset", sample_ids[3])
pdf(file = paste0(output_dir, "/", "smooth_non_zero_scatter_pearson_PAbeta.pdf"))
for( i in seq_len(nrow(non_zero_groups))){
    gene_id = non_zero_groups$gene_id[i]

    plot(res_list[[gene_id]]$x, res_list[[gene_id]]$y,
         xlab = gene_id, ylab = 'PAbeta',
         main= paste0('Absolute pearson corr = ',non_zero_groups$abs_corr[i])
    )

}
dev.off()


smoothened_log_counts <- readRDS(here('input_data',
                                    'sample3_smoothened_log_counts.RDS'))
informative_genes <- c("ENSG00000175899", "ENSG00000133392", "ENSG00000107796",
                       "ENSG00000159251", "ENSG00000198467")


ix_known <- which(rownames(smoothened_log_counts) %in% informative_genes)

predictor_vars <- t(smoothened_log_counts[ix_known, ])

df <- merge(predictor_vars, mat2, by = 0)

colnames(df)[7] <- "PAbeta"

df <- df[rowSums(df[2:6])>0,]
df <- df |> filter(PAbeta > 0)
nrow(df)
# > nrow(df)
# [1] 763
df <- as.data.frame(df)

model2 <- lm(PAbeta ~ ENSG00000175899 + ENSG00000133392 +ENSG00000107796+
                 ENSG00000159251 +ENSG00000198467 ,
             data = df)
summary(model2)

df$predicted_PAbeta <- predict(model2, newdata = df[,2:6])

ggplot(df, aes(x =PAbeta, y = predicted_PAbeta)) + geom_point() +
    xlab('observed percentage Abeta') + ylab('predicted percentage Abeta') +
    ggtitle('Linear regression results')



df <-df |> mutate(binary_pabeta = case_when(
    PAbeta <= 0.2 ~ 0,
    PAbeta > 0.2 ~ 1
) )

table(df$binary_pabeta)
# 0   1
# 535 128




model.log = glm(binary_pabeta ~ ENSG00000175899 + ENSG00000133392 +ENSG00000107796+
                    ENSG00000159251 +ENSG00000198467,
                data=df,
                family = binomial(link="logit")
)

summary(model.log)

df$predicted_prob <- predict(model.log,
                             type="response")


df<- df |> mutate(predicted_binary = case_when(
    predicted_prob <= 0.2 ~ 0,
    predicted_prob > 0.2 ~ 1 ))

table(df$predicted_binary)

#install.packages('caret')
library(caret)
cf <- confusionMatrix(data= binary_pabeta, reference = predicted_binary)

#Display results
cf

####Spearman correlation#####
non_zero_genes <- c()
non_zero_corrs <- c()
for( i in seq_len(length(res_list))){
    if( nrows_df$nrows[i]>=10) {
        print(i)
        non_zero_genes <- c(non_zero_genes,names(res_list)[i] )
        non_zero_corrs <- c( non_zero_corrs,
                             cor(res_list[[i]][,'x'],res_list[[i]][,'y'], method = 'spearman'))
    }
}



non_zero_pabeta_df <- as.data.frame(cbind(non_zero_genes, non_zero_corrs))
colnames(non_zero_pabeta_df) <- c("gene_id", "corr_value")

non_zero_pabeta_df <- na.omit(non_zero_pabeta_df)
non_zero_pabeta_df$corr_value <- as.numeric(non_zero_pabeta_df$corr_value)

non_zero_pabeta_df$abs_corr <- abs(non_zero_pabeta_df$corr_value)
non_zero_pabeta_df <- non_zero_pabeta_df |> arrange(corr_value, desc = TRUE)

non_zero_pabeta_df <- non_zero_pabeta_df |> mutate(bin = case_when(
    abs_corr >= 0 & abs_corr < 0.1 ~ '1',
    abs_corr >= 0.1 & abs_corr < 0.2 ~ '2',
    abs_corr >= 0.2 & abs_corr < 0.3 ~ '3',
    abs_corr >= 0.3 & abs_corr < 0.4 ~ '4',
    abs_corr >= 0.4 & abs_corr < 0.5 ~ '5',
    abs_corr >= 0.5 & abs_corr < 0.6 ~ '6',
    abs_corr >= 0.6 & abs_corr < 0.7 ~ '7',
    abs_corr >= 0.7 & abs_corr < 0.8 ~ '8',
    abs_corr >= 0.8 & abs_corr < 0.9 ~ '9',
    abs_corr >= 0.9 & abs_corr <= 1.0  ~ '10',

))


non_zero_pabeta_df <- merge(non_zero_pabeta_df, nrows_df, by = "gene_id")

non_zero_groups <- non_zero_pabeta_df  |>
    group_by(bin)  |> top_n(n = 50, nrows) |> arrange(by = abs_corr, desc = TRUE)
# arrange(dep_delay, .by_group = TRUE)

output_dir <- here("plots","02_sample_subset", sample_ids[3])
pdf(file = paste0(output_dir, "/", "smooth_non_zero_scatter_PAbeta_spearman.pdf"))
for( i in seq_len(nrow(non_zero_groups))){
    gene_id = non_zero_groups$gene_id[i]

    plot(res_list[[gene_id]]$x, res_list[[gene_id]]$y,
         xlab = gene_id, ylab = 'PAbeta',
         main= paste0('Absolute spearman corr = ',non_zero_groups$abs_corr[i])
    )

}
dev.off()
