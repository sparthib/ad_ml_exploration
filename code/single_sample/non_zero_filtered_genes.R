library(ggplot2)

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



spe <-readRDS(here::here("input_data",
                         paste0("spe_wholegenome_postqc.rds")))



ix <- colData(spe)$sample_id_short == sample_ids[s]
spe<- spe[, ix]

spe <- nnSVG::filter_genes(spe,
                           filter_mito = FALSE)
spe<- scuttle::logNormCounts(spe)



# ENSG00000133056 + ENSG00000130193 + ENSG00000198887 +
#     ENSG00000165028 + ENSG00000168781 + ENSG00000161298 + ENSG00000273492



####PAbeta spearman#######

PAbeta <- as.matrix(colData(spe)$PAbeta)
names(PAbeta) <- rownames(colData(spe))

mat2 = PAbeta[PAbeta !=0]

res_list <- list()
for( s in seq_len(nrow(logcounts(spe)))){


    mat_nonzero <- which(logcounts(spe)[s,] != 0, arr.ind = T)

    res_list[[s]]<- merge( logcounts(spe)[s,][mat_nonzero], mat2, by = 0)


}


names(res_list) <- rownames(logcounts(spe))



output_dir <- here("corr_outputs",
                   sample_ids[3])
saveRDS(res_list, paste0(output_dir, "/filtered_non_zero_PAbeta.RDS"))


non_zero_genes <- c()
non_zero_corrs <- c()

nrows <- c()

for(i in seq_len(length(res_list))) {
    nrows[i] <- nrow(res_list[[i]])
}

nrows_df <- as.data.frame(cbind(names(res_list), nrows))
colnames(nrows_df) <- c("gene_id", "nrows")
nrows_df$nrows <- as.numeric(nrows_df$nrows)

for( i in seq_len(length(res_list))){
    non_zero_genes <- c(non_zero_genes,names(res_list)[i] )
    non_zero_corrs <- c( non_zero_corrs,
                         cor(res_list[[i]][,'x'],res_list[[i]][,'y'], method = 'pearson'))
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
# arrange(dep_delay, .by_group = TRUE)

output_dir <- here("plots","02_sample_subset", sample_ids[3])
pdf(file = paste0(output_dir, "/", "filtered_non_zero_scatter_pearson_PAbeta.pdf"))
for( i in seq_len(nrow(non_zero_groups))){
    gene_id = non_zero_groups$gene_id[i]

    plot(res_list[[gene_id]]$x, res_list[[gene_id]]$y,
         xlab = gene_id, ylab = 'PAbeta',
         main= paste0('Absolute pearson corr = ',non_zero_groups$abs_corr[i])
    )

}
dev.off()



informative_genes <- c("ENSG00000175899", "ENSG00000133392", "ENSG00000107796",
                       "ENSG00000159251", "ENSG00000198467")


ix_known <- which(rowData(spe)$gene_id %in% informative_genes)

spe_sub <- spe[ix_known, ]

predictor_vars <- t(as.matrix(logcounts(spe_sub)))
target_var <- as.matrix(colData(spe_sub)$PAbeta)

df <- cbind(predictor_vars, target_var)
colnames(df) <- c(colnames(df)[1:5], "PAbeta")

df <- as.data.frame(df)
df <- df[rowSums(df[])>0,]
df <- df |> dplyr::filter(PAbeta > 0)
nrow(df)
#763

df <-df |> dplyr::mutate(binary_pabeta = case_when(
    PAbeta <= 0.2 ~ 0,
    PAbeta > 0.2 ~ 1
) )

table(df$binary_pabeta)
# 0   1
# 616 147
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
# 0   1
# 684  79

#install.packages('caret')
library(caret)
cf <- confusionMatrix(data= as.factor(df$predicted_binary), reference = as.factor(df$binary_pabeta))

#Display results
cf


install.packages('pROC')
library(pROC)
roc_score=roc(test_data[,13], logit_P) #AUC score
plot(roc_score ,main ="ROC curve -- Logistic Regression ")


##########

model2 <- lm(PAbeta ~ ENSG00000175899 + ENSG00000133392 +ENSG00000107796+
             ENSG00000159251 +ENSG00000198467 ,
             data = df)
summary(model2)
df$predicted_PAbeta <- predict(model2, newdata = df[,1:5])

par(mar = c(1, 1, 1, 1))

ggplot(df, aes(x =PAbeta, y = predicted_PAbeta)) + geom_point() +
    xlab('observed percentage Abeta') + ylab('predicted percentage Abeta') +
    ggtitle('Linear regression results')
