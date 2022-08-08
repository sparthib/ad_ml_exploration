library(spatialLIBD)
library(patchwork)
install.packages("gridGraphics")
library(gridGraphics)


################

top_15_genes_PAbeta<- c("MYL9",  "CNN1", "TPM1", "CAVIN1", "CALD1", "MYH11",
                        "THBS1" , "ACTA2" , "PPP1R12B", "FN1", "TAGLN", "MGP",
                        "ACTC1" ,  "FLNA",  "TPM2")
top_15_genes_PpTau <- c("ARHGEF7" , "IL6ST",    "ADCY2" ,   "UNC13A",   "GDA" ,
                        "CCDC91" ,  "SERPINE1", "IK",   "CNN1",
                        "PRELP" ,   "NEXN" ,    "MGP" ,     "CARNS1",  "MYH11")

intersect(top_15_genes_PAbeta, top_15_genes_PpTau)
# "CNN1"  "MYH11" "MGP"

ad_related_genes <- c("CLU", "MEF2C", "SORL1", "BIN1", "APOE")




###### spe with rowname = gene name and logcounts available ####

ix_known <- which(rowData(spe)$gene_name %in% top_15_genes_PAbeta)

spe_sub <- spe[ix_known, ]

predictor_vars <- t(as.matrix(logcounts(spe_sub)))
target_var <- as.matrix(colData(spe_sub)$PAbeta)

df <- cbind(predictor_vars, target_var)
colnames(df) <- c(colnames(df)[1:15], "PAbeta")

df <- as.data.frame(df)
df <- df[rowSums(df[])>0,]
# nrow(df)
# [1] 2004
df <- df |> dplyr::filter(PAbeta > 0)
nrow(df)
#763

df <-df |> dplyr::mutate(binary_pabeta = case_when(
    PAbeta <= 0.2 ~ 0,
    PAbeta > 0.2 ~ 1
) )

####Linear regression

models_summary <- list()
models_res <- list()
models_fitted <- list()

for (i in seq_len(length(top_15_genes_PAbeta))){

    model_lr = lm(as.formula(paste0("PAbeta ~", top_15_genes_PAbeta[i])), data = df)
    models_summary[[i]] <- summary(model_lr)
    models_res[[i]] <- resid(model_lr)
    models_fitted[[i]] <- fitted(model_lr)


}

names(models_summary) <- top_15_genes_PAbeta
names(models_res) <- top_15_genes_PAbeta
names(models_fitted) <- top_15_genes_PAbeta


#####residual plots
output_dir <- here("plots","02_sample_subset", sample_ids[3])
pdf(file = paste0(output_dir, "/", "PAbeta_Residual_plots_simple_lin_reg_top_corr_genes.pdf"))

for(i in seq_len(length(top_15_genes_PAbeta))){
    par(mar = c(4,4,4,4))
    p <- plot(models_fitted[[i]], models_res[[i]], font.axis=1, cex.axis=1,
              xlab = 'Fitted', ylab = 'Residual', main = paste0("Fitted vs Residual, ",
                                                                top_15_genes_PAbeta[i]))
}
dev.off()
#####Rsq for all plots
rsq <- list()
for(i in seq_len(length(top_15_genes_PAbeta))){
    rsq[[i]] <- models_summary[[i]]$r.squared
}


#plot pathology and gene expression
colData(spe_sub)$sample_name <- colData(spe_sub)$sample_id

output_dir <- here("plots","02_sample_subset", sample_ids[3])
pdf(file = paste0(output_dir, "/", "spot_plots_gene_vs_abeta.pdf"),
    width = 10)
for(i in seq_len(length(top_15_genes_PAbeta))){
    p <- vis_gene(
        spe_sub,
        sampleid = "V10A27106_C1_Br3873",
        geneid = top_15_genes_PAbeta[i],
        spatial = TRUE,
        assayname = "logcounts",
        minCount = 0,
        viridis = TRUE,
        image_id = "lowres",
        alpha = 1,
        point_size = 1.25
    ) +ggtitle(top_15_genes_PAbeta[i])

    q <- vis_gene(
        spe_sub,
        sampleid = "V10A27106_C1_Br3873",
        geneid = "PAbeta",
        spatial = TRUE,
        assayname = "logcounts",
        minCount = 0,
        viridis = TRUE,
        image_id = "lowres",
        alpha = 1,
        point_size = 1.25
    )+ggtitle("PAbeta")
    par(mar = c(4,4,4,4))
    print(p | q)
}


dev.off()


ad_related_genes <- c("CLU", "MEF2C", "SORL1", "BIN1", "APOE")
colData(spe)$sample_name <- colData(spe)$sample_id
spe<- logNormCounts(spe)

output_dir <- here("plots","02_sample_subset", sample_ids[3])
pdf(file = paste0(output_dir, "/", "spot_plots_ad_gene_vs_abeta.pdf"),
    width = 10)
for(i in seq_len(length(ad_related_genes))){
    p <- vis_gene(
        spe,
        sampleid = "V10A27106_C1_Br3873",
        geneid = ad_related_genes[i],
        spatial = TRUE,
        assayname = "logcounts",
        minCount = 0,
        viridis = TRUE,
        image_id = "lowres",
        alpha = 1,
        point_size = 1.25
    ) +ggtitle(ad_related_genes[i])

    q <- vis_gene(
        spe,
        sampleid = "V10A27106_C1_Br3873",
        geneid = "PAbeta",
        spatial = TRUE,
        assayname = "logcounts",
        minCount = 0,
        viridis = TRUE,
        image_id = "lowres",
        alpha = 1,
        point_size = 1.25
    )+ggtitle("PAbeta")
    par(mar = c(4,4,4,4))
    print(p | q)
}


dev.off()

rsq_df <- as.data.frame(cbind(top_15_genes_PAbeta, rsq))

ggplot(rsq_df , aes(x=factor(0), y=rsq)) + #+
    geom_boxplot() # +
    geom_text(data = rsq_df,
              aes(x = factor(0), y = rsq, label = top_15_genes_PAbeta),
              nudge_x = .5)


















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



