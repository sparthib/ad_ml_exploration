library(spatialLIBD)
library(patchwork)
library(ggplot2)
library(dplyr)
library(here)
library(scuttle)
library(nnSVG)
#install.packages("gridGraphics")
#library(gridGraphics)
install.packages("nnet")
library("nnet")


################

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

#s = as.numeric(Sys.getenv("SGE_TASK_ID"))
s = 2
print(sample_ids[s])



spe <-readRDS(here::here("input_data",
                         paste0("spe_wholegenome_postqc.rds")))



ix <- colData(spe)$sample_id_short == sample_ids[s]
spe<- spe[, ix]

rownames(spe) <- rowData(spe)$gene_name
spe <- nnSVG::filter_genes(spe,filter_genes_ncounts = 3,
                           filter_mito = FALSE)

spe<- scuttle::logNormCounts(spe)


###########

top_genes_NAbeta<- c("PHACTR1", "IL6ST", "DNAJB61", "SST")
# top_15_genes_PpTau <- c("ARHGEF7" , "IL6ST",    "ADCY2" ,   "UNC13A",   "GDA" ,
#                         "CCDC91" ,  "SERPINE1", "IK",   "CNN1",
#                         "PRELP" ,   "NEXN" ,    "MGP" ,     "CARNS1",  "MYH11")
#
# intersect(top_15_genes_PAbeta, top_15_genes_PpTau)
# "CNN1"  "MYH11" "MGP"

ad_related_genes <- c("CLU", "MEF2C", "SORL1", "BIN1", "APOE")




###### spe with rowname = gene name and logcounts available ####

ix_known <- which(rowData(spe)$gene_name %in% top_genes_NAbeta)

spe_sub <- spe[ix_known, ]


predictor_vars <- t(as.matrix(logcounts(spe_sub)))
target_var <- as.matrix(colData(spe_sub)$NAbeta)

df <- cbind(predictor_vars, target_var)

colnames(df) <- c("SST",  "IL6ST", "PHACTR1", "NAbeta")

df <- as.data.frame(df)
df <- df[rowSums(df[])>0,]
nrow(df) #[1] 1790
# [1] 2004
df <- df |> dplyr::filter(NAbeta > 0)
nrow(df)
#327

#### Multinom  Regression ####

models_summary <- list()
models_res <- list()
models_fitted <- list()


for (i in seq_len(3)){

    model_lr = nnet::multinom(as.formula(paste0("NAbeta ~", colnames(df)[i])), data = df)
    models_summary[[i]] <- summary(model_lr)
    models_res[[i]] <- resid(model_lr)
    models_fitted[[i]] <- fitted(model_lr)


}

model <- nnet::multinom(as.formula( "NAbeta ~ SST + IL6ST + PHACTR1"), data = df)
fitted(model)
summary(model)

names(models_summary) <- colnames(df)[1:3]
names(models_res) <- colnames(df)[1:3]
names(models_fitted) <- colnames(df)[1:3]


#####residual plots
output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "PAbeta_Residual_plots_simple_lin_reg_top_corr_genes.pdf"))

gene_names <- c("SST", "IL6ST", "PHACTR1")
for(i in seq_len(3)){
    par(mar = c(4,4,4,4))
    p <- plot(models_fitted[[i]], models_res[[i]], font.axis=1, cex.axis=1,
              xlab = 'Fitted', ylab = 'Residual', main = paste0("Fitted vs Residual, ",
                                                                gene_names))
}
dev.off()

#####Rsq for all plots
rsq <- list()
for(i in seq_len(3)){
    rsq[[i]] <- models_summary[[i]]$r.squared
}


#plot pathology and gene expression
colData(spe_sub)$sample_name <- colData(spe_sub)$sample_id_short

output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "spot_plots_gene_vs_abeta.pdf"),
    width = 10)

for(i in seq_len(gene_names[1:3])){
    p <- spatialLIBD::vis_gene(
        spe,
        sampleid = sample_ids[s],
        geneid = "SST",
        spatial = TRUE,
        assayname = "logcounts",
        minCount = 0,
        viridis = TRUE,
        image_id = "lowres",
        alpha = 1,
        point_size = 1.25
    ) +ggtitle("SST")

    q <- vis_gene(
        spe_sub,
        sampleid = sample_ids[s],
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

output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "spot_plots_ad_gene_vs_abeta.pdf"),
    width = 10)
for(i in seq_len(length(top_15_genes_PAbeta))){
    p <- vis_gene(
        spe,
        sampleid = sample_ids[s],
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
        sampleid = sample_ids[s],
        geneid = "PAbeta",
        spatial = TRUE,
        assayname = "logcounts",
        minCount = 0,
        viridis = TRUE,
        image_id = "lowres",
        alpha = 1,
        point_size = 1.25
    )+ggtitle("PAbeta")
    print(p | q)
}


dev.off()

rsq_df <- as.data.frame(cbind(top_15_genes_PAbeta, rsq))

ggplot(rsq_df , aes(x=factor(0), y=rsq)) + #+
    geom_boxplot() # +
    geom_text(data = rsq_df,
              aes(x = factor(0), y = rsq, label = top_15_genes_PAbeta),
              nudge_x = .5)


output_dir <- here("plots","02_sample_subset", sample_ids[s])
p <- spatialLIBD::vis_grid_gene(
    spe,
    geneid = "SST",
    spatial = TRUE,
    assayname = "logcounts",
    minCount = 0,
    image_id = "lowres",
    alpha = 1,
    viridis = TRUE,
    point_size = 3,
    pdf = paste0(output_dir, "/", "ad_gene_vs_Nabeta.pdf")
    )
dev.off()














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



