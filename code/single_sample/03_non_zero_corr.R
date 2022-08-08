#library(sgejobs)
# sgejobs::job_single('non_zero_corr',
#                     create_shell = TRUE,
#                     create_logdir = TRUE, task_num = 10,
#                     queue = "bluejay",
#                     cores = 5L,
#                     command = "Rscript non_zero_corr.R",
#                     memory = "20G")
#



##This script computes all combos of correlation values for each of the 10 samples.
library(SpatialExperiment)
library(rtracklayer)
library(GenomicRanges)
library(doParallel)
library(BiocParallel)
library(here)
library(sessioninfo)
library(dplyr)

#### cores ####


#numCores <- makeCluster(detectCores(), type='PSOCK') # grabs max available
#switch to 4 for now
# numCores <- 4
#
# options('mc.cores' = numCores)
# registerDoParallel(numCores)

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
names(PAbeta) <- rownames(colData(spe_sub))

mat2 = PAbeta[PAbeta !=0]
i = 1
res_list <- list()


# for (i in seq_len(nrow(logcounts(spe_sub)))){
#     #mat1 : for a given gene, all spots with non-zero expression
#     mat1 <- as.matrix(logcounts(spe_sub))[i,][as.matrix(logcounts(spe_sub))[i,] != 0]
#     if(length(mat1)!=0){
#         #if there remain some spots with non-zero expression
#         #samecols <- intersect(names(mat1),names(mat2))
#         res_list[[length(res_list)+1]] <- merge(mat1,mat2, by = 0)
#         gene_names <- c(gene_names, rownames(logcounts(spe_sub))[i])
#
#
#     }
# }



fun <- function(i) {
    mat1 <- as.matrix(logcounts(spe_sub))[i,][as.matrix(logcounts(spe_sub))[i,] != 0]
    merge(mat1,mat2, by = 0)
}

t0 = Sys.time()
res_list <- bplapply(1:nrow(logcounts(spe_sub)), fun, BPPARAM = MulticoreParam(5))

t1 = Sys.time()
difftime(t1, t0, unit = "secs")
names(res_list) <- rownames(logcounts(spe_sub))



output_dir <- here("corr_outputs",
                   sample_ids[s])
saveRDS(res_list, paste0(output_dir, "/non_zero_PAbeta.RDS"))

#
# PAbeta <- as.matrix(colData(spe_sub)$PAbeta)
# PpTau <- as.matrix(colData(spe_sub)$PpTau)
# NAbeta <- as.matrix(colData(spe_sub)$NAbeta)
# NpTau <- as.matrix(colData(spe_sub)$NpTau)
#
# n_genes <- nrow(logcounts(spe_sub))
# spe_counts <- logcounts(spe_sub)
#


res_list <- readRDS(here('corr_outputs',
                         sample_ids[s],
                         'non_zero_PAbeta.RDS'))


###### Calculate Pearson ######
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
    if( nrows_df$nrows[i]>=10) {
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

output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "non_zero_scatter_PAbeta.pdf"))
for( i in seq_len(nrow(non_zero_groups))){
    gene_id = non_zero_groups$gene_id[i]

    plot(res_list[[gene_id]]$x, res_list[[gene_id]]$y,
         xlab = gene_id, ylab = 'PAbeta',
         main= paste0('Absolute pearson corr = ',non_zero_groups$abs_corr[i])
    )

}
dev.off()




#bin them by 0.1 corr values (Absolute)
#within bin, rank by most spots.
#plot ones with most spots
# filter ones with largest range of y-axis



informative_genes <- c("ENSG00000130193", "ENSG00000165028", "ENSG00000273492",
                       "ENSG00000198887", "ENSG00000168781", "ENSG00000133056",
                       "ENSG00000161298")


ix_known <- which(rowData(spe_sub)$gene_id %in% informative_genes)

spe_sub <- spe_sub[ix_known, ]

predictor_vars <- t(as.matrix(logcounts(spe_sub)))
target_var <- as.matrix(colData(spe_sub)$PAbeta)

df <- cbind(predictor_vars, target_var)
colnames(df) <- c(colnames(df)[1:7], "PAbeta")



nrow(df)
df <- as.data.frame(df)
df <- df[rowSums(df[])>0,]
df <- df |> filter(PAbeta > 0)
#1068

model2 <- lm(PAbeta ~ ENSG00000133056 + ENSG00000130193 + ENSG00000198887 +
                 ENSG00000165028 + ENSG00000168781 + ENSG00000161298 + ENSG00000273492,
             data = df
)
summary(model2)
df$predicted_PAbeta <- predict(model2, newdata = df[,1:7])

ggplot(df, aes(x =PAbeta, y = predicted_PAbeta)) + geom_point() +
    xlab('observed percentage Abeta') + ylab('predicted percentage Abeta') +
    ggtitle('Linear regression results')




df <-df |> dplyr::mutate(binary_pabeta = case_when(
    PAbeta <= 0.2 ~ 0,
    PAbeta > 0.2 ~ 1
) )

table(df$binary_pabeta)
# 0   1
# 616 147
model.log = glm(binary_pabeta ~ ENSG00000133056 + ENSG00000130193 + ENSG00000198887 +
                    ENSG00000165028 + ENSG00000168781 + ENSG00000161298 + ENSG00000273492,
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





###### Calculate Spearman ######
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

output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "non_zero_scatter_PAbeta_spearman.pdf"))
for( i in seq_len(nrow(non_zero_groups))){
    gene_id = non_zero_groups$gene_id[i]

    plot(res_list[[gene_id]]$x, res_list[[gene_id]]$y,
         xlab = gene_id, ylab = 'PAbeta',
         main= paste0('Absolute spearman corr = ',non_zero_groups$abs_corr[i])
    )

}
dev.off()



informative_genes <- c("ENSG00000130193", "ENSG00000165028", "ENSG00000273492",
                       "ENSG00000198887", "ENSG00000168781", "ENSG00000133056",
                       "ENSG00000161298")


ix_known <- which(rowData(spe_sub)$gene_id %in% informative_genes)

spe_sub <- spe_sub[ix_known, ]

predictor_vars <- t(as.matrix(logcounts(spe_sub)))
target_var <- as.matrix(colData(spe_sub)$PAbeta)

df <- cbind(predictor_vars, target_var)
colnames(df) <- c(colnames(df)[1:7], "PAbeta")



nrow(df)
df <- as.data.frame(df)
df <- df[rowSums(df[])>0,]
df <- df |> filter(PAbeta > 0)
#1068

model2 <- lm(PAbeta ~ ENSG00000133056 + ENSG00000130193 + ENSG00000198887 +
                 ENSG00000165028 + ENSG00000168781 + ENSG00000161298 + ENSG00000273492,
             data = df
)
summary(model2)
df$predicted_PAbeta <- predict(model2, newdata = df[,1:7])



par(mar = c(1, 1, 1, 1))
plot(df$PAbeta, df$predicted_PAbeta)

svmfit = svm(PAbeta ~ ENSG00000133056 + ENSG00000130193 + ENSG00000198887 +
                 ENSG00000165028 + ENSG00000168781 + ENSG00000161298 + ENSG00000273492,
             data = df, kernel = "linear", cost = 10, scale = FALSE)


#### Logistic Regression
df <-df |> mutate(binary_pabeta = case_when(
    PAbeta <= 0.2 ~ 0,
    PAbeta > 0.2 ~ 1
) )

table(df$binary_pabeta)
# 0   1
# 921 147


model.log = glm(binary_pabeta ~ ENSG00000133056 + ENSG00000130193 + ENSG00000198887 +
                    ENSG00000165028 + ENSG00000168781 + ENSG00000161298 + ENSG00000273492,
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
# 0
# 1068


set.seed(5555)

df_0 <- df |> filter(binary_pabeta == 0) |>  sample_n(400)
df_1 <- df |> filter(binary_pabeta == 1)


df_balanced <- rbind(df_0, df_1)

model.log = glm(binary_pabeta ~  ENSG00000198887 +
                    ENSG00000165028,

                # ENSG00000133056 + ENSG00000130193 + ENSG00000198887 +
                #     ENSG00000165028 + ENSG00000168781 + ENSG00000161298 + ENSG00000273492,
                data=df_balanced,
                family = binomial(link="logit")
)

summary(model.log)

df_balanced$predicted_prob <- predict(model.log,
                                      type="response")


df_balanced<- df_balanced |> mutate(predicted_binary = case_when(
    predicted_prob <= 0.2 ~ 0,
    predicted_prob > 0.2 ~ 1 ))

table(df_balanced$predicted_binary)
# 0
# 1068





