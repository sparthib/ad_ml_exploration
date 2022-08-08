# sgejobs::job_single('non_zero_filtered_genes_PpTau',
#                     create_shell = TRUE,
#                     create_logdir = TRUE, task_num = 10,
#                     queue = "bluejay",
#                     cores = 5L,
#                     command = "Rscript 04b_non_zero_filtered_genes_PpTau.R",
#                     memory = "20G")
#

library(SpatialExperiment)
library(here)
library(sessioninfo)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)


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

s = as.numeric(Sys.getenv("SGE_TASK_ID"))
print(s)



spe <-readRDS(here::here("input_data",
                         paste0("spe_wholegenome_postqc.rds")))


ix <- colData(spe)$sample_id_short == sample_ids[s]
spe<- spe[, ix]

rownames(spe) <- rowData(spe)$gene_name
spe <- nnSVG::filter_genes(spe,filter_genes_ncounts = 3,
                           filter_mito = FALSE)
spe<- scuttle::logNormCounts(spe)


dim(spe)
# [1] 1990 2665


spe<- scuttle::logNormCounts(spe)




####PpTau spearman#######

PpTau <- as.matrix(colData(spe)$PpTau)
names(PpTau) <- rownames(colData(spe))

mat2 = PpTau[PpTau !=0]

res_list <- list()

## loop for removing zeros from gene as well as PpTau
for( s in seq_len(nrow(logcounts(spe)))){


    mat_nonzero <- which(logcounts(spe)[s,] != 0, arr.ind = T)

    res_list[[s]]<- merge( logcounts(spe)[s,][mat_nonzero], mat2, by = 0)


}


names(res_list) <- rownames(logcounts(spe))



output_dir <- here("corr_outputs",
                   sample_ids[s])
saveRDS(res_list, paste0(output_dir, "/filtered_non_zero_PpTau.RDS"))


non_zero_genes <- c()
non_zero_corrs <- c()


nrows <- c()

## loop stores the number of observations for correlation between PpTau
## and gene for each gene
for(i in seq_len(length(res_list))) {
    nrows[i] <- nrow(res_list[[i]])
}

nrows_df <- as.data.frame(cbind(names(res_list), nrows))
colnames(nrows_df) <- c("gene_id", "nrows")
nrows_df$nrows <- as.numeric(nrows_df$nrows)



######calculate Pearson ########

print("calculate pearson")

##calculate correlation for each filtered gene
for( i in seq_len(length(res_list))){
    non_zero_genes <- c(non_zero_genes,names(res_list)[i] )
    non_zero_corrs <- c( non_zero_corrs,
                         cor(res_list[[i]][,'x'],res_list[[i]][,'y'], method = 'pearson'))
}



######Corr df
non_zero_PpTau_df <- as.data.frame(cbind(non_zero_genes, non_zero_corrs))
colnames(non_zero_PpTau_df) <- c("gene_id", "corr_value")

non_zero_PpTau_df <- na.omit(non_zero_PpTau_df)
non_zero_PpTau_df$corr_value <- as.numeric(non_zero_PpTau_df$corr_value)

non_zero_PpTau_df$abs_corr <- abs(non_zero_PpTau_df$corr_value)
non_zero_PpTau_df <- non_zero_PpTau_df |> arrange(corr_value, desc = TRUE)

non_zero_PpTau_df <- non_zero_PpTau_df |> mutate(bin = case_when(
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


non_zero_PpTau_df <- merge(non_zero_PpTau_df, nrows_df, by = "gene_id")

non_zero_groups <- non_zero_PpTau_df  |>
    group_by(bin)  |> top_n(n = 50, nrows) |> arrange(by = abs_corr, desc = TRUE)
# arrange(dep_delay, .by_group = TRUE)

tail(non_zero_groups, n = 15)$gene_id




output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "filtered_non_zero_scatter_pearson_PpTau.pdf"))
for( i in seq_len(nrow(non_zero_groups))){
    gene_name = non_zero_groups$gene_id[i]

    p <- ggplot(data=as.data.frame(res_list[[gene_name]]), aes(x=x, y= y)) +
        geom_point() +
        stat_density(aes(x=x, y=((..scaled..))), position="identity", geom="line") +
        xlab(gene_name) + ylab('Percentage Abeta') +
        ggtitle(paste0('Absolute Pearson Correlation = ', as.character(non_zero_groups$abs_corr[i])))

    print(p)
    #print(ggMarginal(p,type = "density", margins = "y", size = 4))


    #stuff <- ggplot_build(p)
    #xrange <- stuff[[2]]$ranges[[1]]$x.range  # extract the x range, to make the new densities align with y-axis

    ## Get densities of dim2
    #dens <- with(res_list[["RPS24"]], density(y))
    #ds = data.frame(x=dens$y+xrange[1], y=dens$x)
    #ds = data.frame(x= -dens$y , y=  dens$x)

    #p + geom_path(data=ds, aes(x=x, y=y))

}

dev.off()

######## Calculate Spearman ##########
print("calculate spearman")


##calculate correlation for each filtered gene
for( i in seq_len(length(res_list))){
    non_zero_genes <- c(non_zero_genes,names(res_list)[i] )
    non_zero_corrs <- c( non_zero_corrs,
                         cor(res_list[[i]][,'x'],res_list[[i]][,'y'], method = 'spearman'))
}



######Corr df
non_zero_PpTau_df <- as.data.frame(cbind(non_zero_genes, non_zero_corrs))
colnames(non_zero_PpTau_df) <- c("gene_id", "corr_value")

non_zero_PpTau_df <- na.omit(non_zero_PpTau_df)
non_zero_PpTau_df$corr_value <- as.numeric(non_zero_PpTau_df$corr_value)

non_zero_PpTau_df$abs_corr <- abs(non_zero_PpTau_df$corr_value)
non_zero_PpTau_df <- non_zero_PpTau_df |> arrange(corr_value, desc = TRUE)

non_zero_PpTau_df <- non_zero_PpTau_df |> mutate(bin = case_when(
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


non_zero_PpTau_df <- merge(non_zero_PpTau_df, nrows_df, by = "gene_id")

non_zero_groups <- non_zero_PpTau_df  |>
    group_by(bin)  |> top_n(n = 50, nrows) |> arrange(by = abs_corr, desc = TRUE)
# arrange(dep_delay, .by_group = TRUE)

###### make spot plots ####
print("make spotplots")

output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "filtered_non_zero_scatter_spearman_PpTau.pdf"))
for( i in seq_len(nrow(non_zero_groups))){
    gene_name = non_zero_groups$gene_id[i]

    p <- ggplot(data=as.data.frame(res_list[[gene_name]]), aes(x=x, y= y)) +
        geom_point() +
        stat_density(aes(x=x, y=((..scaled..))), position="identity", geom="line") +
        xlab(gene_name) + ylab('Percentage Abeta') +
        ggtitle(paste0('Absolute Spearman Correlation = ', as.character(non_zero_groups$abs_corr[i])))

    print(p)
    #print(ggMarginal(p,type = "density", margins = "y", size = 4))


    #stuff <- ggplot_build(p)
    #xrange <- stuff[[2]]$ranges[[1]]$x.range  # extract the x range, to make the new densities align with y-axis

    ## Get densities of dim2
    #dens <- with(res_list[["RPS24"]], density(y))
    #ds = data.frame(x=dens$y+xrange[1], y=dens$x)
    #ds = data.frame(x= -dens$y , y=  dens$x)

    #p + geom_path(data=ds, aes(x=x, y=y))

}

dev.off()

# gene_id  corr_value abs_corr bin   nrows
# <chr>         <dbl>    <dbl> <chr> <dbl>
#     1 SLC9A6        0.283    0.283 3       107
# 2 NEXN         -0.3      0.3   4         5
# 3 GDA           0.309    0.309 4        84
# 4 EVI2A        -0.309    0.309 4        38
# 5 THBS1         0.310    0.310 4        10
# 6 CNN1          0.328    0.328 4         5
# 7 MGP           0.330    0.330 4        20
# 8 IK           -0.339    0.339 4        83
# 9 CNTN2         0.342    0.342 4        75
# 10 PRELP         0.347    0.347 4        26
# 11 PRELP         0.363    0.363 4        26
# 12 SERPINE1     -0.365    0.365 4         9
# 13 NEXN         -0.400    0.400 4         5
# 14 FLNC          0.4      0.4   5         5
# 15 MGP           0.407    0.407 5        20
# 16 CARNS1        0.458    0.458 5        34
# 17 MYH11         0.5      0.5   6         3
# 18 MYH11         0.795    0.795 8         3
# 19 DES          -1        1     10        2
# 20 DES          -1        1     10        2
