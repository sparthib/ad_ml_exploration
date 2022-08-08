library(SpatialExperiment)
library(here)
library(sessioninfo)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(scran)


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

set.seed(123)

spe <- nnSVG::filter_genes(spe,filter_genes_ncounts = 3,
                           filter_mito = FALSE)
spe <- logNormCounts(spe, size.factors = NULL)

assayNames(spe)


#ad_related_genes <- c("CLU", "ABCA7", "CD33", "CR1", "TREM2", "EPHA1",
# "MEF2C", "INPP5D", "HLA-DRB5", "HLA-DRB1",
# "CD2AP", "PICALM", "SORL1", "BIN1", "ABCA7", "APOE",
# "DSG2")
ad_related_genes <- c("CLU", "MEF2C", "SORL1", "BIN1", "APOE")

ix_known <- which(rowData(spe)$gene_name %in% ad_related_genes)
spe <- spe[ix_known, ]
dim(spe)

rownames(spe) <- rowData(spe)$gene_name

input_dir <- here("corr_outputs",
                  sample_ids[3])
res_list <- readRDS( paste0(input_dir, "/filtered_non_zero_PpTau.RDS"))



#### PpTau pearson #######

PpTau <- as.matrix(colData(spe)$PpTau)
names(PpTau) <- rownames(colData(spe))

mat2 = PpTau[PpTau !=0]

res_list <- list()

## loop for removing zeros from gene as well as PpTau
for( s in seq_len(nrow(logcounts(spe)))){


    mat_nonzero <- which(logcounts(spe)[s,] != 0, arr.ind = T)

    res_list[[s]]<- merge( logcounts(spe)[s,][mat_nonzero], mat2, by = 0)


}

rownames(spe) <- rowData(spe)$gene_name
names(res_list) <- rownames(logcounts(spe))



output_dir <- here("corr_outputs",
                   sample_ids[3])
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

non_zero_PpTau_df <- merge(non_zero_PpTau_df, nrows_df, by = "gene_id")



output_dir <- here("plots","02_sample_subset", sample_ids[3])
pdf(file = paste0(output_dir, "/", "filtered_ad_pearson_PpTau.pdf"))
for( i in seq_len(nrow(non_zero_PpTau_df))){
    gene_name = non_zero_PpTau_df$gene_id[i]

    p <- ggplot(data=as.data.frame(res_list[[gene_name]]), aes(x=x, y= y)) +
        geom_point() +
        stat_density(aes(x=x, y=((..scaled..))), position="identity", geom="line") +
        xlab(gene_name) + ylab('Percentage Abeta') +
        ggtitle(paste0('Absolute Pearson Correlation = ', as.character(non_zero_PpTau_df$abs_corr[i])))

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


###### Spearman PpTau ######

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

non_zero_PpTau_df <- merge(non_zero_PpTau_df, nrows_df, by = "gene_id")



output_dir <- here("plots","02_sample_subset", sample_ids[3])
pdf(file = paste0(output_dir, "/", "filtered_ad_spearman_PpTau.pdf"))
for( i in seq_len(nrow(non_zero_PpTau_df))){
    gene_name = non_zero_PpTau_df$gene_id[i]

    p <- ggplot(data=as.data.frame(res_list[[gene_name]]), aes(x=x, y= y)) +
        geom_point() +
        stat_density(aes(x=x, y=((..scaled..))), position="identity", geom="line") +
        xlab(gene_name) + ylab('Percentage Abeta') +
        ggtitle(paste0('Absolute Spearman Correlation = ', as.character(non_zero_PpTau_df$abs_corr[i])))

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
