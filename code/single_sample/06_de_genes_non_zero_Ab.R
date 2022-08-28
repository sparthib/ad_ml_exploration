library(here)
library(dplyr)
library(sessioninfo)
library(ggplot2)


library(SpatialExperiment)
library(spatialLIBD)
library(scuttle)
library(nnSVG)
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

s = 4


#load wholegenome modeling results
load(here('/Users/sparthib/Documents/Visium_IF_AD/processed-data/11_grey_matter_only/wholegenome/Visium_IF_AD_modeling_results.Rdata'))


output_dir <- here("corr_outputs",
                   sample_ids[s])
filtered_non_zero_PAbeta <- readRDS( paste0(output_dir, "/filtered_non_zero_PAbeta.RDS"))

length(filtered_non_zero_PAbeta)
#1566

colnames(modeling_results$enrichment)
# [1] "t_stat_none"       "t_stat_Ab+"        "t_stat_next_Ab+"
# [4] "t_stat_pT+"        "t_stat_next_pT+"   "t_stat_both"
# [7] "t_stat_next_both"  "p_value_none"      "p_value_Ab+"
# [10] "p_value_next_Ab+"  "p_value_pT+"       "p_value_next_pT+"
# [13] "p_value_both"      "p_value_next_both" "fdr_none"
# [16] "fdr_Ab+"           "fdr_next_Ab+"      "fdr_pT+"
# [19] "fdr_next_pT+"      "fdr_both"          "fdr_next_both"
# [22] "ensembl"           "gene"

nextAb_enriched_genes <- modeling_results$enrichment |> filter(`fdr_next_Ab+` < 0.1) |> select(gene)
nrow(nextAb_enriched_genes)
# [1] 2380

Ab_enriched_genes <- modeling_results$enrichment |> filter(`fdr_Ab+` < 0.1) |> select(gene)
# [1] 696
PAbeta_genes <- names(filtered_non_zero_PAbeta)


nextAb_enriched_genes <- intersect(nextAb_enriched_genes$gene, PAbeta_genes)
# [1] 431

Ab_enriched_genes <-intersect(Ab_enriched_genes$gene, PAbeta_genes)
# [1] 47


#### load spe ####
s = 2




spe <-readRDS(here::here("input_data",
                         paste0("spe_wholegenome_postqc.rds")))


ix <- colData(spe)$sample_id_short == sample_ids[s]
spe<- spe[, ix]

set.seed(123)

spe <- nnSVG::filter_genes(spe,filter_genes_ncounts = 3,
                           filter_mito = FALSE)
spe <- logNormCounts(spe, size.factors = NULL)

assayNames(spe)


ix_known <- which(rowData(spe)$gene_name %in% Ab_enriched_genes)
spe <- spe[ix_known, ]
dim(spe)

rownames(spe) <- rowData(spe)$gene_name

input_dir <- here("corr_outputs",
                  sample_ids[2])
res_list <- readRDS( paste0(input_dir, "/filtered_non_zero_PAbeta.RDS"))


#### PAbeta pearson #######
non_zero_genes <- c()
non_zero_corrs <- c()


nrows <- c()


## loop stores the number of observations for correlation between PAbeta
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
non_zero_pabeta_df <- as.data.frame(cbind(non_zero_genes, non_zero_corrs))
colnames(non_zero_pabeta_df) <- c("gene_id", "corr_value")

non_zero_pabeta_df <- na.omit(non_zero_pabeta_df)
non_zero_pabeta_df$corr_value <- as.numeric(non_zero_pabeta_df$corr_value)

non_zero_pabeta_df$abs_corr <- abs(non_zero_pabeta_df$corr_value)
non_zero_pabeta_df <- non_zero_pabeta_df |> arrange(corr_value, desc = TRUE)

non_zero_pabeta_df <- merge(non_zero_pabeta_df, nrows_df, by = "gene_id")

# 1       NPY -0.348291659 0.348291659
# 2    CACNG8 -0.161535698 0.161535698
# 3     MXRA7 -0.157835091 0.157835091
# 4     SHTN1 -0.153103967 0.153103967
# 5     VDAC3 -0.148084497 0.148084497
# 6   NEUROD2 -0.132480817 0.132480817
# 7      PKIA -0.099258264 0.099258264
# 8      GLRB -0.089132237 0.089132237
top_Ab_corr_genes <- c("NPY", "CACNG8", "MXRA7", "SHTN1",
                       "VDAC3", "NEUROD2", "PKIA", "GLRB")
non_zero_Ab_enriched_df <- non_zero_pabeta_df |> filter(gene_id %in% Ab_enriched_genes)
non_zero_Ab_enriched_df

non_zero_Ab_enriched_df <- non_zero_Ab_enriched_df |> arrange(corr_value, desc = TRUE)

non_zero_nextAb_enriched_df <- non_zero_pabeta_df |> filter(gene_id %in% nextAb_enriched_genes)
non_zero_nextAb_enriched_df <- non_zero_nextAb_enriched_df |> arrange(corr_value, desc = TRUE)


#Abenriched_res_list
ix_known <- which(names(res_list) %in% Ab_enriched_genes)
Abenriched_res_list  <- res_list[ix_known]



non_zero_genes <- c()
non_zero_corrs <- c()


nrows <- c()

## loop stores the number of observations for correlation between PAbeta
## and gene for each gene
for(i in seq_len(length(Abenriched_res_list))) {
    nrows[i] <- nrow(Abenriched_res_list[[i]])
}

nrows_df <- as.data.frame(cbind(names(Abenriched_res_list), nrows))
colnames(nrows_df) <- c("gene_id", "nrows")
nrows_df$nrows <- as.numeric(nrows_df$nrows)



##calculate correlation for each filtered gene
for( i in seq_len(length(Abenriched_res_list))){
    non_zero_genes <- c(non_zero_genes,names(Abenriched_res_list)[i] )
    non_zero_corrs <- c( non_zero_corrs,
                         cor(Abenriched_res_list[[i]][,'x'],Abenriched_res_list[[i]][,'y'], method = 'pearson'))
}





output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "Abenriched_pearson_PAbeta.pdf"))
for( i in seq_len(nrow(non_zero_Ab_enriched_df))){
    gene_name = non_zero_Ab_enriched_df$gene_id[i]

    p <- ggplot(data=as.data.frame(res_list[[gene_name]]), aes(x=x, y= y)) +
        geom_point() +
        stat_density(aes(x=x, y=((..scaled..))), position="identity", geom="line") +
        xlab(gene_name) + ylab('Percentage Abeta') +
        ggtitle(paste0('Absolute Pearson Correlation = ', as.character(non_zero_Ab_enriched_df$abs_corr[i])))

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




spe <-readRDS(here::here("input_data",
                         paste0("spe_wholegenome_postqc.rds")))
spe<- logNormCounts(spe)
colData(spe)$sample_name <- colData(spe)$sample_id
rownames(spe) <- rowData(spe)$gene_name


output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "spot_plots_Abenrichedgenes_Pabeta.pdf"),
    width = 10)
for(i in seq_len(length(Ab_enriched_genes))){
    p <- spatialLIBD::vis_gene(
        spe,
        sampleid =sample_ids[s],
        geneid = Ab_enriched_genes[i],
        spatial = TRUE,
        assayname = "logcounts",
        minCount = 0,
        viridis = TRUE,
        image_id = "lowres",
        alpha = 1,
        point_size = 1.25
    ) +ggtitle(Ab_enriched_genes[i])


    # q <- spatialLIBD::vis_gene(
    #     spe,
    #     sampleid =sample_ids[s],
    #     geneid = "PAbeta",
    #     spatial = TRUE,
    #     assayname = "logcounts",
    #     minCount = 0,
    #     viridis = TRUE,
    #     image_id = "lowres",
    #     alpha = 1,
    #     point_size = 1.25
    # )+ggtitle("PAbeta")
    # print(p | q)
}


dev.off()



##### nextAbenriched_res_list ####
ix_known <- which(names(res_list) %in% nextAb_enriched_genes)
nextAbenriched_res_list  <- res_list[ix_known]




non_zero_genes <- c()
non_zero_corrs <- c()


nrows <- c()

## loop stores the number of observations for correlation between PAbeta
## and gene for each gene
for(i in seq_len(length(nextAbenriched_res_list))) {
    nrows[i] <- nrow(nextAbenriched_res_list[[i]])
}

nrows_df <- as.data.frame(cbind(names(nextAbenriched_res_list), nrows))
colnames(nrows_df) <- c("gene_id", "nrows")
nrows_df$nrows <- as.numeric(nrows_df$nrows)



##calculate correlation for each filtered gene
for( i in seq_len(length(nextAbenriched_res_list))){
    non_zero_genes <- c(non_zero_genes,names(nextAbenriched_res_list)[i] )
    non_zero_corrs <- c( non_zero_corrs,
                         cor(nextAbenriched_res_list[[i]][,'x'],
                             nextAbenriched_res_list[[i]][,'y'], method = 'pearson'))
}





output_dir <- here("plots","02_sample_subset", sample_ids[s])
pdf(file = paste0(output_dir, "/", "next_Abenriched_pearson_PAbeta.pdf"))
pdf(file = paste0(output_dir, "/", "top_Ab_corr_genes.pdf"))
for( i in seq_len(8)){
    gene_name = non_zero_Ab_enriched_df$gene_id[i]

    p <- ggplot(data=as.data.frame(res_list[[gene_name]]), aes(x=x, y= y)) +
        geom_point() +
        stat_density(aes(x=x, y=((..scaled..))), position="identity", geom="line") +
        xlab(gene_name) + ylab('Percentage Abeta') +
        ggtitle(paste0('Absolute Pearson Correlation = ',
                       as.character(non_zero_Ab_enriched_df$abs_corr[i])))

    print(p)
    #print(ggMarginal(p,type = "density", margins = "y", size = 4))


    #stuff <- ggplot_build(p)
    #xrange <- stuff[[2]]$ranges[[1]]$x.range  # extract the x range, to make the
    #new densities align with y-axis

    ## Get densities of dim2
    #dens <- with(res_list[["RPS24"]], density(y))
    #ds = data.frame(x=dens$y+xrange[1], y=dens$x)
    #ds = data.frame(x= -dens$y , y=  dens$x)

    #p + geom_path(data=ds, aes(x=x, y=y))

}

dev.off()


    #print(ggMarginal(p,type = "density", margins = "y", size = 4))


    #stuff <- ggplot_build(p)
    #xrange <- stuff[[2]]$ranges[[1]]$x.range  # extract the x range, to make the
    #new densities align with y-axis

    ## Get densities of dim2
    #dens <- with(res_list[["RPS24"]], density(y))
    #ds = data.frame(x=dens$y+xrange[1], y=dens$x)
    #ds = data.frame(x= -dens$y , y=  dens$x)

    #p + geom_path(data=ds, aes(x=x, y=y))


