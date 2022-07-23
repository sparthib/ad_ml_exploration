library(SpatialExperiment)
library(here)
library(ggplot2)
library(dplyr)
library(sessioninfo)
# install.packages("DescTools") #was used for converting rho to fisher's Z
# library(DescTools)

#### read in pearson pAbeta correlations ####
pearson_PAbeta <- as.data.frame(readRDS(here("corr_outputs",
                              "pearson",
                              "pearson_pAbeta.RDS") ))

pearson_PpTau<- readRDS(here("corr_outputs",
                               "pearson",
                               "pearson_PpTau.RDS") )

spearman_PAbeta <- readRDS(here("corr_outputs",
                                "spearman",
                                "PAbeta.RDS") )

spearman_PpTau<- readRDS(here("corr_outputs",
                                "spearman",
                                "PpTau.RDS") )


##### convert RDS to dataframe ####


corr_vals<- as.data.frame(cbind(pearson_PAbeta[,1:2],
                                      pearson_PpTau[,2],
                                      spearman_PAbeta[,2],
                                      spearman_PpTau[,2])
                          )


#rename pAbeta
colnames(corr_vals) <- c("gene_id", "pearson_PAbeta",
                         "pearson_PpTau",
                         "spearman_PAbeta",
                         "spearman_PpTau")


#in descending order

# pearson_pAbeta
### make correlation plots for the top 50 genes.
names(pearson_PAbeta) <- c("gene_id", "pAbeta_value")
pearson_PAbeta[,2] <- as.numeric(pearson_PAbeta[,2])
pearson_PAbeta <- pearson_PAbeta |> arrange(desc(pAbeta_value))



spe_postqc <-
    readRDS(
        here::here(
            "input_data",
            paste0("spe_wholegenome_postqc.rds")
        )
    )

spe_postqc <- logNormCounts(spe_postqc)

known_gene_ids <- pearson_PAbeta$gene_id[1:100]
ix_known <- which(rowData(spe_postqc)$gene_id %in% known_gene_ids)

spe <- spe_postqc[ix_known, ]

output_dir <- here("plots",
                   "01_correlation_exploration",
                   "02_explore_corr_values")


pdf(file = paste0(output_dir, "/top_100_genes_pAbeta_pearson.pdf"))


for (i in seq_len(100)){
    par(mfrow=c(1,1))
    plot(as.vector(logcounts(spe)[i,]), as.vector(colData(spe)$PAbeta),
         xlab = paste0(rowData(spe)$gene_name[i], " ",
                       rowData(spe)$gene_id[i]," ", "expression"),
         ylab = "PAbeta_value")


}

dev.off()



corr_vals_pAbeta <
corr_vals$pearson_PAbeta <- as.numeric(corr_vals$pearson_PAbeta)
pearson_pAbeta <-  corr_vals |> arrange(desc(corr_vals$pearson_PAbeta))

##check last x values
tail(pearson_pAbeta$pAbeta_value, 600)
#lot of NAs
#NAs can appear if there are attributes with zero variance
#(with all elements equal)

#number of values that are not NA
sum(is.na(corr_vals$pearson_PAbeta))
sum(is.na(corr_vals$pearson_PpTau))

#27349
#num NAs is 504.

sum(is.na(corr_vals$spearman_PAbeta))
sum(is.na(corr_vals$spearman_PpTau))


corr_vals <- na.omit(corr_vals)
#they all seem to have only 504 <- this means that NA is caused by zero
#variation in gene expression in these 500 genes.



# dir.create(here("plots",
#            "01_correlation_exploration",
#            "02_explore_corr_values",
#            recursive = TRUE))

output_dir <- here("plots",
                   "01_correlation_exploration",
                   "02_explore_corr_values")


##### make plots  #####


# Give the chart file a name.
pdf(file = paste0(output_dir, "/across_all_spots.pdf"))

# Plot the bar chart.
ggplot(corr_vals, aes(x = corr_and_path_type, y = correlation, fill = corr_and_path_type)) +
    geom_boxplot() +
    stat_summary(fun = "mean", geom = "point", shape = 8,
                 size = 2, color = "white")+
    theme(legend.position="none")

# Save the file.
dev.off()

##### List top 100 expressed and depleted genes

top_100 expressed and depleted genes.










