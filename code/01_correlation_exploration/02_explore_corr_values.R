library(SpatialExperiment)
library(here)
install.packages('ggvis')
library(ggvis)
library(dplyr)
library(sessioninfo)
# install.packages("DescTools") #was used for converting rho to fisher's Z
# library(DescTools)

#read in pearson pAbeta correlations
pearson_PAbeta <- readRDS(here("corr_outputs",
                              "pearson",
                              "pearson_pAbeta.RDS") )

pearson_PpTau<- readRDS(here("corr_outputs",
                               "pearson",
                               "pearson_PpTau.RDS") )

spearman_PAbeta <- readRDS(here("corr_outputs",
                                "spearman",
                                "PAbeta.RDS") )

spearman_PpTau<- readRDS(here("corr_outputs",
                                "spearman",
                                "PpTau.RDS") )


#convert RDS to dataframe


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

corr_vals <- corr_vals |> pivot_longer(!gene_id, names_to = "corr_and_path_type",
                          values_to = "correlation")
corr_vals$correlation <- as.numeric(corr_vals$correlation)

#
# #in descending order
# pearson_pAbeta$pAbeta_value <- as.numeric(pearson_pAbeta$pAbeta_value)
# pearson_pAbeta <- pearson_pAbeta |> arrange(desc(pAbeta_value))
# pearson_pAbeta

##make plots to explore correlation values trend pearson pAbeta


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


#make plots and patch them


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

###


