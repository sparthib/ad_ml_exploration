library(rtracklayer)
library(GenomicRanges)
library(SpatialExperiment)
library(here)
library(dplyr)
library(Matrix)
# install.packages("DescTools") #was used for converting rho to fisher's Z
# library(DescTools)

#read in pearson pAbeta correlations
pearson_pAbeta <- readRDS(here("corr_outputs",
                              "pearson",
                              "pearson_pAbeta.RDS") )


#convert RDS to dataframe
pearson_pAbeta <- as.data.frame(pearson_pAbeta)

#rename pAbeta
colnames(pearson_pAbeta) <- c("gene_id", "pAbeta_value")


#in descending order
pearson_pAbeta <- pearson_pAbeta |> arrange(desc(pAbeta_value))
pearson_pAbeta

##make plots to explore correlation values trend pearson pAbeta

dir.create(here("plots",
           "01_correlation_exploration",
           "02_explore_corr_values",
           recursive = TRUE))
output_dir <- here("plots",
                   "01_correlation_exploration",
                   "02_explore_corr_values")
# Give the chart file a name.
pdf(file = output_dir, "pAbeta_pearson_values.pdf")

# Plot the bar chart.
plot(v,type = "o")

# Save the file.
dev.off()

