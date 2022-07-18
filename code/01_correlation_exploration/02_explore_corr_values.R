library(rtracklayer)
library(GenomicRanges)
library(SpatialExperiment)
library(here)
library(ggplot2)
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
pearson_pAbeta$pAbeta_value <- as.numeric(pearson_pAbeta$pAbeta_value)
pearson_pAbeta <- pearson_pAbeta |> arrange(desc(pAbeta_value))
pearson_pAbeta

##make plots to explore correlation values trend pearson pAbeta

# dir.create(here("plots",
#            "01_correlation_exploration",
#            "02_explore_corr_values",
#            recursive = TRUE))
output_dir <- here("plots",
                   "01_correlation_exploration",
                   "02_explore_corr_values")

##check last x values
tail(pearson_pAbeta$pAbeta_value, 600)
#lot of NAs
#NAs can appear if there are attributes with zero variance
#(with all elements equal)

#number of values that are not NA
sum(!is.na(pearson_pAbeta$pAbeta_value))
#27349
#num NAs is 504.


# Give the chart file a name.
pdf(file = paste0(output_dir, "/pAbeta_pearson_values.pdf"))
# Plot the bar chart.
ggplot(data = pearson_pAbeta,
       aes(x = "", y =pAbeta_value)) +
    geom_boxplot(fill="slateblue", alpha=0.2) +
    xlab("") +
    ylab("Range of pAbeta correlation for each gene")
# Save the file.
dev.off()

###


