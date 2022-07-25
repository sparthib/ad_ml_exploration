library(SpatialExperiment)
library(spatialLIBD)
library(doParallel)
library(here)
library(dplyr)
library(scuttle)
library(sessioninfo)
library(nnSVG)
library(BiocParallel)
here()

####sample id names ####
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
print("sample id names listed")


#### get task id and subset spe ####
#s = as.numeric(Sys.getenv("SGE_TASK_ID"))
s = 1


compare_hvg_and_svg  <- function(s){
    #read in top hvgs for sample
    hvg <- readRDS(here("corr_outputs", sample_ids[s], "top_hvgs.RDS"))
    #read in top SVGs for sample
    svg <- as.data.frame(readRDS(here("corr_outputs", sample_ids[s], "top_svgs.RDS")))

    #arrange SVG by rank
    svg_rank <- svg |> select(c('gene_name', 'gene_id', 'rank')) |>
        arrange(rank) |> slice(1:length(hvg))

    #arrange SVG by padj
    svg_pval <- svg |> select(c('gene_name', 'gene_id', 'padj')) |>
        arrange(padj) |> slice(1:length(hvg))

    print(sample_ids[s])
    print(paste0("common genes b/w rank and p-adj SVGs: ",
                 length(intersect(svg_pval$gene_id, svg_rank$gene_id))))
    print(paste0("Number of common genes between highly ranked SVGs and top HVGS: ",
                 length(intersect(hvg, svg_rank$gene_id))))

    print(paste0("Number of common genes between smallest p-adj SVGs and top HVGS: ",
                       length(intersect(hvg, svg_pval$gene_id))))
    print(paste0("number of top 10% HVGs: ", length(hvg)))
    intersect(hvg, svg_rank$gene_id)


}




res_list <- list()
for(i in seq_len(10)){
    res_list[[i]] <- compare_hvg_and_svg(i)
}

head(res_list)
names(res_list) <- sample_ids

#Control
donor_1 <- intersect(res_list$S1_A1_Br3874, intersect(res_list$S2_A1_Br3874,
                     res_list$S3_A1_Br3874))
length(donor_1)
#97

#AD
donor_2 <- intersect(res_list$S1_B1_Br3854, res_list$S2_B1_Br3854)
length(donor_2)
#83

donor_3 <- intersect(res_list$S1_C1_Br3873, res_list$S2_C1_Br3873)
length(donor_3)
#283

donor_4 <- intersect(res_list$S1_D1_Br3880,
                     intersect(res_list$S2_D1_Br3880,
                               res_list$S3_D1_Br3880))
length(donor_4)
#231

#AD cases
ad_cases <- intersect(donor_2, intersect(donor_3, donor_4))
length(ad_cases)
#51

#All
all <- intersect(ad_cases, donor_1)
length(all)

#42
spe_postqc <-
    readRDS(
        here::here(
            "input_data",
            paste0("spe_wholegenome_postqc.rds")
        )
    )

spe_postqc <- logNormCounts(spe_postqc)
ix <- which(rowData(spe_postqc)$gene_id %in% all)
spe_postqc <- spe_postqc[ix, ]


dir.create(here("plots", "plotspots"))


for(i in seq_len(nrow(rowData(spe_postqc)))) {
vis_grid_gene(
    spe = spe_postqc,
    geneid = rowData(spe_postqc)$gene_id[i],
    assayname = "logcounts",
    spatial = TRUE,
    pdf_file = here("plots", "plotspots",
                    paste0(rowData(spe_postqc)$gene_id[i], "_",
                           rowData(spe_postqc)$gene_name[i],".pdf")),
    return_plots = FALSE,
    viridis = TRUE,
    height = 24,
    width = 36,
    image_id = "lowres",
    alpha = 1,
    cont_colors =  c("aquamarine4","springgreen", "goldenrod", "red"),
    sample_order = unique(spe_postqc$sample_id),

)
}

#### All samples and 10X panels  ####
ad_genes <- readxl::read_xlsx(here("10x files",
"10X_NS targeted gene_AD genes.xlsx"))
ad_genes <- ad_genes |> select(c(`Gene name`, `Ensemble ID`))

length(intersect(ad_genes$`Ensemble ID`, all))
ad_genes <- ad_genes |> filter (`Ensemble ID`
                    %in% c("ENSG00000130203", "ENSG00000154277"))

# `Gene name` `Ensemble ID`
# <chr>       <chr>
# 1 APOE        ENSG00000130203
# 2 UCHL1       ENSG00000154277

ns_genes <- readxl::read_xlsx(here("10x files",
                                   "10x_Annotated_Human_Neuroscience_Panel.xlsx"
                                   ),
                              sheet = "Gene Information")

head(ns_genes)
#ensemble_id     gene_name


length(intersect(ns_genes$ensemble_id, all)) #14

ns_genes <- ns_genes |> filter (ensemble_id %in%
                                    intersect(ns_genes$ensemble_id, all))
# ensemble_id     gene_name
# <chr>           <chr>
# 1 ENSG00000130203 APOE
# 2 ENSG00000171885 AQP4    SCZ
# 3 ENSG00000109846 CRYAB   multiple sclerosis
# 4 ENSG00000131095 GFAP
# 5 ENSG00000115756 HPCAL1
# 6 ENSG00000197971 MBP     SCZ
# 7 ENSG00000123560 PLP1    multiple sclerosis
# 8 ENSG00000107317 PTGDS   SCZ
# 9 ENSG00000160307 S100B   SCZ, psychiatric disorders
# 10 ENSG00000110436 SLC1A2 scz
# 11 ENSG00000132639 SNAP25 scz, psych, weight gain
# 12 ENSG00000118785 SPP1   multiple sclerosis
# 13 ENSG00000067715 SYT1
# 14 ENSG00000154277 UCHL1  AD, Parkinson's


#### AD cases and 10X panels  ####
length(intersect(ad_genes$`Ensemble ID`, ad_cases))
ad_genes <- ad_genes |> filter (`Ensemble ID`
                                %in% c("ENSG00000130203", "ENSG00000154277"))

# `Gene name` `Ensemble ID`
# <chr>       <chr>
# 1 APOE        ENSG00000130203
# 2 UCHL1       ENSG00000154277


head(ns_genes)
#ensemble_id     gene_name


length(intersect(ns_genes$ensemble_id, ad_cases)) #14

ns_genes <- ns_genes |> filter (ensemble_id %in%
                                    intersect(ns_genes$ensemble_id, all))
# ensemble_id     gene_name
# <chr>           <chr>
# 1 ENSG00000130203 APOE
# 2 ENSG00000171885 AQP4    SCZ
# 3 ENSG00000109846 CRYAB   multiple sclerosis
# 4 ENSG00000131095 GFAP
# 5 ENSG00000115756 HPCAL1
# 6 ENSG00000197971 MBP     SCZ
# 7 ENSG00000123560 PLP1    multiple sclerosis
# 8 ENSG00000107317 PTGDS   SCZ
# 9 ENSG00000160307 S100B   SCZ, psychiatric disorders
# 10 ENSG00000110436 SLC1A2 scz
# 11 ENSG00000132639 SNAP25 scz, psych, weight gain
# 12 ENSG00000118785 SPP1   multiple sclerosis
# 13 ENSG00000067715 SYT1
# 14 ENSG00000154277 UCHL1  AD, Parkinson's
# [1] "S1_A1_Br3874"
# [1] "common genes b/w rank and p-val SVGs: 1465"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 416"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 416"
# [1] "number of top 10% HVGs: 1465"
# [1] "S1_B1_Br3854"
# [1] "common genes b/w rank and p-val SVGs: 1092"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 159"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 159"
# [1] "number of top 10% HVGs: 1092"
# [1] "S1_C1_Br3873"
# [1] "common genes b/w rank and p-val SVGs: 1410"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 516"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 516"
# [1] "number of top 10% HVGs: 1410"
# [1] "S1_D1_Br3880"
# [1] "common genes b/w rank and p-val SVGs: 1143"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 473"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 473"
# [1] "number of top 10% HVGs: 1143"
# [1] "S2_A1_Br3874"
# [1] "common genes b/w rank and p-val SVGs: 1299"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 291"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 291"
# [1] "number of top 10% HVGs: 1299"
# [1] "S2_B1_Br3854"
# [1] "common genes b/w rank and p-val SVGs: 1125"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 313"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 313"
# [1] "number of top 10% HVGs: 1150"
# [1] "S2_C1_Br3873"
# [1] "common genes b/w rank and p-val SVGs: 1245"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 370"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 370"
# [1] "number of top 10% HVGs: 1245"
# [1] "S2_D1_Br3880"
# [1] "common genes b/w rank and p-val SVGs: 1289"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 485"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 485"
# [1] "number of top 10% HVGs: 1289"
# [1] "S3_A1_Br3874"
# [1] "common genes b/w rank and p-val SVGs: 1263"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 456"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 467"
# [1] "number of top 10% HVGs: 1688"
# [1] "S3_D1_Br3880"
# [1] "common genes b/w rank and p-val SVGs: 1567"
# [1] "Number of common genes between highly ranked SVGs and top HVGS: 622"
# [1] "Number of common genes between smallest p-adj SVGs and top HVGS: 622"
# [1] "number of top 10% HVGs: 1567"


