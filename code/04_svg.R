
library(SpatialExperiment)
library(here)
library(sgejobs)
sgejobs::job_single(
    "calculate_svg",
    create_shell = TRUE,
    queue = "bluejay",
    cores = 10L,
    memory = "10G",
    command = "Rscript 04_svg.R",
    create_logdir = TRUE
)

BiocManager::install('nnSVG')
library(nnSVG)
library(scuttle)



spe_postqc <-
    readRDS(
        here::here(
            "input_data",
            paste0("spe_wholegenome_postqc.rds")
        )
    )



#adapted multi-sample nnSVG code from LC
#https://github.com/lmweber/locus-c/blob/main/code/analyses/08a_LC_nnSVG.R#L71-L99
sample_ids <- c(
    "S1_A1_Br3874" , "S1_B1_Br3854", "S1_C1_Br3873" ,
    "S1_D1_Br3880" , "S2_A1_Br3874" , "S2_B1_Br3854",
    "S2_C1_Br3873",  "S2_D1_Br3880",
    "S3_A1_Br3874" , "S3_D1_Br3880")
sample_ids



res_list <- as.list(rep(NA, length(sample_ids)))

names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {

    s = 1
    ix <- colData(spe_postqc)$sample_id_short == sample_ids[s]
    spe_sub <- spe_postqc[, ix]

    #run nnSVG filtering for mitochondrial gene and low-expressed genes
    spe_sub <- nnSVG::filter_genes(spe_sub)

    # Gene filtering: retaining genes with at least 2 counts
    # in at least 1% (n = 40) of spatial locations
    # removed 23215 out of 27840 genes due to low expression
    # re-calculate logcounts after filtering
    set.seed(123)
    qclus <- quickCluster(spe_sub)
    # calculate size factors
    spe_sub <- computeSumFactors(spe_sub, cluster = qclus)
    spe_sub <- logNormCounts(spe_sub)

    # run nnSVG
    set.seed(123)
    spe_sub <- nnSVG(spe_sub)

    # store results
    res_list[[s]] <- rowData(spe_sub)
}


##### nnSVG function
#BiocManager::install('BRISC', force = TRUE)
library('BRISC')

library(BiocParallel)
spatial_coords = NULL
X = NULL
assay_name = "logcounts"
n_neighbors = 10
order = "AMMD"
n_threads = 4
BPPARAM = 4
verbose = FALSE



# -----------------------
# run BRISC for each gene
# -----------------------

if (is(input, "SpatialExperiment")) {
    y <- assays(spe)[[assay_name]]
    coords <- spatialCoords(spe)


# scale coordinates proportionally
range_all <- max(apply(coords, 2, function(col) diff(range(col))))
coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)

# calculate ordering of coordinates
order_brisc <- BRISC_order(coords, order = order, verbose = verbose)

# calculate nearest neighbors
nn_brisc <- BRISC_neighbor(coords, n.neighbors = n_neighbors, n_omp = 1,
                           search.type = "tree", ordering = order_brisc,
                           verbose = verbose)

# run BRISC using parallelization
ix <- seq_len(nrow(y))
out_brisc <- bplapply(1, function(i) {
    # fit model (intercept-only model if x is NULL)
    y_i <- y[i, ]
    suppressWarnings({
        runtime <- system.time({
            out_i <- BRISC_estimation(coords = coords, y = y_i, x = X,
                                      cov.model = "exponential",
                                      ordering = order_brisc, neighbor = nn_brisc,
                                      verbose = verbose)
        })
    })
    res_i <- c(
        out_i$Theta,
        loglik = out_i$log_likelihood,
        runtime = runtime[["elapsed"]]
    )
    res_i
}, BPPARAM = BPPARAM)

# collapse output list into matrix
mat_brisc <- do.call("rbind", out_brisc)

i= 1:2
y_i <- y[i, ]
suppressWarnings({
    runtime <- system.time({
        out_i <- BRISC_estimation(coords = coords, y = y_i, x = X,
                                  cov.model = "exponential",
                                  ordering = order_brisc, neighbor = nn_brisc,
                                  verbose = verbose)
    })
})


BiocManager::install('STexampleData')
library(STexampleData)
library(scran)
spe <- Visium_humanDLPFC()

dim(spe)

spe <- nnSVG::filter_genes(spe)

set.seed(123)
qclus <- quickCluster(spe)
spe <- computeSumFactors(spe, cluster = qclus)
spe <- logNormCounts(spe)

assayNames(spe)

set.seed(123)
ix_random <- sample(seq_len(nrow(spe)), 10)

known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
ix_known <- which(rowData(spe)$gene_name %in% known_genes)

ix <- c(ix_known, ix_random)

spe <- spe[ix, ]
set.seed(123)
ix_random <- sample(seq_len(nrow(spe)), 10)

known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
ix_known <- which(rowData(spe)$gene_name %in% known_genes)

ix <- c(ix_known, ix_random)

spe <- spe[ix, ]
spe <- nnSVG(spe)



