##---------------------------##------------------##
## Process ESBM results for binary SBM simulation
##- Reads all post_*.rds in results folder
##- For each:
##     * estimate the partition using salso 
##     * computes ARI vs true partition
##     * records scenario, alpha, seed, K_hat
##---------------------------##------------------##
library(here)
library(mclust)     # ARI
library(tools)      # file_path_sans_ext
library(salso)      # point estimate for the clustering
library(cluster)    # to compute silhouette
##---------------------------##
## Settings
##---------------------------##

## Paths to results directory and simulation files
results_dir <- here("simulation", "results", "oneCov")
sim_path_none <- here("simulation", "data", "binarySBM_none.rds")
sim_path_one <- here("simulation", "data", "binarySBM_one.rds")

## MCMC
burn_frac <- 0.2 # brunin
thin_step <- 1   # eventual thinning


files_res   <- list.files(results_dir, pattern = "^postOneCov_.*\\.rds$", full.names = TRUE)

##---------------------------------------------------------##
## List files and split by scenario keyword
##---------------------------------------------------------##
files_all <- list.files(
  results_dir,
  pattern = "^postOneCov_.*\\.rds$",
  full.names = TRUE
)

files_none <- files_all[grepl("_none", basename(files_all), ignore.case = TRUE)]
files_one  <- files_all[grepl("_one",  basename(files_all), ignore.case = TRUE)]

cat("Found", length(files_none), "files for scenario NONE\n")
cat("Found", length(files_one),  "files for scenario ONE\n\n")

res_for_scenario <- function(sim_path, scenario_files) {

  ## load simulation object
  sim_obj <- readRDS(sim_path)
  z_true  <- sim_obj$partition

  ## prepare combined structure — NO meta
  combined <- list(
    sim = list(
      Y              = sim_obj$Y,
      x              = sim_obj$x,
      partition_true = sim_obj$partition,
      block_probs    = sim_obj$block_probs,
      clust_sizes    = sim_obj$clust_sizes,
      cov_dep        = sim_obj$cov_dep
    ),
    results = list()   # ONLY results grouped by alpha
  )

  ## LOOP OVER POSTERIOR FILES
  for (f in scenario_files) {

    cat("Processing:", basename(f), "\n")

    res <- readRDS(f)

    Z_post <- res$Z_post
    alpha  <- res$alpha
    seed   <- res$seed
    N_iter <- res$N_iter

    ## burn-in & thinning
    T_total <- ncol(Z_post)
    burn    <- floor(burn_frac * T_total)
    keep_i  <- seq(from = burn + 1L, to = T_total, by = thin_step)

    Z_keep <- Z_post[, keep_i, drop = FALSE]

    ## PSM
    psm_mat <- psm(t(Z_keep))

    ## VI representative clustering
    z_hat <- salso(t(Z_keep), loss = VI())

    ## ARI
    ARI_val <- mclust::adjustedRandIndex(z_hat, z_true)

    ## silhouette (distance = 1 − PSM)
    diss_mat <- 1 - psm_mat
    sil_obj  <- silhouette(z_hat, dist(diss_mat))
    sil_mean <- mean(sil_obj[, 3])

    ## tag (clean list name)
    tag <- sprintf("alpha_%0.2f", alpha)

    ## store entry
    combined$results[[tag]] <- list(
      alpha      = alpha,
      Z_post     = Z_post,
      z_hat      = z_hat,
      ARI        = ARI_val,
      silhouette = list(object = sil_obj, mean = sil_mean),
      psm        = psm_mat,
      seed       = seed,
      N_iter     = N_iter
    )
  }

  combined
}

##---------------------------------------------------------##
## Build & save for scenario NONE
##---------------------------------------------------------##
combined_none <- res_for_scenario(sim_path_none, files_none)

saveRDS(
  combined_none,
  here("simulation", "results", "scenario_none_results.rds")
)


##---------------------------------------------------------##
## Build & save for scenario ONE
##---------------------------------------------------------##
combined_one <- res_for_scenario(sim_path_one, files_one)
 
 
saveRDS(
  combined_one,
  here("simulation", "results", "scenario_one_results.rds")
)


