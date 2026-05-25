##---------------------------##------------------##
## Process ESBM results for binary SBM simulation (two covariates)
## - Reads all postTwoCov_*.rds in results folder
## - For each:
##     * estimate the partition using salso
##     * computes ARI vs true partition
##     * records scenario, alpha1, alpha2, seed, K_hat
##---------------------------##------------------##
library(here)
library(mclust)     # ARI
library(tools)      # file_path_sans_ext
library(salso)      # point estimate for the clustering
library(cluster)    # to compute silhouette

##---------------------------##
## Settings
##---------------------------##
results_dir <- here("simulation", "fixed_alpha", "results", "twoCov")

sim_path_NN <- here("simulation", "data", "binarySBM_2cov_NN.rds")
sim_path_IN <- here("simulation", "data", "binarySBM_2cov_IN.rds")
sim_path_NI <- here("simulation", "data", "binarySBM_2cov_NI.rds")
sim_path_II <- here("simulation", "data", "binarySBM_2cov_II.rds")

burn_frac <- 0.2
thin_step <- 1

##---------------------------------------------------------##
## List files and split by scenario keyword
##---------------------------------------------------------##
files_all <- list.files(
  results_dir,
  pattern = "^postTwoCov_.*\\.rds$",
  full.names = TRUE
)

files_NN <- files_all[grepl("2cov_NN", basename(files_all), ignore.case = TRUE)]
files_IN <- files_all[grepl("2cov_IN", basename(files_all), ignore.case = TRUE)]
files_NI <- files_all[grepl("2cov_NI", basename(files_all), ignore.case = TRUE)]
files_II <- files_all[grepl("2cov_II", basename(files_all), ignore.case = TRUE)]

cat("Found", length(files_NN), "files for scenario NN\n")
cat("Found", length(files_IN), "files for scenario IN\n")
cat("Found", length(files_NI), "files for scenario NI\n")
cat("Found", length(files_II), "files for scenario II\n\n")

res_for_scenario <- function(sim_path, scenario_files) {

  ## load simulation object
  sim_obj <- readRDS(sim_path)
  z_true  <- sim_obj$partition

  ## prepare combined structure - info from simulation
  combined <- list(
    sim = list(
      Y              = sim_obj$Y,
      x              = sim_obj$x,
      partition_true = sim_obj$partition,
      block_probs    = sim_obj$block_probs,
      clust_sizes    = sim_obj$clust_sizes,
      nCov           = sim_obj$nCov,
      cov1Dep        = sim_obj$cov1Dep,
      cov2Dep        = sim_obj$cov2Dep
    ),
    results = list()   # ONLY results grouped by (alpha1, alpha2)
  )

  ## loop over results files
  for (f in scenario_files) {

    cat("Processing:", basename(f), "\n")

    res <- readRDS(f)

    Z_post <- res$mcmcpost$z_post
    seed   <- res$seed
    N_iter <- res$N_iter

    ## alpha saved as vector c(alpha1, alpha2)
    alpha1 <- if (!is.null(res$alpha) && length(res$alpha) >= 1) res$alpha[1] else NA_real_
    alpha2 <- if (!is.null(res$alpha) && length(res$alpha) >= 2) res$alpha[2] else NA_real_

    ## burn-in & thinning
    T_total <- ncol(Z_post)
    burn    <- floor(burn_frac * T_total)
    keep_i  <- seq(from = burn + 1, to = T_total, by = thin_step)

    Z_keep <- Z_post[, keep_i, drop = FALSE]

    ## PSM
    psm_mat <- psm(t(Z_keep))

    ## VI estimate
    z_hat <- salso(t(Z_keep), loss = VI())

    ## ARI for VI representative
    ARI_vi <- adjustedRandIndex(z_hat, z_true)

    ## ARI across iterations
    ARI_iter <- apply(
      Z_keep, 2,
      function(z_t) adjustedRandIndex(z_t, z_true)
    )
    ARI_mean <- mean(ARI_iter)
    ARI_sd   <- sd(ARI_iter)

    ## silhouette (distance = 1 − PSM)
    diss_mat <- 1 - psm_mat
    sil_obj  <- silhouette(z_hat, dist(diss_mat))
    sil_mean <- mean(sil_obj[, 3])

    tag <- sprintf("alpha1_%0.2f__alpha2_%0.2f", alpha1, alpha2)

    combined$results[[tag]] <- list(
      alpha1      = alpha1,
      alpha2      = alpha2,
      Z_post      = Z_post,
      z_hat       = as.integer(z_hat),
      ARI_vi      = ARI_vi,
      ARI_mean    = ARI_mean,
      ARI_sd      = ARI_sd,
      silhouette  = list(object = sil_obj, mean = sil_mean),
      psm         = psm_mat,
      seed        = seed,
      N_iter      = N_iter
    )
  }

  combined
}

##---------------------------------------------------------##
## Build res & save for each scenario
##---------------------------------------------------------##
out_dir <- here("simulation", "fixed_alpha", "results", "processed", "twoCov")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

combined_NN <- res_for_scenario(sim_path_NN, files_NN)
combined_IN <- res_for_scenario(sim_path_IN, files_IN)
combined_NI <- res_for_scenario(sim_path_NI, files_NI)
combined_II <- res_for_scenario(sim_path_II, files_II)

saveRDS(combined_NN, file.path(out_dir, "scenario_2cov_NN_results.rds"))
saveRDS(combined_IN, file.path(out_dir, "scenario_2cov_IN_results.rds"))
saveRDS(combined_NI, file.path(out_dir, "scenario_2cov_NI_results.rds"))
saveRDS(combined_II, file.path(out_dir, "scenario_2cov_II_results.rds"))
