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
results_dir <- here("simulation", "learn_alpha", "results", "oneCov")
sim_path_neutral   <- here("simulation", "data", "binarySBM_1cov_neutral.rds")
sim_path_info      <- here("simulation", "data", "binarySBM_1cov_informative.rds")
sim_path_mis_rand  <- here("simulation", "data", "binarySBM_1cov_mislead_random.rds")
sim_path_mis_shift <- here("simulation", "data", "binarySBM_1cov_mislead_shifted.rds")

## MCMC
burn_frac <- 0.4 # brunin
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

files_neutral   <- files_all[grepl("1cov_neutral",   basename(files_all), ignore.case = TRUE)]
files_info      <- files_all[grepl("1cov_informative", basename(files_all), ignore.case = TRUE)]
files_mis_rand  <- files_all[grepl("1cov_mislead_random", basename(files_all), ignore.case = TRUE)]
files_mis_shift <- files_all[grepl("1cov_mislead_shifted", basename(files_all), ignore.case = TRUE)]

cat("Found", length(files_neutral), "files for scenario NEUTRAL\n")
cat("Found", length(files_info), "files for scenario INFO \n")
cat("Found", length(files_mis_rand), "files for scenario mis_rand\n")
cat("Found", length(files_mis_shift), "files for scenario mis_shift\n")

res_for_scenario <- function(sim_path, scenario_files) {

  ## load simulation object
  sim_obj <- readRDS(sim_path)
  z_true  <- sim_obj$partition

  ## prepare combined structure - info form simulation
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
    results = list()   # ONLY results grouped by alpha
  )

  ## loop over results files
  for (f in scenario_files) {

    tag <- file_path_sans_ext(basename(f))

    cat("Processing:", basename(f), "\n")

    res <- readRDS(f)
    alpha_g  <- res$alpha_g

    Z_post <- res$mcmcpost$z_post
    seed   <- res$seed
    N_iter <- res$N_iter

    ## burn-in & thinning
    T_total <- ncol(Z_post)
    burn    <- floor(burn_frac * T_total)
    keep_i  <- seq(from = burn + 1L, to = T_total, by = thin_step)

    Z_keep <- Z_post[, keep_i, drop = FALSE]

    alpha_post <- res$mcmcpost$alpha_g_post

    if (!is.null(alpha_post)) {
      alpha_keep <- alpha_post[, keep_i, drop = FALSE]

    }
    ## PSM
    psm_mat <- psm(t(Z_keep))

    ## estimate of the clusering using VI 
    z_hat <- salso(t(Z_keep), loss = VI())
    
    ## ARI for VI representative
    ARI_vi <- adjustedRandIndex(z_hat, z_true)

    ## ARI across iterations
    ARI_iter <- apply(
      Z_keep, 2L,
      function(z_t) adjustedRandIndex(z_t, z_true)
    )
    ARI_mean <- mean(ARI_iter)
    ARI_sd   <- sd(ARI_iter)

    ## silhouette (distance = 1 − PSM)
    diss_mat <- 1 - psm_mat
    sil_obj  <- silhouette(z_hat, dist(diss_mat))
    sil_mean <- mean(sil_obj[, 3])

    combined$results[[tag]] <- list(
      alpha_g     = alpha_g,
      alpha_post  = alpha_keep, 
      Z_post      = Z_keep,
      z_hat       = as.integer(z_hat),
      ARI_vi      = ARI_vi,        # ARI of VI estimate
      ARI_mean    = ARI_mean,      # posterior mean ARI
      ARI_sd      = ARI_sd,        # posterior sd ARI
      silhouette  = list(object = sil_obj, mean = sil_mean),
      psm         = psm_mat,
      seed        = seed,
      N_iter      = N_iter
    )
  }


  combined
}

##---------------------------------------------------------##
## Build res & save for seach scenarion
##---------------------------------------------------------##
out_dir <- here("simulation", "learn_alpha", "results", "processed", "oneCov")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##---------------------------------------------------------##
## Build combined objects (one-cov scenarios)
##---------------------------------------------------------##
combined_neutral   <- res_for_scenario(sim_path_neutral,   files_neutral)
combined_info      <- res_for_scenario(sim_path_info,      files_info)
combined_mis_rand  <- res_for_scenario(sim_path_mis_rand,  files_mis_rand)
combined_mis_shift <- res_for_scenario(sim_path_mis_shift, files_mis_shift)

saveRDS(combined_neutral,   file.path(out_dir, "scenario_1cov_neutral_results.rds"))
saveRDS(combined_info,      file.path(out_dir, "scenario_1cov_informative_results.rds"))
saveRDS(combined_mis_rand,  file.path(out_dir, "scenario_1cov_mislead_random_results.rds"))
saveRDS(combined_mis_shift, file.path(out_dir, "scenario_1cov_mislead_shifted_results.rds"))