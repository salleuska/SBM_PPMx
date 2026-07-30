##---------------------------##------------------##
## Process ESBM results for binary SBM simulation (one covariate)
##- Reads all postOneCov_*.rds in results folder
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

args <- R.utils::commandArgs(asValue = TRUE)

##---------------------------##
## Settings
##---------------------------##

## Only results run with this calibration are processed. Runs produced before
## the calibration tag existed used the sampler default ("normalized") and are
## therefore ignored here.
similarity_calibration <- if (is.null(args$similarity_calibration)) {
  "raw"
} else {
  as.character(args$similarity_calibration)
}

## Paths to results directory and simulation files
results_dir <- here("simulation", "fixed_alpha", "results", "oneCov")
data_dir    <- here("simulation", "data")

## 1-cov scenarios: N = neutral, I = informative, M = misleading
scenarios <- c("N", "I", "M")

## MCMC
burn_frac <- 0.4 # burnin
thin_step <- 1   # eventual thinning

cat("Results dir :", results_dir, "\n")
cat("Calibration :", similarity_calibration, "\n\n")

##---------------------------------------------------------##
## List files and split by scenario code
##---------------------------------------------------------##
files_all <- list.files(
  results_dir,
  pattern = "^postOneCov_.*\\.rds$",
  full.names = TRUE
)

files_for_scenario <- function(code) {
  pat <- sprintf(
    "^postOneCov_binarySBM_1cov_%s_cal-%s_alpha-",
    code, similarity_calibration
  )
  files_all[grepl(pat, basename(files_all))]
}

scenario_files <- lapply(scenarios, files_for_scenario)
names(scenario_files) <- scenarios

for (code in scenarios) {
  cat("Found", length(scenario_files[[code]]), "files for scenario", code, "\n")
}
cat("\n")

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
      covDep         = sim_obj$covDep
    ),
    similarity_calibration = similarity_calibration,
    results = list()   # ONLY results grouped by alpha
  )

  ## loop over results files
  for (f in scenario_files) {

    cat("Processing:", basename(f), "\n")

    res <- readRDS(f)

    Z_post <- res$mcmcpost$z_post
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

    ## estimate of the clustering using VI
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

    ## silhouette (distance = 1 - PSM)
    diss_mat <- 1 - psm_mat
    sil_obj  <- silhouette(z_hat, dist(diss_mat))
    sil_mean <- mean(sil_obj[, 3])

    tag <- sprintf("alpha_%0.2f", alpha)

    combined$results[[tag]] <- list(
      alpha       = alpha,
      Z_post      = Z_post,
      z_hat       = as.integer(z_hat),
      ARI_vi      = ARI_vi,        # ARI of VI estimate
      ARI_mean    = ARI_mean,      # posterior mean ARI
      ARI_sd      = ARI_sd,        # posterior sd ARI
      silhouette  = list(object = sil_obj, mean = sil_mean),
      psm         = psm_mat,
      seed        = seed,
      N_iter      = N_iter,
      similarity_calibration = res$similarity_calibration
    )
  }


  combined
}

##---------------------------------------------------------##
## Build res & save for each scenario
##---------------------------------------------------------##
out_dir <- here("simulation", "fixed_alpha", "results", "processed", "oneCov")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (code in scenarios) {

  if (length(scenario_files[[code]]) == 0) {
    warning("No results found for scenario ", code,
            " with calibration ", similarity_calibration, ". Skipping.")
    next
  }

  sim_path <- file.path(data_dir, sprintf("binarySBM_1cov_%s.rds", code))

  combined <- res_for_scenario(sim_path, scenario_files[[code]])

  out_file <- file.path(out_dir, sprintf("binarySBM_1cov_%s_results.rds", code))
  saveRDS(combined, out_file)

  cat("Saved:", out_file, "\n\n")
}
