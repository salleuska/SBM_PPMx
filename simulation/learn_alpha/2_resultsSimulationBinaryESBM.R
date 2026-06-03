##---------------------------##------------------##
## Process ESBM results for binary SBM simulation
##---------------------------##------------------##

library(here)
library(mclust)
library(tools)
library(salso)
library(cluster)

##---------------------------##
## Settings
##---------------------------##

similarity_calibration <- "raw"

burn_frac <- 0.4
thin_step <- 1

results_root <- here(
  "simulation",
  "learn_alpha",
  similarity_calibration,
  "results"
)

sim_dir <- here("simulation", "data")

out_dir <- here(
  "simulation",
  "learn_alpha",
  similarity_calibration,
  "results",
  "processed"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##---------------------------##
## Helper extractors
##---------------------------##

get_sim_name <- function(path) {
  sub(
    "^post_(binarySBM_[0-9]+cov_[A-Z]+)_alphaLearn-.*$",
    "\\1",
    tools::file_path_sans_ext(basename(path))
  )
}

get_nCov <- function(path) {
  as.integer(
    sub(
      "^.*binarySBM_([0-9]+)cov_[A-Z]+.*$",
      "\\1",
      basename(path)
    )
  )
}

get_scenario <- function(path) {
  sub(
    "^.*binarySBM_[0-9]+cov_([A-Z]+).*$",
    "\\1",
    basename(path)
  )
}

##---------------------------##
## Find files
##---------------------------##

files_all <- list.files(
  results_root,
  pattern = "^post_binarySBM_.*\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files_all) == 0) {
  stop("No posterior files found in: ", results_root)
}

sim_names <- get_sim_name(files_all)

cat("Posterior files found:\n")
print(table(sim_names))

##---------------------------##
## Processing function
##---------------------------##

res_for_sim <- function(sim_name, scenario_files) {

  sim_path <- here("simulation", "data", paste0(sim_name, ".rds"))

  if (!file.exists(sim_path)) {
    stop("Missing simulation file: ", sim_path)
  }

  sim_obj <- readRDS(sim_path)
  z_true  <- sim_obj$partition

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
    results = list()
  )

  for (f in scenario_files) {

    tag <- tools::file_path_sans_ext(basename(f))

    cat("Processing:", basename(f), "\n")

    res <- readRDS(f)

    Z_post <- res$mcmcpost$z_post
    seed   <- res$seed
    N_iter <- res$N_iter

    T_total <- ncol(Z_post)
    burn    <- floor(burn_frac * T_total)
    keep_i  <- seq(from = burn + 1L, to = T_total, by = thin_step)

    Z_keep <- Z_post[, keep_i, drop = FALSE]

    alpha_post <- res$mcmcpost$alpha_g_post
    alpha_keep <- if (!is.null(alpha_post)) alpha_post[keep_i] else NULL

    alpha_rate_post <- res$mcmcpost$alpha_rate_post
    alpha_rate_keep <- if (!is.null(alpha_rate_post)) alpha_rate_post[keep_i] else NULL

    psm_mat <- psm(t(Z_keep))

    z_hat <- salso(t(Z_keep), loss = VI())

    ARI_vi <- adjustedRandIndex(z_hat, z_true)

    ARI_iter <- apply(
      Z_keep,
      2L,
      function(z_t) adjustedRandIndex(z_t, z_true)
    )

    ARI_mean <- mean(ARI_iter)
    ARI_sd   <- sd(ARI_iter)

    diss_mat <- 1 - psm_mat
    sil_obj  <- silhouette(z_hat, dist(diss_mat))
    sil_mean <- mean(sil_obj[, 3])

    combined$results[[tag]] <- list(
      alpha_g         = res$alpha_g,
      alpha_post      = alpha_keep,
      alpha_rate_post = alpha_rate_keep,
      Z_post          = Z_keep,
      z_hat           = as.integer(z_hat),
      ARI_vi          = ARI_vi,
      ARI_mean        = ARI_mean,
      ARI_sd          = ARI_sd,
      silhouette      = list(object = sil_obj, mean = sil_mean),
      psm             = psm_mat,
      seed            = seed,
      scenario        = res$scenario,
      nCov            = res$nCov,
      N_iter          = N_iter
    )
  }

  combined
}

##---------------------------##
## Process all simulations
##---------------------------##

for (sim_name in unique(sim_names)) {

  scenario_files <- files_all[sim_names == sim_name]

  cat("\n========================================\n")
  cat("Simulation:", sim_name, "\n")
  cat("n files:", length(scenario_files), "\n")
  cat("========================================\n")

  combined <- res_for_sim(sim_name, scenario_files)

  nCov <- get_nCov(sim_name)

  out_subdir <- file.path(out_dir, paste0(nCov, "Cov"))
  dir.create(out_subdir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(
    combined,
    file.path(out_subdir, paste0(sim_name, "_results.rds"))
  )
}

