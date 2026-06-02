##---------------------------##
## Process learned-alpha ESBM results: two covariates
##---------------------------##
library(here)
library(mclust)
library(tools)
library(salso)
library(cluster)

##---------------------------##
## Settings
##---------------------------##

results_dir <- here("simulation", "learn_alpha", "results", "twoCov")

sim_paths <- list(
  NN = here("simulation", "data", "binarySBM_2cov_NN.rds"),
  IN = here("simulation", "data", "binarySBM_2cov_IN.rds"),
  NI = here("simulation", "data", "binarySBM_2cov_NI.rds"),
  II = here("simulation", "data", "binarySBM_2cov_II.rds")
)

burn_frac <- 0.4
thin_step <- 1

files_all <- list.files(
  results_dir,
  pattern = "^postTwoCov_.*\\.rds$",
  full.names = TRUE
)

files_by_scenario <- list(
  NN = files_all[grepl("2cov_NN", basename(files_all), ignore.case = TRUE)],
  IN = files_all[grepl("2cov_IN", basename(files_all), ignore.case = TRUE)],
  NI = files_all[grepl("2cov_NI", basename(files_all), ignore.case = TRUE)],
  II = files_all[grepl("2cov_II", basename(files_all), ignore.case = TRUE)]
)

for (sc in names(files_by_scenario)) {
  cat("Found", length(files_by_scenario[[sc]]), "files for scenario", sc, "\n")
}

res_for_scenario <- function(sim_path, scenario_files) {

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
      cov1Dep        = sim_obj$cov1Dep,
      cov2Dep        = sim_obj$cov2Dep
    ),
    results = list()
  )

  for (f in scenario_files) {

    tag <- file_path_sans_ext(basename(f))
    cat("Processing:", basename(f), "\n")

    res <- readRDS(f)

    alpha_g <- res$alpha_g
    Z_post  <- res$mcmcpost$z_post
    seed    <- res$seed
    N_iter  <- res$N_iter

    T_total <- ncol(Z_post)
    burn    <- floor(burn_frac * T_total)
    keep_i  <- seq(from = burn + 1L, to = T_total, by = thin_step)

    Z_keep <- Z_post[, keep_i, drop = FALSE]

    alpha_post <- res$mcmcpost$alpha_g_post
    alpha_keep <- NULL
    if (!is.null(alpha_post)) {
      alpha_keep <- alpha_post[, keep_i, drop = FALSE]
    }

    alpha_rate_post <- res$mcmcpost$alpha_rate_post
    alpha_rate_keep <- NULL
    if (!is.null(alpha_rate_post)) {
      alpha_rate_keep <- alpha_rate_post[, keep_i, drop = FALSE]
    }

    psm_mat <- psm(t(Z_keep))
    z_hat   <- salso(t(Z_keep), loss = VI())

    ARI_vi <- adjustedRandIndex(z_hat, z_true)

    ARI_iter <- apply(
      Z_keep, 2L,
      function(z_t) adjustedRandIndex(z_t, z_true)
    )
    ARI_mean <- mean(ARI_iter)
    ARI_sd   <- sd(ARI_iter)

    diss_mat <- 1 - psm_mat
    sil_obj  <- silhouette(z_hat, dist(diss_mat))
    sil_mean <- mean(sil_obj[, 3])

    combined$results[[tag]] <- list(
      alpha_g         = alpha_g,
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
      N_iter          = N_iter
    )
  }

  combined
}

##---------------------------##
## Build and save
##---------------------------##

out_dir <- here("simulation", "learn_alpha", "results", "processed", "twoCov")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

combined_NN <- res_for_scenario(sim_paths$NN, files_by_scenario$NN)
combined_IN <- res_for_scenario(sim_paths$IN, files_by_scenario$IN)
combined_NI <- res_for_scenario(sim_paths$NI, files_by_scenario$NI)
combined_II <- res_for_scenario(sim_paths$II, files_by_scenario$II)

saveRDS(combined_NN, file.path(out_dir, "scenario_2cov_NN_results.rds"))
saveRDS(combined_IN, file.path(out_dir, "scenario_2cov_IN_results.rds"))
saveRDS(combined_NI, file.path(out_dir, "scenario_2cov_NI_results.rds"))
saveRDS(combined_II, file.path(out_dir, "scenario_2cov_II_results.rds"))# p_rate

