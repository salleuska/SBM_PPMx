##---------------------------##------------------##
## Process network-only SBM baseline results
##---------------------------##------------------##

args <- R.utils::commandArgs(asValue = TRUE)

library(here)
library(mclust)
library(tools)
library(salso)
library(cluster)

##---------------------------##
## Settings
##---------------------------##

burn_frac <- if (is.null(args$burn_frac)) 0.4 else as.numeric(args$burn_frac)
thin_step <- if (is.null(args$thin_step)) 1 else as.integer(args$thin_step)

results_root <- here(
  "simulation",
  "learn_alpha",
  "SBM",
  "results"
)

out_dir <- here(
  "simulation",
  "learn_alpha",
  "SBM",
  "results",
  "processed"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

##---------------------------##
## Helper extractors
##---------------------------##

get_sim_name <- function(path) {
  sub(
    "^postSBM_(binarySBM_[0-9]+cov_[A-Z0-9]+)_seed-[0-9]+$",
    "\\1",
    tools::file_path_sans_ext(basename(path))
  )
}

get_nCov <- function(path) {
  as.integer(
    sub(
      "^.*binarySBM_([0-9]+)cov_[A-Z0-9]+.*$",
      "\\1",
      basename(path)
    )
  )
}

##---------------------------##
## Find files
##---------------------------##

files_all <- list.files(
  results_root,
  pattern = "^postSBM_binarySBM_.*\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files_all) == 0) {
  stop("No SBM posterior files found in: ", results_root)
}

sim_names <- get_sim_name(files_all)

if (any(sim_names == basename(files_all))) {
  bad_files <- basename(files_all)[sim_names == basename(files_all)]
  stop(
    "Could not extract simulation names from: ",
    paste(bad_files, collapse = ", ")
  )
}

cat("SBM posterior files found:\n")
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
    ARI_lwr  <- quantile(ARI_iter, 0.025)
    ARI_upr  <- quantile(ARI_iter, 0.975)

    diss_mat <- 1 - psm_mat
    sil_obj  <- silhouette(z_hat, dist(diss_mat))
    sil_mean <- mean(sil_obj[, 3])

    combined$results[[tag]] <- list(
      alpha_g         = NULL,
      alpha_post      = NULL,
      alpha_rate_post = NULL,
      Z_post          = Z_keep,
      z_hat           = as.integer(z_hat),
      ARI_vi          = ARI_vi,
      ARI_mean        = ARI_mean,
      ARI_sd          = ARI_sd,
      ARI_lwr         = ARI_lwr,
      ARI_upr         = ARI_upr,
      silhouette      = list(object = sil_obj, mean = sil_mean),
      psm             = psm_mat,
      seed            = seed,
      scenario        = res$scenario,
      nCov            = res$nCov,
      N_iter          = N_iter,
      model           = "SBM"
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
