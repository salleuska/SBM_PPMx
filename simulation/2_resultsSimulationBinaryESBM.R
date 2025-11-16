##---------------------------------------------##
## Process ESBM results for binary SBM simulation
## - Reads all post_*.rds in results folder
## - For each:
##     * computes Dahl clustering from Z_DM
##     * computes ARI vs true partition
##     * records scenario, alpha, seed, K_hat
## - Saves a summary CSV and an optional ARI vs alpha plot
## ---------------------------------------------##

suppressPackageStartupMessages({
  library(here)
  library(mclust)     # ARI
  library(ggplot2)    # plotting
  library(tools)      # file_path_sans_ext
  library(salso)      # Dahl via salso::dahl
})

## ---------------------------
## Settings
## ---------------------------

results_dir <- here("simulation", "results")

burn_frac <- 0.5   # drop first 50% iterations
thin_step <- 1     # keep every 

## ---------------------------
## List result files
## ---------------------------

files_res <- list.files(
  results_dir,
  pattern = "^post_.*\\.rds$",
  full.names = TRUE
)

if (length(files_res) == 0L) {
  stop("No result files found in ", results_dir)
}

cat("Found", length(files_res), "result files.\n")

## ---------------------------
## Process each result
## ---------------------------

summ_list <- vector("list", length(files_res))

for (i in seq_along(files_res)) {
  f_res <- files_res[i]
  cat("Processing:", basename(f_res), "...\n")

  res <- readRDS(f_res)

  # Extract posterior samples and metadata
  Z_post   <- res$Z_DM      # V x N_iter
  alpha    <- res$alpha
  seed     <- res$seed
  scen_res <- res$scenario
  sim_path <- res$file      # path to original simulation .rds

  sim_obj <- readRDS(here(sim_path))
  z_true  <- sim_obj$partition
  V_true  <- length(z_true)

  ##---------------------------##
  ## estimate clustering via salso
  ##---------------------------##
  ## keep iterations after burnin
  T_total <- ncol(Z_post)
  burn    <- floor(T_total * burn_frac)
  keep_i  <- seq(from = burn + 1L, to = T_total, by = thin_step)

  Z_keep <- Z_post[, keep_i, drop = FALSE]   # V x T_keep

  if (ncol(Z_keep) < 2L) {
    stop("Too few post-burn-in samples after thinning in ", f_res)
  }

  coClust_prob <- psm(t(Z_keep))
  ## use variation of information to estimate clustering
  z_hat <- salso(t(Z_keep), loss=VI())
  ## ---------------------------
  ## ARI vs true clustering
  ## ---------------------------

  ari <- mclust::adjustedRandIndex(z_hat, z_true)

  # Number of inferred clusters
  K_hat <- length(unique(z_hat))

  # Extract alpha tag & base name from filename (optional)
  base_name <- file_path_sans_ext(basename(f_res))

  summ_list[[i]] <- data.frame(
    file_res   = basename(f_res),
    sim_file   = sim_path,
    base_name  = base_name,
    scenario   = scen_res,
    alpha      = alpha,
    seed       = seed,
    coClust_prob = coClust_prob,
    K_hat      = K_hat,
    ARI        = ari,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summ_list)

## ---------------------------
## Save summary
## ---------------------------

out_csv <- file.path(results_dir, "summary_binaryESBM.csv")
write.csv(summary_df, out_csv, row.names = FALSE)

cat("Saved summary to:\n  ", out_csv, "\n")

## ---------------------------
## Optional: quick ARI vs alpha plot
## ---------------------------

p <- ggplot(summary_df,
            aes(x = alpha, y = ARI,
                color = scenario,
                group = interaction(scenario, seed))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(title = "ARI vs alpha",
       x = expression(alpha),
       y = "Adjusted Rand Index")

print(p)

