##---------------------------##------------------##
## Process ESBM results for binary SBM simulation
##- Reads all post_*.rds in results folder
##- For each:
##     * computes Dahl clustering from Z_DM
##     * computes ARI vs true partition
##     * records scenario, alpha, seed, K_hat
##- Saves a summary CSV and an optional ARI vs alpha plot
##---------------------------##------------------##
library(here)
library(mclust)     # ARI
library(ggplot2)    # plotting
library(tools)      # file_path_sans_ext
library(salso)      # point estimate for the clustering

library(ggplot2)
library(RColorBrewer)

plot_psm_gg <- function(psm, title = "") {
  V <- nrow(psm)
  
  # long data.frame in base R
  grid_idx <- expand.grid(
    row = seq_len(V),
    col = seq_len(V)
  )
  grid_idx$value <- as.vector(psm)
  
  # flip rows so origin is bottom-left like a matrix
  grid_idx$row_f <- factor(grid_idx$row, levels = rev(seq_len(V)))
  grid_idx$col_f <- factor(grid_idx$col, levels = seq_len(V))
  
  pal <- colorRampPalette(brewer.pal(9, "Blues"))(50)
  
  ggplot(grid_idx, aes(x = col_f, y = row_f, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = pal, limits = c(0, 1)) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      axis.title   = element_blank(),
      panel.grid   = element_blank(),
      plot.title   = element_text(hjust = 0.5)
    ) +
    ggtitle(title)
}
##---------------------------##
## Settings
##---------------------------##

## results_dir <- here("simulation", "results", "oneCov")

## if two covariates
results_dir <- here("simulation", "results", "twoCov")

burn_frac <- 0.2 # brunin
thin_step <- 1   # eventual thinning

##---------------------------##
## List result files
##---------------------------##
## ---------------------------
## Process each result
## ---------------------------

summ_list <- vector("list", length(files_res))
psm_list  <- vector("list", length(files_res))

for (i in seq_along(files_res)) {
  f_res <- files_res[i]
  cat("Processing:", basename(f_res), "...\n")

  res <- readRDS(f_res)

  ## Extract posterior samples and metadata
  Z_post   <- res$Z_post      # V x N_iter
  scen_res <- res$scenario
  seed     <- res$seed
  sim_path <- res$file

  # Two alphas: be slightly defensive
  alpha1 <- if (!is.null(res$alpha1)) res$alpha1 else {
    if (!is.null(res$alpha)) res$alpha else NA_real_
  }
  alpha2 <- if (!is.null(res$alpha2)) res$alpha2 else 0

  # Load simulation to get true partition
  sim_obj <- readRDS(here(sim_path))
  z_true  <- sim_obj$partition
  V_true  <- length(z_true)

  if (nrow(Z_post) != V_true) {
    stop("Dimension mismatch: nrow(Z_post) != length(z_true) in ", f_res)
  }

  ## ---------------------------
  ## Posterior similarity + VI clustering
  ## ---------------------------
  T_total <- ncol(Z_post)
  burn    <- floor(T_total * burn_frac)
  keep_i  <- seq(from = burn + 1L, to = T_total, by = thin_step)

  Z_keep <- Z_post[, keep_i, drop = FALSE]   # V x T_keep

  if (ncol(Z_keep) < 2L) {
    stop("Too few post-burn-in samples after thinning in ", f_res)
  }

  # PSM (iterations x items → transpose)
  coClust_prob <- psm(t(Z_keep))

  # VI representative clustering
  z_hat <- salso(t(Z_keep), loss = VI())

  ## ---------------------------
  ## ARI vs true clustering
  ## ---------------------------
  ari   <- mclust::adjustedRandIndex(z_hat, z_true)
  K_hat <- length(unique(z_hat))

  # file base name
  base_name <- file_path_sans_ext(basename(f_res))

  summ_list[[i]] <- data.frame(
    file_res   = basename(f_res),
    sim_file   = sim_path,
    base_name  = base_name,
    scenario   = scen_res,
    alpha1     = alpha1,
    alpha2     = alpha2,
    seed       = seed,
    K_hat      = K_hat,
    ARI        = ari,
    stringsAsFactors = FALSE
  )

  psm_list[[i]] <- coClust_prob
}

summary_df <- do.call(rbind, summ_list)
names(psm_list) <- summary_df$file_res

## ---------------------------
## Save summary + PSMs
## ---------------------------

out_csv <- file.path(results_dir, "summary_twoAlpha_binaryESBM.csv")
write.csv(summary_df, out_csv, row.names = FALSE)
cat("Saved summary to:\n  ", out_csv, "\n")

psm_out <- file.path(results_dir, "psm_twoAlpha_binaryESBM.rds")
saveRDS(psm_list, psm_out)
cat("Saved PSM list to:\n  ", psm_out, "\n")

## ---------------------------
## ARI heatmaps over (alpha1, alpha2)
## ---------------------------

# Simple ARI heatmap per scenario
p_heat <- ggplot(summary_df,
                 aes(x = alpha1, y = alpha2, fill = ARI)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  facet_wrap(~ scenario) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  labs(
    title = "ARI over (alpha1, alpha2) by scenario",
    x = expression(alpha[1]),
    y = expression(alpha[2]),
    fill = "ARI"
  )

print(p_heat)

