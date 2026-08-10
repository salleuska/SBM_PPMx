##------------------------------------------------------------##
## Process ESBM results for brain application — learned alpha
## - Reads postBrain_*alphaLearn*.rds
## - For each run:
##     * burn-in + thinning of z_post
##     * VI representative via salso
##     * PSM
##     * 3 silhouettes (post, net, cov) using z_hat
##     * alpha_g posterior trace (posterior mean + credible interval)
##------------------------------------------------------------##

library(here)
library(salso)    # VI(), salso(), psm()
library(cluster)  # silhouette()

##---------------------------##
## Settings
##---------------------------##
results_dir <- here("application_brain", "results")
data_path   <- here("application_brain", "consensus_scale33_tau50.RData")

burn_frac <- 0.2
thin_step <- 1

##---------------------------##
## Load application data + preprocess identically to run script
##---------------------------##
load(data_path)  # expects: A_bin_tau, coords

stopifnot(exists("A_bin_tau"), exists("coords"))
diag(A_bin_tau) <- 0

if (!is.null(rownames(A_bin_tau))) {
  m <- match(rownames(A_bin_tau), coords$node)
  stopifnot(!anyNA(m))
  coords <- coords[m, , drop = FALSE]
}

deg  <- rowSums(A_bin_tau)
keep <- deg > 0

Y_app       <- A_bin_tau[keep, keep, drop = FALSE]
coords_kept <- coords[keep, , drop = FALSE]

X_raw <- as.matrix(coords_kept[, c("x","y","z")])
X_std <- scale(X_raw, center = TRUE, scale = TRUE)

D_Y <- as.matrix(dist(Y_app))
D_X <- as.matrix(dist(X_std))

##---------------------------##
## List learned-alpha result files only
##---------------------------##
files_all <- list.files(
  results_dir,
  pattern = "^postBrain_.*alphaLearn.*\\.rds$",
  full.names = TRUE
)

cat("Found", length(files_all), "learned-alpha result files.\n\n")

if (length(files_all) == 0) {
  stop("No learned-alpha result files found in: ", results_dir)
}

##---------------------------##
## Combined structure
##---------------------------##
combined <- list(
  app = list(
    data_path        = data_path,
    n_nodes          = nrow(Y_app),
    n_edges          = sum(Y_app) / 2,
    dropped_isolates = sum(!keep),
    covariates       = c("x","y","z"),
    burn_frac        = burn_frac,
    thin_step        = thin_step
  ),
  results = list()
)

##---------------------------##
## Process each result file
##---------------------------##
for (f in files_all) {

  cat("Processing:", basename(f), "\n")
  res <- readRDS(f)

  # handle both old (matrix) and new (list) esbm return structures
  if (is.list(res$Z_post) && !is.null(res$Z_post$z_post)) {
    z_post       <- res$Z_post$z_post
    alpha_g_post <- res$Z_post$alpha_g_post   # J x N_iter or NULL
  } else {
    z_post       <- res$Z_post
    alpha_g_post <- NULL
  }

  if (is.null(z_post)) stop("Missing z_post in: ", f)

  if (nrow(z_post) != nrow(Y_app)) {
    stop("Dimension mismatch for file: ", basename(f),
         "\n  nrow(z_post)=", nrow(z_post), " nrow(Y_app)=", nrow(Y_app))
  }

  # alpha_g spec (list with init/a_alpha/b_alpha for learned runs)
  alpha_spec <- res$alpha_g

  seed   <- if (!is.null(res$seed))   res$seed   else NA_integer_
  N_iter <- if (!is.null(res$N_iter)) res$N_iter else NA_integer_

  # burn-in & thinning
  T_total <- ncol(z_post)
  burn    <- floor(burn_frac * T_total)
  keep_i  <- seq(from = burn + 1, to = T_total, by = thin_step)
  Z_keep  <- z_post[, keep_i, drop = FALSE]

  # PSM + VI representative
  psm_mat <- psm(t(Z_keep))
  z_hat   <- as.integer(salso(t(Z_keep), loss = VI()))
  K_hat   <- length(unique(z_hat))

  # silhouettes
  sil_mean <- function(diss_mat) {
    if (length(unique(z_hat)) < 2) return(NA_real_)
    tryCatch({
      mean(silhouette(z_hat, dist(diss_mat))[, 3])
    }, error = function(e) NA_real_)
  }

  sil_post <- sil_mean(1 - psm_mat)
  sil_net  <- sil_mean(D_Y)
  sil_cov  <- sil_mean(D_X)

  # alpha_g posterior summary (after burn-in)
  alpha_summary <- NULL
  if (!is.null(alpha_g_post)) {
    alpha_kept <- alpha_g_post[, keep_i, drop = FALSE]
    alpha_summary <- data.frame(
      covariate = c("x", "y", "z")[seq_len(nrow(alpha_kept))],
      mean      = rowMeans(alpha_kept),
      sd        = apply(alpha_kept, 1, sd),
      q025      = apply(alpha_kept, 1, function(x) bayestestR::hdi(x, ci = 0.95)$CI_low),
      q975      = apply(alpha_kept, 1, function(x) bayestestR::hdi(x, ci = 0.95)$CI_high)
    )
    cat(sprintf("  K_hat = %d | alpha means: %s\n",
                K_hat,
                paste(sprintf("%.3f", alpha_summary$mean), collapse = ", ")))
  } else {
    cat(sprintf("  K_hat = %d | alpha_g_post not available\n", K_hat))
  }

  tag <- basename(f)
  tag <- sub("\\.rds$", "", tag)

  combined$results[[tag]] <- list(
    file          = f,
    alpha_spec    = alpha_spec,
    alpha_summary = alpha_summary,
    alpha_g_post  = if (!is.null(alpha_g_post)) alpha_g_post[, keep_i, drop = FALSE] else NULL,
    seed          = seed,
    N_iter        = N_iter,
    z_hat         = z_hat,
    K_hat         = K_hat,
    silhouettes   = c(post = sil_post, net = sil_net, cov = sil_cov),
    psm           = psm_mat
  )
}

##---------------------------##
## Print summary table
##---------------------------##
cat("\n--- Summary ---\n")
for (tag in names(combined$results)) {
  r <- combined$results[[tag]]
  sil <- r$silhouettes
  alpha_str <- if (!is.null(r$alpha_summary)) {
    paste(sprintf("%.3f", r$alpha_summary$mean), collapse = "/")
  } else {
    "NA"
  }
  cat(sprintf("  K=%d | sil_post=%.3f | alpha_post_mean=[%s] | %s\n",
              r$K_hat, sil["post"], alpha_str, basename(r$file)))
}

##---------------------------##
## Save combined object
##---------------------------##
out_dir <- here("application_brain", "results", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "brain_xyz_learnAlpha_processed_results.rds")
saveRDS(combined, out_file)

cat("\nSaved processed results to:\n  ", out_file, "\n")
