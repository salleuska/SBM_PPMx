##------------------------------------------------------------##
## Process ESBM results for brain application (three covariates)
## - Reads postBrain_*.rds
## - For each run:
##     * burn-in + thinning of Z_post
##     * VI representative via salso
##     * PSM
##     * 3 silhouettes (post, net, cov) using z_hat
##------------------------------------------------------------##

library(here)
library(salso)      # VI(), salso(), psm()
library(cluster)    # silhouette()

##---------------------------##
## Settings
##---------------------------##
results_dir <- here("application_brain", "results")
data_path   <- here("application_brain", "consensus_scale33_tau50.RData")

burn_frac <- 0.2
thin_step <- 1

##---------------------------##
## Load application data + preprocess identically
##---------------------------##
load(data_path)  # expects: A_bin_tau, coords

stopifnot(exists("A_bin_tau"), exists("coords"))
diag(A_bin_tau) <- 0

# align coords to adjacency order if rownames exist
if (!is.null(rownames(A_bin_tau))) {
  m <- match(rownames(A_bin_tau), coords$node)
  stopifnot(!anyNA(m))
  coords <- coords[m, , drop = FALSE]
}

# drop isolates
deg  <- rowSums(A_bin_tau)
keep <- deg > 0

Y_app <- A_bin_tau[keep, keep, drop = FALSE]
coords_kept <- coords[keep, , drop = FALSE]

# standardized covariates
X_raw <- as.matrix(coords_kept[, c("x","y","z")])
X_std <- scale(X_raw, center = TRUE, scale = TRUE)

# fixed distances (do not depend on alpha)
D_Y <- as.matrix(dist(Y_app))   # network: adjacency-profile distance
D_X <- as.matrix(dist(X_std))   # covariates: Euclidean on standardized xyz

##---------------------------##
## List result files
##---------------------------##
files_all <- list.files(
  results_dir,
  pattern = "^postBrain_.*\\.rds$",
  full.names = TRUE
)

cat("Found", length(files_all), "application result files.\n\n")

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

  # handle both old (bare matrix) and new (list with $z_post) esbm return structures
  if (is.list(res$Z_post) && !is.null(res$Z_post$z_post)) {
    Z_post <- res$Z_post$z_post
  } else {
    Z_post <- res$Z_post
  }
  if (is.null(Z_post)) stop("Missing Z_post in: ", f)

  # must match n after isolate removal
  if (nrow(Z_post) != nrow(Y_app)) {
    stop("Dimension mismatch for file: ", basename(f),
         "\n  nrow(Z_post)=", nrow(Z_post), " nrow(Y_app)=", nrow(Y_app))
  }

  # alphas: support both old (res$alpha) and new (res$alpha_g numeric vec) field names
  alpha_vec <- if (!is.null(res$alpha)) res$alpha else if (is.numeric(res$alpha_g)) res$alpha_g else NULL
  if (is.null(alpha_vec) || length(alpha_vec) < 3) {
    stop("Missing/invalid alpha (need length 3) in: ", basename(f))
  }
  alpha_x <- as.numeric(alpha_vec[1])
  alpha_y <- as.numeric(alpha_vec[2])
  alpha_z <- as.numeric(alpha_vec[3])

  seed   <- if (!is.null(res$seed))   res$seed   else NA_integer_
  N_iter <- if (!is.null(res$N_iter)) res$N_iter else NA_integer_

  # burn-in & thinning
  T_total <- ncol(Z_post)
  burn    <- floor(burn_frac * T_total)
  keep_i  <- seq(from = burn + 1, to = T_total, by = thin_step)
  Z_keep  <- Z_post[, keep_i, drop = FALSE]

  # PSM + VI representative
  psm_mat <- psm(t(Z_keep))
  z_hat   <- as.integer(salso(t(Z_keep), loss = VI()))
  K_hat   <- length(unique(z_hat))

  # silhouettes (return NA if only 1 cluster or if silhouette fails)
  sil_mean <- function(diss_mat) {
    if (length(unique(z_hat)) < 2) return(NA_real_)
    out <- tryCatch({
      sil_obj <- silhouette(z_hat, dist(diss_mat))
      mean(sil_obj[, 3])
    }, error = function(e) NA_real_)
    out
  }

  sil_post <- sil_mean(1 - psm_mat)
  sil_net  <- sil_mean(D_Y)
  sil_cov  <- sil_mean(D_X)

  tag <- sprintf("ax_%0.2f__ay_%0.2f__az_%0.2f", alpha_x, alpha_y, alpha_z)

  combined$results[[tag]] <- list(
    file        = f,
    alpha       = c(alpha_x, alpha_y, alpha_z),   # fixed alpha values
    seed        = seed,
    N_iter      = N_iter,
    z_hat       = z_hat,
    K_hat       = K_hat,
    silhouettes = c(post = sil_post, net = sil_net, cov = sil_cov),
    psm         = psm_mat
  )
}

##---------------------------##
## Save combined object
##---------------------------##
out_dir <- here("application_brain", "results", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "brain_xyz_processed_results.rds")
saveRDS(combined, out_file)

cat("\nSaved processed results to:\n  ", out_file, "\n")
