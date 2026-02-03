library(Rcpp)
library(here)
library(tools)
library(R.utils)

args <- R.utils::commandArgs(asValue = TRUE)

## ------------------------------- ##
## Parse args
## ------------------------------- ##
if (is.null(args$data_path)) {
  stop("Please provide data_path=path/to/brain_data.RData")
} else {
  data_path <- args$data_path
}

seed   <- if (is.null(args$seed)) 123 else as.integer(args$seed)
N_iter <- if (is.null(args$N_iter)) 10000 else as.integer(args$N_iter)

alpha_x <- if (is.null(args$alpha_x)) 1 else as.numeric(args$alpha_x)
alpha_y <- if (is.null(args$alpha_y)) 1 else as.numeric(args$alpha_y)
alpha_z <- if (is.null(args$alpha_z)) 1 else as.numeric(args$alpha_z)

cat("  alpha_x =", alpha_x, "\n")
cat("  alpha_y =", alpha_y, "\n")
cat("  alpha_z =", alpha_z, "\n\n")

# Always drop isolated nodes
drop_isolates <- TRUE

cat("Running ESBM (brain application) with:\n")
cat("  file          =", data_path, "\n")
cat("  alpha         =", paste(alpha_vec, collapse = ", "), "\n")
cat("  seed          =", seed, "\n")
cat("  N_iter        =", N_iter, "\n")
cat("  covariates    = x, y, z (all)\n")
cat("  drop_isolates = TRUE\n\n")

set.seed(seed)

## ------------------------------- ##
## Load code
## ------------------------------- ##
source(here("R_utilties", "esbm.R"))
source(here("R_utilties", "similarity_functions.R"))
Rcpp::sourceCpp(here("R_utilties", "stirling.cpp"))

## ------------------------------- ##
## Hyperparameters
## ------------------------------- ##
a <- 1
b <- 1

## use a gnedin prior
gamma_GN <- 0.3

## ------------------------------- ##
## Load brain data
## ------------------------------- ##
# expects:
#   A_bin_tau : symmetric binary adjacency (N x N)
#   coords    : data.frame with columns node, x, y, z (N rows)
load(data_path)

## ------------------------------- ##
## Drop isolated nodes (degree 0)
## ------------------------------- ##
deg  <- rowSums(A_bin_tau)
keep <- deg > 0

cat(sprintf("Dropping %d isolated nodes (degree 0) out of %d total.\n\n",
            sum(!keep), length(keep)))

Y <- A_bin_tau[keep, keep, drop = FALSE]
coords_kept <- coords[keep, , drop = FALSE]

## ------------------------------- ##
## Covariates: use ALL coordinates (x,y,z), standardized
## ------------------------------- ##
X_raw <- as.matrix(coords_kept[, c("x", "y", "z")])
colnames(X_raw) <- c("x", "y", "z")

# center + scale each column to mean 0, sd 1
X_std <- scale(X_raw, center = TRUE, scale = TRUE)

# quick sanity
mu <- colMeans(X_std)
sdv <- apply(X_std, 2, sd)
cat(sprintf("X_std means: x=%.3f y=%.3f z=%.3f\n", mu[1], mu[2], mu[3]))
cat(sprintf("X_std sds:   x=%.3f y=%.3f z=%.3f\n\n", sdv[1], sdv[2], sdv[3]))

## ------------------------------- ##
## Similarity setup (multivariate)
## ------------------------------- ##
## define data frame with three covariates (standardized)
x_df <- data.frame(
  x1 = X_std[, 1],
  x2 = X_std[, 2],
  x3 = X_std[, 3]
)

#  provide three similarity components
similarity_fun <- list(
  similarity_ppmx_gaussian_mean,
  similarity_ppmx_gaussian_mean,
  similarity_ppmx_gaussian_mean
)

#  provide args for all three components
sim_args <- list(
  list(m0 = 0, s0 = sqrt(1)),
  list(m0 = 0, s0 = sqrt(1)),
  list(m0 = 0, s0 = sqrt(1))
)

#  provide three alphas
alpha_g_vec <- c(alpha_x, alpha_y, alpha_z)
## ------------------------------- ##
## Run ESBM
## ------------------------------- ##
cat("Starting ESBM...\n")

Z_post <- esbm(
  Y              = Y,
  seed           = seed,
  N_iter         = N_iter,
  prior          = "GN",
  z_init         = rep(1, nrow(Y)),   # safe default for GN
  a              = a,
  b              = b,
  gamma_GN       = gamma_GN,
  x              = X_std,
  similarity_fun = similarity_fun,
  sim_args       = sim_args,
  alpha_g        = alpha_vec
)

cat("Done.\n\n")

## ------------------------------- ##
## Save output
## ------------------------------- ##
out_dir <- here("application_brain", "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

alpha_tag <- paste0(
  "ax-", gsub("\\.", "p", sprintf("%.2f", alpha_x)),
  "_ay-", gsub("\\.", "p", sprintf("%.2f", alpha_y)),
  "_az-", gsub("\\.", "p", sprintf("%.2f", alpha_z))
)

base_name <- file_path_sans_ext(basename(data_path))

out_file <- file.path(
  out_dir,
  sprintf("postBrain_xyz_%s_%s_seed-%d.rds",
          base_name, alpha_tag, seed)
)

saveRDS(
  list(
    Z_post        = Z_post,
    file          = data_path,
    alpha         = alpha_vec,
    seed          = seed,
    N_iter        = N_iter,
    gamma_GN      = gamma_GN,
    covariates    = c("x", "y", "z"),
    drop_isolates = TRUE,
    n_nodes       = nrow(Y),
    n_edges       = sum(Y) / 2
  ),
  file = out_file
)

cat("Saved posterior to:\n  ", out_file, "\n")