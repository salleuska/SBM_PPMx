# ---- Load packages ----
suppressPackageStartupMessages({
  library(Rcpp)
  library(here)
  library(igraph)
  library(tools)
})
## ------------------------------- ##
## Parse command-line args
## example
# args <- list()
# args$data_path= "simulation/data/binarySBM_one.rds"
# args$alpha1=1
# args$alpha2=2
# args$seed=13
## ------------------------------- ##
args <- R.utils::commandArgs(asValue = TRUE)

# data_path  (REQUIRED)
if (is.null(args$data_path)) {
  stop("Please provide data_path=path/to/simulated.rds")
} else {
  data_path <- args$data_path
}

# alpha1, alpha2 for two covariates
alpha1_in <- if (is.null(args$alpha1)) 0 else as.numeric(args$alpha1)
alpha2_in <- if (is.null(args$alpha2)) 0 else as.numeric(args$alpha2)

cat("  alpha1  =", alpha1_in, "\n")
cat("  alpha2  =", alpha2_in, "\n\n")

# seed (default 123)
if (is.null(args$seed)) {
  seed <- 123
} else {
  seed <- as.integer(args$seed)
}

# N_iter (default 10000)
if (is.null(args$N_iter)) {
  N_iter <- 10000
} else {
  N_iter <- as.integer(args$N_iter)
}

## ------------------------------- ##
## Source ESBM code
## ------------------------------- ##
source(here("R_utilties", "esbm.R"))
source(here("R_utilties", "similarity_functions.R"))
Rcpp::sourceCpp(here("R_utilties", "stirling.cpp"))

## ------------------------------- ##
## Hyperparameters
## ------------------------------- ##
a       <- 1
b       <- 1
beta_dm <- 1
H_dm    <- 4

set.seed(seed)

## ------------------------------- ##
## Load the simulated object
## ------------------------------- ##
sim_obj <- readRDS(data_path)

Y  <- sim_obj$Y
x  <- sim_obj$x
z0 <- sim_obj$partition

## ------------------------------- ##
## Infer scenario from file name
## ------------------------------- ##
fname <- basename(data_path)

if (grepl("2cov_NN", fname, ignore.case = TRUE)) {
  scenario <- "NN"
} else if (grepl("2cov_IN", fname, ignore.case = TRUE)) {
  scenario <- "IN"
} else if (grepl("2cov_NI", fname, ignore.case = TRUE)) {
  scenario <- "NI"
} else if (grepl("2cov_II", fname, ignore.case = TRUE)) {
  scenario <- "II"
} else {
  stop("Cannot infer two-cov scenario from file name: ", fname,
       "\nExpected it to include one of: 2cov_NN / 2cov_IN / 2cov_NI / 2cov_II.")
}

cat("Simulation scenario =", scenario, "\n\n")

## ------------------------------- ##
## Prepare x, similarity_fun, sim_args, alpha
## ------------------------------- ##

## define data frame with two covariates
x_df <- data.frame(
  x1 = x[, 1],
  x2 = x[, 2]
)

# always provide two similarity components
similarity_fun <- list(
  similarity_ppmx_gaussian_mean,
  similarity_ppmx_gaussian_mean
)

# always provide args for both components
sim_args <- list(
  list(m0 = 0, s0 = sqrt(1)),
  list(m0 = 0, s0 = sqrt(1))
)

# always provide two alphas
alpha_g_vec <- c(alpha1_in, alpha2_in)

## ------------------------------- ##
## Run ESBM
## ------------------------------- ##
cat("Starting ESBM sampler...\n")

Z_post <- esbm(
  Y         = Y,
  seed      = seed,
  N_iter    = N_iter,
  prior  = list(
    name    = "DM",
    beta_DM = beta_dm,
    H_DM    = H_dm
  ),
  z_init    = sample(1:4, nrow(Y), replace = TRUE),
  a         = a,
  b         = b,
  x         = x_df,
  similarity_fun = similarity_fun,
  sim_args       = sim_args,
  alpha_g        = alpha_g_vec
)

cat("Sampler finished.\n\n")

## ------------------------------- ##
## Save results
## ------------------------------- ##
out_dir <- here("simulation", "results", "twoCov")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

alpha1_tag <- gsub("\\.", "p", sprintf("%.2f", alpha1_in))
alpha2_tag <- gsub("\\.", "p", sprintf("%.2f", alpha2_in))
base_name <- file_path_sans_ext(basename(data_path))

out_file <- file.path(
  out_dir,
  sprintf("postTwoCov_%s_a1-%s_a2-%s_seed-%d.rds",
          base_name, alpha1_tag, alpha2_tag, seed)
)


saveRDS(
  list(
    Z_post    = Z_post,
    file    = data_path,
    alpha   = alpha_g_vec,
    N_iter  = N_iter,
    scenario= scenario,
    seed = seed
  ),
  file = out_file
)

cat("Saved posterior to:\n  ", out_file, "\n")