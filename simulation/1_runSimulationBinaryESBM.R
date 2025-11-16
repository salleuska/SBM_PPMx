# ---- Load packages ----
suppressPackageStartupMessages({
  library(Rcpp)
  library(here)
  library(igraph)
  library(tools)
})
## -------------------------------
## Parse command-line args
##   file=path/to/file.rds
##   alpha=1
##   seed=123
##   N_iter=10000
## -------------------------------
args <- R.utils::commandArgs(asValue = TRUE)

# data_path  (REQUIRED)
if (is.null(args$data_path)) {
  stop("Please provide data_path=path/to/simulated.rds")
} else {
  data_path <- args$data_path
}

# alpha (default 1)
if (is.null(args$alpha)) {
  alpha_in <- 1
} else {
  alpha_in <- as.numeric(args$alpha)
}

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

cat("Running ESBM with:\n")
cat("  file     =", data_path, "\n")
cat("  alpha    =", alpha_in, "\n")
cat("  seed     =", seed, "\n")
cat("  N_iter   =", N_iter, "\n\n")

## -------------------------------
## Source ESBM code
## -------------------------------
source(here("R_utilties", "esbm.R"))
source(here("R_utilties", "similarity_functions.R"))
Rcpp::sourceCpp(here("R_utilties", "stirling.cpp"))

## -------------------------------
## Hyperparameters
## -------------------------------
a       <- 1
b       <- 1
beta_dm <- 1
H_dm    <- 4

set.seed(seed)

## -------------------------------
## Load the simulated object
## -------------------------------
sim_obj <- readRDS(data_path)

Y  <- sim_obj$Y
x  <- sim_obj$x
z0 <- sim_obj$partition

## -------------------------------
## Infer scenario from file name
## -------------------------------
fname <- basename(data_path)

if (grepl("none", fname, ignore.case = TRUE)) {
  scenario <- "none"
} else if (grepl("one", fname, ignore.case = TRUE)) {
  scenario <- "one"
} else if (grepl("both", fname, ignore.case = TRUE)) {
  scenario <- "both"
} else {
  stop("Cannot infer scenario from file name: ", fname,
       "\nExpected it to include 'none', 'one', or 'both'.")
}

cat("Simulation scenario =", scenario, "\n\n")

## -------------------------------
## Prepare x, similarity_fun, sim_args, alpha
## -------------------------------
if (scenario == "none") {

  x_df           <- NULL
  similarity_fun <- NULL
  sim_args       <- list()
  alpha_g_vec    <- 1

} else if (scenario == "one") {

  x_df <- data.frame(x1 = x[, 1])

  similarity_fun <- similarity_ppmx_gaussian_mean
  sim_args       <- list(list(m0 = 0, s0 = sqrt(1)))
  alpha_g_vec    <- alpha_in

} else if (scenario == "both") {

  x_df <- data.frame(
    x1 = x[, 1],
    x2 = x[, 2]
  )

  similarity_fun <- list(
    similarity_ppmx_gaussian_mean,
    similarity_ppmx_gaussian_mean
  )

  sim_args <- list(
    list(m0 = 0, s0 = sqrt(1)),
    list(m0 = 0, s0 = sqrt(1))
  )

  alpha_g_vec <- rep(alpha_in, 2)
}

## -------------------------------
## Run ESBM
## -------------------------------
cat("Starting ESBM sampler...\n")

Z_post <- esbm(
  Y         = Y,
  seed      = seed,
  N_iter    = N_iter,
  prior     = "DM",
  z_init    = sample(1:4, nrow(Y), replace = TRUE),
  a         = a,
  b         = b,
  beta_DM   = beta_dm,
  H_DM      = H_dm,
  x         = x_df,
  similarity_fun = similarity_fun,
  sim_args       = sim_args,
  alpha_g        = alpha_g_vec
)

cat("Sampler finished.\n\n")

## -------------------------------
## Save results
## -------------------------------
out_dir <- here("simulation", "results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

alpha_tag <- gsub("\\.", "p", sprintf("%.2f", alpha_in))

base_name <- file_path_sans_ext(basename(data_path))

out_file <- file.path(
  out_dir,
  sprintf("post_%s_alpha-%s_seed-%d.rds",
          base_name, alpha_tag, seed)
)

saveRDS(
  list(
    Z_post    = Z_post,
    file    = data_path,
    alpha   = alpha_in,
    N_iter  = N_iter,
    scenario= scenario,
    seed = seed
  ),
  file = out_file
)

cat("Saved posterior to:\n  ", out_file, "\n")