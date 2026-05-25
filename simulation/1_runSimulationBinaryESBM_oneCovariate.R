library(Rcpp)
library(here)
library(igraph)
library(tools)

args <- R.utils::commandArgs(asValue = TRUE)

## ------------------------------- ##
## Parse args
## ------------------------------- ##

if (is.null(args$data_path)) {
  stop("Please provide data_path=path/to/simulated.rds")
} else {
  data_path <- args$data_path
}

alpha_in <- if (is.null(args$alpha)) 1 else as.numeric(args$alpha)
seed     <- if (is.null(args$seed)) 123 else as.integer(args$seed)
N_iter   <- if (is.null(args$N_iter)) 10000 else as.integer(args$N_iter)

cat("Running ESBM with:\n")
cat("  file     =", data_path, "\n")
cat("  alpha    =", alpha_in, "\n")
cat("  seed     =", seed, "\n")
cat("  N_iter   =", N_iter, "\n\n")

## ------------------------------- ##
## Load code
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
## Load simulation
## ------------------------------- ##
sim_obj <- readRDS(data_path)

Y  <- sim_obj$Y
x  <- sim_obj$x              # matrix n x 2; use only first column for one-cov runs
z0 <- sim_obj$partition

## ------------------------------- ##
## Infer scenario 
## ------------------------------- ##
## Infer scenario from filename
fname <- basename(data_path)

if (grepl("1cov_neutral", fname, ignore.case = TRUE)) {
  scenario <- "neutral"
} else if (grepl("1cov_informative", fname, ignore.case = TRUE)) {
  scenario <- "informative"
} else if (grepl("1cov_mislead_random", fname, ignore.case = TRUE)) {
  scenario <- "mislead_random"
} else if (grepl("1cov_mislead_shifted", fname, ignore.case = TRUE)) {
  scenario <- "mislead_shifted"
} else {
  stop("Filename must contain one of: 1cov_neutral / 1cov_informative / 1cov_mislead_random / 1cov_mislead_shifted. Found: ", fname)
}

cat("Scenario:", scenario, "\n\n")

## ------------------------------- ##
## Prepare covariate and similarity
## ------------------------------- ##

x_df <- data.frame(x1 = x[, 1])

similarity_fun <- similarity_ppmx_gaussian_mean
sim_args       <- list(list(m0 = 0, s0 = sqrt(1)))
alpha_vec      <- alpha_in

## ------------------------------- ##
## Run ESBM
## ------------------------------- ##
cat("Starting ESBM...\n")

Z_post <- esbm(
  Y      = Y,
  seed   = seed,
  N_iter = N_iter,
  prior  = list(
    name    = "DM",
    beta_DM = beta_dm,
    H_DM    = H_dm
  ),
  z_init = sample(1:4, nrow(Y), replace = TRUE),
  a = a,
  b = b,
  x = x_df,
  similarity_fun = similarity_fun,
  sim_args = sim_args,
  alpha_g = alpha_vec
)


cat("Done.\n\n")

## ------------------------------- ##
## Save output
## ------------------------------- ##
out_dir <- here("simulation", "results", "oneCov")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

alpha_tag  <- gsub("\\.", "p", sprintf("%.2f", alpha_in))
base_name  <- file_path_sans_ext(basename(data_path))

out_file <- file.path(
  out_dir,
  sprintf("postOneCov_%s_alpha-%s_seed-%d.rds",
          base_name, alpha_tag, seed)
)

saveRDS(
  list(
    Z_post   = Z_post,
    file     = data_path,
    alpha    = alpha_in,
    seed     = seed,
    scenario = scenario,
    N_iter   = N_iter
  ),
  file = out_file
)

cat("Saved posterior to:\n  ", out_file, "\n")