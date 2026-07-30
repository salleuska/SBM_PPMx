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

# data_path  
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

## Similarity calibration: "raw" leaves log g() unstandardized (default here).
## "normalized" subtracts a log-sum-exp from log g() at every node update.
similarity_calibration <- if (is.null(args$similarity_calibration)) {
  "raw"
} else {
  as.character(args$similarity_calibration)
}

cat("  calibration =", similarity_calibration, "\n")

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
## Expected: binarySBM_2cov_<CODE>.rds, one letter per covariate
##   N = neutral, I = informative, M = misleading
## ------------------------------- ##
fname <- file_path_sans_ext(basename(data_path))

scenario <- sub("^binarySBM_2cov_([A-Z]{2})$", "\\1", fname)

valid_scenarios <- c("NN", "IN", "II", "MN", "MI", "MM")

if (scenario == fname || !scenario %in% valid_scenarios) {
  stop("Cannot infer two-cov scenario from file name: ", basename(data_path),
       "\nExpected binarySBM_2cov_<CODE>.rds with CODE in: ",
       paste(valid_scenarios, collapse = " / "), ".")
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

mcmcpost <- esbm(
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
  alpha_g        = alpha_g_vec,
  similarity_calibration = similarity_calibration
)

cat("Sampler finished.\n\n")

## ------------------------------- ##
## Save results
## ------------------------------- ##
out_dir <- here("simulation", "fixed_alpha", "results", "twoCov")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

alpha1_tag <- gsub("\\.", "p", sprintf("%.2f", alpha1_in))
alpha2_tag <- gsub("\\.", "p", sprintf("%.2f", alpha2_in))
base_name <- file_path_sans_ext(basename(data_path))

out_file <- file.path(
  out_dir,
  sprintf("postTwoCov_%s_cal-%s_a1-%s_a2-%s_seed-%d.rds",
          base_name, similarity_calibration, alpha1_tag, alpha2_tag, seed)
)


saveRDS(
  list(
    mcmcpost = mcmcpost,
    file    = data_path,
    alpha   = alpha_g_vec,
    N_iter  = N_iter,
    scenario= scenario,
    nCov    = 2L,
    seed = seed,
    similarity_calibration = similarity_calibration
  ),
  file = out_file
)

cat("Saved posterior to:\n  ", out_file, "\n")