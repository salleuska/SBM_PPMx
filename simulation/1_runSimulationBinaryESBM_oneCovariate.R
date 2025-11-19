SuppressPackageStartupMessages({
  library(Rcpp)
  library(here)
  library(igraph)
  library(tools)
})

args <- R.utils::commandArgs(asValue = TRUE)

## -------------------------------
## Parse args
## -------------------------------

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

## -------------------------------
## Load code
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
## Load simulation
## -------------------------------
sim_obj <- readRDS(data_path)

Y  <- sim_obj$Y
x  <- sim_obj$x              # matrix n x 1 or n x 2 but we use only first
z0 <- sim_obj$partition

## -------------------------------
## Infer scenario (none / one)
## -------------------------------
fname <- basename(data_path)

if (grepl("none", fname, ignore.case = TRUE)) {
  scenario <- "none"
} else if (grepl("one", fname, ignore.case = TRUE)) {
  scenario <- "one"
} else {
  stop("Filename must contain 'none' or 'one'. Found: ", fname)
}

cat("Scenario:", scenario, "\n\n")

## -------------------------------
## Prepare covariate and similarity
## -------------------------------
if (scenario == "none") {

  x_df           <- NULL
  similarity_fun <- NULL
  sim_args       <- list()
  alpha_vec      <- 1

} else {  # scenario == "one"

  x_df <- data.frame(x1 = x[, 1])

  similarity_fun <- similarity_ppmx_gaussian_mean
  sim_args       <- list(list(m0 = 0, s0 = sqrt(1)))
  alpha_vec      <- alpha_in
}

## -------------------------------
## Run ESBM
## -------------------------------
cat("Starting ESBM...\n")

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
  alpha_g        = alpha_vec
)

cat("Done.\n\n")

## -------------------------------
## Save output
## -------------------------------
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