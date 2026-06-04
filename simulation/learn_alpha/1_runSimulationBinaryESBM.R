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

alpha_init <- if (is.null(args$alpha_init)) 2 else as.numeric(args$alpha_init)
a_alpha    <- if (is.null(args$a_alpha)) 2 else as.numeric(args$a_alpha)
b_alpha    <- if (is.null(args$b_alpha)) 1 else as.numeric(args$b_alpha)

seed     <- if (is.null(args$seed)) 123 else as.integer(args$seed)
N_iter   <- if (is.null(args$N_iter)) 10000 else as.integer(args$N_iter)

similarity_calibration <- if (is.null(args$similarity_calibration)) {
  "normalized"
} else {
  as.character(args$similarity_calibration)
}

model <- if (is.null(args$model)) {"ESBM"} else {as.character(args$model)}

cat("Running the following:\n")
cat("  model = ", model, "\n")
cat("  file       =", data_path, "\n")
cat("  alpha init =", alpha_init, "\n")
cat("  a_alpha    =", a_alpha, "\n")
cat("  b_alpha    =", b_alpha, "\n")
cat("  calibration =", similarity_calibration, "\n")
cat("  seed       =", seed, "\n")
cat("  N_iter     =", N_iter, "\n\n")

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
x <- sim_obj[["x"]]            # matrix n x 2; use only first column for one-cov runs
z0 <- sim_obj$partition

## ------------------------------- ##
## Infer scenario from file name
## ------------------------------- ##
fname <- tools::file_path_sans_ext(basename(data_path))

scenario <- sub(
  "^binarySBM_[0-9]+cov_([A-Z]+)$",
  "\\1",
  fname
)

if (scenario == fname) {
  stop("Could not infer scenario from filename: ", fname)
}

nCov <- as.integer(sub(
  "^binarySBM_([0-9]+)cov_[A-Z]+$",
  "\\1",
  fname
))

cat("Scenario:", scenario, "\n")
cat("nCov:", nCov, "\n\n")


## ------------------------------- ##
## Prepare covariates and similarity
## ------------------------------- ##

if (model == "SBM") {

  x_df <- NULL
  similarity_fun <- NULL
  sim_args <- NULL
  alpha_spec <- NULL

} else {

  x_df <- as.data.frame(x)
  colnames(x_df) <- paste0("x", seq_len(nCov))

  similarity_fun <- replicate(
    nCov,
    similarity_ppmx_gaussian_mean,
    simplify = FALSE
  )

  sim_args <- replicate(
    nCov,
    list(m0 = 0, s0 = sqrt(1)),
    simplify = FALSE
  )

  alpha_spec <- list(
    init    = alpha_init,
    a_alpha = a_alpha,
    b_alpha = b_alpha
  )
}


## ------------------------------- ##
## Run ESBM
## ------------------------------- ##
cat("Starting sampling...\n")

if (model == "SBM") {

  mcmcpost <- esbm(
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
    x = NULL
  )

} else {

  mcmcpost <- esbm(
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
    alpha_g = alpha_spec,
    similarity_calibration = similarity_calibration
  )

}
cat("Done.\n\n")

## ------------------------------- ##
## Save output
## ------------------------------- ##
base_name <- tools::file_path_sans_ext(basename(data_path))

model_tag <- ifelse(model == "SBM", "SBM", similarity_calibration)

out_dir <- here(
  "simulation",
  "learn_alpha",
  model_tag,
  "results",
  paste0(nCov, "Cov")
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (model == "SBM") {
  out_name <- sprintf(
    "postSBM_%s_seed-%d.rds",
    base_name,
    seed
  )
} else {
  alpha_tag <- sprintf(
    "cal-%s_init-%s_a-%s_b-%s",
    similarity_calibration,
    gsub("\\.", "p", sprintf("%.2f", alpha_init)),
    gsub("\\.", "p", sprintf("%.2f", a_alpha)),
    gsub("\\.", "p", sprintf("%.2f", b_alpha))
  )

  out_name <- sprintf(
    "post_%s_alphaLearn-%s_seed-%d.rds",
    base_name,
    alpha_tag,
    seed
  )
}

out_file <- file.path(out_dir, out_name)

saveRDS(
  list(
    mcmcpost = mcmcpost,
    file     = data_path,
    alpha_g  = alpha_spec,
    seed     = seed,
    scenario = scenario,
    nCov     = nCov,
    N_iter   = N_iter,
    model    = model
  ),
  file = out_file
)

cat("Saved posterior to:\n  ", out_file, "\n")

