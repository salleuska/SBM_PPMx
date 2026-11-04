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
## a_alpha and b_alpha accept one value per covariate, comma separated
## (e.g. b_alpha=1,75); a single value is reused for all covariates.
a_alpha    <- if (is.null(args$a_alpha)) 2 else as.numeric(strsplit(args$a_alpha, ",")[[1]])
b_alpha    <- if (is.null(args$b_alpha)) 1 else as.numeric(strsplit(args$b_alpha, ",")[[1]])

## alpha_prior:
##   "gamma"      -> Gamma(a_alpha, b_alpha) as given
##   "recentered" -> a_j = c, b_j = c + log G0_j, with log G0_j the log
##                   similarity of all nodes in one cluster; a_alpha and
##                   b_alpha are ignored. c_alpha defaults to n/2.
alpha_prior <- if (is.null(args$alpha_prior)) "gamma" else as.character(args$alpha_prior)
c_alpha     <- if (is.null(args$c_alpha)) NULL else as.numeric(args$c_alpha)

if (!alpha_prior %in% c("gamma", "recentered")) {
  stop("alpha_prior must be either 'gamma' or 'recentered'.")
}

seed     <- if (is.null(args$seed)) 123 else as.integer(args$seed)
N_iter   <- if (is.null(args$N_iter)) 10000 else as.integer(args$N_iter)

similarity_calibration <- if (is.null(args$similarity_calibration)) {
  "raw"
} else {
  as.character(args$similarity_calibration)
}

model <- if (is.null(args$model)) {"ESBM"} else {as.character(args$model)}

cat("Running the following:\n")
cat("  model = ", model, "\n")
cat("  file       =", data_path, "\n")
cat("  alpha init =", alpha_init, "\n")
cat("  alpha prior =", alpha_prior, "\n")
if (alpha_prior == "gamma") {
  cat("  a_alpha    =", a_alpha, "\n")
  cat("  b_alpha    =", b_alpha, "\n")
} else {
  cat("  c_alpha    =", if (is.null(c_alpha)) "n/2" else c_alpha, "\n")
}
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
  "^binarySBM_[0-9]+cov_([A-Z0-9]+)$",
  "\\1",
  fname
)

if (scenario == fname) {
  stop("Could not infer scenario from filename: ", fname)
}

nCov <- as.integer(sub(
  "^binarySBM_([0-9]+)cov_[A-Z0-9]+$",
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

  log_G0 <- NULL

  if (alpha_prior == "recentered") {

    if (similarity_calibration != "raw") {
      stop("alpha_prior=recentered requires similarity_calibration=raw.")
    }

    n <- nrow(x_df)
    if (is.null(c_alpha)) c_alpha <- n / 2

    ## log G0_j on the same standardized covariates the sampler uses
    x_std <- standardize_covariates(x_df)
    Z_one <- matrix(1, n, 1)

    log_G0 <- vapply(seq_len(nCov), function(j) {
      sum(similarity_fun[[j]](
        mode = "partition",
        Z    = Z_one,
        x    = x_std[[j]],
        args = sim_args[[j]]
      ))
    }, numeric(1))

    a_alpha <- rep(c_alpha, nCov)
    b_alpha <- c_alpha + log_G0

    cat("Recentered prior:\n")
    cat("  log G0     =", round(log_G0, 3), "\n")
    cat("  a_alpha    =", round(a_alpha, 3), "\n")
    cat("  b_alpha    =", round(b_alpha, 3), "\n\n")
  }

  alpha_spec <- list(
    init    = alpha_init,
    a_alpha = a_alpha,
    b_alpha = b_alpha,
    prior   = alpha_prior,
    c_alpha = c_alpha,
    log_G0  = log_G0
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
  num_tag <- function(v) paste(gsub("\\.", "p", sprintf("%.2f", v)), collapse = "-")

  if (alpha_prior == "recentered") {
    ## a and b are per covariate and b can be negative: tag by c instead
    alpha_tag <- sprintf(
      "cal-%s_init-%s_recentered_c-%s",
      similarity_calibration,
      num_tag(alpha_init),
      num_tag(c_alpha)
    )
  } else {
    alpha_tag <- sprintf(
      "cal-%s_init-%s_a-%s_b-%s",
      similarity_calibration,
      num_tag(alpha_init),
      num_tag(a_alpha),
      num_tag(b_alpha)
    )
  }

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

