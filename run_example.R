library(here)

## load data
dat <- readRDS(here("simulation/data/binarySBM_one.rds"))
Y   <- dat$Y
x   <- dat$x

source(here("R_utilities/similarity_functions.R"))
source(here("R_utilities/esbm.R"))

fit <- esbm(
  Y       = Y,
  seed    = 1,
  N_iter  = 2000,
  prior   = "DM",
  x       = x,
  similarity_fun = sim_normal_aux,
  sim_args = list(alpha = 1)
)

str(fit)