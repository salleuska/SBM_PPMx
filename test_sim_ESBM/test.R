## this script is used to check changes in the code

## ---- Setup ----
rm(list = ls())
suppressPackageStartupMessages({
  library(Rcpp)
  library(here)
  library(reshape)
  library(gdata)
  library(igraph)
  library(mcclust.ext)
  library(RColorBrewer)
  library(pheatmap)
  library(gridExtra)
  library(grid)
  library(cowplot)
  library(ggplot2)
  library(coda)
##   library(dummies) - not available
  library(randnet)
  library(greed)
  library(LaplacesDemon)
})


## ---- Source ----
source(here("R_utilties", "esbm.R"))
Rcpp::sourceCpsp(here("R_utilties", "stirling.cpp"))

## ---- Data ----
load(here("test_sim_ESBM", "network_1.RData"))  # expects: Y, z_0 (true labels)
V <- nrow(Y); diag(Y) <- 0

## ---- Hyperparameters (priors) ----
# Beta(a,b) for block probabilities
a <- 1; b <- 1

# Choose Gibbs-type priors to target E[H] ~ 10 for V=80
# (uses helper funcs in esbm.R: expected_cl_py, HGnedin)
sigma_dm <- 0; H_dm <- 50; beta_dm <- 3.5 / H_dm

expected_cl_py(V, sigma = sigma_dm, theta = beta_dm * H_dm, H = H_dm)

## ---- Collapsed Gibbs sampler (unsupervised) ----
N_iter <- 2000
my_seed <- 1
my_z <- seq_len(V)

# DM
Z_DM <- esbm(Y, my_seed, N_iter, "DM", my_z, a = a, b = b,
             beta_DM = beta_dm, H_DM = H_dm)

## ---- Collapsed Gibbs sampler (supervised with node attributes) ----
# Here we use the *true labels* as attributes for the supervised version (as in paper)
N_iter <- 2000
my_seed <- 1
my_z <- seq_len(V)
my_x <- as.factor(sample(1:4, V, replace = T))                        # node attribute (supervised)
params <- rep(1, length(unique(my_x)))

# DM + x
Z_DM_x <- esbm(Y, my_seed, N_iter, "DM", my_z, a = a, b = b,
               beta_DM = beta_dm, H_DM = H_dm, 
               x = my_x, similarity_fun = similarity_dirichlet, 
               sim_args = list(params = params)
              )

# DP + x
Z_DP_x <- esbm(Y, my_seed, N_iter, "PY", my_z, a = a, b = b,
               alpha_PY = alpha_dp, sigma_PY = sigma_dp, x = my_x, alpha_xi = my_alpha_xi)

