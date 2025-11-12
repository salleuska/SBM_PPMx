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
source(here("R_utilties", "esbm_original.R"))
source(here("R_utilties", "esbm.R"))
Rcpp::sourceCpp(here("R_utilties", "stirling.cpp"))

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

colnames(Z_DM)

# DM
Z_DM2 <- esbm_original(Y, my_seed, N_iter, "DM", my_z, a = a, b = b,
             beta_DM = beta_dm, H_DM = H_dm)
head(Z_DM) == head(Z_DM2)

## ---- Collapsed Gibbs sampler (supervised with node attributes) ----
# Here we use the *true labels* as attributes for the supervised version (as in paper)

## SP - sanity check -- should get the same results
N_iter <- 2000
my_seed <- 1
my_z <- seq_len(V)
# define the vector with node attributes
my_x <- c(rep(1,20),rep(2,20),rep(3,15),rep(4,15),rep(5,10))
params <- rep(1, length(unique(my_x)))

Z_DM_x <- esbm(Y, my_seed, N_iter, "DM", my_z, a = a, b = b,
               beta_DM = beta_dm, H_DM = H_dm, 
               x = as.factor(my_x), similarity_fun = similarity_dirichlet, 
               sim_args = list(params = params)
              )

# DP + x
my_alpha_xi <- rep(1, length(unique(my_x)))

Z_DP_2 <- esbm_original(Y, my_seed, N_iter,"DM", my_z, a = a, b = b,
               beta_DM = beta_dm, H_DM = H_dm, 
               x = my_x,  alpha_xi = my_alpha_xi)

tail(Z_DM_x) == tail(Z_DP_2)
