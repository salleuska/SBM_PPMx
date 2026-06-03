#------------------------------------------------------------------#
# This script generate synthetic data from a binary SBM 
#------------------------------------------------------------------#

library(MASS)
library(here)
library(igraph)
library(ggplot2)
library(RColorBrewer)

#------------------------------------------------------------------#
# Part 1: simBinarySBM simulates a network Y of size n x n from 
# stochastic block model with 4 clusters
#         - Clusters have sizes ns = (n1, n2, n3, n4).
#         - Edge probability within clusters = p_within.
#         - Edge probability between clusters = p_between.
#
# Part 2: Covariate generator (the network is fixed)
#   (2) A covariate matrix x (n x 2), with dependence on clusters
#       controlled by cov_dep:
#         - "both": both covariates informative (cluster-specific means).
#         - "one" : only first covariate informative.
#         - "neutral": covariates independent of cluster structure.
#
# Returns:
#       Y               adjacency matrix
#       x               covariates
#       partition       true memberships
#       block_probs     block-probability matrix
#       clust_sizes     cluster sizes
#
#------------------------------------------------------------------#
# simulate a binary network from the SBM
simBinarySBM <- function(
  clust_sizes = c(50, 40, 30, 30),   # cluster sizes
  p_within    = 0.6,                 # within-cluster edge prob
  p_between   = 0.3                  # between-cluster edge prob
) {

  n <- sum(clust_sizes)                 # sample size
  K <- length(clust_sizes)              # number of clusters
  z <- rep(1:K, times = clust_sizes)    # cluster memberships

  # Block matrix Theta
  Theta <- matrix(p_between, K, K)
  diag(Theta) <- p_within

  #--- SAMPLE ADJACENCY MATRIX Y ---
  # (this step does NOT vary with cov_dep)
  # NOTE: you can even use the same Y for all the 3 scenarios
  Y <- matrix(0L, n, n)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      Y[i, j] <- Y[j, i] <- rbinom(1, 1, Theta[z[i],z[j]])
    }
  }
  
  # Return network and information
  list(
    Y             = Y,
    block_probs   = Theta,
    clust_sizes   = clust_sizes,
    partition     = z
  )
}

add_covariates_SBMsim <- function(
  baseObj,
  covDep = c("informative"),
  seed = NULL
) {

  if (!is.null(seed)) set.seed(seed)

  z <- as.integer(factor(baseObj$partition))
  n <- length(z)
  K <- length(unique(z))
  nCov <- length(covDep)

  covDep <- match.arg(
    covDep,
    choices = c("informative", "neutral", "mislead_random"),
    several.ok = TRUE
  )

  mx <- c(-1.5, -0.5, 0.5, 1.5)
  sx <- 0.5

  x <- matrix(NA_real_, nrow = n, ncol = nCov)

  for (j in seq_len(nCov)) {

    if (covDep[j] == "neutral") {

      x[, j] <- rnorm(n)

    } else {

      lab <- if (covDep[j] == "informative") {
        z
      } else {
        sample(rep(seq_len(K), length.out = n))
      }

      x[, j] <- rnorm(n, mean = mx[lab], sd = sx)
    }
  }

  colnames(x) <- paste0("x", seq_len(nCov))

  baseObj$x      <- x
  baseObj$nCov   <- nCov
  baseObj$covDep <- covDep

  baseObj
}

#------------------------------------------------------------------#
# Plotting functions
#------------------------------------------------------------------#
plot_network_geo <- function(sim_obj, main = NULL) {
  Y <- sim_obj$Y
  x <- sim_obj$x
  z <- sim_obj$partition

  g <- igraph::graph_from_adjacency_matrix(
    Y,
    mode = "undirected",
    diag = FALSE
  )

  layout_xy <- as.matrix(x)   # use (x1, x2) as coordinates

  K    <- length(unique(z))
  pal  <- c(1, 2, 3,4)
  vcol <- pal[z]

  if (is.null(main)) {
    main <- if (sim_obj$nCov == 2) {
      paste0("Network (x1: ", sim_obj$cov1Dep,
             ", x2: ", sim_obj$cov2Dep, ")")
    } else {
      paste0("Network (x1: ", sim_obj$cov1Dep, ")")
    }
  }

  plot(
    g,
    layout         = layout_xy,
    vertex.color   = vcol,
    vertex.size    = 6,
    vertex.label   = NA,
    edge.arrow.mode = 0,
    main           = main
  )
}


plot_adj_matrix_gg <- function(sim_obj) {
  Y  <- sim_obj$Y
  z0 <- sim_obj$partition
  V  <- length(z0)
  
  # remove self-loops
  diag(Y) <- 0
  
  #------------------------------------------#
  # Order nodes by true cluster to show blocks
  #------------------------------------------#
  ord   <- order(z0)
  Y_ord <- Y[ord, ord]
  z_ord <- z0[ord]
  
  #------------------------------------------#
  # Build long data.frame in base R
  #------------------------------------------#
  # indices
  grid_idx <- expand.grid(
    row = seq_len(V),
    col = seq_len(V)
  )
  
  # adjacency values (row-major)
  grid_idx$value <- as.vector(Y_ord)
  
  # factors for plotting:
  # flip row so cluster 1 is at bottom (heatmap style)
  grid_idx$row_f <- factor(grid_idx$row,
                           levels = rev(seq_len(V)))
  grid_idx$col_f <- factor(grid_idx$col,
                           levels = seq_len(V))
  
  # cluster label per row (after ordering)
  grid_idx$cluster <- factor(z_ord[grid_idx$row])
  
  #------------------------------------------#
  # Colors
  #------------------------------------------#
  adj_pal <- colorRampPalette(
    brewer.pal(9, "Greys")[c(1, 9)]
  )(50)
  
  #------------------------------------------#
  # Plot
  #------------------------------------------#
  ggplot(grid_idx, aes(x = col_f, y = row_f, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = adj_pal, guide = "none") +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      axis.title   = element_blank(),
      panel.grid   = element_blank()
    ) +
    ggtitle("Adjacency matrix")
}



#------------------------------------------------------------------#
# Simulations
#------------------------------------------------------------------#
set.seed(1116)

# ------------------------------------------------------------
# 1) Generate one network (fixed Y and true partition z)
# ------------------------------------------------------------
base <- simBinarySBM(p_within = 0.4, p_between = 0.2)

# ------------------------------------------------------------
# 2) One-covariate scenarios
# ------------------------------------------------------------
sim_1_N <- add_covariates_SBMsim(base, covDep = "neutral")
sim_1_I <- add_covariates_SBMsim(base, covDep = "informative")
sim_1_M <- add_covariates_SBMsim(base, covDep = "mislead_random")

# ------------------------------------------------------------
# 3) Two-covariate scenarios
#    Skip symmetric cases.
# ------------------------------------------------------------
sim_2_NN <- add_covariates_SBMsim(base, covDep = c("neutral",        "neutral"))
sim_2_IN <- add_covariates_SBMsim(base, covDep = c("informative",    "neutral"))
sim_2_II <- add_covariates_SBMsim(base, covDep = c("informative",    "informative"))
sim_2_MN <- add_covariates_SBMsim(base, covDep = c("mislead_random", "neutral"))
sim_2_MI <- add_covariates_SBMsim(base, covDep = c("mislead_random", "informative"))
sim_2_MM <- add_covariates_SBMsim(base, covDep = c("mislead_random", "mislead_random"))

# ------------------------------------------------------------
# 4) More-covariate stress-test scenarios
# ------------------------------------------------------------
sim_5_INNNN <- add_covariates_SBMsim(
  base,
  covDep = c("informative", rep("neutral", 4))
)

sim_5_IMMMM <- add_covariates_SBMsim(
  base,
  covDep = c("informative", rep("mislead_random", 4))
)

sim_5_MIIII <- add_covariates_SBMsim(
  base,
  covDep = c("mislead_random", rep("informative", 4))
)

# ------------------------------------------------------------
# 5) Save
# ------------------------------------------------------------
dir.create(here("simulation/data"), recursive = TRUE, showWarnings = FALSE)

saveRDS(sim_1_N, here("simulation/data/binarySBM_1cov_N.rds"))
saveRDS(sim_1_I, here("simulation/data/binarySBM_1cov_I.rds"))
saveRDS(sim_1_M, here("simulation/data/binarySBM_1cov_M.rds"))

saveRDS(sim_2_NN, here("simulation/data/binarySBM_2cov_NN.rds"))
saveRDS(sim_2_IN, here("simulation/data/binarySBM_2cov_IN.rds"))
saveRDS(sim_2_II, here("simulation/data/binarySBM_2cov_II.rds"))
saveRDS(sim_2_MN, here("simulation/data/binarySBM_2cov_MN.rds"))
saveRDS(sim_2_MI, here("simulation/data/binarySBM_2cov_MI.rds"))
saveRDS(sim_2_MM, here("simulation/data/binarySBM_2cov_MM.rds"))

saveRDS(sim_5_INNNN, here("simulation/data/binarySBM_5cov_INNNN.rds"))
saveRDS(sim_5_IMMMM, here("simulation/data/binarySBM_5cov_IMMMM.rds"))
saveRDS(sim_5_MIIII, here("simulation/data/binarySBM_5cov_MIIII.rds"))

# ------------------------------------------------------------
# 6) Quick plots
# ------------------------------------------------------------
plot_adj_matrix_gg(sim_1_N)

plot_network_geo(sim_1_N)
plot_network_geo(sim_1_I)
plot_network_geo(sim_1_M)

plot_network_geo(sim_2_NN)
plot_network_geo(sim_2_IN)
plot_network_geo(sim_2_II)
plot_network_geo(sim_2_MN)
plot_network_geo(sim_2_MI)
plot_network_geo(sim_2_MM)

plot_network_geo(sim_5_INNNN)
plot_network_geo(sim_5_IMMMM)
plot_network_geo(sim_5_MIIII)

