# ------------------------------------------------------------------#
# This script generate synthetic data from a binary SBM 
#
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
#         - "none": covariates independent of cluster structure.
#
# Returns:
#       Y               adjacency matrix
#       x               covariates
#       partition       true memberships
#       block_probs     block-probability matrix
#       clust_sizes     cluster sizes
#
# ------------------------------------------------------------------ #
library(MASS)
library(here)

## simulate a binary network from the SBM
simBinarySBM <- function(
  clust_sizes = c(50, 40, 30, 30),   # cluster sizes
  p_within    = 0.8,                 # within-cluster edge prob
  p_between   = 0.1                  # between-cluster edge prob
) {

  n <- sum(clust_sizes)                 # sample size
  K <- length(clust_sizes)              # number of clusters
  z <- rep(1:K, times = clust_sizes)    # cluster memberships

  # Block matrix Theta
  Theta <- matrix(p_between, K, K)
  diag(Theta) <- p_within

  # --- SAMPLE ADJACENCY MATRIX Y ---
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


add_covariates_SBMsim <- function(base_obj, cov_dep = c("none", "one", "both")) {
  cov_dep <- match.arg(cov_dep)

  Y <- base_obj$Y
  z <- base_obj$partition

  n <- length(z)
  K <- length(unique(z))
  if (K != 4L) stop("Current covariate design assumes exactly 4 clusters.")

  x <- matrix(NA_real_, n, 2)

  if (cov_dep == "both") {

    mx <- matrix(c( 1,  1,
                    1, -1,
                   -1,  1,
                   -1, -1),
                 K, 2, byrow = TRUE)

    s  <- 0.5
    Sx <- s * diag(2)

    for (i in 1:n) {
      x[i, ] <- MASS::mvrnorm(1, mx[z[i], ], Sx)
    }

  } else if (cov_dep == "one") {

    mx1 <- c(-1.5, -0.5, +0.5, +1.5)
    sx1 <- 0.5

    mx2 <- 0
    sx2 <- 1

    for (i in 1:n) {
      x[i, ] <- c(
        rnorm(1, mx1[z[i]], sx1),
        rnorm(1, mx2, sx2)
      )
    }

  } else { # none

    x <- MASS::mvrnorm(n, mu = c(0,0), Sigma = diag(2))
  }

  base_obj$x       <- x
  base_obj$cov_dep <- cov_dep
  base_obj
}

add_covariates_SBMsim <- function(base_obj, 
  cov_dep = c("none", "one", "both")) {
  cov_dep <- match.arg(cov_dep)

  Y <- base_obj$Y
  z <- base_obj$partition

  n <- length(z)
  K <- length(unique(z))
  if (K != 4L) stop("Current covariate design assumes exactly 4 clusters.")

  x <- matrix(NA_real_, n, 2)

  if (cov_dep == "both") {

    mx <- matrix(c( 1,  1,
                    1, -1,
                   -1,  1,
                   -1, -1),
                 K, 2, byrow = TRUE)

    s  <- 0.5
    Sx <- s * diag(2)

    for (i in 1:n) {
      x[i, ] <- MASS::mvrnorm(1, mx[z[i], ], Sx)
    }

  } else if (cov_dep == "one") {

    mx1 <- c(-1.5, -0.5, +0.5, +1.5)
    sx1 <- 0.5

    mx2 <- 0
    sx2 <- 1

    for (i in 1:n) {
      x[i, ] <- c(
        rnorm(1, mx1[z[i]], sx1),
        rnorm(1, mx2, sx2)
      )
    }

  } else { # none

    x <- MASS::mvrnorm(n, mu = c(0,0), Sigma = diag(2))
  }

  base_obj$x       <- x
  base_obj$cov_dep <- cov_dep
  base_obj
}

#######################
## simulations
#######################

set.seed(1116)

# Generate one network
base <- simBinarySBM()

# Generate three covariate scenarios using same Y and same z
sim_none  <- add_covariates_SBMsim(base, "none")
sim_one   <- add_covariates_SBMsim(base, "one")
sim_both  <- add_covariates_SBMsim(base, "both")

## create directory if it does not exists
dir.create(here("simulation/data"), 
  recursive = TRUE, showWarnings = FALSE)

saveRDS(sim_none, here("simulation/data/binarySBM_none.rds"))
saveRDS(sim_one, here("simulation/data/binarySBM_one.rds"))
saveRDS(sim_both, here("simulation/data/binarySBM_both.rds"))
#######################
## Plotting
#######################
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
  pal  <- grDevices::rainbow(K)
  vcol <- pal[z]

  if (is.null(main)) {
    main <- paste0("Network with geographic layout (cov_dep = '",
                   sim_obj$cov_dep, "')")
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

## plot_network_geo(sim_none)
#######################################################
#### BACKUP: one function
# simBinarySBMx <- function(
#   clust_sizes = c(50, 40, 30, 30),        # cluster sizes
#   p_within    = 0.8,                      # within-cluster edge prob
#   p_between   = 0.1,                      # between-cluster edge prob
#   cov_dep     = c("none", "one", "both")  # covariate–cluster dependence
# ) {
#   cov_dep <- match.arg(cov_dep)

#   n <- sum(clust_sizes)       # sample size
#   K <- length(clust_sizes)    # number of clusters
#   z <- rep(1:K, times = clust_sizes)   # cluster memberships

#   if (K != 4L) stop("Current covariate design assumes exactly 4 clusters.")

#   # Block matrix Theta
#   Theta <- matrix(p_between, K, K)
#   diag(Theta) <- p_within

#   # --- SAMPLE ADJACENCY MATRIX Y ---
#   Y <- matrix(0L, n, n)

#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       pij      <- Theta[z[i], z[j]]
#       Y_ij     <- rbinom(1, 1, pij)
#       Y[i, j]  <- Y_ij
#       Y[j, i]  <- Y_ij
#     }
#   }

#   # --- SAMPLE COVARIATES x (depends on cov_dep) ---
#   x <- matrix(NA_real_, n, 2)

#   if (cov_dep == "both") {

#     mx <- matrix(c( 1,  1,
#                     1, -1,
#                    -1,  1,
#                    -1, -1),
#                  K, 2, byrow = TRUE)

#     s  <- 0.5
#     Sx <- s * diag(2)

#     for (i in 1:n) {
#       x[i, ] <- mvrnorm(1, mx[z[i], ], Sx)
#     }

#   } else if (cov_dep == "one") {

#     mx1 <- c(-1.5, -0.5, +0.5, +1.5)
#     sx1 <- 0.5

#     mx2 <- 0
#     sx2 <- 1

#     for (i in 1:n) {
#       x[i, ] <- c(
#         rnorm(1, mx1[z[i]], sx1),
#         rnorm(1, mx2, sx2)
#       )
#     }

#   } else {  # cov_dep == "none"

#     mx <- c(0, 0)
#     Sx <- diag(2)

#     x <- mvrnorm(n, mx, Sx)
#   }

#   list(
#     Y           = Y,
#     x           = x,
#     partition   = z,
#     block_probs = Theta,
#     clust_sizes = clust_sizes,
#     cov_dep     = cov_dep
#   )
# }