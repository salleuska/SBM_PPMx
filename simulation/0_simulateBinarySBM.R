# ------------------------------------------------------------------#
# This script generate synthetic data from a binary SBM 
#
# The function generates:
#   (1) A network Y of size n x n from a 4-cluster stochastic block model.
#         - Clusters have sizes ns = (n1, n2, n3, n4).
#         - Edge probability within clusters = p_within.
#         - Edge probability between clusters = p_between.
#
#   (2) A covariate matrix x (n x 2), with dependence on clusters
#       controlled by cov_dep:
#         - "both": both covariates informative (cluster-specific means).
#         - "one" : only first covariate informative.
#         - "none": covariates independent of cluster structure.
#
# Returns:
#       Y       adjacency matrix
#       x       covariates
#       z       true memberships
#       Theta   block-probability matrix
#       ns      cluster sizes
#
# ------------------------------------------------------------------ #

library(MASS)

simBinarySBMx <- function(
  ns         = c(50, 40, 30, 30),   # cluster sizes
  p_within   = 0.8,                 # within-cluster edge prob
  p_between  = 0.1,                 # between-cluster edge prob
  cov_dep    = c("none", "one", "both")  # covariate–cluster dependence
) {
  cov_dep <- match.arg(cov_dep)

  n <- sum(ns)                 # sample size
  K <- length(ns)              # number of clusters
  z <- rep(1:K, times = ns)    # cluster memberships

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

  # --- SAMPLE COVARIATES x (depends on cov_dep) ---
  x <- matrix(NA, n, 2)

  if (cov_dep == "both") {

    mx <- matrix(c( 1,  1,
                    1, -1,
                   -1,  1,
                   -1, -1),
                 K, 2, byrow = TRUE)

    s  <- 0.5
    Sx <- s * diag(2)

    for (i in 1:n) {
      x[i, ] <- mvrnorm(1, mx[z[i], ], Sx)
    }

  } else if (cov_dep == "one") {

    ## NOTE: this also currently assumes K = 4
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

  } else {  # cov_dep == "none"

    mx <- c(0, 0)
    Sx <- diag(2)

    x <- mvrnorm(n, mx, Sx)
  }

  # Return what you'll likely need downstream
  list(
    Y     = Y,
    x     = x,
    z     = z,
    Theta = Theta,
    ns    = ns
  )
}


############
plot_covariates <- function(sim_obj, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting.")
  }

  x <- sim_obj$x
  z <- sim_obj$z

  df <- data.frame(
    x1 = x[, 1],
    x2 = x[, 2],
    z  = factor(z)
  )

  if (is.null(title)) {
    title <- paste0("Covariates vs clusters (cov_dep = '",
                    sim_obj$cov_dep, "')")
  }

  ggplot2::ggplot(df, ggplot2::aes(x = x1, y = x2, color = z)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "x1",
      y = "x2",
      color = "Cluster",
      title = title
    )
}

set.seed(123)

sim_none  <- my_sim(cov_dep = "none")
sim_one   <- my_sim(cov_dep = "one")
sim_both  <- my_sim(cov_dep = "both")

plot_covariates(sim_none)
plot_covariates(sim_one)
plot_covariates(sim_both)