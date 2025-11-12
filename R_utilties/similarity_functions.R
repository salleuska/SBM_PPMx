################################################################################
# similarity_functions.R
#
# This file defines similarity (cohesion) functions g(X_h*) used in the PPMx
# component of the Gibbs sampler for the extended stochastic block model (ESBM).
#
# Each similarity function must:
#   1. Accept the following arguments:
#        - Z_minus_v     : (V-1) × H binary membership matrix excluding node v
#        - cluster_sizes : numeric vector of length H with current cluster sizes
#        - v_index       : index of the focal unit v (integer scalar)
#        - H             : number of existing clusters
#        - x             : vector or matrix of covariates (user-provided format)
#        - args          : list of hyperparameters specific to the similarity
#   2. Return a numeric vector of log-weights of length (H + 1):
#        one log-weight for each existing cluster and one for a potential new one.
#   3. Perform basic checks (type or coavariates, length, etc.) and avoid side effects.
#
# The returned log-weights are combined with the urn (partition) term and the
# network likelihood to sample a new cluster assignment for node v.
################################################################################


# Similarity g() for PPMx with Gaussian auxiliary:
#   x_i | mu_h ~ N(mu_h, 1),   mu_h ~ N(m0, s0^2)
# Returns log-weights for H existing clusters + 1 "new" cluster.
similarity_ppmx_gaussian_mean <- function(Z_minus_v, cluster_sizes, v_index, H, x, args) {
  # args: list(m0 = ..., s0 = ...)  with s0 > 0
  m0 <- args$m0
  s0 <- args$s0
  if (!is.numeric(m0) || !is.numeric(s0) || s0 <= 0) stop("args must contain m0 (numeric) and s0 > 0.")
  if (!is.numeric(x)) stop("x must be a numeric vector.")

  # Sufficient stats over clusters (excluding v):
  x_minus <- x[-v_index]                          # length V-1
  n_h     <- as.numeric(cluster_sizes)            # length H
  s1_h    <- as.numeric(crossprod(Z_minus_v, x_minus))     # sum x in cluster h
  # We don't actually need SS_h for the predictive; only n_h and bar x_h.
  bar_h   <- ifelse(n_h > 0, s1_h / n_h, 0)

  # Posterior variance and mean of mu_h given x_{-v,h}
  inv_s02 <- 1 / (s0 * s0)
  s2_h    <- 1 / (inv_s02 + n_h)                 # s_h^2
  m_h     <- s2_h * (inv_s02 * m0 + n_h * bar_h)

  # Predictive variance for x_v in existing clusters: 1 + s_h^2
  var_old <- 1 + s2_h
  # Prior predictive for new cluster: N(m0, 1 + s0^2)
  var_new <- 1 + s0 * s0

  # Log-densities
  xv <- x[v_index]
  log_old <- -0.5 * (log(2*pi*var_old) + (xv - m_h)^2 / var_old)  # length H
  log_new <- -0.5 * (log(2*pi*var_new) + (xv - m0)^2 / var_new)   # scalar

  c(log_old, log_new)  # length H + 1
}

similarity_dirichlet <- function(Z_minus_v, cluster_sizes, v_index, H, 
                                x, args) {
  # x: factor with C categories
  stopifnot(is.factor(x))
  
  params <- as.numeric(args$params)
  if (length(params) != nlevels(x))
    stop("args$params must have length equal to nlevels(x).")
  params_sum <- sum(params)

  # One-hot representation computed on the fly (fast for one covariate)
  X <- stats::model.matrix(~ x - 1)               # V x C
  Vx <- crossprod(Z_minus_v, X[-v_index, , drop = FALSE])  # H x C
  cat_idx <- as.integer(x[v_index])               # 1..C

  add_old <- (Vx[, cat_idx] + params[cat_idx]) / (cluster_sizes + params_sum)
  add_new <-  params[cat_idx] / params_sum

  log(c(add_old, add_new))                        # numeric length H+1, log-scale
}

