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
# 2. Gaussian auxiliary similarity with random mean and fixed variance
#
#    Model:
#        x_i | μ_h ~ N(μ_h, 1)
#        μ_h ~ N(m0, s0^2)
#
#    Integrated predictive:
#        Existing cluster h:  N(x_v | m_h, 1 + s_h^2)
#          s_h^2 = (1/s0^2 + n_h)^(-1)
#          m_h   = s_h^2 * (m0/s0^2 + n_h * x̄_h)
#        New cluster:          N(x_v | m0, 1 + s0^2)
#
#    Inputs:
#        args = list(m0 = ..., s0 = ...)  with s0 > 0
################################################################################
similarity_ppmx_gaussian_mean <- function(Z_minus_v, cluster_sizes, v_index, H, x, args) {
  
#   browser()
  # --- hyperparameters
  m0   <- args$m0
  var0 <- args$s0^2          # prior variance
  if (!is.numeric(m0) || !is.numeric(var0) || var0 <= 0)
    stop("args must contain m0 (numeric) and s0 > 0.")
  if (!is.numeric(x))
    stop("x must be a numeric vector.")

  # --- cluster sufficient statistics (excluding v)
  x_minus <- x[-v_index]
  n_h  <- as.numeric(cluster_sizes)                 # cluster sizes
  s1_h <- as.numeric(crossprod(Z_minus_v, x_minus)) # sum of x in each cluster
  bar_h <- ifelse(n_h > 0, s1_h / n_h, 0)

  # --- posterior mean and variance of μ_h | x_{S_h^{(-v)}}
  inv_var0 <- 1 / var0
  var_h <- 1 / (inv_var0 + n_h)
  m_h   <- var_h * (inv_var0 * m0 + n_h * bar_h)

  # --- predictive variances: τ_h^2 = 1 + var_h
  var_pred_old <- 1 + var_h    # existing clusters
  var_pred_new <- 1 + var0     # new cluster

  # --- log predictive densities (up to additive constant)
  xv <- x[v_index]
  log_old <- -0.5 * (log(2 * pi * var_pred_old) + (xv - m_h)^2 / var_pred_old)
  log_new <- -0.5 * (log(2 * pi * var_pred_new) + (xv - m0)^2 / var_pred_new)

  c(log_old, log_new)          # vector of length H + 1
}


similarity_dirichlet <- function(Z_minus_v, cluster_sizes, v_index, H, 
                                x, args) {
  # --- check inputs
  # x: factor with C categories
  stopifnot(is.factor(x))
  
  params <- as.numeric(args$params)
  if (length(params) != nlevels(x))
    stop("args$params must have length equal to nlevels(x).")
  params_sum <- sum(params)

  # --- useful quantities
  # One-hot representation computed on the fly (fast for one covariate)
  X <- stats::model.matrix(~ x - 1)               # V x C
  Vx <- crossprod(Z_minus_v, X[-v_index, , drop = FALSE])  # H x C
  cat_idx <- as.integer(x[v_index])               # 1..C

  # --- log densities
  log_old <- (Vx[, cat_idx] + params[cat_idx]) / (cluster_sizes + params_sum)
  log_new <-  params[cat_idx] / params_sum

  log(c(log_old, log_new))                        # numeric length H+1, log-scale
}

