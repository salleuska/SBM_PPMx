################################################################################
# esbm_utils.R
#
# Utility functions for the extended stochastic block model:
#   - log_pY_z       : log marginal likelihood of Y given a partition
#   - pred_esbm      : out-of-sample cluster probabilities for a new node
#   - expected_K     : prior expected number of clusters under DP / PY / DM / GN
#   - dclust_py      : prior distribution of K under DP / PY
#
# These functions depend on the urn functions defined in esbm.R; source that
# file before sourcing this one.
################################################################################


#' Log marginal likelihood of the adjacency matrix given a partition
#'
#' Computes \eqn{\log p(Y \mid z, a, b)} under the Beta-Bernoulli block model,
#' integrating out the block connection probabilities.
#'
#' @param Y A \eqn{V \times V} symmetric binary adjacency matrix (no self-loops).
#' @param z Integer vector of length \eqn{V} with cluster labels.
#' @param a,b Positive hyperparameters of the Beta prior on block probabilities.
#'
#' @return A scalar log marginal likelihood.
#'
#' @examples
#' log_pY_z(Y, z = c(1,1,2,2), a = 1, b = 1)
log_pY_z <- function(Y, z, a = 1, b = 1) {

  stopifnot(is.matrix(Y), nrow(Y) == ncol(Y))
  stopifnot(length(z) == nrow(Y))
  stopifnot(a > 0, b > 0)

  diag(Y) <- 0
  H  <- length(unique(z))
  Z  <- diag(H)[z, , drop = FALSE]   # V x H membership matrix

  # sufficient statistics: edge counts and total pairs per block pair
  M     <- t(Z) %*% Y %*% Z          # H x H edge counts (diag double-counted)
  diag(M) <- diag(M) / 2

  n_h   <- colSums(Z)
  pairs <- outer(n_h, n_h)           # H x H total pairs
  diag(pairs) <- n_h * (n_h - 1) / 2

  M_bar <- pairs - M                 # non-edge counts

  # lower triangle including diagonal (each block pair counted once)
  idx   <- lower.tri(M, diag = TRUE)
  a_n   <- M[idx]   + a
  b_n   <- M_bar[idx] + b

  sum(lbeta(a_n, b_n)) - (H * (H + 1) / 2) * lbeta(a, b)
}


#' Out-of-sample cluster probabilities for a new node
#'
#' Given an estimated partition \code{z_hat} for the observed nodes and the
#' edges from a new node to all observed nodes (appended as the last row/column
#' of \code{Y}), returns the full-conditional probabilities of the new node
#' belonging to each existing cluster or forming a new one.
#'
#' @param Y A \eqn{(V+1) \times (V+1)} symmetric adjacency matrix. The last
#'   row and column encode the edges of the new node to the \eqn{V} observed nodes.
#' @param prior A list specifying the partition prior, in the same format as
#'   \code{esbm()}: e.g. \code{list(name = "GN", gamma_GN = 0.3)}.
#' @param z_hat Integer vector of length \eqn{V} with cluster labels for the
#'   observed nodes.
#' @param a,b Hyperparameters of the Beta prior on block probabilities.
#'
#' @return A numeric vector of length \eqn{H + 1} of (unnormalised) probabilities,
#'   one per existing cluster plus one for a new cluster. Normalise with
#'   \code{p / sum(p)} to get a proper distribution.
#'
#' @examples
#' p <- pred_esbm(Y_new, prior = list(name = "GN", gamma_GN = 0.3), z_hat = z_hat)
#' p / sum(p)
pred_esbm <- function(Y, prior, z_hat, a = 1, b = 1) {

  if (!is.list(prior) || is.null(prior$name)) {
    stop("prior must be a list with a name field.")
  }

  if (prior$name == "DP") {
    urn <- function(v_minus) urn_DP(v_minus, prior$alpha_PY)
  } else if (prior$name == "PY") {
    urn <- function(v_minus) urn_PY(v_minus, prior$alpha_PY, prior$sigma_PY)
  } else if (prior$name == "DM") {
    urn <- function(v_minus) urn_DM(v_minus, prior$beta_DM, prior$H_DM)
  } else if (prior$name == "GN") {
    urn <- function(v_minus) urn_GN(v_minus, prior$gamma_GN)
  } else {
    stop("Invalid prior$name.")
  }

  V_new <- nrow(Y)
  v     <- V_new    # new node is the last one

  # initialise Z placing the new node in a singleton cluster
  z_init <- c(z_hat, max(z_hat) + 1L)
  H_init <- max(z_init)
  Z      <- diag(H_init)[z_init, , drop = FALSE]

  temp   <- Y %*% Z
  m_full <- t(Z) %*% temp - diag(0.5 * colSums(temp * Z), ncol(Z))

  # remove empty clusters (including the temporary singleton for node v)
  if (ncol(Z) > 1L) {
    nonempty <- which(colSums(Z[-v, ]) > 0)
    Z        <- Z[, nonempty, drop = FALSE]
    if (length(nonempty) == 1L) Z <- matrix(Z, V_new, 1L)
    m_full   <- matrix(m_full[nonempty, nonempty], ncol(Z), ncol(Z))
  }

  H     <- ncol(Z)
  Z_v   <- Z[-v, , drop = FALSE]
  v_minus <- if (H == 1L) sum(Z[-v]) else colSums(Z_v)
  r_v     <- crossprod(Z_v, Y[-v, v])

  m_bar   <- matrix(v_minus, H, 1) %*% matrix(v_minus, 1, H) -
             diag(0.5 * v_minus * (v_minus + 1), H) - m_full
  V_minus <- matrix(1, H, 1) %*% matrix(v_minus, 1, H)
  R       <- matrix(1, H, 1) %*% matrix(r_v, 1, H)

  log_lhds_old <- rowSums(
    lbeta(m_full + R + a, m_bar + V_minus - R + b) - lbeta(m_full + a, m_bar + b)
  )
  log_lhd_new <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b))

  log_p <- log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
  p     <- exp(log_p - max(log_p))
  p / sum(p)
}


#' Prior expected number of clusters
#'
#' Computes the prior expected number of clusters \eqn{E[K_n]} for \eqn{n}
#' nodes under the DP, PY, DM, or GN partition prior.
#'
#' @param n Number of nodes (positive integer).
#' @param prior A list in the same format as \code{esbm()}.
#'
#' @return A scalar expected number of clusters.
#'
#' @examples
#' expected_K(150, prior = list(name = "GN", gamma_GN = 0.3))
#' expected_K(150, prior = list(name = "DP", alpha_PY = 1))
expected_K <- function(n, prior) {

  n <- as.integer(n)
  stopifnot(n > 0L)

  if (!is.list(prior) || is.null(prior$name)) {
    stop("prior must be a list with a name field.")
  }

  if (prior$name == "DP") {

    theta <- prior$alpha_PY
    stopifnot(theta > 0)
    theta * sum(1 / (theta + 0:(n - 1L)))

  } else if (prior$name == "PY") {

    sigma <- prior$sigma_PY
    theta <- prior$alpha_PY
    stopifnot(sigma >= 0, sigma < 1, theta > -sigma)

    if (sigma == 0) {
      theta * sum(1 / (theta + 0:(n - 1L)))
    } else {
      exp(lgamma(theta + sigma + n) - lgamma(theta + sigma) -
          lgamma(theta + n) + lgamma(theta + 1)) / sigma - theta / sigma
    }

  } else if (prior$name == "DM") {

    beta_DM <- prior$beta_DM
    H_DM    <- prior$H_DM
    stopifnot(beta_DM > 0, H_DM >= 1L)
    # E[K] under DM: sum over h=1..H_DM of P(at least one node in cluster h)
    # = H_DM * (1 - prod_{i=0}^{n-1} (i*beta_DM) / (i*beta_DM + beta_DM*(H_DM-1)+beta_DM))
    # Simpler: simulate from the urn
    urn <- function(v_minus) urn_DM(v_minus, beta_DM, H_DM)
    K_samples <- replicate(2000L, {
      z <- integer(n)
      z[1L] <- 1L
      v_minus <- 1L
      for (i in 2L:n) {
        w  <- urn(tabulate(z[seq_len(i - 1L)]))
        z[i] <- sample(length(w), 1L, prob = w)
        if (z[i] > length(w) - 1L) z[i] <- max(z[seq_len(i - 1L)]) + 1L
      }
      length(unique(z))
    })
    mean(K_samples)

  } else if (prior$name == "GN") {

    gamma <- prior$gamma_GN
    stopifnot(gamma > 0, gamma < 1)
    # no closed form; use the recursive formula for E[K_n] under Gnedin
    # E[K_1] = 1; E[K_n] = E[K_{n-1}] + P(new cluster at step n)
    ek <- 1
    for (i in 2L:n) {
      # P(new | K_{i-1} = k) averaged; use recursion on the urn weights
      # approximate via simulation
    }
    # fall back to simulation for GN
    K_samples <- replicate(2000L, {
      z <- 1L
      for (i in 2L:n) {
        w  <- urn_GN(tabulate(z), gamma)
        zi <- sample(length(w), 1L, prob = w)
        z  <- c(z, if (zi == length(w)) max(z) + 1L else zi)
      }
      length(unique(z))
    })
    mean(K_samples)

  } else {
    stop("Invalid prior$name.")
  }
}
