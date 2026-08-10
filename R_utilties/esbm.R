####################################################################################
# FUNCTIONS TO COMPUTE THE DIFFERENT URN SCHEMES FOR THE PRIORS IN THE ARTICLE #####
####################################################################################

urn_DP <- function(v_minus,alpha_PY){
  return(c(v_minus,alpha_PY))
}

urn_PY <- function(v_minus,alpha_PY,sigma_PY){
  H<-length(v_minus)
  return(c(v_minus-sigma_PY,alpha_PY+H*sigma_PY))
}

urn_DM <- function(v_minus,beta_DM,H_DM){
  H<-length(v_minus)
  return(c(v_minus+beta_DM,beta_DM*(H_DM-H)*(H_DM>H)))
}

urn_GN <- function(v_minus,gamma_GN){
  H<-length(v_minus)
  return(c((v_minus+1)*(sum(v_minus)-H+gamma_GN),H^2-H*gamma_GN))
}

####################################################################################
# GIBBS SAMPLER FOR THE EXTENDED STOCHASTIC BLOCK MODEL  ###########################
####################################################################################
#' Gibbs sampler for the extended stochastic block model
#'
#' Runs a Gibbs sampler for the extended stochastic block model with a Gibbs-type
#' partition prior and optional covariate-dependent PPMx similarity weights.
#'
#' @param Y A \eqn{V \times V} symmetric adjacency matrix.
#' @param seed Integer seed for reproducibility.
#' @param N_iter Number of MCMC iterations.
#' @param prior A list specifying the partition prior. The entry \code{name} must
#'   be one of \code{"DP"}, \code{"PY"}, \code{"DM"}, or \code{"GN"}. Additional
#'   entries define the prior-specific hyperparameters; for example,
#'   \code{list(name = "DP", alpha_PY = 1)},
#'   \code{list(name = "PY", alpha_PY = 1, sigma_PY = 0.2)},
#'   \code{list(name = "DM", beta_DM = 1, H_DM = 10)}, or
#'   \code{list(name = "GN", gamma_GN = 0.3)}.
#' @param z_init Integer vector of length \eqn{V} giving the initial cluster
#'   assignment for each node. Defaults to one cluster per node.
#' @param a,b Hyperparameters of the Beta prior on the block connection
#'   probabilities.
#' @param x Optional data frame of node covariates. If provided, \code{nrow(x)}
#'   must equal \code{nrow(Y)}.
#' @param similarity_fun Optional similarity function or list of similarity
#'   functions. If a list is provided, its length must equal \code{ncol(x)}.
#'   Each function must support two modes: \code{mode = "node"}, returning log
#'   similarity values for assigning the current node to each existing cluster
#'   and to a new cluster; and \code{mode = "partition"}, returning log
#'   similarity values for the clusters in the current partition.
#' @param sim_args Optional list of additional arguments passed to the similarity
#'   functions. If multiple similarity functions are used, this should be a list
#'   of the same length as \code{similarity_fun}.
#' @param alpha_g Scalar, vector, or list specifying the covariate weights.
#'   If numeric, \code{alpha_g} is treated as fixed. It may be a scalar, reused
#'   for all similarity functions, or a vector with one entry per similarity
#'   function. If a list, \code{alpha_g} is learned within the Gibbs sampler via
#'   a conjugate Gamma update and should contain \code{init} (positive starting
#'   value), \code{a_alpha} (Gamma shape), and \code{b_alpha} (Gamma rate). For
#'   example, \code{list(init = 1, a_alpha = 2, b_alpha = 1)} or
#'   \code{list(init = c(1,2), a_alpha = c(2,2), b_alpha = c(1,1))}.
#' @param similarity_calibration Character, either \code{"raw"} (default) or
#'   \code{"normalized"}. Controls whether the log similarity \eqn{\log g()} is
#'   used as-is (\code{"raw"}) or normalised by its log-sum-exp across clusters
#'   (\code{"normalized"}) before being weighted by \code{alpha_g}. The
#'   \code{"raw"} form corresponds to the PPMx model as defined in the paper;
#'   \code{"normalized"} follows an alternative formulation from the literature
#'   and is provided for comparison.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{z_post}}{A \eqn{V \times N_{\mathrm{iter}}} integer matrix of
#'   posterior samples of cluster labels. Entry \code{z_post[v, t]} gives the
#'   cluster label of node \eqn{v} at iteration \eqn{t}.}
#'   \item{\code{alpha_g_post}}{A \eqn{J \times N_{\mathrm{iter}}} matrix of
#'   posterior samples of the covariate weights. Returned only when \code{alpha_g}
#'   is specified as a list; \code{NULL} otherwise.}
#'   \item{\code{alpha_rate_post}}{A \eqn{J \times N_{\mathrm{iter}}} matrix of
#'   the Gamma rate parameter used at each update of \code{alpha_g}. Useful for
#'   diagnosing the sampler. Returned only when \code{alpha_g} is a list;
#'   \code{NULL} otherwise.}
#' }
#'
#' @details
#' The Gibbs update for \eqn{\alpha_g^{(j)}} exploits the conjugacy of the
#' Gamma prior with the likelihood contribution of the partition-level similarity
#' \eqn{\log g()}: the posterior is
#' \eqn{\alpha_g^{(j)} \mid \text{rest} \sim
#'   \mathrm{Gamma}(a_\alpha^{(j)},\; b_\alpha^{(j)} - \sum_h \log g(X_h^*))}.
#' The rate \eqn{b_\alpha - \sum_h \log g_h} must be positive for the Gamma to
#' be well-defined. If it is not (which can occur when log similarities are
#' large and positive), the current value of \code{alpha_g} is kept and a
#' warning is issued.
#'
#' @examples
#' prior <- list(name = "GN", gamma_GN = 0.3)
#' fit   <- esbm(Y = Y, seed = 1, N_iter = 1000, prior = prior)
#'
#' # with learned alpha
#' fit2  <- esbm(Y = Y, seed = 1, N_iter = 1000, prior = prior,
#'               x = x_df, similarity_fun = sim_fun, sim_args = sim_args,
#'               alpha_g = list(init = 1, a_alpha = 2, b_alpha = 1))

 
esbm <- function(Y, seed, N_iter, prior, 
                 z_init = seq_len(nrow(Y)), 
                 a = 1, b = 1,
                 x = NULL, similarity_fun = NULL, sim_args = list(),
                 alpha_g = 1,
                 similarity_calibration = c("raw", "normalized")) {
 
  # ----------------------------------------------
  # Selection of the prior distribution to be used
  # ----------------------------------------------
  
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

  
  # ----------------------------------------------
  # SP: checks if covariates have been provided
  # and if so we check on whether there is a prior
  # for alpha or not
  # ----------------------------------------------  
  ## Default: no covariates, no learned alpha

  ## ----------------------------------------------
  ## Covariate and alpha checks
  ## ----------------------------------------------

  has_covariates <- !is.null(x)
  learn_alpha_g  <- FALSE

  J <- 0L
  a_alpha <- NULL
  b_alpha <- NULL

  if (!has_covariates) {

    similarity_fun <- NULL
    sim_args <- list()
    alpha_g <- NULL

  } else {

    message("Note: covariates have been provided.")

    if (!is.data.frame(x)) {
      stop("x must be a data.frame when covariates are provided.")
    }

    if (ncol(x) == 0L) {
      stop("x must have at least one column.")
    }

    if (is.null(similarity_fun)) {
      stop("similarity_fun must be provided when x is provided.")
    }

    if (is.function(similarity_fun)) {
      similarity_fun <- list(similarity_fun)
    }

    if (!is.list(similarity_fun) || !all(vapply(similarity_fun, is.function, logical(1)))) {
      stop("similarity_fun must be a function or a list of functions.")
    }

    J <- length(similarity_fun)

    if (ncol(x) != J) {
      stop("ncol(x) must equal length(similarity_fun).")
    }

    if (!is.list(sim_args)) {
      stop("sim_args must be a list.")
    }

    if (length(sim_args) == 0L) {
      sim_args <- replicate(J, list(), simplify = FALSE)
    }

    if (length(sim_args) == 1L && J > 1L) {
      sim_args <- replicate(J, sim_args[[1L]], simplify = FALSE)
    }

    if (length(sim_args) != J) {
      stop("sim_args must have length 0, 1, or J.")
    }

    similarity_calibration <- match.arg(similarity_calibration)

    learn_alpha_g <- is.list(alpha_g)

    if (learn_alpha_g) {

      required_names <- c("init", "a_alpha", "b_alpha")

      if (!all(required_names %in% names(alpha_g))) {
        stop("When alpha_g is a list, it must contain init, a_alpha, and b_alpha.")
      }

      alpha_g_init <- alpha_g$init
      a_alpha <- alpha_g$a_alpha
      b_alpha <- alpha_g$b_alpha

      if (length(alpha_g_init) == 1L) alpha_g_init <- rep(alpha_g_init, J)
      if (length(a_alpha) == 1L) a_alpha <- rep(a_alpha, J)
      if (length(b_alpha) == 1L) b_alpha <- rep(b_alpha, J)

      if (length(alpha_g_init) != J) stop("alpha_g$init must have length 1 or J.")
      if (length(a_alpha) != J) stop("alpha_g$a_alpha must have length 1 or J.")
      if (length(b_alpha) != J) stop("alpha_g$b_alpha must have length 1 or J.")

      if (any(alpha_g_init <= 0)) stop("alpha_g$init must be positive.")
      if (any(a_alpha <= 0)) stop("alpha_g$a_alpha must be positive.")
      if (any(b_alpha <= 0)) stop("alpha_g$b_alpha must be positive.")

      alpha_g <- alpha_g_init

    } else {

      if (length(alpha_g) == 1L) {
        alpha_g <- rep(alpha_g, J)
      }

      if (length(alpha_g) != J) {
        stop("alpha_g must have length 1 or J.")
      }

      if (any(alpha_g < 0)) {
        stop("Fixed alpha_g values must be nonnegative.")
      }
    }
  }

  

  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  
  set.seed(seed)
  
  V <- nrow(Y)
  
  # cluster assignments are encoded in two equivalent ways:
  # [i] a VxH matrix Z, s.t. Z[v,h]=1{node v in cluster h}, faster to use within each iteration
  Z <- diag(max(z_init))[z_init, , drop = FALSE]
  # Z <- vec2mat() - remove unnecessary functions
  
  
  # [ii] a vector of length V containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[v,t]=h if node v is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  z_post <- matrix(NA,V,N_iter)
  
  # set up speace to save the posterior for alpha if present
  alpha_g_post <- NULL
  alpha_rate_post <- NULL
  if (learn_alpha_g) {
    alpha_g_post <- matrix(NA, nrow = J, ncol = N_iter)

    ### TMP for diagnostic
    alpha_rate_post <- matrix(NA, nrow = J, ncol = N_iter)
  }
  
  # Create the matrix with block connections
  temp   <- Y%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  # ----------------------------------------------
  # Beginning of the Gibbs sampler
  # ----------------------------------------------
  
  for (t in 1:N_iter){
    for (v in 1:V){
      
      # remove empty clusters and
      # if the cluster containing node v has no other node, discard it as well
      if(ncol(Z) > 1){
        nonempty_v <- which(colSums(Z[-v,]) > 0)  
        Z <- Z[, nonempty_v]
        if (length(nonempty_v)==1){Z <- matrix(Z,V,1)}
        
        # Reduce the dimensions of the m_full matrix
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y[-v,v])
      h_v         <- which(Z[v,] > 0)
      
      # Compute the m matrix by difference
      if(length(h_v) == 1){
        resid1       <- matrix(0,H,H)
        resid1[,h_v] <- r_v; resid1[h_v,] <- r_v
        m            <- m_full - resid1
      } else {m <- m_full} # No need to update m in this case
      
      # m_bar = number of non-edges between cluster h and cluster k, excluding node v 
      m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
      
      # ----------------------------------------------
      # Computing the probabilities
      # ----------------------------------------------
      
      log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # vector of length H
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
      

          
      # ----------------------------------------------
      # SP - Covariates part: Similarity g()
      # log-weights for H existing + 1 new cluster
      # covariates weighted with alpha
      # ----------------------------------------------
 
      if (is.null(similarity_fun)) {

        log_similarity_g <- rep(0, H + 1)

      } else {

        # initialize combined log g
        log_similarity_g <- rep(0, H + 1)

        # loop over similarities
        for (j in seq_len(J)) {

          log_g_j <- similarity_fun[[j]](
            mode = "node",
            Z_minus_v     = Z_v,
            cluster_sizes = v_minus,
            v_index       = v,
            H             = H,
            x             = x[, j],
            args          = sim_args[[j]]
          )

          if (similarity_calibration == "normalized") {

            lse_j <- max(log_g_j) +
              log(sum(exp(log_g_j - max(log_g_j))))

            log_g_j <- log_g_j - lse_j

            log_similarity_g <- log_similarity_g + alpha_g[j] * log_g_j

          } else if (similarity_calibration == "raw") {

            # log g() left unstandardized: weighted sum across covariates
            log_similarity_g <- log_similarity_g + alpha_g[j] * log_g_j
          }
        }
      }
     
      # ----------------------------------------------
      # SP: clustering probabilities
      # ----------------------------------------------
      # # Probabilities
      log_p <- log_similarity_g + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
      p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
      
      # ----------------------------------------------
      # Sampling the indicator
      # ----------------------------------------------
      
      new_sample <- which(rmultinom(1,1,p) > 0)
      
      # ----------------------------------------------
      # Adjusting Z, H, r_v and m
      # ----------------------------------------------
      
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      # If a new value is sampled, increase the dimension of m_full
      if(new_sample== H+1){
        Z               <- cbind(Z,rep(0,V))
        Z[v,new_sample] <- 1
        m               <- rbind(cbind(m,0),0)
        H               <- H + 1
        r_v             <- crossprod(Z[-v,], Y[-v,v])
      } else {Z[v, new_sample] <- 1}
      
      # Updating m_full
      resid2              <- matrix(0,H,H)
      resid2[,new_sample] <- r_v; resid2[new_sample,] <- r_v
      m_full              <- m + resid2
    }
    
 
    # update and store alpha_g if learned
    if (learn_alpha_g) {

      for (j in seq_len(J)) {

        log_g_h <- similarity_fun[[j]](
          mode = "partition",
          Z    = Z,
          x    = x[, j],
          args = sim_args[[j]]
        )

        if (similarity_calibration == "normalized") {

          lse <- max(log_g_h) +
            log(sum(exp(log_g_h - max(log_g_h))))

          log_g_h <- log_g_h - lse

          rate_alpha <- b_alpha[j] - sum(log_g_h)

        } else if (similarity_calibration == "raw") {

          rate_alpha <- b_alpha[j] - sum(log_g_h)
        }

        alpha_rate_post[j, t] <- rate_alpha

        if (rate_alpha <= 0) {
          warning(sprintf(
            "Iteration %d, covariate %d: Gamma rate = %.4f <= 0; keeping current alpha_g.",
            t, j, rate_alpha
          ))
        } else {
          alpha_g[j] <- rgamma(1, shape = a_alpha[j], rate = rate_alpha)
        }
      }

      alpha_g_post[, t] <- alpha_g
    }


    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[, t] <- Z %*% seq_len(ncol(Z))

    
    #print(table(z_post[,t])) 
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }

  return(list(
    z_post = z_post,
    alpha_g_post = alpha_g_post, 
    alpha_rate_post = alpha_rate_post
  ))

}


# Utility functions (log_pY_z, pred_esbm, expected_K) live in esbm_utils.R.

