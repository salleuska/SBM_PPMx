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

# Inputs:
# Y = VxV symmetric adjacency matrix
# z_init = V-vector of initialization assignment for each node (default = one cluster for each node)
# a,b = parameters of the Beta prior on the thetas
# prior = string in c("DP","PY","DM","GN")
# if prior=="DP", alpha_PY must be set
# if prior=="PY", alpha_PY and sigma_PY must be set
# if prior=="DM", beta_DM and H_DM must be set
# if prior=="GN", gamma_GN must be set

## SP 
# x  is a data.frame 
# similarity_fun - either one R function or a list of R functions 
# assume the user passes a well-formed x (vector, matrix, or factor) and the function knows how to handle it.
# sim_args = list of arguments for the similiarity function. Note that the similarity_fun needs to handle checks    # V x C


# Output:
# Posterior samples of the community labels for each node v=1,...,V

 
  ## possibile controllo
  ## if x is a dataframe
  ## ncol(Y) = nrow(Y) = nrow(x)
  ## ncol(x) = length(similarity_fun)
  ## similarity_fun is list 
  ## da aggiungere ora un vettore alpha - length(alpha) = length(similarity_fun)


esbm <- function(Y, seed, N_iter, prior, z_init=c(1:nrow(Y)), a=1, b=1,
                 alpha_PY=NA, sigma_PY=NA, beta_DM=NA, H_DM=NA, gamma_GN=NA, 
                 x=NULL,
                 similarity_fun = NULL, 
                 sim_args = list(), 
                 alpha_g = 1){
 
  # ----------------------------------------------
  # Selection of the prior distribution to be used
  # ----------------------------------------------
  
  if (prior=="DP"){
    urn<-function(v_minus){return(urn_DP(v_minus,alpha_PY))}
  } else if (prior=="PY"){
    urn<-function(v_minus){return(urn_PY(v_minus,alpha_PY,sigma_PY))}
  } else if (prior=="DM"){
    urn<-function(v_minus){return(urn_DM(v_minus,beta_DM,H_DM))}
  } else if (prior=="GN"){
    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  } else { 
    stop("Invalid value for prior")  
  }
  
  # ----------------------------------------------
  # SP: checks if covariates have been provided
  # ----------------------------------------------  
  if (!is.null(x)){
    print("Note : covariates have been provided")
    # x must be a data.frame
    if (!is.data.frame(x)) {
      stop("When similarity_fun is used, x must be a data.frame.")
    }

    # coerce similarity_fun to list so Gibbs loop is unified
    if (is.function(similarity_fun)) {
      similarity_fun <- list(similarity_fun)
    } else if (!is.list(similarity_fun)) {
      stop("similarity_fun must be either a function or a list of functions.")
    }

    J <- length(similarity_fun)

    # x must have exactly one column per similarity function
    if (ncol(x) != J) {
      stop("ncol(x) must equal length(similarity_fun).")
    }

    ## sim_args must also be a list of length J
    if (!is.list(sim_args)) {
      stop("sim_args must be a list when similarity_fun is provided.")
    }
    if (length(sim_args) == 1L && J > 1L) {
      sim_args <- replicate(J, sim_args, simplify = FALSE)
    }
    if (length(sim_args) != J) {
      stop("sim_args must have length equal to length(similarity_fun).")
    }

    # alpha_g must be scalar or length J
    if (length(alpha_g) == 1L) {
      alpha_g <- rep(alpha_g, J)
    } else if (length(alpha_g) != J) {
      stop("alpha_g must be scalar or a vector of length equal to similarity_fun.")
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
            Z_minus_v     = Z_v,
            cluster_sizes = v_minus,
            v_index       = v,
            H             = H,
            x             = x[, j],        # one column per similarity
            args          = sim_args[[j]]
          )
          # --- normalize g_j across clusters ---
          lse_j     <- log(sum(exp(log_g_j)))          # log(sum exp(log g_j))
          log_g_j_n <- log_g_j - lse_j                 # log normalized g_j

          # add covariate similarity, weighted by alpha_g[j]
          log_similarity_g <- log_similarity_g + (alpha_g[j] * log_g_j_n)
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
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z %*% c(1:ncol(Z))
    
    #print(table(z_post[,t])) 
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }
  return(z_post)
}


# ----------------------------------------------
# INTERNAL USAGE ONLY
# ----------------------------------------------

lV_py <- function(n, k, sigma, theta){
  if(k==1) return(lgamma(theta + 1) - lgamma(theta + n))
  index <- 1:(k-1)
  sum(log(theta + index*sigma)) + lgamma(theta + 1) - lgamma(theta + n)
}

# ----------------------------------------------
# INTERNAL USAGE ONLY
# ----------------------------------------------

lV_py <- Vectorize(lV_py, vectorize.args = "k")

dclust_py <- function(n, sigma, theta, H, log_scale=FALSE){
  
  probs <- numeric(n) - Inf
  if(H == Inf){
    index <- 1:n
    if(sigma > 0) {
      probs  <- c(lgfactorial(n, sigma)) - index*log(sigma) + lV_py(n,index,sigma,theta)
    } else {
      probs  <- c(lastirling1(n)) + lV_py(n,index,sigma,theta)
    }
  }
  if(H < Inf){
    index <- 1:min(n,H)
    if(sigma == 0){
      probs[index]  <- lfactorial(H) - lfactorial(H-index) + c(lgfactorial_ns(n,-theta/H))[index] + lgamma(theta) - lgamma(theta + n)
    }
  }
  if(log_scale) return(probs)
  return(exp(probs))
}


# ----------------------------------------------
# GNEDIN
# ----------------------------------------------

HGnedin <- function(V, h, gamma=0.5){
  exp(lchoose(V, h) + lgamma(h-gamma) - lgamma(1-gamma) + log(gamma) + lgamma(V+ gamma - h) - lgamma(V +gamma))
}

