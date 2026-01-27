library(MASS)

my_sim(
  ns=c(50,40,30,30),                  # cluster sizes (they imply n and K, see below)
  p_within=.8, p_between=.1,          # edge probabilities (default = clear assortativity) 
  cov_dep=c("neutral","one","both")){    # 3 scenarios of covariates' dependence on clusters
  
  cov_dep <- match.arg(cov_dep)
  n = sum(ns)                             # sample size
  K = length(ns)                          # number of clusters
  z = rep(1:length(ns), times = ns)       # cluster memberships

  Theta = matrix(p_between, K, K)         
  diag(Theta) = p_within                  # block matrix 
  
  # SAMPLE ADJACENCY MATRIX Y 
  # (this step does NOT vary with cov_dep)
  # NOTE: you can even use the same Y for all the 3 scenarios
  
  Y <- matrix(0L, n, n)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      Y[i, j] <- Y[j, i] <- rbinom(1, 1, Theta[z[i],z[j]])
    }
  }
  
  # SAMPLE COVARIATES (this step varies with cov_dep)
  
  x<- matrix(NA, n, 2)
  
  if (cov_dep="both"){
    mx = matrix(c(1,1,
                  1,-1,
                  -1,1,
                  -1,-1),
                K,2,byrow = T)
    s = 0.5
    Sx = s*diag(2)
    for (i in 1:n) x[i,] = mvrnorm(1,mx[z[i],],Sx) 
  } else if (cov_dep=="one"){
    mx1 = c(-1.5, -0.5, +0.5, +1.5)
    sx1 = 0.5
    mx2 = 0
    sx2 = 1
    for (i in 1:n) x[i,] = c(rnorm(1,mx1[z[i]],sx1), rnorm(1,mx2,sx2)) 
  } else { # that is, if cov_dep="neutral"
    mx = c(0,0)
    Sx = diag(2)
    x = mvrnorm(n,mx,Sx) 
  }
}