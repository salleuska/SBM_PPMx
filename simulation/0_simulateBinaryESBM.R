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


add_covariates_SBMsim <- function(baseObj,
  nCov = 1,
  cov1Dep = c("informative", "neutral", "mislead_random", "mislead_shifted"),
  cov2Dep = c("neutral", "informative"),
  seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  zRaw <- baseObj$partition
  z    <- as.integer(factor(zRaw))   # ensure labels are 1:4
  n    <- length(z)
  K    <- length(unique(z))
  if (K != 4) stop("Current covariate design assumes exactly 4 clusters.")
  if (!(nCov %in% c(1, 2))) stop("nCov must be 1 or 2.")

  cov1Dep <- match.arg(cov1Dep)
  cov2Dep <- match.arg(cov2Dep)

  x <- matrix(NA_real_, n, 2)

  ## ------------------------------------------------------------##
  ## Parameters
  ## ------------------------------------------------------------##
  mx1 <- c(-1.5, -0.5, 0.5, 1.5)   # informative means
  sx1 <- 0.5

  mxN <- 0                         # noise
  sxN <- 1

  ## ------------------------------------------------------------##
  ## Build labels driving x1
  ## ------------------------------------------------------------##
  if (cov1Dep == "informative") {

    lab1 <- z

  } else if (cov1Dep == "neutral") {

    lab1 <- NULL

  } else if (cov1Dep == "mislead_random") {

    ## independent, balanced wrong partition
    lab1 <- sample(rep(seq_len(K), length.out = n))

  } else if (cov1Dep == "mislead_shifted") {

    ## shift half of each true cluster to the next (cyclic)
    lab1 <- z
    for (k in 1:K) {
      idxK     <- which(z == k)
      idxShift <- sample(idxK, size = floor(length(idxK) / 2), replace = FALSE)
      kNext    <- if (k == K) 1 else (k + 1)
      lab1[idxShift] <- kNext
    }
  }

  ## ------------------------------------------------------------##
  ## Generate x1
  ## ------------------------------------------------------------##
  if (is.null(lab1)) {
    x[, 1] <- rnorm(n, mxN, sxN)
  } else {
    for (i in 1:n) {
      x[i, 1] <- rnorm(1, mx1[lab1[i]], sx1)
    }
  }

  ## ------------------------------------------------------------##
  ## Generate x2
  ## ------------------------------------------------------------##
  if (nCov == 1) {

    ## plotting-only coordinate
    x[, 2] <- rnorm(n, mxN, sxN)

  } else {

    ## second covariate may be informative or neutral
    if (cov2Dep == "neutral") {
      x[, 2] <- rnorm(n, mxN, sxN)
    } else {
      for (i in 1:n) {
        x[i, 2] <- rnorm(1, mx1[z[i]], sx1)
      }
    }
  }

  ## ------------------------------------------------------------##
  ## Attach and return
  ## ------------------------------------------------------------##
  baseObj$x       <- x
  baseObj$nCov    <- nCov
  baseObj$cov1Dep <- cov1Dep
  baseObj$cov2Dep <- if (nCov == 2) cov2Dep else "neutral"

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
# 2) One-covariate scenarios (x1 varies; x2 is plotting noise)
# ------------------------------------------------------------
sim_1_neutral   <- add_covariates_SBMsim(base, nCov = 1, cov1Dep = "neutral")
sim_1_info      <- add_covariates_SBMsim(base, nCov = 1, cov1Dep = "informative")
sim_1_mis_rand  <- add_covariates_SBMsim(base, nCov = 1, cov1Dep = "mislead_random")
sim_1_mis_shift <- add_covariates_SBMsim(base, nCov = 1, cov1Dep = "mislead_shifted")

# ------------------------------------------------------------
# 3) Two-covariate scenarios (4 combos: {neutral, informative}^2)
# ------------------------------------------------------------
sim_2_NN <- add_covariates_SBMsim(base, nCov = 2, cov1Dep = "neutral",     cov2Dep = "neutral")
sim_2_IN <- add_covariates_SBMsim(base, nCov = 2, cov1Dep = "informative", cov2Dep = "neutral")
sim_2_NI <- add_covariates_SBMsim(base, nCov = 2, cov1Dep = "neutral",     cov2Dep = "informative")
sim_2_II <- add_covariates_SBMsim(base, nCov = 2, cov1Dep = "informative", cov2Dep = "informative")

# ------------------------------------------------------------
# 4) Save
# ------------------------------------------------------------
dir.create(here("simulation/data"), recursive = TRUE, showWarnings = FALSE)

saveRDS(sim_1_neutral,   here("simulation/data/binarySBM_1cov_neutral.rds"))
saveRDS(sim_1_info,      here("simulation/data/binarySBM_1cov_informative.rds"))
saveRDS(sim_1_mis_rand,  here("simulation/data/binarySBM_1cov_mislead_random.rds"))
saveRDS(sim_1_mis_shift, here("simulation/data/binarySBM_1cov_mislead_shifted.rds"))

saveRDS(sim_2_NN, here("simulation/data/binarySBM_2cov_NN.rds"))
saveRDS(sim_2_IN, here("simulation/data/binarySBM_2cov_IN.rds"))
saveRDS(sim_2_NI, here("simulation/data/binarySBM_2cov_NI.rds"))
saveRDS(sim_2_II, here("simulation/data/binarySBM_2cov_II.rds"))

# ------------------------------------------------------------
# 5) Quick plots (optional)
# ------------------------------------------------------------
plot_adj_matrix_gg(sim_1_neutral)

plot_network_geo(sim_1_neutral)
plot_network_geo(sim_1_info)
plot_network_geo(sim_1_mis_rand)
plot_network_geo(sim_1_mis_shift)

plot_network_geo(sim_2_NN)
plot_network_geo(sim_2_IN)
plot_network_geo(sim_2_NI)
plot_network_geo(sim_2_II)



# set.seed(1116)

# # Generate one network - sparse
# base <- simBinarySBM(p_within = 0.4, p_between = 0.2)

# # Generate three covariate scenarios using same Y and same z
# sim_none  <- add_covariates_SBMsim(base, "neutral")
# sim_one   <- add_covariates_SBMsim(base, "one")
# sim_both  <- add_covariates_SBMsim(base, "both")

# ## create directory if it does not exists
# dir.create(here("simulation/data"), 
#   recursive = TRUE, showWarnings = FALSE)

# saveRDS(sim_none, here("simulation/data/binarySBM_neutral.rds"))
# saveRDS(sim_one, here("simulation/data/binarySBM_one.rds"))
# saveRDS(sim_both, here("simulation/data/binarySBM_both.rds"))


# plot_adj_matrix_gg(sim_none)
# plot_network_geo(sim_none)
# plot_network_geo(sim_one)
# plot_network_geo(sim_both)

#######################################################
#######################################################
#### BACKUP: one function
# simBinarySBMx <- function(
#   clust_sizes = c(50, 40, 30, 30),        # cluster sizes
#   p_within    = 0.8,                      # within-cluster edge prob
#   p_between   = 0.1,                      # between-cluster edge prob
#   cov_dep     = c("neutral", "one", "both")  # covariate–cluster dependence
# ) {
#   cov_dep <- match.arg(cov_dep)

#   n <- sum(clust_sizes)       # sample size
#   K <- length(clust_sizes)    # number of clusters
#   z <- rep(1:K, times = clust_sizes)   # cluster memberships

#   if (K != 4L) stop("Current covariate design assumes exactly 4 clusters.")

#   # Block matrix Theta
#   Theta <- matrix(p_between, K, K)
#   diag(Theta) <- p_within

#   #--- SAMPLE ADJACENCY MATRIX Y ---
#   Y <- matrix(0L, n, n)

#   for (i in 1:(n - 1)) {
#     for (j in (i + 1):n) {
#       pij      <- Theta[z[i], z[j]]
#       Y_ij     <- rbinom(1, 1, pij)
#       Y[i, j]  <- Y_ij
#       Y[j, i]  <- Y_ij
#     }
#   }

#   #--- SAMPLE COVARIATES x (depends on cov_dep) ---
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

#   } else {  # cov_dep == "neutral"

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