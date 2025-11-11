## simulate_spatial_esbm.R
## Spatial ESBM-style simulator: 4 clusters, 100 nodes
## - Nodes have (lon, lat) in [0,1]^2
## - Clusters = quadrants (NW, NE, SW, SE)
## - Edge logit: block baseline  −  gamma * distance
## - Undirected, binary, no self-loops

## simulate_spatial_three_scenarios.R
## Generates 3 spatial network scenarios (n=100, K=4)
## - Clusters depend on location as specified by scenario
## - Undirected, binary, no self-loops
## - Edge logit = block_baseline - gamma * distance

logit_inv <- function(eta) 1 / (1 + exp(-eta))

simulate_spatial_esbm <- function(
    n = 100,
    seed = 1,
    cluster_rule = c("lat", "lon", "both"),  # how clusters depend on location
    theta_within = 1.5,      # within-cluster baseline (logit)
    theta_between = -1.5,    # between-cluster baseline (logit)
    gamma = 4.0,             # distance penalty on logit
    scale_dist = TRUE
) {
  set.seed(seed)
  cluster_rule <- match.arg(cluster_rule)
  
  # Coordinates on unit square
  lon <- runif(n)
  lat <- runif(n)
  coords <- cbind(lon = lon, lat = lat)
  
  # Cluster labels (K=4) from location
  if (cluster_rule == "lat") {
    # 4 horizontal bands by latitude quartiles
    q <- quantile(lat, probs = c(.25, .5, .75))
    z_0 <- cut(lat, breaks = c(-Inf, q, Inf), labels = FALSE)
  } else if (cluster_rule == "lon") {
    # 4 vertical bands by longitude quartiles
    q <- quantile(lon, probs = c(.25, .5, .75))
    z_0 <- cut(lon, breaks = c(-Inf, q, Inf), labels = FALSE)
  } else {
    # "both": quadrants (NW=1, NE=2, SW=3, SE=4)
    z_0 <- ifelse(lon < .5 & lat >= .5, 1,
                  ifelse(lon >= .5 & lat >= .5, 2,
                         ifelse(lon < .5 & lat < .5, 3, 4)))
  }
  
  # Distances
  D <- as.matrix(dist(coords, method = "euclidean"))
  if (scale_dist) {
    mD <- mean(D[upper.tri(D)])
    if (mD > 0) D <- D / mD
  }
  
  # Block baselines (logit)
  K <- 4
  Theta <- matrix(theta_between, K, K)
  diag(Theta) <- theta_within
  
  # Edges with spatial decay
  Y <- matrix(0L, n, n)
  P <- matrix(0,  n, n)
  for (i in 1:(n-1)) {
    gi <- z_0[i]
    for (j in (i+1):n) {
      gj <- z_0[j]
      eta <- Theta[gi, gj] - gamma * D[i, j]
      p <- logit_inv(eta)
      P[i, j] <- P[j, i] <- p
      Y[i, j] <- Y[j, i] <- rbinom(1, 1, p)
    }
  }
  diag(Y) <- 0L
  
  list(Y = Y, z_0 = z_0, coords = coords, P = P, D = D, Theta = Theta,
       params = list(theta_within = theta_within, theta_between = theta_between,
                     gamma = gamma, rule = cluster_rule, seed = seed))
}


# 1) Latitude-only
s1 <- simulate_spatial_esbm(n = 100, seed = 101, cluster_rule = "lat")
with(s1, save(Y, z_0, coords, P, D, Theta, file = "network_spatial_scen1_lat.RData"))

# 2) Longitude-only
s2 <- simulate_spatial_esbm(n = 100, seed = 202, cluster_rule = "lon")
with(s2, save(Y, z_0, coords, P, D, Theta, file = "network_spatial_scen2_lon.RData"))

# 3) Both latitude & longitude (quadrants)
s3 <- simulate_spatial_esbm(n = 100, seed = 303, cluster_rule = "both")
with(s3, save(Y, z_0, coords, P, D, Theta, file = "network_spatial_scen3_both.RData"))

cat("Saved:\n",
    "  Simulation/network_spatial_scen1_lat.RData\n",
    "  Simulation/network_spatial_scen2_lon.RData\n",
    "  Simulation/network_spatial_scen3_both.RData\n")

###
plot_spatial_network <- function(sim) {
  stopifnot(all(c("Y", "z_0", "coords") %in% names(sim)))
  
  Y <- sim$Y
  z <- sim$z_0
  coords <- sim$coords
  K <- length(unique(z))
  cols <- rainbow(K)
  ord <- order(z)
  
  op <- par(mfrow = c(1, 2), mar = c(3, 3, 2, 2))
  
  ## 1) Spatial scatterplot
  plot(coords, col = cols[z], pch = 19, cex = 1.2, asp = 1,
       xlab = "Longitude", ylab = "Latitude",
       main = paste("Spatial layout:", sim$params$rule))
  legend("topright", legend = paste("Cluster", 1:K),
         col = cols, pch = 19, cex = 0.8, bty = "n")
  
  ## 2) Adjacency matrix (nodes ordered by cluster)
  image(t(Y[ord, ord][, nrow(Y):1]),
        col = gray.colors(20, start = 1, end = 0),
        axes = FALSE, main = "Adjacency matrix")
  box()
  
  par(op)
}

plot_spatial_network(s1)
plot_spatial_network(s2)
plot_spatial_network(s3)
