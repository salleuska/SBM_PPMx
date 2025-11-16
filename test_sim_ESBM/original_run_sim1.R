## scenario1_run.R — ESBM: Scenario 1 end-to-end
## Reproduces simulation analyses for network_1.RData
## ---- Setup ----
rm(list = ls())
suppressPackageStartupMessages({
  library(Rcpp)
  library(here)
  library(reshape)
  library(gdata)
  library(igraph)
  library(mcclust.ext)
  library(RColorBrewer)
  library(pheatmap)
  library(gridExtra)
  library(grid)
  library(cowplot)
  library(ggplot2)
  library(coda)
##   library(dummies) - not available
  library(randnet)
  library(greed)
  library(LaplacesDemon)
})


## ---- Source ----
source(here("R_utilties", "esbm.R"))
Rcpp::sourceCpsp(here("R_utilties", "stirling.cpp"))

## ---- Data ----
load(here("test_sim_ESBM", "network_1.RData"))  # expects: Y, z_0 (true labels)
V <- nrow(Y); diag(Y) <- 0

## ---- Hyperparameters (priors) ----
# Beta(a,b) for block probabilities
a <- 1; b <- 1

# Choose Gibbs-type priors to target E[H] ~ 10 for V=80
# (uses helper funcs in esbm.R: expected_cl_py, HGnedin)
sigma_dm <- 0; H_dm <- 50; beta_dm <- 3.5 / H_dm
sigma_dp <- 0; alpha_dp <- 3
sigma_py <- 0.6; alpha_py <- -0.3
gamma_gn <- 0.45

invisible(expected_cl_py(V, sigma = sigma_dm, theta = beta_dm * H_dm, H = H_dm))
invisible(expected_cl_py(V, sigma = sigma_dp, theta = alpha_dp, H = Inf))
invisible(expected_cl_py(V, sigma = sigma_py, theta = alpha_py, H = Inf))
invisible(sum(1:V * HGnedin(V, 1:V, gamma = gamma_gn)))

## ---- Collapsed Gibbs sampler (unsupervised) ----
N_iter <- 2000
my_seed <- 1
my_z <- seq_len(V)

# DM
Z_DM <- esbm(Y, my_seed, N_iter, "DM", my_z, a = a, b = b,
             beta_DM = beta_dm, H_DM = H_dm)

# DP (CRP)
Z_DP <- esbm(Y, my_seed, N_iter, "PY", my_z, a = a, b = b,
             alpha_PY = alpha_dp, sigma_PY = sigma_dp)

# PY
Z_PY <- esbm(Y, my_seed, N_iter, "PY", my_z, a = a, b = b,
             alpha_PY = alpha_py, sigma_PY = sigma_py)

# GN
Z_GN <- esbm(Y, my_seed, N_iter, "GN", my_z, a = a, b = b,
             gamma_GN = gamma_gn)

save(Z_DP, Z_PY, Z_GN, Z_DM, file = here("Simulation", "Posterior_No_Attributes1.RData"))
rm(Z_DP, Z_PY, Z_GN, Z_DM)

## ---- Collapsed Gibbs sampler (supervised with node attributes) ----
# Here we use the *true labels* as attributes for the supervised version (as in paper)
N_iter <- 10000
my_seed <- 1
my_z <- seq_len(V)
my_x <- z_0                          # node attribute (supervised)
my_alpha_xi <- rep(1, length(unique(z_0)))

# DM + x
Z_DM_x <- esbm(Y, my_seed, N_iter, "DM", my_z, a = a, b = b,
               beta_DM = beta_dm, H_DM = H_dm, x = my_x, alpha_xi = my_alpha_xi)

# DP + x
Z_DP_x <- esbm(Y, my_seed, N_iter, "PY", my_z, a = a, b = b,
               alpha_PY = alpha_dp, sigma_PY = sigma_dp, x = my_x, alpha_xi = my_alpha_xi)

# PY + x
Z_PY_x <- esbm(Y, my_seed, N_iter, "PY", my_z, a = a, b = b,
               alpha_PY = alpha_py, sigma_PY = sigma_py, x = my_x, alpha_xi = my_alpha_xi)

# GN + x
Z_GN_x <- esbm(Y, my_seed, N_iter, "GN", my_z, a = a, b = b,
               gamma_GN = gamma_gn, x = my_x, alpha_xi = my_alpha_xi)

save(Z_DP_x, Z_PY_x, Z_GN_x, Z_DM_x, file = here("Simulation", "Posterior_Attributes1.RData"))
rm(Z_DP_x, Z_PY_x, Z_GN_x, Z_DM_x)

## ---- Posterior inference (WAIC, trace sanity) ----
burn_in <- 10000
set.seed(1)
index_traceplot <- sample.int(V * (V - 1) / 2, 1)

LL <- matrix(nrow = V * (V - 1) / 2, ncol = N_iter - burn_in)

load(here("Simulation", "Posterior_No_Attributes1.RData"))
load(here("Simulation", "Posterior_Attributes1.RData"))

waic_model <- function(Zmat) {
  Zs <- Zmat[, (burn_in + 1):N_iter, drop = FALSE]
  for (t in seq_len(ncol(Zs))) {
    LL[, t] <- sampleLL(Zs[, t], Y, a, b)
    if (t %% 10000 == 0) message("Iteration: ", t)
  }
  list(WAIC = WAIC(LL)$WAIC, trace = LL[index_traceplot, ])
}

# Unsupervised
Z_DM_WAIC <- waic_model(Z_DM)
Z_DP_WAIC <- waic_model(Z_DP)
Z_PY_WAIC <- waic_model(Z_PY)
Z_GN_WAIC <- waic_model(Z_GN)

# Supervised
Z_DM_WAIC_x <- waic_model(Z_DM_x)
Z_DP_WAIC_x <- waic_model(Z_DP_x)
Z_PY_WAIC_x <- waic_model(Z_PY_x)
Z_GN_WAIC_x <- waic_model(Z_GN_x)

## ---- VI vs truth, posterior H quantiles ----
vi_pair <- function(Zmat) c(
  VI_unsup = VI(z_0, t(Zmat[, (burn_in + 1):N_iter, drop = FALSE]))
)

vi_unsup <- list(
  DM = vi_pair(Z_DM),
  DP = vi_pair(Z_DP),
  PY = vi_pair(Z_PY),
  GN = vi_pair(Z_GN)
)

vi_sup <- list(
  DMx = VI(z_0, t(Z_DM_x[, (burn_in + 1):N_iter, drop = FALSE])),
  DPx = VI(z_0, t(Z_DP_x[, (burn_in + 1):N_iter, drop = FALSE])),
  PYx = VI(z_0, t(Z_PY_x[, (burn_in + 1):N_iter, drop = FALSE])),
  GNx = VI(z_0, t(Z_GN_x[, (burn_in + 1):N_iter, drop = FALSE]))
)

H_q <- function(Zmat) {
  quantile(apply(Zmat[, (burn_in + 1):N_iter, drop = FALSE], 2L, max))[2:4]
}

Hq_unsup <- rbind(
  DM = H_q(Z_DM), DP = H_q(Z_DP), PY = H_q(Z_PY), GN = H_q(Z_GN)
)
Hq_sup <- rbind(
  DMx = H_q(Z_DM_x), DPx = H_q(Z_DP_x), PYx = H_q(Z_PY_x), GNx = H_q(Z_GN_x)
)

## ---- Point estimates, credible balls, misclassification ----
analyze_model <- function(Zmat) {
  Zs <- Zmat[, (burn_in + 1):N_iter, drop = FALSE]
  C  <- pr_cc(Zs)
  mv <- minVI(C, method = "avg", max.k = 20)
  cl <- mv$cl
  list(
    cl = cl,
    credball = credibleball(mv$cl, t(Zs))[[5]],
    misclass = misclass(cl, Y, a = a, b = b)
  )
}

res_DM   <- analyze_model(Z_DM)
res_DM_x <- analyze_model(Z_DM_x)
res_DP   <- analyze_model(Z_DP)
res_DP_x <- analyze_model(Z_DP_x)
res_PY   <- analyze_model(Z_PY)
res_PY_x <- analyze_model(Z_PY_x)
res_GN   <- analyze_model(Z_GN)
res_GN_x <- analyze_model(Z_GN_x)

## ---- Competitors (Table 3) ----
# True edge probs used to simulate the scenario (block SBM with 5 groups)
pi_true <- matrix(0, V, V)
block_matrix <- matrix(0.25, 5, 5); diag(block_matrix) <- 0.75
for (v in 2:V) for (u in 1:(v - 1)) pi_true[v, u] <- block_matrix[z_0[v], z_0[u]]
pi_true <- pi_true + t(pi_true)

lowerTriangle <- function(M) M[lower.tri(M, diag = FALSE)]

edge_err <- function(cl) {
  mean(abs(lowerTriangle(edge_est(cl, Y, a = a, b = b)) - lowerTriangle(pi_true)))
}

# Louvain
net <- graph.adjacency(Y, mode = "undirected", weighted = NULL, diag = FALSE)
Louv <- cluster_louvain(net)$membership
louvain_stats <- c(H = length(table(Louv)), VI = VI(z_0, t(Louv)), ERR = edge_err(Louv))

# H selection for spectral (randnet)
set.seed(1)
H_select <- integer(8)
H_select[1] <- BHMC.estimate(Y, K.max = 10)$K
lr <- LRBIC(Y, Kmax = 10); H_select[2] <- lr$SBM.K
ncv <- NCV.select(Y, max.K = 10); H_select[3] <- which.min(ncv$l2); H_select[4] <- which.min(ncv$dev)
ecv <- ECV.block(Y, max.K = 10); H_select[5] <- which.min(ecv$l2); H_select[6] <- which.min(ecv$dev)
er  <- ECV.Rank(Y, 10, weighted = FALSE, mode = "undirected")
H_select[7] <- er$sse.rank; H_select[8] <- er$auc.rank
sel_H <- round(median(H_select))

# Spectral + Regularized Spectral
sc   <- reg.SP(Y, K = sel_H, lap = TRUE, tau = 0)$cluster
r_sc <- reg.SP(Y, K = sel_H, lap = TRUE, tau = 1)$cluster

spectral_stats    <- c(H = length(table(sc)),   VI = VI(z_0, t(sc)),   ERR = edge_err(sc))
regspectral_stats <- c(H = length(table(r_sc)), VI = VI(z_0, t(r_sc)), ERR = edge_err(r_sc))

# GREED: SBM and DC-SBM
set.seed(1)
g_sbm   <- greed(Y, K = sel_H, model = new("sbm",   alpha = beta_dm, type = "undirected"),
                 alg = methods::new("hybrid"), verbose = FALSE)@cl
g_dcsbm <- greed(Y, K = sel_H, model = new("dcsbm", alpha = beta_dm, type = "undirected"),
                 alg = methods::new("hybrid"), verbose = FALSE)@cl

greed_sbm_stats   <- c(H = length(table(g_sbm)),   VI = VI(z_0, t(g_sbm)),   ERR = edge_err(g_sbm))
greed_dcsbm_stats <- c(H = length(table(g_dcsbm)), VI = VI(z_0, t(g_dcsbm)), ERR = edge_err(g_dcsbm))

# JCDC (attribute-assisted)
A <- Y; K <- sel_H; N <- V
phi <- outer(z_0, z_0, function(a, b) as.integer(a == b))
Rcpp::sourceCpp(here("Source", "JCDC.cpp"))

D.inv <- diag(1 / (sqrt(rowSums(A)) + 1e-7))
Lapl  <- D.inv %*% A %*% D.inv
U.K   <- svd(Lapl)$u[, 1:K, drop = FALSE]
spec.cluster <- kmeans(U.K, K, nstart = 10)$cluster
G.fit <- array(0, c(N, K)); for (k in 1:K) G.fit[spec.cluster == k, k] <- 1

run_jcdc <- function(W_max) {
  p <- 1L
  res <- JCDC(A, phi, p, G.fit, 1, K, 20, 30, 20, W_max, 2, 0, 1)
  cl  <- as.vector(res$G.fit %*% (1:K))
  c(H = length(table(cl)), VI = VI(z_0, t(cl)), ERR = edge_err(cl))
}
jcdc_5_stats   <- run_jcdc(5)
jcdc_1_5_stats <- run_jcdc(1.5)

## ---- Predictive performance (GN supervised) ----
set.seed(1)
V_new <- 300
memb_new <- rep(1:6, each = 50)

Y_new <- matrix(0, V_new, V)
block_matrix_pred <- matrix(0, 6, 6); block_matrix_pred[1:5, 1:5] <- block_matrix
block_matrix_pred[6, ] <- block_matrix_pred[, 6] <- 0.05

for (v in 1:V_new) for (u in 1:V) {
  Y_new[v, u] <- rbinom(1, 1, prob = block_matrix_pred[memb_new[v], z_0[u]])
}

# Point estimate for GN (unsupervised) to plug into prediction
c_Z_GN   <- pr_cc(Z_GN[, (burn_in + 1):N_iter, drop = FALSE])
memb_GN  <- minVI(c_Z_GN, method = "avg", max.k = 20)$cl

Y_aug <- matrix(0, V + 1, V + 1); Y_aug[1:V, 1:V] <- Y
Post_Prob <- matrix(0, V_new, 6)

for (v in 1:V_new) {
  Y_aug[V + 1, 1:V] <- Y_aug[1:V, V + 1] <- Y_new[v, ]
  Post_Prob[v, ] <- pred_esbm(Y_aug, prior = "GN", z_hat = memb_GN, a = a, b = b, gamma_GN = gamma_gn)
}
pred_err <- (V_new - sum(diag(table(memb_new, apply(Post_Prob, 1, which.max))))) / V_new

## ---- Plots (Figure-style quick reproduction) ----
# Figure 3 (matrix)
row_plot_Y <- data.frame(z_0 = as.factor(z_0))
rownames(Y) <- rownames(row_plot_Y) <- paste0("v", seq_len(V))
mycolors <- c(brewer.pal(10, "RdBu")[c(4, 7)],
              brewer.pal(10, "PRGn")[c(7, 4)],
              brewer.pal(9, "YlOrBr")[3])
names(mycolors) <- levels(row_plot_Y$z_0)
Network <- pheatmap(
  Y,
  color = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30),
  cluster_cols = FALSE, cluster_rows = FALSE,
  annotation_row = row_plot_Y, annotation_names_row = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  legend = FALSE, border_color = FALSE,
  annotation_legend = FALSE,
  annotation_colors = list(z_0 = mycolors)
)
g <- grid.arrange(Network[[4]], nrow = 1, vp = viewport(width = 1, height = 1))
print(cowplot::ggdraw(g) + theme(plot.background = element_rect(fill = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30)[8])))

# Figure 4 (first column)
c_Z_GN <- pr_cc(Z_GN[, (burn_in + 1):N_iter, drop = FALSE])
memb_Z_GN <- minVI(c_Z_GN, method = "avg", max.k = 20)$cl
row_plot_GN <- data.frame(memb_Z_GN = as.factor(memb_Z_GN))
rownames(c_Z_GN) <- rownames(row_plot_GN) <- paste0("v", seq_len(V))
mycolors2 <- mycolors[seq_along(levels(row_plot_GN$memb_Z_GN))]
Marg <- pheatmap(
  c_Z_GN,
  color = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30),
  cluster_cols = FALSE, cluster_rows = FALSE,
  annotation_row = row_plot_GN, annotation_names_row = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  legend = FALSE, border_color = FALSE,
  annotation_legend = FALSE,
  annotation_colors = list(memb_Z_GN = mycolors2)
)
g <- grid.arrange(Marg[[4]], ncol = 1, vp = viewport(width = 1, height = 1))
print(cowplot::ggdraw(g) + theme(plot.background = element_rect(fill = colorRampPalette(brewer.pal(9, "Greys")[c(1, 8)])(30)[8])))

## ---- Console summary (quick glance) ----
cat("\n=== WAIC (unsupervised) ===\n")
print(sapply(list(DM = Z_DM_WAIC, DP = Z_DP_WAIC, PY = Z_PY_WAIC, GN = Z_GN_WAIC), `[[`, "WAIC"))
cat("\n=== WAIC (supervised) ===\n")
print(sapply(list(DMx = Z_DM_WAIC_x, DPx = Z_DP_WAIC_x, PYx = Z_PY_WAIC_x, GNx = Z_GN_WAIC_x), `[[`, "WAIC"))
cat("\n=== H quantiles (unsup) ===\n"); print(Hq_unsup)
cat("\n=== H quantiles (sup) ===\n");  print(Hq_sup)
cat("\n=== Competitors ===\n")
print(rbind(
  Louvain = louvain_stats,
  Spectral = spectral_stats,
  RegSpectral = regspectral_stats,
  Greed_SBM = greed_sbm_stats,
  Greed_DCSBM = greed_dcsbm_stats,
  JCDC_W5 = jcdc_5_stats,
  JCDC_W1_5 = jcdc_1_5_stats
))
cat("\nPredictive misclassification (GN supervised framework): ", round(pred_err, 4), "\n")