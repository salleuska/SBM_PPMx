#------------------------------------------------------------------#
# Figures describing the simulated data.
#
# Reads the saved .rds files: it does NOT regenerate them.
# (Sourcing 0_simulateBinaryESBM.R would redraw the data and change
#  every downstream result, so keep the two scripts separate.)
#
# Produces, in simulation/data/plots/ :
#   - network_adjacency.pdf   adjacency matrix ordered by true cluster
#   - covariate_by_cluster.pdf covariate distribution per true cluster,
#                              for the neutral / informative / misleading
#                              one-covariate scenarios
#------------------------------------------------------------------#
library(here)
library(ggplot2)
library(RColorBrewer)

outDir <- here("simulation", "data", "plots")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

## The network Y and the true partition are shared across every scenario
## (0_simulateBinaryESBM.R draws one `base` object and only varies the
## covariates), so any of the files gives the same network.
sim_N <- readRDS(here("simulation", "data", "binarySBM_1cov_N.rds"))
sim_I <- readRDS(here("simulation", "data", "binarySBM_1cov_I.rds"))
sim_M <- readRDS(here("simulation", "data", "binarySBM_1cov_M.rds"))

z0 <- sim_N$partition
V  <- length(z0)

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#D55E00",
               "#CC79A7", "#0072B2", "#F0E442", "#999999")

#------------------------------------------------------------------#
# 1) Adjacency matrix, ordered by true cluster
#------------------------------------------------------------------#
Y <- sim_N$Y
diag(Y) <- 0

ord   <- order(z0)
Y_ord <- Y[ord, ord]
z_ord <- z0[ord]

grid_idx <- expand.grid(row = seq_len(V), col = seq_len(V))
grid_idx$value <- as.vector(Y_ord)
grid_idx$row_f <- factor(grid_idx$row, levels = rev(seq_len(V)))
grid_idx$col_f <- factor(grid_idx$col, levels = seq_len(V))

adj_pal <- colorRampPalette(brewer.pal(9, "Greys")[c(1, 9)])(50)

## cluster boundaries, for guide lines
brk <- cumsum(table(z_ord))
brk <- brk[-length(brk)]

p_adj <- ggplot(grid_idx, aes(x = col_f, y = row_f, fill = value)) +
  geom_tile() +
  geom_vline(xintercept = brk + 0.5, colour = "#D55E00", linewidth = 0.4) +
  geom_hline(yintercept = V - brk + 0.5, colour = "#D55E00", linewidth = 0.4) +
  scale_fill_gradientn(colors = adj_pal, guide = "none") +
  coord_fixed() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text  = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title    = "Simulated network, nodes ordered by true cluster",
    subtitle = sprintf("V = %d,  K = %d,  block sizes %s",
                       V, length(unique(z0)),
                       paste(as.integer(table(z0)), collapse = "/"))
  )

ggsave(file.path(outDir, "network_adjacency.pdf"), p_adj, width = 6, height = 6.4)

#------------------------------------------------------------------#
# 2) Covariate distribution by true cluster, per scenario
#------------------------------------------------------------------#
cov_df <- rbind(
  data.frame(scenario = "neutral",     x = sim_N$x[, 1], cluster = factor(z0)),
  data.frame(scenario = "informative", x = sim_I$x[, 1], cluster = factor(z0)),
  data.frame(scenario = "misleading",  x = sim_M$x[, 1], cluster = factor(z0))
)

cov_df$scenario <- factor(
  cov_df$scenario,
  levels = c("neutral", "informative", "misleading")
)

p_cov <- ggplot(cov_df, aes(x = x, fill = cluster, colour = cluster)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  facet_wrap(~ scenario, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = okabe_ito) +
  scale_colour_manual(values = okabe_ito) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    strip.text   = element_text(face = "bold")
  ) +
  labs(
    title    = "Covariate distribution within each true cluster",
    x        = expression(x[1]),
    y        = "Density",
    fill     = "True cluster",
    colour   = "True cluster"
  )

ggsave(file.path(outDir, "covariate_by_cluster.pdf"), p_cov, width = 7, height = 7)

#------------------------------------------------------------------#
# 3) Node-link view of the same network, nodes coloured by true cluster
#    Dense graph, so edges are drawn very light.
#------------------------------------------------------------------#
library(igraph)

g <- igraph::graph_from_adjacency_matrix(Y, mode = "undirected", diag = FALSE)

cat("Network has", igraph::gorder(g), "nodes and",
    igraph::gsize(g), "edges\n")

set.seed(1)
lay <- igraph::layout_with_fr(g)

pdf(file.path(outDir, "network_nodelink.pdf"), width = 7, height = 7)
par(mar = c(0, 0, 3, 0))
plot(
  g,
  layout          = lay,
  vertex.color    = okabe_ito[z0],
  vertex.frame.color = "grey20",
  vertex.size     = 5,
  vertex.label    = NA,
  edge.color      = adjustcolor("grey50", alpha.f = 0.06),
  edge.width      = 0.4,
  main            = "Simulated network, nodes coloured by true cluster"
)
legend(
  "bottomleft",
  legend = paste("cluster", sort(unique(z0))),
  col    = okabe_ito[sort(unique(z0))],
  pch    = 19,
  bty    = "n",
  cex    = 0.9
)
invisible(dev.off())

#------------------------------------------------------------------#
# 4) Network drawn with the COVARIATE as the layout coordinate.
#    Same network and same true clusters in all three panels; only the
#    horizontal position changes, because that is the covariate.
#    This is the figure that shows network and covariate together:
#      informative -> connected nodes are also co-located
#      neutral     -> position carries no cluster information
#      misleading  -> position is well separated but wrong
#------------------------------------------------------------------#

## edge list, taken once (the network is shared across scenarios)
ends <- which(Y == 1 & upper.tri(Y), arr.ind = TRUE)

## one fixed vertical jitter, reused in every panel so that the only
## thing differing between panels is the covariate on the x axis
set.seed(7)
y_jit <- runif(V, -1, 1)

geo_panel <- function(sim_obj, label) {
  xv <- sim_obj$x[, 1]

  nodes <- data.frame(
    x = xv, y = y_jit, cluster = factor(z0), scenario = label
  )

  edges <- data.frame(
    x        = xv[ends[, 1]],
    y        = y_jit[ends[, 1]],
    xend     = xv[ends[, 2]],
    yend     = y_jit[ends[, 2]],
    scenario = label
  )

  list(nodes = nodes, edges = edges)
}

panels <- list(
  geo_panel(sim_N, "neutral"),
  geo_panel(sim_I, "informative"),
  geo_panel(sim_M, "misleading")
)

nodes_df <- do.call(rbind, lapply(panels, `[[`, "nodes"))
edges_df <- do.call(rbind, lapply(panels, `[[`, "edges"))

lev <- c("neutral", "informative", "misleading")
nodes_df$scenario <- factor(nodes_df$scenario, levels = lev)
edges_df$scenario <- factor(edges_df$scenario, levels = lev)

p_geo <- ggplot() +
  geom_segment(
    data = edges_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    colour = "grey50", alpha = 0.05, linewidth = 0.2
  ) +
  geom_point(
    data = nodes_df,
    aes(x = x, y = y, fill = cluster),
    shape = 21, size = 2.6, stroke = 0.3, colour = "grey20"
  ) +
  facet_wrap(~ scenario, ncol = 1) +
  scale_fill_manual(values = okabe_ito) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y  = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title   = element_text(face = "bold"),
    strip.text   = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title    = "Network positioned by the covariate",
    subtitle = "Same network and same true clusters in every panel; only x changes",
    x        = expression(x[1]),
    y        = NULL,
    fill     = "True cluster"
  )

ggsave(file.path(outDir, "network_by_covariate.pdf"), p_geo,
       width = 7.5, height = 8.5)

cat("Saved to:", outDir, "\n")
