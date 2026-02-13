suppressPackageStartupMessages({
  library(here)
  library(R.utils)
  library(ggplot2)
  library(cowplot)
  library(scales)   # for hue_pal()
  library(RColorBrewer)
})

args <- R.utils::commandArgs(asValue = TRUE)

# args <-list(alpha_x= 0,
#   alpha_y= 0,
#   alpha_z= 0, 
#   show_edges = TRUE)

## ---------------------------
## Parse args
## ---------------------------
alpha_x <- if (is.null(args$alpha_x)) stop("Provide alpha_x=") else as.numeric(args$alpha_x)
alpha_y <- if (is.null(args$alpha_y)) stop("Provide alpha_y=") else as.numeric(args$alpha_y)
alpha_z <- if (is.null(args$alpha_z)) stop("Provide alpha_z=") else as.numeric(args$alpha_z)

label <- if (is.null(args$label)) "run" else as.character(args$label)
edges_flag <- if (is.null(args$show_edges)) FALSE else as.logical(as.integer(args$show_edges))

## ---------------------------
## Paths
## ---------------------------
data_path <- here("application_brain", "consensus_scale33_tau50.RData")
proc_path <- here("application_brain", "results", "processed",
                  "brain_xyz_processed_results.rds")

outDir <- here("application_brain", "figures",
               sprintf("%s_ax%0.2f_ay%0.2f_az%0.2f", label, alpha_x, alpha_y, alpha_z))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

cat("Posterior plotting:\n")
cat("  alpha_x =", alpha_x, "\n")
cat("  alpha_y =", alpha_y, "\n")
cat("  alpha_z =", alpha_z, "\n")
cat("  label   =", label, "\n")
cat("  outDir  =", outDir, "\n")
cat("  show_edges =", edges_flag, "\n\n")

## ---------------------------
## Load utility plotting function
## ---------------------------
source(here("R_utilties", "plottingFunctions.R"))

## ---------------------------
## Load processed results
## ---------------------------
res_app <- readRDS(proc_path)

find_tag_by_alpha <- function(obj, ax, ay, az, tol = 1e-8) {
  a_mat <- t(sapply(obj$results, function(r) r$alpha[1:3]))
  idx <- which(abs(a_mat[,1]-ax) < tol & abs(a_mat[,2]-ay) < tol & abs(a_mat[,3]-az) < tol)
  if (length(idx) == 0) stop(sprintf("No run found for alpha=(%.2f,%.2f,%.2f).", ax, ay, az))
  names(obj$results)[idx[1]]
}

tag <- find_tag_by_alpha(res_app, alpha_x, alpha_y, alpha_z)
rr  <- res_app$results[[tag]]

psm_mat <- rr$psm
z_hat   <- as.integer(rr$z_hat)
K_hat   <- rr$K_hat

cat("Matched tag:", tag, "\n")
cat("K_hat =", K_hat, "\n\n")

## ---------------------------
## Load data + drop isolates (must match fitting)
## ---------------------------
load(data_path)  # A_bin_tau, coords
diag(A_bin_tau) <- 0
deg  <- rowSums(A_bin_tau)
keep <- deg > 0

A <- A_bin_tau[keep, keep, drop = FALSE]
diag(A) <- 0
coords_kept <- coords[keep, , drop = FALSE]

stopifnot(nrow(A) == length(z_hat))

## ---------------------------
## Hemisphere labels (for annotation bar)
## ---------------------------
hemi <- ifelse(grepl("^lh\\.", coords_kept$node), "Left",
        ifelse(grepl("^rh\\.", coords_kept$node), "Right",
        ifelse(grepl("^Left-",  coords_kept$node), "Left",
        ifelse(grepl("^Right-", coords_kept$node), "Right", "Other"))))

## ---------------------------
## ORDERING FOR MATRICES: cluster first, within cluster keep hemisphere separated
## ---------------------------
hemi_fac <- factor(hemi, levels = c("Left","Right","Other"))

ord_clust <- order(
  factor(z_hat),      # cluster blocks
  hemi_fac,           # hemisphere within cluster
  coords_kept$node     # stable tie-breaker
)

A_ord2   <- A[ord_clust, ord_clust, drop = FALSE]
psm_ord  <- psm_mat[ord_clust, ord_clust, drop = FALSE]
z_ord    <- z_hat[ord_clust]
hemi_ord <- hemi[ord_clust]

## ---------------------------
## Colors: keep "original" discrete palette (ggplot2 hue),
## and reuse it BOTH for matrix est bar and 3-view scatter.
## ---------------------------
lev_est <- sort(unique(as.integer(z_ord)))
n_cl <- length(lev_est)

library(colorspace)

est_colors <- setNames(
  colorspace::qualitative_hcl(12, palette = "Dark 3"),
  as.character(sort(unique(z_hat)))
)

## ---------------------------
## Save PSM (cluster-ordered)
## ---------------------------
title_psm <- sprintf("PSM, alpha=(%.2f,%.2f,%.2f), K_hat=%d",
                     alpha_x, alpha_y, alpha_z, K_hat)

plot_matrix_with_annotations(
  M           = psm_ord,
  type        = "psm",
  mat_palette = "Blues",
  n_mat_cols  = 120,
  ann_df      = data.frame(
    estimate        = factor(z_ord, levels = lev_est),
    hemisphere = factor(hemi_ord, levels = c("Left","Right","Other"))
  ),
  show_cols   = c("estimate", "hemisphere"),
  ann_colors  = list(
    estimate        = est_colors,
    hemisphere = c(Left = "#4D4D4D", Right = "#A6A6A6", Other = "#D9D9D9")
  ),
  gaps_by     = z_ord,
  main        = title_psm,
  filename    = file.path(outDir, "psm_clusterOrder.pdf")
)

## ---------------------------
## Save adjacency (same order; hemisphere bar first, cluster second)
## ---------------------------
title_adj <- sprintf("Adjacency, alpha=(%.2f,%.2f,%.2f)",
                     alpha_x, alpha_y, alpha_z)

plot_matrix_with_annotations(
  M          = A_ord2,
  type       = "adj",
  ann_df     = data.frame(
    hemisphere = factor(hemi_ord, levels = c("Left","Right","Other")),
    estimate        = factor(z_ord, levels = lev_est)
  ),
  show_cols   = c("estimate", "hemisphere"),
  ann_colors = list(
    hemisphere = c(Left = "#4D4D4D", Right = "#A6A6A6", Other = "#D9D9D9"),
    estimate        = est_colors
  ),
  gaps_by    = z_ord,
  main       = title_adj,
  filename   = file.path(outDir, "adjacency_clusterOrder.pdf")
)

## ---------------------------
## 3-view cluster plot helper
## ---------------------------
# Edge list for optional edges in 3-view plot (in original kept-node order)
edges <- which(A == 1, arr.ind = TRUE)
edges <- edges[edges[, 1] < edges[, 2], , drop = FALSE]

edges_df <- data.frame(
  x1 = coords_kept$x[edges[, 1]],
  y1 = coords_kept$y[edges[, 1]],
  z1 = coords_kept$z[edges[, 1]],
  x2 = coords_kept$x[edges[, 2]],
  y2 = coords_kept$y[edges[, 2]],
  z2 = coords_kept$z[edges[, 2]]
)

make_view_plot <- function(nodes_df, edges_df, view = c("xy","xz","yz"),
                           show_edges = TRUE,
                           edge_alpha = 0.06, edge_lwd = 0.25,
                           node_alpha = 0.95, node_size = 2.0) {

  view <- match.arg(view)

  if (view == "xy") {
    xvar <- "x"; yvar <- "y"
    ex1 <- "x1"; ey1 <- "y1"; ex2 <- "x2"; ey2 <- "y2"
    ttl <- "Axial"
    xlab <- "x1"; ylab <- "x2"
  } else if (view == "xz") {
    xvar <- "x"; yvar <- "z"
    ex1 <- "x1"; ey1 <- "z1"; ex2 <- "x2"; ey2 <- "z2"
    ttl <- "Coronal"
    xlab <- "x1"; ylab <- "x3"
  } else {
    xvar <- "y"; yvar <- "z"
    ex1 <- "y1"; ey1 <- "z1"; ex2 <- "y2"; ey2 <- "z2"
    ttl <- "Sagittal"
    xlab <- "x2"; ylab <- "x3"
  }

  p <- ggplot()

  if (show_edges && nrow(edges_df) > 0) {
    p <- p + geom_segment(
      data = edges_df,
      aes(x = .data[[ex1]], y = .data[[ey1]],
          xend = .data[[ex2]], yend = .data[[ey2]]),
      alpha = edge_alpha,
      linewidth = edge_lwd
    )
  }

  p +
    geom_point(
      data = nodes_df,
      aes(x = .data[[xvar]], y = .data[[yvar]], color = cluster),
      alpha = node_alpha,
      size = node_size
    ) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    labs(title = ttl, x = xlab, y = ylab, color = "Cluster") +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

## ---------------------------
## Save 3-view cluster plot (kept-node geometry; colored by z_hat)
## One shared legend across the 3 plots.
## ---------------------------
nodes_df <- coords_kept
nodes_df$cluster <- factor(z_hat, levels = lev_est)

lev_est_chr <- levels(nodes_df$cluster)
# ensure est_colors is NAMED with those exact level strings
est_colors <- est_colors[lev_est_chr]
names(est_colors) <- lev_est_chr

p_xy <- make_view_plot(nodes_df, edges_df, "xy", show_edges = edges_flag) + scale_color_manual(values = est_colors)

p_xz <- make_view_plot(nodes_df, edges_df, "xz", show_edges = edges_flag) +
  scale_color_manual(values = est_colors)

p_yz <- make_view_plot(nodes_df, edges_df, "yz", show_edges = edges_flag) +
  scale_color_manual(values = est_colors)

# Extract one legend and remove legends from panels
leg <- cowplot::get_legend(p_xy + theme(legend.position = "right"))

p_xy_nl <- p_xy + theme(legend.position = "none")
p_xz_nl <- p_xz + theme(legend.position = "none")
p_yz_nl <- p_yz + theme(legend.position = "none")

# legend + panels
leg <- cowplot::get_legend(p_xy + theme(legend.position = "right"))

p_final <- cowplot::plot_grid(
  cowplot::plot_grid(p_xy_nl, p_xz_nl, p_yz_nl, nrow = 1, align = "none"),
  leg,
  nrow = 1,
  rel_widths = c(1, 0.12)
)

ggsave(file.path(outDir, "clusters_3views.pdf"), p_final, width = 12.8, height = 4.4, dpi = 300)