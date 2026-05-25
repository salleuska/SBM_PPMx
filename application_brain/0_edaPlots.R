# ============================================================ 
# Explorative data analysis plots
# ============================================================
library(Rcpp)
library(here)
library(igraph)
library(tools)
library(ggplot2)
library(cowplot)   # for arranging plots

# Load utility plotting function
source(here("R_utilties", "plottingFunctions.R"))
# ---------------------------
# 1) load data, remove nodes with no connections
# ---------------------------
# load data
load(here("application_brain", "consensus_scale33_tau50.RData"))

deg <- rowSums(A_bin_tau)
keep <- deg > 0

message(sprintf("Dropping %d isolated nodes out of %d total.",
                sum(!keep), length(keep)))

A <- A_bin_tau[keep, keep, drop = FALSE]
coords_kept <- coords[keep, , drop = FALSE]

# ---------------------------
# Decide node order for network plots
# ---------------------------
# hemisfere labels
hemi <- ifelse(grepl("^lh\\.", coords_kept$node), "Left",
        ifelse(grepl("^rh\\.", coords_kept$node), "Right",
        ifelse(grepl("^Left-",  coords_kept$node), "Left",
        ifelse(grepl("^Right-", coords_kept$node), "Right", "Other"))))

## if you want to order using labels
# ord_fixed <- order(hemi, coords_kept$node)
# hemi_ord <- hemi[ord_fixed]

## ordered using the x coordinate
ord_fixed <- order(coords_kept$x)

A_ord <- A[ord_fixed, ord_fixed, drop = FALSE]
coords_ord <- coords_kept[ord_fixed, ]

# ---------------------------
# 2) Edge list (upper triangle)
# ---------------------------

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


# ---------------------------
# 3) Helper to make a view plot
# ---------------------------
make_view_plot <- function(nodes_df, edges_df, view = c("xy", "xz", "yz"),
                           show_edges = TRUE,
                           edge_alpha = 0.10, edge_lwd = 0.3,
                           node_alpha = 0.90, node_size = 1.8) {

  view <- match.arg(view)
 
  if (view == "xy") {
    xvar <- "x"; yvar <- "y"
    ex1 <- "x1"; ey1 <- "y1"; ex2 <- "x2"; ey2 <- "y2"
    ttl <- "Axial";    xlab <- "x1"; ylab <- "x2"
  } else if (view == "xz") {
    xvar <- "x"; yvar <- "z"
    ex1 <- "x1"; ey1 <- "z1"; ex2 <- "x2"; ey2 <- "z2"
    ttl <- "Coronal";  xlab <- "x1"; ylab <- "x3"
  } else { # "yz"
    xvar <- "y"; yvar <- "z"
    ex1 <- "y1"; ey1 <- "z1"; ex2 <- "y2"; ey2 <- "z2"
    ttl <- "Sagittal"; xlab <- "x2"; ylab <- "x3"
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

  p <- p +
    geom_point(
      data = nodes_df,
      aes(x = .data[[xvar]], y = .data[[yvar]]),
      alpha = node_alpha,
      size = node_size
    ) +
    coord_equal() +
    labs(title = ttl, x = xlab, y = ylab) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  p
}

# ---------------------------
# 4) Make  plots
# ---------------------------

p_axial <- make_view_plot(coords_kept, edges_df, "xy")
p_coronal <- make_view_plot(coords_kept, edges_df, "xz")
p_sagittal <- make_view_plot(coords_kept, edges_df, "yz")


zlim <- range(coords_kept$z, na.rm = TRUE)

p_coronal  <- p_coronal  + coord_equal(ylim = zlim)
p_sagittal <- p_sagittal + coord_equal(ylim = zlim)

xlim <- range(coords_kept$x, na.rm = TRUE)

p_axial   <- p_axial   + coord_equal(xlim = xlim)
p_coronal <- p_coronal + coord_equal(xlim = xlim)  # keeps the zlim you set above too


p_all <- cowplot::plot_grid(p_axial, p_coronal, p_sagittal, 
  nrow = 1, align = "none",   axis = "tblr")

# p_axial   <- make_view_plot(coords_kept, edges_df, view = "xy")
# p_coronal <- make_view_plot(coords_kept, edges_df, view = "xz")
# p_sagittal <- make_view_plot(coords_kept, edges_df, view = "yz")
# p_all <- cowplot::plot_grid(
#   p_axial,no
#   p_coronal,
#   p_sagittal,
#   nrow = 1
# )

# ---------------------------
# 5) Optional: save to file
# ---------------------------
# Create figures directory if needed

fig_dir <- here("application_brain", "figures")
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

ggsave(here(fig_dir, "raw_sagittal.pdf"), p_sagittal, width=5, height=5, dpi=300)
ggsave(here(fig_dir, "raw_coronal.pdf"), p_coronal, width=5, height=5, dpi=300)
ggsave(here(fig_dir, "raw_axial.pdf"), p_axial, width=5, height=5, dpi=300)


# Save the 3-view plot
ggsave(
  filename = here("application_brain", "figures", "brain_network_3views.pdf"),
  plot = p_all,
  width = 12,
  height = 4,
  dpi = 300
)

# ---------------------------
# 6) Adjacency matrix visualization (ggplot)
# ---------------------------

plot_matrix_with_annotations(
  M        = A_ord,
  type     = "adj",
  filename = here("application_brain", "figures", "raw_adjacency.pdf")
)

########
hemi_ord <- hemi[ord_fixed]

ann_df <- data.frame(
  hemisphere = factor(hemi_ord, levels = c("Right","Other","Left"))
)


plot_matrix_with_annotations(
  M        = A_ord,
  type     = "adj",
  ann_df   = ann_df, 
  show_cols = "hemisphere",
  ann_colors = list(
    hemisphere = c(Left = "#4D4D4D",  Other = "#D9D9D9", Right = "#A6A6A6")),
  gaps_by  = hemi_ord,
  filename = here("application_brain", "figures", "adjacency_left_right_order.pdf")
)
