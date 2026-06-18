suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(cowplot)
  library(colorspace)
  library(RColorBrewer)
  library(tidyr)
  library(dplyr)
})

## ---------------------------
## Paths
## ---------------------------
data_path <- here("application_brain", "consensus_scale33_tau50.RData")
proc_path <- here("application_brain", "results", "processed",
                  "brain_xyz_learnAlpha_processed_results.rds")

outDir <- here("application_brain", "figures", "learnAlpha")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------
## Load utilities
## ---------------------------
source(here("R_utilties", "plottingFunctions.R"))

## ---------------------------
## Load processed results
## ---------------------------
res_app <- readRDS(proc_path)
results  <- res_app$results

cat("Found", length(results), "learned-alpha run(s) to plot.\n\n")

## ---------------------------
## Load data + drop isolates (must match fitting)
## ---------------------------
load(data_path)   # A_bin_tau, coords
diag(A_bin_tau) <- 0

if (!is.null(rownames(A_bin_tau))) {
  m <- match(rownames(A_bin_tau), coords$node)
  stopifnot(!anyNA(m))
  coords <- coords[m, , drop = FALSE]
}

deg  <- rowSums(A_bin_tau)
keep <- deg > 0

A           <- A_bin_tau[keep, keep, drop = FALSE]
diag(A)     <- 0
coords_kept <- coords[keep, , drop = FALSE]

hemi <- ifelse(grepl("^lh\\.",    coords_kept$node), "Left",
        ifelse(grepl("^rh\\.",    coords_kept$node), "Right",
        ifelse(grepl("^Left-",    coords_kept$node), "Left",
        ifelse(grepl("^Right-",   coords_kept$node), "Right", "Other"))))

edges     <- which(A == 1, arr.ind = TRUE)
edges     <- edges[edges[, 1] < edges[, 2], , drop = FALSE]
edges_df  <- data.frame(
  x1 = coords_kept$x[edges[, 1]], y1 = coords_kept$y[edges[, 1]],
  z1 = coords_kept$z[edges[, 1]],
  x2 = coords_kept$x[edges[, 2]], y2 = coords_kept$y[edges[, 2]],
  z2 = coords_kept$z[edges[, 2]]
)

## ---------------------------
## Helper: 3-view brain scatter
## ---------------------------
make_view_plot <- function(nodes_df, edges_df, view = c("xy", "xz", "yz"),
                           show_edges = FALSE,
                           edge_alpha = 0.06, edge_lwd = 0.25,
                           node_alpha = 0.95, node_size = 2.0) {
  view <- match.arg(view)
  if (view == "xy") {
    xvar <- "x"; yvar <- "y"; ex1 <- "x1"; ey1 <- "y1"; ex2 <- "x2"; ey2 <- "y2"
    ttl <- "Axial"; xlab <- "x1"; ylab <- "x2"
  } else if (view == "xz") {
    xvar <- "x"; yvar <- "z"; ex1 <- "x1"; ey1 <- "z1"; ex2 <- "x2"; ey2 <- "z2"
    ttl <- "Coronal"; xlab <- "x1"; ylab <- "x3"
  } else {
    xvar <- "y"; yvar <- "z"; ex1 <- "y1"; ey1 <- "z1"; ex2 <- "y2"; ey2 <- "z2"
    ttl <- "Sagittal"; xlab <- "x2"; ylab <- "x3"
  }

  p <- ggplot()
  if (show_edges && nrow(edges_df) > 0) {
    p <- p + geom_segment(
      data = edges_df,
      aes(x = .data[[ex1]], y = .data[[ey1]],
          xend = .data[[ex2]], yend = .data[[ey2]]),
      alpha = edge_alpha, linewidth = edge_lwd
    )
  }
  p +
    geom_point(data = nodes_df,
               aes(x = .data[[xvar]], y = .data[[yvar]], color = cluster),
               alpha = node_alpha, size = node_size) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    labs(title = ttl, x = xlab, y = ylab, color = "Cluster") +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

## ---------------------------
## Loop over runs
## ---------------------------
for (tag in names(results)) {

  rr <- results[[tag]]

  cat("Plotting:", tag, "\n")
  cat("  K_hat =", rr$K_hat, "\n")

  psm_mat <- rr$psm
  z_hat   <- as.integer(rr$z_hat)
  K_hat   <- rr$K_hat

  alpha_g_post  <- rr$alpha_g_post   # J x N_kept (after burn-in), or NULL
  alpha_summary <- rr$alpha_summary  # data.frame or NULL

  run_dir <- file.path(outDir, tag)
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  ## ---- cluster ordering ----
  hemi_fac <- factor(hemi, levels = c("Left", "Right", "Other"))
  ord_clust <- order(factor(z_hat), hemi_fac, coords_kept$node)

  A_ord   <- A[ord_clust, ord_clust, drop = FALSE]
  psm_ord <- psm_mat[ord_clust, ord_clust, drop = FALSE]
  z_ord   <- z_hat[ord_clust]
  hemi_ord <- hemi[ord_clust]

  lev_est    <- sort(unique(as.integer(z_ord)))
  est_colors <- setNames(
    colorspace::qualitative_hcl(max(12, length(lev_est)), palette = "Dark 3"),
    as.character(sort(unique(z_hat)))
  )

  alpha_label <- if (!is.null(alpha_summary)) {
    paste0("alpha_post=[",
           paste(sprintf("%.2f", alpha_summary$mean), collapse = ","),
           "]")
  } else {
    "alpha_post=NA"
  }

  ## ---- PSM ----
  plot_matrix_with_annotations(
    M           = psm_ord,
    type        = "psm",
    mat_palette = "Blues",
    n_mat_cols  = 120,
    ann_df      = data.frame(
      estimate   = factor(z_ord, levels = lev_est),
      hemisphere = factor(hemi_ord, levels = c("Left", "Right", "Other"))
    ),
    show_cols   = c("estimate", "hemisphere"),
    ann_colors  = list(
      estimate   = est_colors,
      hemisphere = c(Left = "#4D4D4D", Right = "#A6A6A6", Other = "#D9D9D9")
    ),
    gaps_by     = z_ord,
    main        = sprintf("PSM | K=%d | %s", K_hat, alpha_label),
    filename    = file.path(run_dir, "psm_clusterOrder.pdf")
  )

  ## ---- Adjacency ----
  plot_matrix_with_annotations(
    M          = A_ord,
    type       = "adj",
    ann_df     = data.frame(
      estimate   = factor(z_ord, levels = lev_est),
      hemisphere = factor(hemi_ord, levels = c("Left", "Right", "Other"))
    ),
    show_cols  = c("estimate", "hemisphere"),
    ann_colors = list(
      estimate   = est_colors,
      hemisphere = c(Left = "#4D4D4D", Right = "#A6A6A6", Other = "#D9D9D9")
    ),
    gaps_by    = z_ord,
    main       = sprintf("Adjacency | K=%d | %s", K_hat, alpha_label),
    filename   = file.path(run_dir, "adjacency_clusterOrder.pdf")
  )

  ## ---- 3-view brain plots ----
  nodes_df <- coords_kept
  nodes_df$cluster <- factor(z_hat, levels = sort(unique(z_hat)))
  lev_chr  <- levels(nodes_df$cluster)
  col_vec  <- est_colors[lev_chr]
  names(col_vec) <- lev_chr

  p_xy <- make_view_plot(nodes_df, edges_df, "xy") + scale_color_manual(values = col_vec)
  p_xz <- make_view_plot(nodes_df, edges_df, "xz") + scale_color_manual(values = col_vec)
  p_yz <- make_view_plot(nodes_df, edges_df, "yz") + scale_color_manual(values = col_vec)

  leg     <- cowplot::get_legend(p_xy + theme(legend.position = "right"))
  p_final <- cowplot::plot_grid(
    cowplot::plot_grid(
      p_xy + theme(legend.position = "none"),
      p_xz + theme(legend.position = "none"),
      p_yz + theme(legend.position = "none"),
      nrow = 1, align = "none"
    ),
    leg,
    nrow = 1, rel_widths = c(1, 0.12)
  )

  ggsave(file.path(run_dir, "brain_3view.pdf"), p_final,
         width = 14, height = 4.5)

  ## ---- Alpha posterior density plots ----
  if (!is.null(alpha_g_post)) {

    J          <- nrow(alpha_g_post)
    cov_names  <- c("x", "y", "z")[seq_len(J)]
    n_kept     <- ncol(alpha_g_post)

    # long-format data frame
    alpha_df <- data.frame(
      t(alpha_g_post),
      check.names = FALSE
    )
    colnames(alpha_df) <- cov_names
    alpha_df$iter <- seq_len(n_kept)

    alpha_long <- tidyr::pivot_longer(
      alpha_df,
      cols      = all_of(cov_names),
      names_to  = "covariate",
      values_to = "alpha"
    )
    alpha_long$covariate <- factor(alpha_long$covariate, levels = cov_names)

    ## density
    p_dens <- ggplot(alpha_long, aes(x = alpha, fill = covariate, color = covariate)) +
      geom_density(alpha = 0.35, linewidth = 0.7) +
      facet_wrap(~ covariate, scales = "free", nrow = 1) +
      theme_minimal(base_size = 13) +
      labs(title = "Posterior density of alpha",
           x = expression(alpha), y = "Density") +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold"))

    ggsave(file.path(run_dir, "alpha_density.pdf"), p_dens,
           width = 4 * J, height = 4)

    ## trace
    p_trace <- ggplot(alpha_long, aes(x = iter, y = alpha, color = covariate)) +
      geom_line(alpha = 0.6, linewidth = 0.4) +
      facet_wrap(~ covariate, scales = "free_y", nrow = J) +
      theme_minimal(base_size = 13) +
      labs(title = "Trace of alpha (post burn-in)",
           x = "Iteration", y = expression(alpha)) +
      theme(legend.position = "none",
            strip.text = element_text(face = "bold"))

    ggsave(file.path(run_dir, "alpha_trace.pdf"), p_trace,
           width = 10, height = 3 * J)

    ## summary table to console
    cat("  Alpha posterior summary:\n")
    print(alpha_summary, row.names = FALSE)
    cat("\n")
  }

  cat("  Saved to:", run_dir, "\n\n")
}

cat("Done. All plots saved under:\n  ", outDir, "\n")
