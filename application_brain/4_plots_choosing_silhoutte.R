library(here)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(Polychrome)

plot_matrix_with_clusters <- function(M, z_est,
                                      type = c("psm", "adj"),
                                      title = NULL) {

  type <- match.arg(type)

  require(pheatmap)
  require(RColorBrewer)

  n <- nrow(M)
  rn <- paste0("i", seq_len(n))
  rownames(M) <- rn
  colnames(M) <- rn

  ann <- data.frame(cluster = factor(z_est))
  rownames(ann) <- rn

  lev  <- levels(ann$cluster)
  cols <- palette36.colors(length(lev))
  names(cols) <- lev

  ord <- order(z_est)
  M_ord   <- M[ord, ord, drop = FALSE]
  ann_ord <- ann[ord, , drop = FALSE]
  gaps <- which(diff(as.integer(ann_ord$cluster)) != 0)

  if (type == "psm") {
    mat_cols <- colorRampPalette(brewer.pal(9, "Greys"))(80)
    breaks   <- NULL
  } else {
    mat_cols <- c("white", "black")
    breaks   <- c(-0.5, 0.5, 1.5)
  }

  pheatmap(
    M_ord,
    color             = mat_cols,
    breaks            = breaks,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    annotation_row    = ann_ord,
    annotation_col    = ann_ord,
    annotation_colors = list(cluster = cols),
    gaps_row          = gaps,
    gaps_col          = gaps,
    border_color      = NA,
    main              = title
  )
}

##--------------------------------------------
## Read processed application results (3-cov)
##--------------------------------------------
proc_path <- here("application_brain", "results", "processed",
                  "brain_xyz_processed_results.rds")
res_app <- readRDS(proc_path)

outDir <- here("application_brain", "figures")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

##--------------------------------------------
## Extract df: one row per alpha triple
##--------------------------------------------
df3 <- data.frame(
  tag      = names(res_app$results),
  alpha_x  = sapply(res_app$results, function(r) r$alpha[1]),
  alpha_y  = sapply(res_app$results, function(r) r$alpha[2]),
  alpha_z  = sapply(res_app$results, function(r) r$alpha[3]),
  K_hat    = sapply(res_app$results, function(r) r$K_hat),
  sil_post = sapply(res_app$results, function(r) unname(r$silhouettes["post"])),
  sil_net  = sapply(res_app$results, function(r) unname(r$silhouettes["net"])),
  sil_cov  = sapply(res_app$results, function(r) unname(r$silhouettes["cov"])),
  stringsAsFactors = FALSE
)

## factors for tiles
df3$alpha_x_f <- factor(df3$alpha_x, levels = sort(unique(df3$alpha_x)))
df3$alpha_y_f <- factor(df3$alpha_y, levels = sort(unique(df3$alpha_y)))
df3$alpha_z_f <- factor(df3$alpha_z, levels = sort(unique(df3$alpha_z)))

## palettes
blue_pal   <- colorRampPalette(brewer.pal(9, "Blues"))(50)
green_pal  <- colorRampPalette(brewer.pal(9, "Greens"))(50)
purple_pal <- colorRampPalette(brewer.pal(9, "Purples"))(50)

##--------------------------------------------
## Pick best alpha by each criterion
##--------------------------------------------
pick_best <- function(df, col) {
  if (all(is.na(df[[col]]))) stop("All ", col, " are NA; cannot pick best alpha.")
  df[which.max(df[[col]]), ]
}

best_post <- pick_best(df3, "sil_post")
best_net  <- pick_best(df3, "sil_net")
best_cov  <- pick_best(df3, "sil_cov")

best_all <- rbind(
  data.frame(criterion = "post", best_post, row.names = NULL),
  data.frame(criterion = "net",  best_net,  row.names = NULL),
  data.frame(criterion = "cov",  best_cov,  row.names = NULL)
)

write.csv(best_all, file.path(outDir, "best_alpha_by_each_silhouette.csv"), row.names = FALSE)

cat("\nBest alpha by silhouette criterion:\n")
print(best_all[, c("criterion","alpha_x","alpha_y","alpha_z","K_hat","sil_post","sil_net","sil_cov")])

##--------------------------------------------
## Heatmap helper
##--------------------------------------------
make_heatmap <- function(df, fill_var, title, fill_lab, pal) {
  ggplot(df, aes(x = alpha_x_f, y = alpha_y_f, fill = .data[[fill_var]])) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", .data[[fill_var]])), size = 3) +
    scale_fill_gradientn(colors = pal) +
    facet_wrap(~ alpha_z_f) +
    coord_equal() +
    theme_minimal(base_size = 13) +
    labs(
      title = title,
      x     = expression(alpha[x]),
      y     = expression(alpha[y]),
      fill  = fill_lab
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
}

##--------------------------------------------
## Heatmaps (3 silhouettes)
##--------------------------------------------
p_post <- make_heatmap(
  df3, "sil_post",
  "Posterior silhouette (D = 1 - PSM) \nover (alpha_x, alpha_y) \nfaceted by alpha_z",
  "Silhouette", blue_pal
)
ggsave(file.path(outDir, "silhouette_post_heatmap_xyz.pdf"), p_post, width = 9, height = 6)

p_net <- make_heatmap(
  df3, "sil_net",
  "Network silhouette (adjacency-profile distance) \nover (alpha_x, alpha_y) \nfaceted by alpha_z",
  "Silhouette", green_pal
)
ggsave(file.path(outDir, "silhouette_net_heatmap_xyz.pdf"), p_net, width = 9, height = 6)

p_cov <- make_heatmap(
  df3, "sil_cov",
  "Covariate silhouette (Euclidean on standardized xyz) \nover (alpha_x, alpha_y) \nfaceted by alpha_z",
  "Silhouette", purple_pal
)
ggsave(file.path(outDir, "silhouette_cov_heatmap_xyz.pdf"), p_cov, width = 9, height = 6)

##--------------------------------------------
## Load data for coordinates (same isolate-drop as model run)
##--------------------------------------------
load(here("application_brain", "consensus_scale33_tau50.RData"))
diag(A_bin_tau) <- 0
deg  <- rowSums(A_bin_tau)
keep <- deg > 0
coords_kept <- coords[keep, , drop = FALSE]

##--------------------------------------------
## Plot helper: make PSM heatmap + 3-view cluster map for a chosen "best"
##--------------------------------------------
make_view_nodes <- function(nodes_df, view = c("xy", "xz", "yz"),
                            node_alpha = 0.95, node_size = 2.0) {
  view <- match.arg(view)

  if (view == "xy") {
    xvar <- "x"; yvar <- "y"; ttl <- "Axial (x–y)"
  } else if (view == "xz") {
    xvar <- "x"; yvar <- "z"; ttl <- "Coronal (x–z)"
  } else {
    xvar <- "y"; yvar <- "z"; ttl <- "Sagittal (y–z)"
  }

  ggplot(nodes_df, aes(x = .data[[xvar]], y = .data[[yvar]], color = cluster)) +
    geom_point(alpha = node_alpha, size = node_size) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    labs(title = ttl, x = xvar, y = yvar, color = "Cluster") +
    theme(plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank())
}

save_best_plots <- function(best_row, criterion_label) {

  best_tag <- best_row$tag
  best_res <- res_app$results[[best_tag]]

  psm_mat <- best_res$psm
  z_hat   <- as.integer(best_res$z_hat)

  ## -------------------------------
  ## 1) PSM heatmap (ggplot)
  ## -------------------------------
  ord <- order(z_hat)
  psm_ord <- psm_mat[ord, ord, drop = FALSE]

  psm_df <- data.frame(
    i = rep(seq_len(nrow(psm_ord)), times = ncol(psm_ord)),
    j = rep(seq_len(ncol(psm_ord)), each  = nrow(psm_ord)),
    p = as.vector(psm_ord)
  )

  title_psm <- sprintf(
    "PSM (ordered by VI partition) — best by %s silhouette: (%.2f, %.2f, %.2f)",
    criterion_label, best_row$alpha_x, best_row$alpha_y, best_row$alpha_z
  )

  p_psm <- ggplot(psm_df, aes(x = j, y = i, fill = p)) +
    geom_tile() +
    scale_fill_gradientn(colors = blue_pal, limits = c(0, 1)) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    labs(title = title_psm, x = "Node index (ordered)", y = "Node index (ordered)", fill = "PSM") +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    file.path(outDir, sprintf("best_%s_psm_heatmap.pdf", criterion_label)),
    p_psm, width = 8, height = 7, dpi = 300
  )

  ## -------------------------------
  ## 2) 3-view spatial plots
  ## -------------------------------
  nodes_df <- coords_kept
  nodes_df$cluster <- factor(z_hat)

  p_xy <- make_view_nodes(nodes_df, "xy")
  p_xz <- make_view_nodes(nodes_df, "xz")
  p_yz <- make_view_nodes(nodes_df, "yz")

  p_views <- cowplot::plot_grid(p_xy, p_xz, p_yz, nrow = 1)

  ggsave(
    filename = file.path(outDir, sprintf("best_%s_clusters_3views.pdf", criterion_label)),
    plot     = p_views,
    width    = 12,
    height   = 4,
    dpi      = 300
  )

  ## -------------------------------
  ## 3) Matrix plots (PSM + adjacency)
  ## -------------------------------
  # adjacency matrix A must already be loaded & isolate-filtered
  pdf(file.path(outDir, sprintf("best_%s_psm_matrix.pdf", criterion_label)),
      width = 7, height = 7)
  plot_matrix_with_clusters(
    psm_mat, z_hat,
    type  = "psm",
    title = sprintf("PSM — best by %s silhouette", criterion_label)
  )
  dev.off()

  pdf(file.path(outDir, sprintf("best_%s_adj_matrix.pdf", criterion_label)),
      width = 7, height = 7)
  plot_matrix_with_clusters(
    A_bin_tau[keep, keep] , z_hat,
    type  = "adj",
    title = sprintf("Adjacency — best by %s silhouette", criterion_label)
  )
  dev.off()

  invisible(NULL)
}

##--------------------------------------------
## Produce posterior plots for all three "best" criteria
##--------------------------------------------
save_best_plots(best_post, "post")
save_best_plots(best_net,  "net")
save_best_plots(best_cov,  "cov")

cat("\nSaved:\n",
    " - silhouette heatmaps (post/net/cov)\n",
    " - best_post/net/cov_psm_heatmap.pdf\n",
    " - best_post/net/cov_clusters_3views.pdf\n",
    "in:\n  ", outDir, "\n")