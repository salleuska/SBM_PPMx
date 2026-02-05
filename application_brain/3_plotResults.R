library(here)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

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
## Pick the best alpha by posterior silhouette
##--------------------------------------------
if (all(is.na(df3$sil_post))) stop("All sil_post are NA; cannot pick best alpha.")
best_i <- which.max(df3$sil_post)
best   <- df3[best_i, ]

cat("\nBest alpha (by posterior silhouette sil_post):\n")
print(best[, c("alpha_x","alpha_y","alpha_z","K_hat","sil_post","sil_net","sil_cov")])

# save a small summary table
write.csv(best, file.path(outDir, "best_alpha_by_sil_post.csv"), row.names = FALSE)

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
  "Posterior silhouette (D = 1 - PSM) over (alpha_x, alpha_y), faceted by alpha_z",
  "Silhouette", blue_pal
)
ggsave(file.path(outDir, "silhouette_post_heatmap_xyz.pdf"), p_post, width = 9, height = 6)

p_net <- make_heatmap(
  df3, "sil_net",
  "Network silhouette (adjacency-profile distance) over (alpha_x, alpha_y), faceted by alpha_z",
  "Silhouette", green_pal
)
ggsave(file.path(outDir, "silhouette_net_heatmap_xyz.pdf"), p_net, width = 9, height = 6)

p_cov <- make_heatmap(
  df3, "sil_cov",
  "Covariate silhouette (Euclidean on standardized xyz) over (alpha_x, alpha_y), faceted by alpha_z",
  "Silhouette", purple_pal
)
ggsave(file.path(outDir, "silhouette_cov_heatmap_xyz.pdf"), p_cov, width = 9, height = 6)

##--------------------------------------------
## Posterior plots for BEST alpha: PSM heatmap + 3-view cluster map
##--------------------------------------------

best_tag <- best$tag
best_res <- res_app$results[[best_tag]]

psm_mat <- best_res$psm
z_hat   <- as.integer(best_res$z_hat)

# order nodes by z_hat, then within cluster by average PSM similarity
ord <- order(z_hat)

psm_ord <- psm_mat[ord, ord, drop = FALSE]
z_ord   <- z_hat[ord]

# build long df for ggplot heatmap
psm_df <- data.frame(
  i = rep(seq_len(nrow(psm_ord)), times = ncol(psm_ord)),
  j = rep(seq_len(ncol(psm_ord)), each  = nrow(psm_ord)),
  p = as.vector(psm_ord)
)

p_psm <- ggplot(psm_df, aes(x = j, y = i, fill = p)) +
  geom_raster() +
  scale_fill_gradientn(colors = blue_pal, limits = c(0, 1)) +
  coord_equal() +
  theme_minimal(base_size = 12) +
  labs(
    title = sprintf("PSM (ordered by VI partition) — best alpha: (%.2f, %.2f, %.2f)",
                    best$alpha_x, best$alpha_y, best$alpha_z),
    x = "Node index (ordered)", y = "Node index (ordered)", fill = "PSM"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(outDir, "best_alpha_psm_heatmap.pdf"), p_psm, width = 7, height = 6)

##--------------------------------------------
## 3-view plot colored by z_hat (best alpha)
##--------------------------------------------
# Load original data to get coordinates and apply the same isolate-drop as the model run
load(here("application_brain", "consensus_scale33_tau50.RData"))

diag(A_bin_tau) <- 0
deg  <- rowSums(A_bin_tau)
keep <- deg > 0

coords_kept <- coords[keep, , drop = FALSE]
coords_kept$cluster <- factor(z_hat)  # z_hat already corresponds to kept nodes

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

p_xy <- make_view_nodes(coords_kept, "xy")
p_xz <- make_view_nodes(coords_kept, "xz")
p_yz <- make_view_nodes(coords_kept, "yz")

p_views <- cowplot::plot_grid(p_xy, p_xz, p_yz, nrow = 1)

ggsave(
  filename = file.path(outDir, "best_alpha_clusters_3views.pdf"),
  plot     = p_views,
  width    = 12,
  height   = 4,
  dpi      = 300
)

cat("\nSaved:\n",
    " - silhouette heatmaps (post/net/cov)\n",
    " - best_alpha_psm_heatmap.pdf\n",
    " - best_alpha_clusters_3views.pdf\n",
    "in:\n  ", outDir, "\n")