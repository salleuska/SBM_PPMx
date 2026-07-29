################################
## script for producing plots (one-cov scenarios)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(mclust)
library(salso)
library(cluster)
library(tools)
library(here)
################################
## set directory for output
outDir <- here("simulation", "fixed_alpha", "results", "processed", "oneCov")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
################################

## read processed results (1-cov)
res_neutral <- readRDS(
  here("simulation", "fixed_alpha", "results", "processed", "oneCov",
       "scenario_1cov_neutral_results.rds")
)

res_info <- readRDS(
  here("simulation", "fixed_alpha", "results", "processed", "oneCov",
       "scenario_1cov_informative_results.rds")
)

res_mis_rand <- readRDS(
  here("simulation", "fixed_alpha", "results", "processed", "oneCov",
       "scenario_1cov_mislead_random_results.rds")
)

# res_mis_shift <- readRDS(
#   here("simulation", "fixed_alpha", "results", "processed", "oneCov",
#        "scenario_1cov_mislead_shifted_results.rds")
# )
################################
extract_df <- function(obj, scen_name) {
  alphas <- names(obj$results)

  data.frame(
    scenario   = scen_name,
    alpha      = as.numeric(sub("alpha_", "", alphas)),
    ARI_mean   = sapply(obj$results, function(r) r$ARI_mean),
    ARI_sd     = sapply(obj$results, function(r) r$ARI_sd),
    ARI_vi     = sapply(obj$results, function(r) r$ARI_vi),
    silhouette = sapply(obj$results, function(r) r$silhouette$mean)
  )
}

#############################################################################
## Overall clustering results as a function of alpha
#############################################################################
df <- rbind(
  extract_df(res_neutral,   "neutral"),
  extract_df(res_info,      "informative"),
  extract_df(res_mis_rand,  "misleading")
#  extract_df(res_mis_shift, "mislead_shifted")
)

## (optional) drop alpha = 4 if present
df <- df[df$alpha != 4, ]

df$scenario <- factor(
  df$scenario,
  levels = c("neutral", "informative", "misleading")
)

df <- df[order(df$scenario, df$alpha), ]

### Plot of average ARI (+/- SD)
p_ari_mean <- ggplot(df, aes(x = alpha, y = ARI_mean, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ARI_mean - ARI_sd,
                    ymax = ARI_mean + ARI_sd),
                width = 0.1, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = sort(unique(df$alpha))) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Posterior mean ARI vs alpha",
    x     = expression(alpha),
    y     = expression(E(ARI(z, z[0]) ~ "|" ~ data))
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.title    = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = file.path(outDir, "average_ARI_withSD_1cov.pdf"),
  plot     = p_ari_mean,
  width    = 7,
  height   = 4
)


### ARI of VI representative clustering
p_ari_vi <- ggplot(df, aes(x = alpha, y = ARI_vi, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = sort(unique(df$alpha))) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "ARI between VI estimate and true clustering",
    x     = expression(alpha),
    y     = expression(ARI(hat(z), z[0]))
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.title    = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(outDir, "ARI_trueVsEstimated_1cov.pdf"), p_ari_vi, width = 7, height = 4)

### Silhouette
p_sil <- ggplot(df, aes(x = alpha, y = silhouette, color = scenario)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = sort(unique(df$alpha))) +
  labs(
    title = "Silhouette vs alpha",
    x = expression(alpha),
    y = "Average silhouette width"
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.title    = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(outDir, "silhouette_1cov.pdf"),         p_sil,    width = 7, height = 4)

#############################################################################
## Details: PSM + estimated partition + true partition
#############################################################################
plot_psm_with_annotations <- function(psm, z_true, z_est) {
  require(pheatmap)
  require(cowplot)

  n <- nrow(psm)

  rn <- paste0("i", seq_len(n))
  rownames(psm) <- rn

  ann <- data.frame(
    est  = factor(z_est),
    true = factor(z_true)
  )
  rownames(ann) <- rn

  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73",
    "#D55E00", "#CC79A7", "#0072B2",
    "#F0E442", "#999999"
  )

  lev_true <- levels(ann$true)
  cols_true <- okabe_ito[seq_len(length(lev_true))]
  names(cols_true) <- lev_true

  lev_est <- levels(ann$est)
  cols_est <- okabe_ito[seq_len(length(lev_est))]
  names(cols_est) <- lev_est

  annotation_colors <- list(true = cols_true, est = cols_est)

  mat_cols <- colorRampPalette(brewer.pal(9, "Greys"))(50)

  pheatmap(
    psm,
    color             = mat_cols,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    annotation_row    = ann,
    annotation_colors = annotation_colors,
    border_color      = NA,
    annotation_legend = FALSE,
    legend            = FALSE,
    gaps_row          = c(which(diff(z_true) != 0)),
    gaps_col          = c(which(diff(z_true) != 0))
  )
}

simple_relabel <- function(z_hat, z_true) {
  sapply(z_hat, function(h) {
    idx <- which(z_hat == h)
    as.integer(names(which.max(table(z_true[idx]))))
  })
}

## Pick which scenario object to showcase in PSM plots
## Example: neutral at alpha=0, informative at alpha=0.5/1/2

alpha0  <- res_neutral$results$alpha_0.00
p0 <- plot_psm_with_annotations(
  psm    = alpha0$psm,
  z_true = res_neutral$sim$partition_true,
  z_est  = alpha0$z_hat
)
p0 <- cowplot::ggdraw() +
  cowplot::draw_label("Neutral covariate (alpha = 0)",
                      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(p0[[4]], x = 0, y = 0, width = 1, height = 0.92)

ggsave(file.path(outDir, "psm_neutral_alpha0.pdf"), p0, width = 6, height = 6)

alpha05 <- res_info$results$alpha_0.50
p05 <- plot_psm_with_annotations(
  psm    = alpha05$psm,
  z_true = res_info$sim$partition_true,
  z_est  = simple_relabel(alpha05$z_hat, z_true = res_info$sim$partition_true)
)
p05 <- cowplot::ggdraw() +
  cowplot::draw_label("Informative covariate (alpha = 0.5)",
                      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(p05[[4]], x = 0, y = 0, width = 1, height = 0.92)

ggsave(file.path(outDir, "psm_neutral_alpha05.pdf"), p05, width = 6, height = 6)

alpha1 <- res_info$results$alpha_1.00
p1 <- plot_psm_with_annotations(
  psm    = alpha1$psm,
  z_true = res_info$sim$partition_true,
  z_est  = simple_relabel(alpha1$z_hat, z_true = res_info$sim$partition_true)
)

p1 <- cowplot::ggdraw() +
  cowplot::draw_label("Informative covariate (alpha = 1)",
                      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(p1[[4]], x = 0, y = 0, width = 1, height = 0.92)

ggsave(file.path(outDir, "psm_informative_alpha1.pdf"), p1, width = 6, height = 6)


alpha2 <- res_info$results$alpha_2.00
p2 <- plot_psm_with_annotations(
  psm    = alpha2$psm,
  z_true = res_info$sim$partition_true,
  z_est  = simple_relabel(alpha2$z_hat, z_true = res_info$sim$partition_true)
)


p2 <- cowplot::ggdraw() +
  cowplot::draw_label("Informative covariate (alpha = 2)",
                      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(p2[[4]], x = 0, y = 0, width = 1, height = 0.92)

ggsave(file.path(outDir, "psm_informative_alpha2.pdf"), p2, width = 6, height = 6)
