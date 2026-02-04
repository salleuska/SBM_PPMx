################################
## script for producing plots 
library(ggplot2)
library(RColorBrewer)
library(cowplot)
################################
## read files with results
res_none <- readRDS("scenario_none_results.rds")
res_one <- readRDS("scenario_one_results.rds")
################################
extract_df <- function(obj, scen_name) {
  ## extract for each alpha
  ## - the ARI between true and estimated partition (salso & VI)
  ## - the average ARI between true and the partition at each iteration
  ## - silhouette between true and estimated partition (salso & VI)
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
## Overall clustering results as a function of the calibration weight
## Extract ARI_VI, averaged ARI and silhouette
#############################################################################

df <- rbind(
  extract_df(res_none, "neutral"),
  extract_df(res_one,  "one")
)

## remove alpha = 4 
df <- df[df$alpha != 4, ]

df$scenario <- factor(df$scenario, levels = c("neutral", "one"), labels = c("neutral", "informative"))
df <- df[order(df$alpha), ]

### Plot of average ari

p_ari_mean <- ggplot(df, aes(x = alpha, y = ARI_mean, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ARI_mean - ARI_sd,
                    ymax = ARI_mean + ARI_sd),
                width = 0.1, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = unique(df$alpha)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "ARI vs alpha",
    x     = expression(alpha),
    y     = expression(E(ARI(z, z[0]) ~ "|" ~ data))
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.title    = element_blank(),
    legend.position = "right"
  )

p_ari_mean

ggsave("average_ARI_withSD.pdf", p_ari_mean,
       width = 6, height = 4)

### Using the VI
p_ari_vi <- ggplot(df, aes(x = alpha, y = ARI_vi, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = unique(df$alpha)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "ARI between the partition estimate and true clustering",
    x     = expression(alpha),
    y     = expression(ARI(hat(z), z[0]))
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.title    = element_blank(),
    legend.position = "right"
  )

ggsave("ARI_trueVsEstimated.pdf", p_ari_vi,
       width = 6, height = 4)

#####################
### Silhouette ###

## Sally - a questa non ci ho pensato molto
p_sil <- ggplot(df, aes(x = alpha, y = silhouette, color = scenario)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Silhouette vs alpha",
    x = expression(alpha),
    y = "Average silhouette width"
  )

print(p_sil)

ggsave("silhouette.pdf", p_sil,
       width = 6, height = 4)

#############################################################################
## Details about the clusering results 
## posterior probability of co-clustering + estimated partition + true partition
#############################################################################

plot_psm_with_annotations <- function(psm, z_true, z_est) {

  require(pheatmap)
  require(cowplot)

  n <- nrow(psm)

  ## rownames for the co-occurence probability matrix
  rn <- paste0("i", seq_len(n))
  rownames(psm) <- rn

  ## dataframe for annotation on the left side (clustering info)
  ann <- data.frame(
    est  = factor(z_est),
    true = factor(z_true)
  )
  rownames(ann) <- rn

  ###############################
  ## colors for the trues clusters
  lev_true <- levels(ann$true)
  K0 <- length(lev_true)

  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73",
    "#D55E00", "#CC79A7", "#0072B2",
    "#F0E442", "#999999"
  )

  ## set colors for true partition
  lev_true <- levels(ann$true)
  cols_true <- okabe_ito[seq_len(length(lev_true))]
  names(cols_true) <- lev_true

  ## set colors for estimated partition
  lev_est <- levels(ann$est)
  cols_est <- okabe_ito[seq_len(length(lev_est))]
  names(cols_est) <- lev_est

  annotation_colors <- list(
    true = cols_true,
    est  = cols_est
  )

  ## Greyscale palette for co-clustering probabilities
  mat_cols <- colorRampPalette(brewer.pal(9, "Greys"))(50)
 
  ## Blue palette
#   mat_cols <- colorRampPalette(brewer.pal(9, "Blues"))(50)

  p <- pheatmap(
    psm,
    color          = mat_cols,
    cluster_rows   = FALSE,
    cluster_cols   = FALSE,
    show_rownames  = FALSE,
    show_colnames  = FALSE,
    ## Annotation row and colors define the columns
    annotation_row = ann,
    annotation_colors = annotation_colors,
    border_color   = NA,
    annotation_legend = FALSE,
    legend         = FALSE,
    ## Introduce gaps between to separate true clusters
    gaps_row=c(which(diff(z_true)!=0)),
    gaps_col=c(which(diff(z_true)!=0))
  )

  p
}

## function to relabel the estimated partition for plotting purposes 
## (the partition is defined in lexicographical order, and the meaning of the 
## last two clusters is swapped)

simple_relabel <- function(z_hat, z_true) {
  sapply(z_hat, function(h) {
    idx <- which(z_hat == h)
    as.integer(names(which.max(table(z_true[idx]))))
  })
}

## No covariate effect
alpha0 <- res_one$results$alpha_0.00

p_alpha0 <- plot_psm_with_annotations(
  psm    = alpha0$psm,
  z_true = res_one$sim$partition_true,
  z_est  = alpha0$z_hat)


p_alpha0 <- cowplot::ggdraw() +
  cowplot::draw_label("No covariate effect (alpha = 0)",
                      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(p_alpha0[[4]], x = 0, y = 0, width = 1, height = 0.92)


ggsave("psm_alpha0.pdf", p_alpha0,
       width = 6, height = 6)
##########

### Scenario with effect of one covariate
alpha05 <- res_one$results$alpha_0.50

## no relabeling
# p_alpha05 <- plot_psm_with_annotations(
#   psm    = alpha05$psm,
#   z_true = res_one$sim$partition_true,
#   z_est  = alpha05$z_hat
# )

p_alpha05 <- plot_psm_with_annotations(
  psm    = alpha05$psm,
  z_true = res_one$sim$partition_true,
  z_est  = simple_relabel(alpha05$z_hat, z_true = res_one$sim$partition_true)
)

p_alpha05 <- cowplot::ggdraw() +
  cowplot::draw_label("Covariate effect (alpha = 0.5)",
                      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(p_alpha05[[4]], x = 0, y = 0, width = 1, height = 0.92)


ggsave("psm_alpha05.pdf", p_alpha05,
       width = 6, height = 6)
##########
alpha1 <- res_one$results$alpha_1.00

## no relabeling
# p_alpha1 <- plot_psm_with_annotations(
#   psm    = alpha1$psm,
#   z_true = res_one$sim$partition_true,
#   z_est  = alpha1$z_hat
# )

p_alpha1 <- plot_psm_with_annotations(
  psm    = alpha1$psm,
  z_true = res_one$sim$partition_true,
  z_est  = simple_relabel(alpha1$z_hat, z_true = res_one$sim$partition_true)
)

p_alpha1 <- cowplot::ggdraw() +
  cowplot::draw_label("Covariate effect (alpha = 1)",
                      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(p_alpha1[[4]], x = 0, y = 0, width = 1, height = 0.92)


ggsave("psm_alpha1.pdf", p_alpha1,
       width = 6, height = 6)
##########
alpha2 <- res_one$results$alpha_2.00

p_alpha2 <- plot_psm_with_annotations(
  psm    = alpha2$psm,
  z_true = res_one$sim$partition_true,
  z_est  = alpha2$z_hat
)

p_alpha2 <- cowplot::ggdraw() +
  cowplot::draw_label("Covariate effect (alpha = 2)",
                      fontface = "bold", x = 0.5, y = 0.98, hjust = 0.5) +
  cowplot::draw_plot(p_alpha2[[4]], x = 0, y = 0, width = 1, height = 0.92)


ggsave("psm_alpha2.pdf", p_alpha2,
       width = 6, height = 6)
