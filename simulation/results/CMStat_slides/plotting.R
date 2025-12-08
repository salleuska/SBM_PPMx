################################
## script for producing plots 
library(ggplot2)
library(RColorBrewer)1
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
    ARI_vi     = sapply(obj$results, function(r) r$ARI_vi),
    silhouette = sapply(obj$results, function(r) r$silhouette$mean)
  )
}
################################
## Extract ARI_VI, averaged ARI and silhouette

df <- rbind(
  extract_df(res_none, "none"),
  extract_df(res_one,  "one")
)

## remove alpha = 4 
df <- df[df$alpha != 4, ]

df$scenario <- factor(df$scenario, levels = c("none", "one"))
df <- df[order(df$alpha), ]

### Plot of average 

p_ari <- ggplot(df, aes(x = alpha, y = ARI_mean, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = unique(df$alpha)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "ARI between posterior partitions and true clustering",
    x     = expression(alpha),
    y     = expression(ARI(hat(z), z[0]))
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.title    = element_blank(),
    legend.position = "right"
  )

p_ari

################################
## Extract ARI_VI, averaged ARI and silhouette

df <- rbind(
  extract_df(res_none, "none"),
  extract_df(res_one,  "one")
)

## exclude alpha = 4
df <- df[!(df$alpha ==4), ]

df$scenario <- factor(df$scenario, levels = c("none", "one"))
df <- df[order(df$alpha), ]

p_ari <- ggplot(df, aes(x = alpha, y = ARI, color = scenario)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = unique(df$alpha)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "ARI between the estimated partition and true clustering",
    x     = expression(alpha),
    y     = expression(ARI(hat(z), z[0]))
  ) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.title    = element_blank(),
    legend.position = "right"
  )
p_ari


### Silhouette ###
p_sil <- ggplot(df, aes(x = alpha, y = silhouette, color = scenario)) +
  geom_point(size = 3) +
  geom_line() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Silhouette vs α",
    x = expression(alpha),
    y = "Average silhouette width"
  )

print(p_sil)

#######
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
  mat_cols <- colorRampPalette(c("white", "black"))(60)

  p <- pheatmap(
    psm,
    color          = mat_cols,
    cluster_rows   = FALSE,
    cluster_cols   = FALSE,
    show_rownames  = FALSE,
    show_colnames  = FALSE,
    annotation_row = ann,
    annotation_colors = annotation_colors,
    border_color   = NA,
    annotation_legend = FALSE,
    legend         = FALSE,
    gaps_row=c(which(diff(z_true)!=0)),
    gaps_col=c(which(diff(z_true)!=0))
  )

  ggdraw(p[[4]])
}


alpha0 <- res_one$results$alpha_0.00

p_alpha0 <- plot_psm_with_annotations(
  psm    = alpha0$psm,
  z_true = res_one$sim$partition_true,
  z_est  = alpha0$z_hat
)

ggsave(alpha0, )
alpha05 <- res_one$results$alpha_0.50

p_alpha05 <- plot_psm_with_annotations(
  psm    = alpha05$psm,
  z_true = res_one$sim$partition_true,
  z_est  = alpha05$z_hat
)

alpha1 <- res_one$results$alpha_1.00

p_alpha1 <- plot_psm_with_annotations(
  psm    = alpha1$psm,
  z_true = res_one$sim$partition_true,
  z_est  = alpha1$z_hat
)

alpha2 <- res_one$results$alpha_2.00

p_alpha2 <- plot_psm_with_annotations(
  psm    = alpha2$psm,
  z_true = res_one$sim$partition_true,
  z_est  = alpha2$z_hat
)

p_alpha0
p_alpha05
p_alpha1
p_alpha2

