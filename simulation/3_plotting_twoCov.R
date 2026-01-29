library(here)
library(ggplot2)
library(RColorBrewer)
library(scales)

## read processed results (2-cov)
res_NN <- readRDS(here("simulation", "results", "processed", "twoCov",
                       "scenario_2cov_NN_results.rds"))
res_IN <- readRDS(here("simulation", "results", "processed", "twoCov",
                       "scenario_2cov_IN_results.rds"))
res_NI <- readRDS(here("simulation", "results", "processed", "twoCov",
                       "scenario_2cov_NI_results.rds"))
res_II <- readRDS(here("simulation", "results", "processed", "twoCov",
                       "scenario_2cov_II_results.rds"))

extract_df_2cov <- function(obj, scen_name) {
  tags <- names(obj$results)

  ## pull alpha1/alpha2 from stored fields (preferred)
  alpha1 <- sapply(obj$results, function(r) r$alpha1)
  alpha2 <- sapply(obj$results, function(r) r$alpha2)

  data.frame(
    scenario   = scen_name,
    alpha1     = alpha1,
    alpha2     = alpha2,
    ARI_vi     = sapply(obj$results, function(r) r$ARI_vi),
    ARI_mean   = sapply(obj$results, function(r) r$ARI_mean),
    ARI_sd     = sapply(obj$results, function(r) r$ARI_sd),
    silhouette = sapply(obj$results, function(r) r$silhouette$mean),
    stringsAsFactors = FALSE
  )
}

df2 <- rbind(
  extract_df_2cov(res_NN, "NN"),
  extract_df_2cov(res_IN, "IN"),
  extract_df_2cov(res_NI, "NI"),
  extract_df_2cov(res_II, "II")
)

## Discrete grid for tiles (ordering only; fill is continuous)
df2$alpha1_f <- factor(df2$alpha1, levels = sort(unique(df2$alpha1)))
df2$alpha2_f <- factor(df2$alpha2, levels = sort(unique(df2$alpha2)))

## helper: single-hue palette
blue_pal <- colorRampPalette(brewer.pal(9, "Blues"))(50)

## Heatmap: ARI (VI estimate)
p_heat_ari <- ggplot(df2, aes(x = alpha1_f, y = alpha2_f, fill = ARI_vi)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", ARI_vi)), size = 3) +
  scale_fill_gradientn(colors = blue_pal, limits = c(0, 1)) +
  facet_wrap(~ scenario) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  labs(
    title = "ARI (VI estimate) over (alpha1, alpha2)",
    x     = expression(alpha[1]),
    y     = expression(alpha[2]),
    fill  = "ARI"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.title = element_text(size = 11)
  )

print(p_heat_ari)

outDir2 <- here("simulation", "results", "processed", "twoCov")
dir.create(outDir2, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(outDir2, "ARI_heatmap_2cov.pdf"),
  plot     = p_heat_ari,
  width    = 8,
  height   = 5
)

## Heatmap: posterior mean ARI
p_heat_mean <- ggplot(df2, aes(x = alpha1_f, y = alpha2_f, fill = ARI_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", ARI_mean)), size = 3) +
  scale_fill_gradientn(colors = blue_pal, limits = c(0, 1)) +
  facet_wrap(~ scenario) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  labs(
    title = "Posterior mean ARI over (alpha1, alpha2)",
    x     = expression(alpha[1]),
    y     = expression(alpha[2]),
    fill  = "ARI"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(outDir2, "ARImean_heatmap_2cov.pdf"),
  plot     = p_heat_mean,
  width    = 8,
  height   = 5
)