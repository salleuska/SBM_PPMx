library(here)
library(ggplot2)
library(RColorBrewer)
library(scales)

## read processed results (2-cov)
## scenario codes, one letter per covariate:
##   N = neutral, I = informative, M = misleading
scenarios <- c("NN", "IN", "II", "MN", "MI", "MM")

read_scenario <- function(code) {
  readRDS(here("simulation", "fixed_alpha", "results", "processed", "twoCov",
               sprintf("binarySBM_2cov_%s_results.rds", code)))
}

res_list <- lapply(scenarios, read_scenario)
names(res_list) <- scenarios

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

df2 <- do.call(
  rbind,
  lapply(scenarios, function(code) extract_df_2cov(res_list[[code]], code))
)

## keep panels in the declared scenario order
df2$scenario <- factor(df2$scenario, levels = scenarios)

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
  facet_wrap(~ scenario, nrow = 2) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  labs(
    title    = "ARI (VI estimate) over (alpha1, alpha2)",
    subtitle = "N = neutral, I = informative, M = misleading",
    x        = expression(alpha[1]),
    y        = expression(alpha[2]),
    fill     = "ARI"
  ) +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.title = element_text(size = 11)
  )

print(p_heat_ari)

outDir2 <- here("simulation", "fixed_alpha", "results", "processed", "twoCov")
dir.create(outDir2, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(outDir2, "ARI_heatmap_2cov.pdf"),
  plot     = p_heat_ari,
  width    = 10,
  height   = 6.5
)

## Heatmap: posterior mean ARI
p_heat_mean <- ggplot(df2, aes(x = alpha1_f, y = alpha2_f, fill = ARI_mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", ARI_mean)), size = 3) +
  scale_fill_gradientn(colors = blue_pal, limits = c(0, 1)) +
  facet_wrap(~ scenario, nrow = 2) +
  coord_equal() +
  theme_minimal(base_size = 13) +
  labs(
    title    = "Posterior mean ARI over (alpha1, alpha2)",
    subtitle = "N = neutral, I = informative, M = misleading",
    x        = expression(alpha[1]),
    y        = expression(alpha[2]),
    fill     = "ARI"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(outDir2, "ARImean_heatmap_2cov.pdf"),
  plot     = p_heat_mean,
  width    = 10,
  height   = 6.5
)