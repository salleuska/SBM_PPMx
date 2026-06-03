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
similarity_calibration <- "raw"

outDir <- here("simulation", "learn_alpha", similarity_calibration, 
  "results", "processed", "oneCov")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
################################

## read processed results (1-cov)
res_neutral <- readRDS(
  here("simulation", "learn_alpha", similarity_calibration, "results", "processed", "oneCov",
       "scenario_1cov_neutral_results.rds")
)

res_info <- readRDS(
  here("simulation", "learn_alpha",similarity_calibration, "results", "processed", "oneCov",
       "scenario_1cov_informative_results.rds")
)

# res_mis_rand <- readRDS(
#   here("simulation", "learn_alpha", "results", "processed", "oneCov",
#        "scenario_1cov_mislead_random_results.rds")
# )

# res_mis_shift <- readRDS(
#   here("simulation", "learn_alpha", "results", "processed", "oneCov",
#        "scenario_1cov_mislead_shifted_results.rds")
# )
################################
extract_df <- function(obj, scen_name) {
  out <- lapply(names(obj$results), function(nm) {
    r <- obj$results[[nm]]
    alpha_vec <- as.numeric(r$alpha_post)
    
    data.frame(
      scenario = scen_name,
      tag = nm,
      seed = r$seed,
      a_alpha = r$alpha_g$a_alpha,
      b_alpha = r$alpha_g$b_alpha,
      prior = paste0("Gamma(", r$alpha_g$a_alpha, ",", r$alpha_g$b_alpha, ")"),
      alpha_mean = mean(alpha_vec),
      alpha_sd = sd(alpha_vec),
      alpha_lwr = quantile(alpha_vec, 0.025),
      alpha_upr = quantile(alpha_vec, 0.975),
      ARI_vi = r$ARI_vi,
      ARI_mean = r$ARI_mean,
      ARI_sd = r$ARI_sd,
      silhouette = r$silhouette$mean
    )
  })

  do.call(rbind, out)
}
#############################################################################
## Overall clustering results as a function of alpha
#############################################################################
df <- rbind(
  extract_df(res_neutral,   "neutral"),
  extract_df(res_info,      "informative")
)

df$scenario <- factor(
  df$scenario,
  levels = c("neutral", "informative")
)


p_alpha <- ggplot(df, aes(x = scenario, y = alpha_mean, color = prior)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(
    aes(ymin = alpha_lwr, ymax = alpha_upr),
    position = position_dodge(width = 0.5),
    width = 0.2,
    alpha = 0.7
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Posterior distribution of learned alpha",
    x = NULL,
    y = expression(alpha[g])
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

p_alpha

ggsave(file.path(outDir, "posterior_alpha_1cov.pdf"), p_alpha, width = 7, height = 4)
p_ari <- ggplot(df, aes(x = scenario, y = ARI_mean, color = prior)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(
    aes(ymin = ARI_mean - ARI_sd, ymax = ARI_mean + ARI_sd),
    position = position_dodge(width = 0.5),
    width = 0.2,
    alpha = 0.7
  ) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Posterior mean ARI under learned alpha",
    x = NULL,
    y = expression(E(ARI(z, z[0]) ~ "|" ~ data))
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

ggsave(file.path(outDir, "ARI_mean_withSD_learnedAlpha_1cov.pdf"), p_ari, width = 7, height = 4)

#############################################################################
## Trace plots: one PDF per scenario/result
#############################################################################
make_trace_df <- function(obj, scenario, result_name) {
  r <- obj$results[[result_name]]

  alpha <- as.numeric(r$alpha_post)
  rate  <- as.numeric(r$alpha_rate_post)
  b_alpha <- r$alpha_g$b_alpha

  data.frame(
    iter = seq_along(alpha),
    alpha = alpha,
    similarity = rate - b_alpha,
    scenario = scenario,
    tag = result_name
  )
}


plot_trace_pair <- function(obj, scenario, result_name, outDir) {

  dat <- make_trace_df(obj, scenario, result_name)

  mu_alpha <- mean(dat$alpha, na.rm = TRUE)
  dat$alpha_rm <- cumsum(dat$alpha) / seq_along(dat$alpha)

  p_alpha <- ggplot(dat, aes(iter)) +
    geom_line(aes(y = alpha), alpha = 0.35, linewidth = 0.3) +
    geom_line(aes(y = alpha_rm), linewidth = 1) +
    geom_hline(yintercept = mu_alpha, linetype = 2) +
    theme_minimal(base_size = 14) +
    labs(title = paste(scenario, "scenario. alpha trace."),
         subtitle = paste("mean =", round(mu_alpha, 2)),
         x = NULL,
         y = expression(alpha[g]))

  mu_sim <- mean(dat$similarity, na.rm = TRUE)
  dat$sim_rm <- cumsum(dat$similarity) / seq_along(dat$similarity)

  p_sim <- ggplot(dat, aes(iter)) +
    geom_line(aes(y = similarity), alpha = 0.35, linewidth = 0.3) +
    geom_line(aes(y = sim_rm), linewidth = 1) +
    geom_hline(yintercept = mu_sim, linetype = 2) +
    theme_minimal(base_size = 14) +
    labs(title = paste(scenario, "scenario. Similarity trace."),
         subtitle = paste("mean =", round(mu_sim, 2)),
         x = "Iteration",
         y = "Similarity")

  p <- cowplot::plot_grid(p_alpha, p_sim, ncol = 1, align = "v")

  ggsave(
    file.path(outDir,
              paste0("trace_alpha_similarity_", scenario, "_", result_name, ".pdf")),
    p,
    width = 7,
    height = 8
  )
}

for (nm in names(res_neutral$results))
  plot_trace_pair(res_neutral, "neutral", nm, outDir)

for (nm in names(res_info$results))
  plot_trace_pair(res_info, "informative", nm, outDir)