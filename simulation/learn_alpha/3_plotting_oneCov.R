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
outDir <- here("simulation", "learn_alpha", "results", "processed", "oneCov")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
################################

## read processed results (1-cov)
res_neutral <- readRDS(
  here("simulation", "learn_alpha", "results", "processed", "oneCov",
       "scenario_1cov_neutral_results.rds")
)

res_info <- readRDS(
  here("simulation", "learn_alpha", "results", "processed", "oneCov",
       "scenario_1cov_informative_results.rds")
)

res_mis_rand <- readRDS(
  here("simulation", "learn_alpha", "results", "processed", "oneCov",
       "scenario_1cov_mislead_random_results.rds")
)

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
  extract_df(res_info,      "informative"),
  extract_df(res_mis_rand,  "misleading")
)

df$scenario <- factor(
  df$scenario,
  levels = c("neutral", "informative", "misleading")
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



ggsave(file.path(outDir, "posterior_alpha_1cov.pdf"), p_alpha, width = 7, height = 4)

p_ari <- ggplot(df, aes(x = scenario, y = ARI_vi, color = prior)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "ARI of VI estimate under learned alpha",
    x = NULL,
    y = expression(ARI(hat(z), z[0]))
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

ggsave(file.path(outDir, "ARI_learnedAlpha_1cov.pdf"), p_ari, width = 7, height = 4)

trace_df <- function(obj, scenario, result_name) {
  r <- obj$results[[result_name]]
  data.frame(
    iter = seq_along(as.numeric(r$alpha_keep)),
    alpha = as.numeric(r$alpha_keep),
    scenario = scenario,
    prior = paste0("Gamma(", r$alpha_g$a_alpha, ",", r$alpha_g$b_alpha, ")"),
    seed = r$seed
  )
}