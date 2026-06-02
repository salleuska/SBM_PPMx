################################
## plots for learned-alpha ESBM: two-covariate scenarios
################################

library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(mclust)
library(salso)
library(cluster)
library(tools)
library(here)

################################
## output directory
################################

outDir <- here("simulation", "learn_alpha", "results", "processed", "twoCov")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

################################
## read processed results
################################

res_NN <- readRDS(here("simulation", "learn_alpha", "results", "processed", "twoCov",
                      "scenario_2cov_NN_results.rds"))

res_IN <- readRDS(here("simulation", "learn_alpha", "results", "processed", "twoCov",
                      "scenario_2cov_IN_results.rds"))

res_NI <- readRDS(here("simulation", "learn_alpha", "results", "processed", "twoCov",
                      "scenario_2cov_NI_results.rds"))

res_II <- readRDS(here("simulation", "learn_alpha", "results", "processed", "twoCov",
                      "scenario_2cov_II_results.rds"))

################################
## extract summaries
################################

extract_df <- function(obj, scen_name) {
  out <- lapply(names(obj$results), function(nm) {
    r <- obj$results[[nm]]

    do.call(rbind, lapply(1:2, function(j) {
      alpha_vec <- as.numeric(r$alpha_post[j, ])

      data.frame(
        scenario = scen_name,
        tag = nm,
        seed = r$seed,
        covariate = paste0("x", j),
        a_alpha = r$alpha_g$a_alpha[j],
        b_alpha = r$alpha_g$b_alpha[j],
        prior = paste0("Gamma(", r$alpha_g$a_alpha[j], ",", r$alpha_g$b_alpha[j], ")"),
        alpha_mean = mean(alpha_vec),
        alpha_sd = sd(alpha_vec),
        alpha_lwr = quantile(alpha_vec, 0.025),
        alpha_upr = quantile(alpha_vec, 0.975),
        ARI_vi = r$ARI_vi,
        ARI_mean = r$ARI_mean,
        ARI_sd = r$ARI_sd,
        silhouette = r$silhouette$mean
      )
    }))
  })

  do.call(rbind, out)
}

df <- rbind(
  extract_df(res_NN, "NN"),
  extract_df(res_IN, "IN"),
  extract_df(res_NI, "NI"),
  extract_df(res_II, "II")
)

df$scenario <- factor(df$scenario, levels = c("NN", "IN", "NI", "II"))
df$covariate <- factor(df$covariate, levels = c("x1", "x2"))

################################
## posterior alpha plot
################################

p_alpha <- ggplot(df, aes(x = scenario, y = alpha_mean, color = covariate)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(
    aes(ymin = alpha_lwr, ymax = alpha_upr),
    position = position_dodge(width = 0.5),
    width = 0.2,
    alpha = 0.7
  ) +
  facet_wrap(~ prior) +
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

ggsave(
  file.path(outDir, "posterior_alpha_2cov.pdf"),
  p_alpha,
  width = 8,
  height = 4
)

################################
## ARI plot
################################

df_ari <- df[!duplicated(df$tag), ]

p_ari <- ggplot(df_ari, aes(x = scenario, y = ARI_mean, color = prior)) +
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

p_ari

ggsave(
  file.path(outDir, "ARI_mean_withSD_learnedAlpha_2cov.pdf"),
  p_ari,
  width = 7,
  height = 4
)

################################
## trace plots
################################

trace_df <- function(obj, scenario, result_name) {
  r <- obj$results[[result_name]]

  rbind(
    data.frame(
      iter = seq_along(r$alpha_post[1, ]),
      alpha = as.numeric(r$alpha_post[1, ]),
      scenario = scenario,
      covariate = "x1"
    ),
    data.frame(
      iter = seq_along(r$alpha_post[2, ]),
      alpha = as.numeric(r$alpha_post[2, ]),
      scenario = scenario,
      covariate = "x2"
    )
  )
}

trace_dat <- rbind(
  trace_df(res_NN, "NN", names(res_NN$results)[1]),
  trace_df(res_IN, "IN", names(res_IN$results)[1]),
  trace_df(res_NI, "NI", names(res_NI$results)[1]),
  trace_df(res_II, "II", names(res_II$results)[1])
)

trace_dat$scenario <- factor(trace_dat$scenario, levels = c("NN", "IN", "NI", "II"))
trace_dat$covariate <- factor(trace_dat$covariate, levels = c("x1", "x2"))

p_trace <- ggplot(trace_dat, aes(x = iter, y = alpha, color = covariate)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ scenario) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Trace plots for learned alpha",
    x = "Iteration",
    y = expression(alpha[g])
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

p_trace

ggsave(
  file.path(outDir, "trace_alpha_2cov.pdf"),
  p_trace,
  width = 8,
  height = 5
)