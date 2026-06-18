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

## read processed results (1-cov) — ESBM with learned alpha
res_neutral <- readRDS(
  here("simulation", "learn_alpha", similarity_calibration, "results", "processed", "oneCov",
       "scenario_1cov_neutral_results.rds")
)

res_info <- readRDS(
  here("simulation", "learn_alpha", similarity_calibration, "results", "processed", "oneCov",
       "scenario_1cov_informative_results.rds")
)

res_mis_rand <- readRDS(
  here("simulation", "learn_alpha", similarity_calibration, "results", "processed", "oneCov",
       "scenario_1cov_mislead_random_results.rds")
)

## read processed results — SBM baseline (no covariate)
sbm_dir <- here("simulation", "learn_alpha", "SBM", "results", "processed", "1Cov")

res_sbm_neutral  <- readRDS(file.path(sbm_dir, "binarySBM_1cov_N_results.rds"))
res_sbm_info     <- readRDS(file.path(sbm_dir, "binarySBM_1cov_I_results.rds"))
res_sbm_mis_rand <- readRDS(file.path(sbm_dir, "binarySBM_1cov_M_results.rds"))

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
      silhouette = r$silhouette$mean,
      model = "ESBM"
    )
  })

  do.call(rbind, out)
}

## SBM baseline extractor (no alpha info)
extract_df_sbm <- function(obj, scen_name) {
  out <- lapply(names(obj$results), function(nm) {
    r <- obj$results[[nm]]
    data.frame(
      scenario   = scen_name,
      tag        = nm,
      seed       = r$seed,
      a_alpha    = NA_real_,
      b_alpha    = NA_real_,
      prior      = "SBM (no covariate)",
      alpha_mean = NA_real_,
      alpha_sd   = NA_real_,
      alpha_lwr  = NA_real_,
      alpha_upr  = NA_real_,
      ARI_vi     = r$ARI_vi,
      ARI_mean   = r$ARI_mean,
      ARI_sd     = r$ARI_sd,
      silhouette = r$silhouette$mean,
      model      = "SBM"
    )
  })
  do.call(rbind, out)
}
#############################################################################
## Overall clustering results as a function of alpha
#############################################################################

## ESBM results (learned alpha, with covariate)
df_esbm <- rbind(
  extract_df(res_neutral,   "neutral"),
  extract_df(res_info,      "informative"),
  extract_df(res_mis_rand,  "misleading")
)

## SBM baseline (no covariate)
df_sbm <- rbind(
  extract_df_sbm(res_sbm_neutral,   "neutral"),
  extract_df_sbm(res_sbm_info,      "informative"),
  extract_df_sbm(res_sbm_mis_rand,  "misleading")
)

## combined (used for ARI plot)
df_all <- rbind(df_esbm, df_sbm)
df_all$scenario <- factor(
  df_all$scenario,
  levels = c("neutral", "informative", "misleading")
)

## ESBM-only (used for alpha plot)
df <- df_esbm
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

p_alpha

ggsave(file.path(outDir, "posterior_alpha_1cov.pdf"), p_alpha, width = 7, height = 4)
p_ari <- ggplot(df_all, aes(x = scenario, y = ARI_mean, color = prior, shape = model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(
    aes(ymin = ARI_mean - ARI_sd, ymax = ARI_mean + ARI_sd),
    position = position_dodge(width = 0.5),
    width = 0.2,
    alpha = 0.7
  ) +
  scale_shape_manual(values = c("ESBM" = 16, "SBM" = 17)) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Posterior mean ARI: ESBM vs SBM baseline",
    x = NULL,
    y = expression(E(ARI(z, z[0]) ~ "|" ~ data)),
    shape = NULL
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