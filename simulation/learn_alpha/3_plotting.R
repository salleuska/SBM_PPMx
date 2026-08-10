################################
## 3_plotting.R
## Plots for processed ESBM results
################################

library(ggplot2)
library(cowplot)
library(tools)
library(here)

args <- R.utils::commandArgs(asValue = TRUE)

## ------------------------------------------------------------
## Arguments
## ------------------------------------------------------------

similarity_calibration <- if (is.null(args$similarity_calibration)) {
  "raw"
} else {
  as.character(args$similarity_calibration)
}

processed_root <- if (is.null(args$processed_root)) {
  here("simulation", "learn_alpha", similarity_calibration, "results", "processed")
} else {
  args$processed_root
}

nCov_plot <- if (is.null(args$nCov)) {
  NULL
} else {
  as.integer(args$nCov)
}

outDir <- file.path(processed_root, "plots")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

cat("Processed root:", processed_root, "\n")
cat("Output dir:", outDir, "\n")
cat("nCov plot:", ifelse(is.null(nCov_plot), "all", nCov_plot), "\n\n")

## ------------------------------------------------------------
## Read processed files
## ------------------------------------------------------------

processed_files <- list.files(
  processed_root,
  pattern = "^binarySBM_[0-9]+cov_[A-Z0-9]+_results\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(processed_files) == 0) {
  stop("No processed result files found in: ", processed_root)
}

## ------------------------------------------------------------
## Helpers
## ------------------------------------------------------------
scenario_caption <- "N = neutral, I = informative, M = misleading"

## ESBM results (learned alpha, with covariate)
extract_df <- function(f) {

  obj <- readRDS(f)

  sim_name <- file_path_sans_ext(basename(f))
  sim_name <- sub("_results$", "", sim_name)

  nCov <- as.integer(sub("^binarySBM_([0-9]+)cov_[A-Z0-9]+$", "\\1", sim_name))
  scenario <- sub("^binarySBM_[0-9]+cov_([A-Z0-9]+)$", "\\1", sim_name)

  out <- lapply(names(obj$results), function(nm) {

    r <- obj$results[[nm]]

    alpha_mat <- r$alpha_post

    if (is.null(dim(alpha_mat))) {
      alpha_mat <- matrix(alpha_mat, nrow = 1)
    }

    ARI_vec <- apply(r$Z_post, 2,
      function(z_t) mclust::adjustedRandIndex(z_t, obj$sim$partition_true)
    )

    covDep <- obj$sim$covDep
    cov_type <- ifelse(covDep == "informative",    "informative",
                ifelse(covDep == "mislead_random", "misleading", "neutral"))

    out_alpha <- lapply(seq_len(nrow(alpha_mat)), function(j) {

      alpha_vec <- as.numeric(alpha_mat[j, ])

      data.frame(
        sim_name = sim_name,
        nCov = nCov,
        scenario_label = scenario,
        covariate = paste0("x", j),
        cov_type = cov_type[j],
        tag = nm,
        seed = r$seed,
        a_alpha = r$alpha_g$a_alpha,
        b_alpha = r$alpha_g$b_alpha,
        prior = paste0("Gamma(", r$alpha_g$a_alpha, ",", r$alpha_g$b_alpha, ")"),
        model = "ESBM",

        alpha_mean = mean(alpha_vec, na.rm = TRUE),
        alpha_lwr = bayestestR::hdi(alpha_vec, ci = 0.95)$CI_low,
        alpha_upr = bayestestR::hdi(alpha_vec, ci = 0.95)$CI_high,

        ARI_vi = r$ARI_vi,
        ARI_mean = mean(ARI_vec, na.rm = TRUE),
        ARI_lwr = quantile(ARI_vec, 0.025, na.rm = TRUE),
        ARI_upr = quantile(ARI_vec, 0.975, na.rm = TRUE),

        silhouette = r$silhouette$mean
      )
    })

    do.call(rbind, out_alpha)
  })

  do.call(rbind, out)
}

## SBM baseline (no covariate)
extract_df_sbm <- function(f) {

  obj <- readRDS(f)

  sim_name <- file_path_sans_ext(basename(f))
  sim_name <- sub("_results$", "", sim_name)

  nCov <- as.integer(sub("^binarySBM_([0-9]+)cov_[A-Z0-9]+$", "\\1", sim_name))
  scenario <- sub("^binarySBM_[0-9]+cov_([A-Z0-9]+)$", "\\1", sim_name)

  out <- lapply(names(obj$results), function(nm) {

    r <- obj$results[[nm]]

    ARI_vec <- apply(r$Z_post, 2,
      function(z_t) mclust::adjustedRandIndex(z_t, obj$sim$partition_true)
    )

    data.frame(
      sim_name       = sim_name,
      nCov           = nCov,
      scenario_label = scenario,
      covariate      = NA_character_,
      tag            = nm,
      seed           = r$seed,
      a_alpha        = NA_real_,
      b_alpha        = NA_real_,
      prior          = "SBM (no covariate)",
      model          = "SBM",

      alpha_mean     = NA_real_,
      alpha_lwr      = NA_real_,
      alpha_upr      = NA_real_,

      ARI_vi         = r$ARI_vi,
      ARI_mean       = mean(ARI_vec, na.rm = TRUE),
      ARI_lwr        = quantile(ARI_vec, 0.025, na.rm = TRUE),
      ARI_upr        = quantile(ARI_vec, 0.975, na.rm = TRUE),

      silhouette     = r$silhouette$mean
    )
  })

  do.call(rbind, out)
}

## ------------------------------------------------------------
## Build plotting data
## ------------------------------------------------------------

## ESBM
df <- do.call(rbind, lapply(processed_files, extract_df))

## SBM baseline
sbm_processed_root <- here("simulation", "learn_alpha", "SBM", "results", "processed")
sbm_files <- list.files(
  sbm_processed_root,
  pattern = "^binarySBM_[0-9]+cov_[A-Z0-9]+_results\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(sbm_files) > 0) {
  df_sbm <- do.call(rbind, lapply(sbm_files, extract_df_sbm))
} else {
  warning("No SBM processed files found in: ", sbm_processed_root,
          "\nRun 2_processSBMresults.R first.")
  df_sbm <- df[0, ]   # empty frame with same columns
}

if (!is.null(nCov_plot)) {
  df      <- subset(df,     nCov == nCov_plot)
  df_sbm  <- subset(df_sbm, nCov == nCov_plot)
}

if (nrow(df) == 0) {
  stop("No ESBM results found for requested nCov.")
}

## combined for ARI plot
df_all <- rbind(df, df_sbm)

## shared factor levels (scenario ordering by nCov)
scen_levels <- unique(df$scenario_label[order(df$nCov, df$scenario_label)])
df$scenario_label     <- factor(df$scenario_label,     levels = scen_levels)
df_all$scenario_label <- factor(df_all$scenario_label, levels = scen_levels)

nCov_tag <- if (is.null(nCov_plot)) "all" else paste0(nCov_plot, "Cov")

## ------------------------------------------------------------
## Posterior alpha plot
## ------------------------------------------------------------

type_colors <- c(informative = "#2166AC", misleading = "#D6604D", neutral = "#878787")

many_cov <- !is.null(nCov_plot) && nCov_plot >= 10

if (many_cov) {

  # one plot per scenario: covariates on x-axis sequentially, colored by type
  df_many <- df[df$nCov == nCov_plot, ]
  df_many$cov_type <- factor(df_many$cov_type, levels = c("informative", "misleading", "neutral"))

  # covariate index for sequential x ordering
  df_many$cov_idx <- as.integer(sub("^x", "", df_many$covariate))

  p_alpha <- ggplot(df_many, aes(x = cov_idx, y = alpha_mean, color = cov_type)) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = alpha_lwr, ymax = alpha_upr), width = 0, alpha = 0.6) +
    scale_color_manual(values = type_colors) +
    facet_wrap(~ scenario_label, ncol = 1) +
    theme_minimal(base_size = 13) +
    labs(
      title = "Posterior of learned alpha (50 covariates)",
      x = "Covariate index",
      y = expression(alpha[g]),
      color = "Covariate type"
    ) +
    theme(
      plot.title    = element_text(face = "bold"),
      legend.title  = element_text(face = "bold"),
      strip.text    = element_text(face = "bold")
    )

} else {

  p_alpha <- ggplot(df, aes(x = scenario_label, y = alpha_mean, color = covariate)) +
    geom_point(position = position_dodge(width = 0.6), size = 2.5) +
    geom_errorbar(
      aes(ymin = alpha_lwr, ymax = alpha_upr),
      position = position_dodge(width = 0.6),
      width = 0.2, alpha = 0.7
    ) +
    theme_minimal(base_size = 14) +
    ylim(c(0, 5)) +
    labs(
      title = "Posterior distribution of learned alpha",
      subtitle = scenario_caption,
      x = "Scenario",
      y = expression(alpha[g]),
      color = "Covariate"
    ) +
    theme(
      plot.title   = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      axis.text.x  = element_text(angle = 0, hjust = 0.5)
    )

  if (is.null(nCov_plot)) {
    p_alpha <- p_alpha + facet_wrap(~ nCov, scales = "free_x", labeller = label_both)
  }

}

ggsave(
  file.path(outDir, paste0("posterior_alpha_", nCov_tag, ".pdf")),
  p_alpha,
  width = if (many_cov) 10 else 8,
  height = if (many_cov) 3 * length(unique(df_many$scenario_label)) else 4.5
)

## ------------------------------------------------------------
## ARI plot with 95% posterior credible interval
## (ESBM points + SBM baseline as line/ribbon)
## ------------------------------------------------------------
df_ari_esbm <- df_all[df_all$model == "ESBM" & !duplicated(df_all$tag), ]
df_ari_sbm  <- df_all[df_all$model == "SBM",  ]

p_ari <- ggplot(mapping = aes(x = scenario_label)) +
  ## SBM baseline: ribbon + dashed line (drawn first, sits behind)
  geom_ribbon(
    data = df_ari_sbm,
    aes(ymin = ARI_lwr, ymax = ARI_upr, group = nCov),
    fill = "grey60",
    alpha = 0.25
  ) +
  geom_line(
    data = df_ari_sbm,
    aes(y = ARI_mean, group = nCov),
    color = "grey40",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  ## ESBM: points + error bars, colored by prior
  geom_errorbar(
    data = df_ari_esbm,
    aes(ymin = ARI_lwr, ymax = ARI_upr, color = prior),
    position = position_dodge(width = 0.5),
    width = 0.2,
    alpha = 0.7
  ) +
  geom_point(
    data = df_ari_esbm,
    aes(y = ARI_mean, color = prior),
    position = position_dodge(width = 0.5),
    size = 2.5
  ) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Posterior ARI: ESBM vs SBM baseline",
    subtitle = scenario_caption,
    x = NULL,
    y = expression(ARI(z, z[0]) ~ "|" ~ data),
    color = NULL
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

if (is.null(nCov_plot)) {
  p_ari <- p_ari + facet_wrap(~ nCov, scales = "free_x", labeller = label_both)
}

ggsave(
  file.path(outDir, paste0("ARI_mean_CI_learnedAlpha_", nCov_tag, ".pdf")),
  p_ari,
  width = 8,
  height = 4.5
)

## ------------------------------------------------------------
## Trace plots
## ------------------------------------------------------------

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

  p_alpha_trace <- ggplot(dat, aes(iter)) +
    geom_line(aes(y = alpha), alpha = 0.35, linewidth = 0.3) +
    geom_line(aes(y = alpha_rm), linewidth = 1) +
    geom_hline(yintercept = mu_alpha, linetype = 2) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(scenario, "scenario: alpha trace"),
      subtitle = paste("mean =", round(mu_alpha, 2)),
      x = NULL,
      y = expression(alpha[g])
    )

  mu_sim <- mean(dat$similarity, na.rm = TRUE)
  dat$sim_rm <- cumsum(dat$similarity) / seq_along(dat$similarity)

  p_sim_trace <- ggplot(dat, aes(iter)) +
    geom_line(aes(y = similarity), alpha = 0.35, linewidth = 0.3) +
    geom_line(aes(y = sim_rm), linewidth = 1) +
    geom_hline(yintercept = mu_sim, linetype = 2) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(scenario, "scenario: similarity trace"),
      subtitle = paste("mean =", round(mu_sim, 2)),
      x = "Iteration",
      y = "Similarity"
    )

  p <- cowplot::plot_grid(p_alpha_trace, p_sim_trace, ncol = 1, align = "v")

  ggsave(
    file.path(
      outDir,
      paste0("trace_alpha_similarity_", scenario, "_", result_name, ".pdf")
    ),
    p,
    width = 7,
    height = 8
  )
}

for (f in processed_files) {

  obj <- readRDS(f)

  sim_name <- file_path_sans_ext(basename(f))
  sim_name <- sub("_results$", "", sim_name)

  this_nCov <- as.integer(sub("^binarySBM_([0-9]+)cov_[A-Z0-9]+$", "\\1", sim_name))

  if (!is.null(nCov_plot) && this_nCov != nCov_plot) next

  trace_dir <- file.path(outDir, "traces", sim_name)
  dir.create(trace_dir, recursive = TRUE, showWarnings = FALSE)

  for (nm in names(obj$results)) {
    plot_trace_pair(obj, sim_name, nm, trace_dir)
  }
}

cat("Done.\n")