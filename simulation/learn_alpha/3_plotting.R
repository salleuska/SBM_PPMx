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
  pattern = "^binarySBM_[0-9]+cov_[A-Z]+_results\\.rds$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(processed_files) == 0) {
  stop("No processed result files found in: ", processed_root)
}

## ------------------------------------------------------------
## Helpers
## ------------------------------------------------------------

scenario_label <- function(sc) {
  d <- c(
    N = "neutral",
    I = "informative",
    M = "misleading"
  )

  paste(d[strsplit(sc, "")[[1]]], collapse = " + ")
}

extract_df <- function(f) {

  obj <- readRDS(f)

  sim_name <- file_path_sans_ext(basename(f))
  sim_name <- sub("_results$", "", sim_name)

  nCov <- as.integer(sub("^binarySBM_([0-9]+)cov_[A-Z]+$", "\\1", sim_name))
  scenario <- sub("^binarySBM_[0-9]+cov_([A-Z]+)$", "\\1", sim_name)

  out <- lapply(names(obj$results), function(nm) {

    r <- obj$results[[nm]]

    alpha_vec <- as.numeric(r$alpha_post)

    ARI_vec <- apply(
      r$Z_post,
      2L,
      function(z_t) mclust::adjustedRandIndex(z_t, obj$sim$partition_true)
    )

    data.frame(
      sim_name = sim_name,
      nCov = nCov,
      scenario = scenario,
      scenario_label = scenario_label(scenario),
      tag = nm,
      seed = r$seed,
      a_alpha = r$alpha_g$a_alpha,
      b_alpha = r$alpha_g$b_alpha,
      prior = paste0("Gamma(", r$alpha_g$a_alpha, ",", r$alpha_g$b_alpha, ")"),

      alpha_mean = mean(alpha_vec, na.rm = TRUE),
      alpha_lwr = quantile(alpha_vec, 0.025, na.rm = TRUE),
      alpha_upr = quantile(alpha_vec, 0.975, na.rm = TRUE),

      ARI_vi = r$ARI_vi,
      ARI_mean = mean(ARI_vec, na.rm = TRUE),
      ARI_lwr = quantile(ARI_vec, 0.025, na.rm = TRUE),
      ARI_upr = quantile(ARI_vec, 0.975, na.rm = TRUE),

      silhouette = r$silhouette$mean
    )
  })

  do.call(rbind, out)
}

## ------------------------------------------------------------
## Build plotting data
## ------------------------------------------------------------

df <- do.call(rbind, lapply(processed_files, extract_df))

if (!is.null(nCov_plot)) {
  df <- subset(df, nCov == nCov_plot)
}

if (nrow(df) == 0) {
  stop("No results found for requested nCov.")
}

df$scenario_label <- factor(
  df$scenario_label,
  levels = unique(df$scenario_label[order(df$nCov, df$scenario)])
)

nCov_tag <- if (is.null(nCov_plot)) "all" else paste0(nCov_plot, "Cov")

## ------------------------------------------------------------
## Posterior alpha plot
## ------------------------------------------------------------

p_alpha <- ggplot(df, aes(x = scenario_label, y = alpha_mean, color = prior)) +
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
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

if (is.null(nCov_plot)) {
  p_alpha <- p_alpha + facet_wrap(~ nCov, scales = "free_x", labeller = label_both)
}

ggsave(
  file.path(outDir, paste0("posterior_alpha_", nCov_tag, ".pdf")),
  p_alpha,
  width = 8,
  height = 4.5
)

## ------------------------------------------------------------
## ARI plot with 95% posterior credible interval
## ------------------------------------------------------------

p_ari <- ggplot(df, aes(x = scenario_label, y = ARI_mean, color = prior)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(
    aes(ymin = ARI_lwr, ymax = ARI_upr),
    position = position_dodge(width = 0.5),
    width = 0.2,
    alpha = 0.7
  ) +
  theme_minimal(base_size = 14) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Posterior ARI under learned alpha",
    x = NULL,
    y = expression(ARI(z, z[0]) ~ "|" ~ data)
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

  this_nCov <- as.integer(sub("^binarySBM_([0-9]+)cov_[A-Z]+$", "\\1", sim_name))

  if (!is.null(nCov_plot) && this_nCov != nCov_plot) next

  trace_dir <- file.path(outDir, "traces", sim_name)
  dir.create(trace_dir, recursive = TRUE, showWarnings = FALSE)

  for (nm in names(obj$results)) {
    plot_trace_pair(obj, sim_name, nm, trace_dir)
  }
}

cat("Done.\n")