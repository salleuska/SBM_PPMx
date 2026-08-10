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

nCov_plot <- if (is.null(args$nCov)) NULL else as.integer(args$nCov)

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
## Shared constants
## ------------------------------------------------------------
TYPE_COLORS    <- c(informative = "#2166AC", misleading = "#D6604D", neutral = "#878787")
TYPE_LEVELS    <- c("informative", "misleading", "neutral")
scenario_caption <- "N = neutral, I = informative, M = misleading"

map_cov_type <- function(covDep) {
  ifelse(covDep == "informative", "informative",
  ifelse(covDep == "mislead_random", "misleading", "neutral"))
}

parse_sim_name <- function(f) {
  nm <- sub("_results$", "", file_path_sans_ext(basename(f)))
  list(
    sim_name = nm,
    nCov     = as.integer(sub("^binarySBM_([0-9]+)cov_[A-Z0-9]+$", "\\1", nm)),
    scenario = sub("^binarySBM_[0-9]+cov_([A-Z0-9]+)$", "\\1", nm)
  )
}

## ------------------------------------------------------------
## Extract functions
## ------------------------------------------------------------

## ESBM results (learned alpha, with covariate)
extract_df <- function(f) {

  obj  <- readRDS(f)
  meta <- parse_sim_name(f)

  do.call(rbind, lapply(names(obj$results), function(nm) {

    r         <- obj$results[[nm]]
    alpha_mat <- if (is.null(dim(r$alpha_post))) matrix(r$alpha_post, nrow = 1) else r$alpha_post
    cov_type  <- map_cov_type(obj$sim$covDep)

    ARI_vec <- apply(r$Z_post, 2,
      function(z_t) mclust::adjustedRandIndex(z_t, obj$sim$partition_true))

    do.call(rbind, lapply(seq_len(nrow(alpha_mat)), function(j) {
      av <- as.numeric(alpha_mat[j, ])
      data.frame(
        sim_name       = meta$sim_name,
        nCov           = meta$nCov,
        scenario_label = meta$scenario,
        covariate      = paste0("x", j),
        cov_type       = cov_type[j],
        tag            = nm,
        seed           = r$seed,
        a_alpha        = r$alpha_g$a_alpha,
        b_alpha        = r$alpha_g$b_alpha,
        prior          = paste0("Gamma(", r$alpha_g$a_alpha, ",", r$alpha_g$b_alpha, ")"),
        model          = "ESBM",
        alpha_mean     = mean(av, na.rm = TRUE),
        alpha_lwr      = bayestestR::hdi(av, ci = 0.95)$CI_low,
        alpha_upr      = bayestestR::hdi(av, ci = 0.95)$CI_high,
        ARI_vi         = r$ARI_vi,
        ARI_mean       = mean(ARI_vec, na.rm = TRUE),
        ARI_lwr        = quantile(ARI_vec, 0.025, na.rm = TRUE),
        ARI_upr        = quantile(ARI_vec, 0.975, na.rm = TRUE),
        silhouette     = r$silhouette$mean
      )
    }))
  }))
}

## SBM baseline (no covariate)
extract_df_sbm <- function(f) {

  obj  <- readRDS(f)
  meta <- parse_sim_name(f)

  do.call(rbind, lapply(names(obj$results), function(nm) {

    r       <- obj$results[[nm]]
    ARI_vec <- apply(r$Z_post, 2,
      function(z_t) mclust::adjustedRandIndex(z_t, obj$sim$partition_true))

    data.frame(
      sim_name       = meta$sim_name,
      nCov           = meta$nCov,
      scenario_label = meta$scenario,
      covariate      = NA_character_,
      cov_type       = NA_character_,
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
  }))
}

## ------------------------------------------------------------
## Build plotting data
## ------------------------------------------------------------
df <- do.call(rbind, lapply(processed_files, extract_df))

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
  df_sbm <- df[0, ]
}

if (!is.null(nCov_plot)) {
  df     <- subset(df,     nCov == nCov_plot)
  df_sbm <- subset(df_sbm, nCov == nCov_plot)
}

if (nrow(df) == 0) stop("No ESBM results found for requested nCov.")

df_all <- rbind(df, df_sbm)

scen_levels           <- unique(df$scenario_label[order(df$nCov, df$scenario_label)])
df$scenario_label     <- factor(df$scenario_label,     levels = scen_levels)
df_all$scenario_label <- factor(df_all$scenario_label, levels = scen_levels)

nCov_tag <- if (is.null(nCov_plot)) "all" else paste0(nCov_plot, "Cov")

## ------------------------------------------------------------
## Posterior alpha plot
## ------------------------------------------------------------
many_cov <- !is.null(nCov_plot) && nCov_plot >= 10

df_alpha <- if (many_cov) df[df$nCov == nCov_plot, ] else df
df_alpha$cov_type <- factor(df_alpha$cov_type, levels = TYPE_LEVELS)

## x-axis: covariate index (many_cov) or scenario label (few_cov)
if (many_cov) {
  df_alpha$x_val <- as.integer(sub("^x", "", df_alpha$covariate))
  x_lab <- "Covariate index"
  title_str <- "Posterior of learned alpha (50 covariates)"
} else {
  df_alpha$x_val <- as.integer(df_alpha$scenario_label)
  x_lab <- "Scenario"
  title_str <- "Posterior distribution of learned alpha"
}

p_alpha <- ggplot(df_alpha, aes(x = x_val, y = alpha_mean, color = cov_type)) +
  geom_point(size = if (many_cov) 1.8 else 2.5) +
  geom_errorbar(aes(ymin = alpha_lwr, ymax = alpha_upr),
                width = 0, alpha = if (many_cov) 0.6 else 0.7) +
  scale_color_manual(values = TYPE_COLORS) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title    = title_str,
    subtitle = if (!many_cov) scenario_caption else NULL,
    x = x_lab,
    y = expression(alpha[g]),
    color = "Covariate type"
  ) +
  theme_minimal(base_size = if (many_cov) 13 else 14) +
  theme(plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

if (many_cov) {
  p_alpha <- p_alpha +
    scale_x_continuous(breaks = scales::breaks_pretty()) +
    facet_wrap(~ scenario_label, ncol = 1) +
    theme(strip.text = element_text(face = "bold"))
} else if (is.null(nCov_plot)) {
  p_alpha <- p_alpha +
    scale_x_continuous(breaks = seq_along(scen_levels), labels = scen_levels) +
    facet_wrap(~ nCov, scales = "free_x", labeller = label_both) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
} else {
  p_alpha <- p_alpha +
    scale_x_continuous(breaks = seq_along(scen_levels), labels = scen_levels) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
}

n_scen <- length(unique(df_alpha$scenario_label))
ggsave(
  file.path(outDir, paste0("posterior_alpha_", nCov_tag, ".pdf")),
  p_alpha,
  width  = if (many_cov) 10 else 8,
  height = if (many_cov) 3 * n_scen else 4.5
)

## ------------------------------------------------------------
## ARI plot
## ------------------------------------------------------------
df_ari_esbm <- df_all[df_all$model == "ESBM" & !duplicated(df_all$tag), ]
df_ari_sbm  <- df_all[df_all$model == "SBM", ]

p_ari <- ggplot(mapping = aes(x = scenario_label)) +
  geom_ribbon(data = df_ari_sbm,
    aes(ymin = ARI_lwr, ymax = ARI_upr, group = nCov),
    fill = "grey60", alpha = 0.25) +
  geom_line(data = df_ari_sbm,
    aes(y = ARI_mean, group = nCov),
    color = "grey40", linewidth = 0.8, linetype = "dashed") +
  geom_errorbar(data = df_ari_esbm,
    aes(ymin = ARI_lwr, ymax = ARI_upr, color = prior),
    position = position_dodge(width = 0.5), width = 0.2, alpha = 0.7) +
  geom_point(data = df_ari_esbm,
    aes(y = ARI_mean, color = prior),
    position = position_dodge(width = 0.5), size = 2.5) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title    = "Posterior ARI: ESBM vs SBM baseline",
    subtitle = scenario_caption,
    x = NULL,
    y = expression(ARI(z, z[0]) ~ "|" ~ data),
    color = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1))

if (is.null(nCov_plot)) {
  p_ari <- p_ari + facet_wrap(~ nCov, scales = "free_x", labeller = label_both)
}

ggsave(
  file.path(outDir, paste0("ARI_mean_CI_learnedAlpha_", nCov_tag, ".pdf")),
  p_ari, width = 8, height = 4.5
)

## ------------------------------------------------------------
## Trace plots
## ------------------------------------------------------------

make_trace_df <- function(obj, scenario, result_name) {

  r        <- obj$results[[result_name]]
  alpha_mat <- if (is.null(dim(r$alpha_post))) matrix(r$alpha_post, nrow = 1) else r$alpha_post
  rate_mat  <- if (is.null(dim(r$alpha_rate_post))) matrix(r$alpha_rate_post, nrow = 1) else r$alpha_rate_post
  b_alpha  <- r$alpha_g$b_alpha
  cov_type <- map_cov_type(obj$sim$covDep[seq_len(nrow(alpha_mat))])

  J      <- nrow(alpha_mat)
  N_iter <- ncol(alpha_mat)

  do.call(rbind, lapply(seq_len(J), function(j) {
    av <- as.numeric(alpha_mat[j, ])
    rv <- as.numeric(rate_mat[j, ])
    data.frame(
      iter       = seq_len(N_iter),
      alpha      = av,
      alpha_rm   = cumsum(av) / seq_len(N_iter),
      similarity = rv - b_alpha,
      sim_rm     = cumsum(rv - b_alpha) / seq_len(N_iter),
      covariate  = paste0("x", j),
      cov_type   = cov_type[j],
      scenario   = scenario,
      tag        = result_name
    )
  }))
}

plot_trace_pair <- function(obj, scenario, result_name, outDir) {

  dat <- make_trace_df(obj, scenario, result_name)
  dat$cov_type <- factor(dat$cov_type, levels = TYPE_LEVELS)

  J      <- length(unique(dat$covariate))
  n_cols <- min(J, 5)
  n_rows <- ceiling(J / n_cols)

  mu_df     <- aggregate(alpha      ~ covariate + cov_type, dat, mean)
  mu_sim_df <- aggregate(similarity ~ covariate + cov_type, dat, mean)

  trace_theme <- list(
    scale_color_manual(values = TYPE_COLORS),
    facet_wrap(~ covariate, ncol = n_cols, scales = "free_y"),
    theme_minimal(base_size = 11),
    theme(strip.text = element_text(face = "bold"), legend.position = "bottom")
  )

  p_alpha_trace <- ggplot(dat, aes(x = iter)) +
    geom_line(aes(y = alpha, color = cov_type), alpha = 0.3, linewidth = 0.3) +
    geom_line(aes(y = alpha_rm, color = cov_type), linewidth = 0.8) +
    geom_hline(data = mu_df, aes(yintercept = alpha, color = cov_type),
               linetype = 2, linewidth = 0.5) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(title = paste0(scenario, ": alpha traces"),
         x = NULL, y = expression(alpha[g]), color = "Covariate type") +
    trace_theme

  p_sim_trace <- ggplot(dat, aes(x = iter)) +
    geom_line(aes(y = similarity, color = cov_type), alpha = 0.3, linewidth = 0.3) +
    geom_line(aes(y = sim_rm, color = cov_type), linewidth = 0.8) +
    geom_hline(data = mu_sim_df, aes(yintercept = similarity, color = cov_type),
               linetype = 2, linewidth = 0.5) +
    labs(title = paste0(scenario, ": similarity traces"),
         x = "Iteration", y = "Similarity", color = "Covariate type") +
    trace_theme

  ggsave(
    file.path(outDir, paste0("trace_alpha_similarity_", scenario, "_", result_name, ".pdf")),
    cowplot::plot_grid(p_alpha_trace, p_sim_trace, ncol = 1, align = "v"),
    width      = n_cols * 3,
    height     = n_rows * 5,
    limitsize  = FALSE
  )
}

for (f in processed_files) {

  meta <- parse_sim_name(f)
  if (!is.null(nCov_plot) && meta$nCov != nCov_plot) next

  obj       <- readRDS(f)
  trace_dir <- file.path(outDir, "traces", meta$sim_name)
  dir.create(trace_dir, recursive = TRUE, showWarnings = FALSE)

  for (nm in names(obj$results)) {
    plot_trace_pair(obj, meta$sim_name, nm, trace_dir)
  }
}

cat("Done.\n")
