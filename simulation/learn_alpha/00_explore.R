library(ggplot2)
library(here)

source(here("R_utilties", "similarity_functions.R"))

set.seed(123)

outDir <- here("simulation", "learn_alpha", "diagnostics")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

make_Z <- function(z) {
  Z <- matrix(0, length(z), max(z))
  for (i in seq_along(z)) Z[i, z[i]] <- 1
  Z
}

compute_alpha_summary <- function(Z, x, calibration, a_alpha, b_alpha, J = 1) {
  x <- as.numeric(scale(x))

  log_g_h <- similarity_ppmx_gaussian_mean(
    mode = "partition",
    Z = Z,
    x = x,
    args = list(m0 = 0, s0 = 1)
  )

  if (calibration == "normalized") {
    lse <- max(log_g_h) + log(sum(exp(log_g_h - max(log_g_h))))
    log_g_h <- log_g_h - lse
    stat <- -sum(log_g_h)
  }

  if (calibration == "geometric") {
    stat <- -sum(log_g_h) / J
  }

  rate <- b_alpha + stat

  data.frame(
    calibration = calibration,
    rate = rate,
    alpha_mean = a_alpha / rate,
    alpha_var = a_alpha / rate^2,
    stat = stat
  )
}

n_per_cluster <- c(50, 40, 30, 30)
K <- length(n_per_cluster)
z_true <- rep(seq_len(K), times = n_per_cluster)
Z_true <- make_Z(z_true)

z_wrong <- sample(z_true)
Z_wrong <- make_Z(z_wrong)

mu <- c(-1.5, -0.5, 0.5, 1.5)
sigma <- 0.5

x_informative <- rnorm(length(z_true), mean = mu[z_true], sd = sigma)
x_neutral <- rnorm(length(z_true), mean = 0, sd = 1)
x_misleading <- rnorm(length(z_true), mean = mu[z_wrong], sd = sigma)

cases <- list(
  informative_trueZ = list(x = x_informative, Z = Z_true),
  neutral_trueZ = list(x = x_neutral, Z = Z_true),
  misleading_trueZ = list(x = x_misleading, Z = Z_true),
  informative_wrongZ = list(x = x_informative, Z = Z_wrong)
)

priors <- data.frame(
  prior = c(
    "Gamma(mean=0.5,var=0.25)",
    "Gamma(mean=1,var=0.25)",
    "Gamma(mean=2,var=1)"
  ),
  m0 = c(0.5, 1, 2),
  v0 = c(0.25, 0.25, 1)
)

priors$a_alpha <- priors$m0^2 / priors$v0
priors$b_alpha <- priors$m0 / priors$v0

df <- do.call(
  rbind,
  lapply(names(cases), function(case_name) {
    do.call(
      rbind,
      lapply(seq_len(nrow(priors)), function(i) {
        do.call(
          rbind,
          lapply(c("normalized", "geometric"), function(cal) {
            out <- compute_alpha_summary(
              Z = cases[[case_name]]$Z,
              x = cases[[case_name]]$x,
              calibration = cal,
              a_alpha = priors$a_alpha[i],
              b_alpha = priors$b_alpha[i],
              J = 1
            )

            out$case <- case_name
            out$prior <- priors$prior[i]
            out$a_alpha <- priors$a_alpha[i]
            out$b_alpha <- priors$b_alpha[i]
            out
          })
        )
      })
    )
  })
)

df$case <- factor(df$case, levels = names(cases))
df$prior <- factor(df$prior, levels = priors$prior)

saveRDS(df, file.path(outDir, "alpha_covariate_diagnostic.rds"))

p_rate <- ggplot(df, aes(x = case, y = rate, fill = calibration)) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_wrap(~ prior, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Posterior rate for alpha under different Gamma priors",
    x = NULL,
    y = "Gamma posterior rate"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

p_rate

ggsave(
  file.path(outDir, "alpha_rate_covariate_diagnostic_by_prior.pdf"),
  p_rate,
  width = 11,
  height = 4
)

p_mean <- ggplot(df, aes(x = case, y = alpha_mean, fill = calibration)) +
  geom_col(position = position_dodge(width = 0.7)) +
  facet_wrap(~ prior, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Implied posterior mean of alpha under different Gamma priors",
    x = NULL,
    y = expression(E(alpha ~ "|" ~ Z * "," * X))
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold"),
    legend.title = element_blank()
  )

p_mean

ggsave(
  file.path(outDir, "alpha_mean_covariate_diagnostic_by_prior.pdf"),
  p_mean,
  width = 11,
  height = 4
)

print(df)