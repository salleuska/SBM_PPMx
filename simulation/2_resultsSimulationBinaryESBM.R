##---------------------------##------------------##
## Process ESBM results for binary SBM simulation
##- Reads all post_*.rds in results folder
##- For each:
##     * computes Dahl clustering from Z_DM
##     * computes ARI vs true partition
##     * records scenario, alpha, seed, K_hat
##- Saves a summary CSV and an optional ARI vs alpha plot
##---------------------------##------------------##

library(here)
library(mclust)     # ARI
library(ggplot2)    # plotting
library(tools)      # file_path_sans_ext
library(salso)      # point estimate for the clustering

##---------------------------##
## Settings
##---------------------------##

results_dir <- here("simulation", "results", "oneCov")

## if two covariates
#results_dir <- here("simulation", "results", "twoCov")

burn_frac <- 0.2 # brunin
thin_step <- 1   # eventual thinning

##---------------------------##
## List result files
##---------------------------##

files_res <- list.files(
  results_dir,
  pattern = "^post.*\\.rds$",
  full.names = TRUE
)

if (length(files_res) == 0L) {
  stop("No result files found in ", results_dir)
}

cat("Found", length(files_res), "result files.\n")

##---------------------------##
## Process each result
##---------------------------##

summ_list <- vector("list", length(files_res))
psm_list  <- vector("list", length(files_res))  # store PSMs here

for (i in seq_along(files_res)) {
  f_res <- files_res[i]
  cat("Processing:", basename(f_res), "...\n")

  res <- readRDS(f_res)

  # Extract posterior samples and metadata
  Z_post   <- res$Z_post      # V x N_iter
  alpha    <- res$alpha
  seed     <- res$seed
  scen_res <- res$scenario
  sim_path <- res$file        # path to original simulation .rds

  sim_obj <- readRDS(here(sim_path))
  z_true  <- sim_obj$partition
  V_true  <- length(z_true)

  ##---------------------------####
  ## estimate clustering via salso
  ##---------------------------####
  ## keep iterations after burnin
  T_total <- ncol(Z_post)
  burn    <- floor(T_total * burn_frac)
  keep_i  <- seq(from = burn + 1L, to = T_total, by = thin_step)

  Z_keep <- Z_post[, keep_i, drop = FALSE]   # V x T_keep

  # Posterior similarity matrix (PSM)
  coClust_prob <- psm(t(Z_keep))    # iterations x items → transpose

  # VI representative clustering
  z_hat <- salso(t(Z_keep), loss = VI())

  ##---------------------------##
  ## ARI vs true clustering
  ##---------------------------##

  ari <- mclust::adjustedRandIndex(z_hat, z_true)

  # Number of inferred clusters
  K_hat <- length(unique(z_hat))

  # Extract alpha tag & base name from filename (optional)
  base_name <- file_path_sans_ext(basename(f_res))

  # store summary row
  summ_list[[i]] <- data.frame(
    file_res   = basename(f_res),
    sim_file   = sim_path,
    base_name  = base_name,
    scenario   = scen_res,
    alpha      = alpha,
    seed       = seed,
    K_hat      = K_hat,
    ARI        = ari,
    stringsAsFactors = FALSE
  )

  # store PSM in list
  psm_list[[i]] <- coClust_prob
}

summary_df <- do.call(rbind, summ_list)
names(psm_list) <- basename(files_res)

##---------------------------##
## Save summary and PSMs
##---------------------------##

out_csv <- file.path(results_dir, "summary_binaryESBM.csv")
write.csv(summary_df, out_csv, row.names = FALSE)
cat("Saved summary to:\n  ", out_csv, "\n")

psm_out <- file.path(results_dir, "psm_binaryESBM.rds")
saveRDS(psm_list, psm_out)
cat("Saved PSM list to:\n  ", psm_out, "\n")

##---------------------------##
## ARI vs alpha plot
##---------------------------##

p <- ggplot(summary_df,
            aes(x = alpha, y = ARI,
                color = scenario,
                group = interaction(scenario, seed))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  labs(title = "ARI vs alpha",
       x = expression(alpha),
       y = "Adjusted Rand Index")

print(p)

psm_list <- readRDS(here(results_dir, "psm_binaryESBM.rds"))
names(psm_list)
# pick one:
psm_both_a1  <- psm_list[[5]]
image(psm_both_a1)  # or use your own heatmap code
