suppressPackageStartupMessages({
  library(here)
  library(mclust)  # ARI
})

proc_path <- here("application_brain", "results", "processed",
                  "brain_xyz_processed_results.rds")
out_csv <- here("application_brain", "results", "processed",
                "brain_distance_to_baseline.csv")

res <- readRDS(proc_path)

df <- data.frame(
  tag     = names(res$results),
  alpha_x = sapply(res$results, \(r) r$alpha[1]),
  alpha_y = sapply(res$results, \(r) r$alpha[2]),
  alpha_z = sapply(res$results, \(r) r$alpha[3]),
  K_hat   = sapply(res$results, \(r) r$K_hat),
  stringsAsFactors = FALSE
)

Z <- lapply(res$results, \(r) as.integer(r$z_hat))

# baseline alpha = (0,0,0)
tol <- 1e-8
i0 <- which(abs(df$alpha_x) < tol & abs(df$alpha_y) < tol & abs(df$alpha_z) < tol)
if (!length(i0)) stop("No baseline run found with alpha=(0,0,0).")
z0 <- Z[[i0[1]]]

# normalized VI: VI / log(n)  (entropy/MI use natural logs)
H  <- \(z) { p <- as.numeric(table(z))/length(z); -sum(p * log(p)) }
MI <- \(a,b){ tab<-table(a,b); p<-tab/sum(tab); p1<-rowSums(p); p2<-colSums(p)
              sum(ifelse(p>0, p*log(p/(p1 %o% p2)), 0)) }
VI <- \(a,b) H(a) + H(b) - 2*MI(a,b)

df$ARI_to_base <- vapply(Z, \(z) adjustedRandIndex(z, z0), numeric(1))
df$nVI_to_base <- vapply(Z, \(z) VI(z, z0)/log(length(z0)), numeric(1))

write.csv(df[order(df$alpha_x, df$alpha_y, df$alpha_z), ], out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n")