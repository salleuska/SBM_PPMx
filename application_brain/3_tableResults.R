suppressPackageStartupMessages({
  library(here)
})

# ---- paths (edit if you want) ----
proc_path <- here("application_brain", "results", "processed",
                  "brain_xyz_processed_results.rds")
out_path  <- here("application_brain", "tables", "alpha_summary_table.tex")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

res_app <- readRDS(proc_path)

# ---- helpers ----
alpha_key <- function(a) paste0(as.integer(a[1]), as.integer(a[2]), as.integer(a[3]))

# Desired column order (8 combos)
order_keys <- c(
  "000", # all 0
  "100", # x only
  "010", # y only
  "001", # z only
  "110", # xy
  "101", # xz
  "011",  # yz
  "111" # all 1
)

col_titles <- c(
  "$(0,0,0)$",
  "$(1,0,0)$",
  "$(0,1,0)$",
  "$(0,0,1)$",
  "$(1,1,0)$",
  "$(1,0,1)$",
  "$(0,1,1)$",
  "$(1,1,1)$"
)

# extract one run per alpha combo (assuming unique)
runs <- res_app$results
a_mat <- t(sapply(runs, function(r) r$alpha[1:3]))
keys  <- apply(a_mat, 1, alpha_key)

# keep only 0/1 alpha runs
keep01 <- keys %in% order_keys
runs01 <- runs[keep01]
keys01 <- keys[keep01]

# map key -> run (first if duplicates)
key_to_run <- setNames(vector("list", length(order_keys)), order_keys)
for (k in order_keys) {
  idx <- which(keys01 == k)
  if (length(idx) > 0) key_to_run[[k]] <- runs01[[idx[1]]] else key_to_run[[k]] <- NULL
}

get_K <- function(r) if (is.null(r)) NA_integer_ else as.integer(r$K_hat)

get_sil_net <- function(r) {
  if (is.null(r)) return(NA_real_)
  # robust to different storage
  if (!is.null(r$silhouettes) && "net" %in% names(r$silhouettes)) return(as.numeric(r$silhouettes[["net"]]))
  if (!is.null(r$sil_net)) return(as.numeric(r$sil_net))
  NA_real_
}

K_row  <- sapply(order_keys, function(k) get_K(key_to_run[[k]]))
silrow <- sapply(order_keys, function(k) get_sil_net(key_to_run[[k]]))

# format
K_str   <- ifelse(is.na(K_row), "--", as.character(K_row))
sil_str <- ifelse(is.na(silrow), "--", sprintf("%.3f", silrow))

# ---- write LaTeX (booktabs) ----


tabular_spec <- paste0("|l|", paste(rep("c|", length(order_keys)), collapse = ""))

tab <- c(
  paste0("\\begin{tabular}{", tabular_spec, "}"),
  "\\hline",
  paste(c("$(\\alpha_1,\\alpha_2,\\alpha_3)$", col_titles), collapse = " & "), " \\\\",
  "\\hline",
  paste(c("$\\widehat{K}$", K_str), collapse = " & "), " \\\\",
  "\\hline",
  paste(c("Silhouette (network)", sil_str), collapse = " & "), " \\\\",
  "\\hline",
  "\\end{tabular}"
)

writeLines(tab, out_path)
cat("Wrote:", out_path, "\n")