library(here)
library(ggplot2)
library(RColorBrewer)

plot_psm_gg <- function(psm, title = "") {
  V <- nrow(psm)
  
  # long data.frame in base R
  grid_idx <- expand.grid(
    row = seq_len(V),
    col = seq_len(V)
  )
  grid_idx$value <- as.vector(psm)
  
  # flip rows so origin is bottom-left like a matrix
  grid_idx$row_f <- factor(grid_idx$row, levels = rev(seq_len(V)))
  grid_idx$col_f <- factor(grid_idx$col, levels = seq_len(V))
  
  pal <- colorRampPalette(brewer.pal(9, "Blues"))(50)
  
  ggplot(grid_idx, aes(x = col_f, y = row_f, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colors = pal, limits = c(0, 1)) +
    coord_fixed() +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      axis.title   = element_blank(),
      panel.grid   = element_blank(),
      plot.title   = element_text(hjust = 0.5)
    ) +
    ggtitle(title)
}


p <- ggplot(summary_df,
            aes(x = alpha, y = ARI,
                color = scenario,
                group = interaction(scenario, seed))) +
  geom_point() +
  geom_line() +
  theme_minimal() +
  ylim(c(0,1)) + 
  labs(title = "ARI vs alpha",
       x = expression(alpha),
       y = "Adjusted Rand Index")


ggsave(p, file = here(results_dir, "alphaVsARI.pdf"))

psm_list <- readRDS(here(results_dir, "psm_binaryESBM.rds"))
psm_pdf <- file.path(results_dir, "psm_heatmaps_binaryESBM.pdf")
pdf(psm_pdf, width = 5, height = 5)

for (i in seq_len(nrow(summary_df))) {
  file_res_i <- summary_df$file_res[i]
  scen_i     <- summary_df$scenario[i]
  alpha_i    <- summary_df$alpha[i]
  seed_i     <- summary_df$seed[i]
  
  psm_i <- psm_list[[file_res_i]]
  
  title_i <- sprintf("Scenario: %s\nalpha = %.2f, seed = %d",
                     scen_i, alpha_i, seed_i)
  
  print(plot_psm_gg(psm_i, title_i))
}