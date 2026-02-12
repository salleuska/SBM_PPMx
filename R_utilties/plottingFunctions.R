#' Plot a matrix (binary adjacency or continuous PSM) with optional annotations
#'
#' Plots either (i) a binary adjacency matrix (0/1) using a white/black palette,
#' or (ii) a posterior similarity matrix (PSM) on [0,1] using a continuous
#' gradient. Ordering is assumed to be already encoded in the input matrix.
#' Optional row annotations (e.g., partitions) can be displayed, and gaps may be
#' introduced according to a supplied grouping vector.
#'
#' @param M Numeric matrix. Either adjacency (0/1) or PSM (values in [0,1]).
#' @param type Character in c("psm","adj"). Controls palette and breaks.
#' @param ann_df Optional data.frame of annotation columns (each column a partition).
#' @param show_cols Optional character vector selecting which columns of ann_df to display.
#' @param ann_colors Optional named list of color vectors for annotation columns.
#' @param gaps_by Optional vector defining cluster boundaries for visual gaps (in plotting order).
#' @param mat_palette Character. RColorBrewer palette name for continuous matrices (default "Greys").
#' @param n_mat_cols Integer. Number of colors for the continuous gradient.
#' @param border_color Cell border color (default "grey90").
#' @param filename Optional file path; if provided, pheatmap writes directly to file.
#'
#' @return A pheatmap object (invisibly).
#' @export

#' Plot a matrix (binary adjacency or continuous PSM) with optional annotations
#'
#' Plots either (i) a binary adjacency matrix (0/1) using a white/black palette,
#' or (ii) a posterior similarity matrix (PSM) on [0,1] using a continuous
#' gradient. Ordering is assumed to be already encoded in the input matrix.
#' Optional row annotations (e.g., partitions) can be displayed, and gaps may be
#' introduced according to a supplied grouping vector.
#'
#' @param M Numeric matrix. Either adjacency (0/1) or PSM (values in [0,1]).
#' @param type Character in c("psm","adj"). Controls palette and breaks.
#' @param ann_df Optional data.frame of annotation columns (each column a partition).
#' @param show_cols Optional character vector selecting which columns of ann_df to display.
#' @param ann_colors Optional named list of color vectors for annotation columns.
#' @param gaps_by Optional vector defining cluster boundaries for visual gaps (in plotting order).
#' @param mat_palette Character. RColorBrewer palette name for continuous matrices (default "Greys").
#' @param n_mat_cols Integer. Number of colors for the continuous gradient (PSM).
#' @param border_color Cell border color (default "grey90").
#' @param main Optional character title passed to \code{pheatmap}.
#' @param filename Optional file path; if provided, pheatmap writes directly to file.
#'
#' @return A pheatmap object (invisibly).
#' @export
plot_matrix_with_annotations <- function(M,
                                         type = c("psm", "adj"),
                                         ann_df = NULL,
                                         show_cols = NULL,
                                         ann_colors = NULL,
                                         gaps_by = NULL,
                                         mat_palette = "Greys",
                                         n_mat_cols = 200,
                                         border_color = "grey90",
                                         main = NULL,
                                         filename = NULL) {

  require(pheatmap)
  require(RColorBrewer)

  if (is.null(main)) main <- NA_character_
  
  type <- match.arg(type)

  n <- nrow(M)
  rn <- paste0("i", seq_len(n))
  rownames(M) <- rn
  colnames(M) <- rn

  ## -------------------------
  ## Annotation handling
  ## -------------------------
  if (!is.null(ann_df)) {
    ann_df <- as.data.frame(ann_df)
    rownames(ann_df) <- rn

    if (!is.null(show_cols)) {
      ann_df <- ann_df[, intersect(show_cols, colnames(ann_df)), drop = FALSE]
    }
    if (ncol(ann_df) == 0) ann_df <- NULL
  }

  ## -------------------------
  ## Default annotation colors
  ## -------------------------
  if (!is.null(ann_df)) {
    if (is.null(ann_colors)) ann_colors <- list()

    for (nm in colnames(ann_df)) {
      if (!nm %in% names(ann_colors)) {
        lev <- levels(factor(ann_df[[nm]]))
        base_cols <- brewer.pal(min(12, max(3, length(lev))), "Set3")
        cols <- rep(base_cols, length.out = length(lev))
        names(cols) <- lev
        ann_colors[[nm]] <- cols
      }
    }
  }

  ## -------------------------
  ## Gaps between groups
  ## -------------------------
  gaps <- NULL
  if (!is.null(gaps_by)) {
    gaps <- which(diff(as.integer(factor(gaps_by))) != 0)
  }

  ## -------------------------
  ## Palette + breaks
  ## -------------------------
  if (type == "adj") {
    mat_cols <- c("white", "black")
    breaks   <- c(-0.5, 0.5, 1.5)
  } else {
    # Clamp, then force a dense 0–1 mapping (appears continuous)
    M[M < 0] <- 0
    M[M > 1] <- 1

    mat_cols <- colorRampPalette(brewer.pal(9, mat_palette))(n_mat_cols)
    breaks   <- seq(0, 1, length.out = n_mat_cols + 1)
  }

  ## -------------------------
  ## Plot
  ## -------------------------
  pheatmap::pheatmap(
    M,
    color             = mat_cols,
    breaks            = breaks,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    show_rownames     = FALSE,
    show_colnames     = FALSE,
    annotation_row    = ann_df,
    annotation_colors = ann_colors,
    border_color      = border_color,
    annotation_legend = FALSE,
    legend            = FALSE,
    gaps_row          = gaps,
    gaps_col          = gaps,
    main              = main,
    filename          = filename
  )
}
