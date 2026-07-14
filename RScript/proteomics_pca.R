# =============================================================================
# proteomics_pca.R
# -----------------------------------------------------------------------------
# PCA quality control for quantitative proteomics data.
# Generates per-axis 2D PCA plots coloured by batch, feeder, and replication.
#
# Output: Figures/PCA_<cell_type>_<n>_axis.pdf
#
# Usage:
#   Rscript proteomics_pca.R
#   Or set parameters below and source() in RStudio.
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/proteomics-quantitative-analysis
# =============================================================================

source("proteomics_utils.R")
.load_pkg("ggplot2")


# =============================================================================
# PARAMETERS — edit here
# =============================================================================

DATA_DIR     <- "."
EXPR_FILE    <- file.path(DATA_DIR, "data", "dataMatrix_wash_groups_newSettings.tsv")
ANNOT_FILE   <- file.path(DATA_DIR, "data", "SampleAnnot.tsv")
FIGURES_DIR  <- file.path(DATA_DIR, "Figures")

N_PCA_AXES   <- 4       # Number of PCA axes to plot pairwise

# Sample annotation columns to colour PCA by
COLOR_VARS   <- list(
  Batch       = "Batch",
  Feeder      = "Feader",
  Replication = "Geminin"
)

# Cell type filter (must match values in sampleAnnot$Tissue)
CELL_TYPES   <- c("Sperm", "Progenitor")


# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading data...")
expr   <- read_tsv(EXPR_FILE)
annot  <- read_tsv(ANNOT_FILE)
check_sample_names(expr, annot)

ordered_expr <- expr[order(rownames(expr)), rownames(annot)]
create_dir(FIGURES_DIR)


# =============================================================================
# PCA FUNCTION
# =============================================================================

#' Run PCA and produce 2D pair plots for a subset of samples
#'
#' @param expr_mat    Matrix (genes x samples). Expression data.
#' @param annot_sub   data.frame. Sample annotation for selected samples.
#' @param color_vars  Named list. Annotation column names to colour by.
#' @param n_axes      Integer. Number of axes to plot pairwise.
#' @param label       Character. Label for output filename.
#' @param out_dir     Character. Output directory.
run_pca_plots <- function(expr_mat, annot_sub, color_vars,
                           n_axes, label, out_dir) {
  message("Running PCA for: ", label)

  # Standard PCA via prcomp (base R — no dependency on home functions)
  pca_res    <- prcomp(t(expr_mat), center = TRUE, scale. = FALSE)
  pct_var    <- (pca_res$sdev^2 / sum(pca_res$sdev^2)) * 100
  scores     <- as.data.frame(pca_res$x[, seq_len(n_axes)])
  scores     <- cbind(scores, annot_sub[rownames(scores), , drop = FALSE])

  pdf_path   <- file.path(out_dir, paste0("PCA_", label, "_", n_axes, "_axes.pdf"))
  pdf(pdf_path, width = 10, height = 10)
  on.exit(dev.off())

  # Scree plot
  barplot(pct_var[seq_len(n_axes)],
          names.arg = paste0("PC", seq_len(n_axes), "\n",
                              round(pct_var[seq_len(n_axes)], 1), "%"),
          main = paste("PCA variance explained -", label),
          ylab = "% Variance", col = "steelblue")

  # Pairwise 2D scatter plots
  for (i in seq_len(n_axes - 1)) {
    for (j in (i + 1):n_axes) {
      pc_x   <- paste0("PC", i)
      pc_y   <- paste0("PC", j)
      x_lab  <- paste0(pc_x, " (", round(pct_var[i], 1), "%)")
      y_lab  <- paste0(pc_y, " (", round(pct_var[j], 1), "%)")

      for (var_label in names(color_vars)) {
        col_col <- color_vars[[var_label]]
        if (!col_col %in% colnames(scores)) next
        grp     <- as.factor(scores[[col_col]])

        p <- ggplot(scores, aes_string(x = pc_x, y = pc_y,
                                        colour = col_col, label = "rownames(scores)")) +
          geom_point(size = 3) +
          geom_text(size = 2.5, vjust = -0.5, show.legend = FALSE) +
          scale_colour_brewer(palette = "Set1", name = var_label) +
          labs(title = paste(label, "-", var_label, ":", pc_x, "vs", pc_y),
               x = x_lab, y = y_lab) +
          theme_bw(base_size = 11) +
          coord_fixed()
        print(p)
      }
    }
  }
  message("  Saved: ", pdf_path)
}


# =============================================================================
# RUN PCA PER CELL TYPE
# =============================================================================

for (ct in CELL_TYPES) {
  ct_samples  <- rownames(annot)[annot$Tissue == ct]
  if (length(ct_samples) == 0) {
    warning("No samples found for cell type: ", ct, " — skipping.")
    next
  }
  ct_expr   <- ordered_expr[, ct_samples, drop = FALSE]
  ct_annot  <- annot[ct_samples, , drop = FALSE]
  run_pca_plots(ct_expr, ct_annot, COLOR_VARS, N_PCA_AXES, ct, FIGURES_DIR)
}

message("PCA complete. Figures in: ", FIGURES_DIR)
