# =============================================================================
# proteomics_utils.R
# -----------------------------------------------------------------------------
# Shared utility functions for quantitative proteomics analysis.
# Sourced by all other scripts in this project.
#
# Provides: read_tsv(), write_tsv(), check_sample_names(),
#           scale_rows(), pdf_height_from_nrow()
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/proteomics-quantitative-analysis
# Project : Quantitative proteomics — sperm vs progenitor cells
# Date    : 2023 (CR2TI, UMR 1064, Nantes Universite)
# =============================================================================

.load_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' not found. Install with install.packages('", pkg, "')")
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

#' Read a tab-separated file with row names
#' @param filepath Character. Path to .tsv file.
#' @return data.frame.
read_tsv <- function(filepath) {
  if (!file.exists(filepath)) stop("File not found: ", filepath)
  read.table(filepath, header = TRUE, sep = "\t",
             row.names = 1, stringsAsFactors = FALSE, quote = "")
}

#' Write a data frame to a tab-separated file
#' @param x data.frame.
#' @param filepath Character. Output path.
write_tsv <- function(x, filepath) {
  write.table(x, file = filepath, sep = "\t",
              col.names = NA, row.names = TRUE, quote = FALSE)
  message("Saved: ", filepath)
}

#' Check and report sample name mismatches between expression matrix and annotation
#' @param expr data.frame. Expression matrix (genes x samples).
#' @param annot data.frame. Sample annotation (samples x variables).
#' @return Invisible NULL. Warns on mismatch.
check_sample_names <- function(expr, annot) {
  missing_in_expr <- rownames(annot)[!rownames(annot) %in% colnames(expr)]
  missing_in_annot <- colnames(expr)[!colnames(expr) %in% rownames(annot)]
  if (length(missing_in_expr) > 0) {
    warning("Samples in annotation but not in expression matrix:\n",
            paste(missing_in_expr, collapse = ", "))
  }
  if (length(missing_in_annot) > 0) {
    warning("Samples in expression matrix but not in annotation (will be removed):\n",
            paste(missing_in_annot, collapse = ", "))
  }
  if (length(missing_in_expr) == 0 && length(missing_in_annot) == 0)
    message("Sample names: OK")
  invisible(NULL)
}

#' Z-score scale rows of a matrix (mean-center + scale by SD)
#' @param mat Matrix or data.frame (genes x samples).
#' @param center Logical. Subtract row mean. Default TRUE.
#' @param scale  Logical. Divide by row SD. Default TRUE.
#' @return Matrix of same dimensions.
scale_rows <- function(mat, center = TRUE, scale = TRUE) {
  t(scale(t(as.matrix(mat)), center = center, scale = scale))
}

#' Compute adaptive PDF height based on number of features
#' @param n Integer. Number of rows (genes/proteins).
#' @return Numeric. PDF height in inches.
pdf_height_from_nrow <- function(n) {
  if      (n < 400)  10
  else if (n < 1000) 15
  else if (n < 2000) 20
  else if (n < 4000) 40
  else if (n < 6000) 50
  else               60
}

#' Parse comparison filename to short label
#' @param filepath Character. Full path to comparison TSV.
#' @return Character. Short label (filename without extension).
parse_comparison_label <- function(filepath) {
  base <- basename(filepath)
  sub("\\.tsv$", "", base, ignore.case = TRUE)
}

#' Create directory if it does not exist
#' @param path Character. Directory path.
create_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}
