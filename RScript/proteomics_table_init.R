# =============================================================================
# proteomics_table_init.R
# -----------------------------------------------------------------------------
# Initialise comparison matrices from differential expression TSV files.
#
# Outputs:
#   results/ComparisonTable.tsv        log2FC per protein per comparison
#   results/BinaryComparisonTable.tsv  +1 / -1 / 0 per protein per comparison
#
# Usage:
#   Rscript proteomics_table_init.R <data_dir>
#   Or set DATA_DIR below and source() in RStudio.
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/proteomics-quantitative-analysis
# =============================================================================

source("proteomics_utils.R")

# =============================================================================
# PARAMETERS
# =============================================================================

# Root data directory (can be overridden by command-line argument)
DATA_DIR <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(DATA_DIR) || DATA_DIR == "") DATA_DIR <- "."

EXPR_FILE       <- file.path(DATA_DIR, "data", "dataMatrix.tsv")
COMPARISONS_DIR <- file.path(DATA_DIR, "Comparaison")
RESULTS_DIR     <- file.path(DATA_DIR, "results")

# Expected column names in comparison files
COMP_COLS <- c("Protein", "AVG_Log2_Ratio", "Absolute_AVG_Log2_Ratio",
               "PValue", "QValue", "PctRatios", "UniProtIds", "Genes",
               "ProteinName", "PctUniqueTotalPeptide", "PctChange", "Ratio")

# =============================================================================
# MAIN
# =============================================================================

message("=== Comparison Table Initialisation ===")

# Discover comparison files
filenames <- list.files(COMPARISONS_DIR, pattern = "\\.tsv$",
                        full.names = TRUE, recursive = FALSE)
if (length(filenames) == 0)
  stop("No .tsv files found in: ", COMPARISONS_DIR)
message("Comparison files found: ", length(filenames))

# Load protein background
message("Loading expression matrix...")
prot_names <- rownames(read_tsv(EXPR_FILE))
message("Total proteins: ", length(prot_names))

# Initialise empty matrix
comp_labels <- sapply(filenames, parse_comparison_label)
comp_table  <- matrix(0, nrow = length(prot_names), ncol = length(filenames),
                      dimnames = list(sort(prot_names), comp_labels))

# Fill matrix with AVG_Log2_Ratio per comparison
message("Filling comparison matrix...")
for (i in seq_along(filenames)) {
  comp      <- read.table(filenames[i], header = FALSE,
                          stringsAsFactors = FALSE)
  colnames(comp) <- COMP_COLS
  rownames(comp) <- comp$Protein
  comp           <- comp[order(rownames(comp)), ]
  shared         <- intersect(rownames(comp), rownames(comp_table))
  comp_table[shared, i] <- comp[shared, "AVG_Log2_Ratio"]
}

# Export
create_dir(RESULTS_DIR)
comp_df        <- as.data.frame(comp_table)
binary_df      <- sign(comp_df)   # +1 / -1 / 0

write_tsv(comp_df,   file.path(RESULTS_DIR, "ComparisonTable.tsv"))
write_tsv(binary_df, file.path(RESULTS_DIR, "BinaryComparisonTable.tsv"))

message("Done. Output in: ", RESULTS_DIR)
