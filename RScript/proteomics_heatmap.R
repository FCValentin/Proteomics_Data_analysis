# =============================================================================
# proteomics_heatmap.R
# -----------------------------------------------------------------------------
# ComplexHeatmap visualisation for quantitative proteomics data.
# Plots a z-score scaled heatmap for a user-defined protein list,
# with hierarchical clustering of proteins (pearson distance).
#
# Output: Figures/(<n> Proteins).pdf
#         results/(<n> Proteins).tsv
#
# Usage:
#   Rscript proteomics_heatmap.R
#   Or set parameters below and source() in RStudio.
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/proteomics-quantitative-analysis
# =============================================================================

source("proteomics_utils.R")
.load_pkg("circlize")
.load_pkg("ComplexHeatmap")


# =============================================================================
# PARAMETERS — edit here
# =============================================================================

DATA_DIR      <- "."
EXPR_FILE     <- file.path(DATA_DIR, "data", "dataMatrix_wash_groups_newSettings.tsv")
ANNOT_FILE    <- file.path(DATA_DIR, "data", "SampleAnnot.tsv")
FIGURES_DIR   <- file.path(DATA_DIR, "Figures")
RESULTS_DIR   <- file.path(DATA_DIR, "results")

# Path to protein list of interest (one protein ID per line, no header)
PROT_LIST_FILE <- file.path(DATA_DIR,
  "ListOfDemandes/list_orthologs_EPI_and_TF_SEGM_DE/list_orthologs_EPI_and_TF_SEGM_DE.txt")

# Comparison file used to retrieve gene names
COMPARISON_FILE <- file.path(DATA_DIR, "Comparaison", "Sperm_prog.tsv")

# Column selection: which samples to include
# 1 = Prog / Sperm / proggem- / spermgem-
# 2 = Prog / spermgem- / proggem- / spermgem+ / proggem+
# 3 = proggem- / spermgem- / proggem+ / spermgem+
# 4 = sperm / sperm DE gem-
# 5 = sperm / sperm DE gem- / sperm DE gem+
# 6 = prog / prog DE gem-
# 7 = prog / prog DE gem- / prog DE gem+
# 8 = Prog / Sperm
# 9 = custom (set CUSTOM_COLUMNS below)
WHAT_TO_COMPARE <- 6

# Custom column indices (used only when WHAT_TO_COMPARE == 9)
CUSTOM_COLUMNS  <- c()

# Annotation colours
ANNOT_COLOURS <- list(
  Batch       = c("1" = "green", "2" = "red",  "3" = "blue", "4" = "yellow"),
  Feeder      = c("DE" = "pink", "ELB" = "brown"),
  Replication = c("WithGeminin" = "grey",   "NoGeminin" = "black"),
  CellType    = c("Sperm" = "lightgreen", "Progenitor" = "cyan")
)


# =============================================================================
# COLUMN SELECTION MAP
# =============================================================================

COLUMN_MAP <- list(
  "1" = c(1:4,  21:24, 5:8,  25:28),
  "2" = c(1:4,  21:24, 5:8,  25:28, 9:12, 29:32),
  "3" = c(5:8,  25:28, 9:12, 29:32),
  "4" = c(21:28),
  "5" = c(21:32),
  "6" = c(1:8),
  "7" = c(1:12),
  "8" = c(1:4,  21:24)
)


# =============================================================================
# LOAD DATA
# =============================================================================

message("Loading data...")
expr   <- read_tsv(EXPR_FILE)
annot  <- read_tsv(ANNOT_FILE)
check_sample_names(expr, annot)

ordered_expr <- expr[order(rownames(expr)), rownames(annot)]

# Protein list of interest
prot_bg <- unique(read.table(PROT_LIST_FILE, stringsAsFactors = FALSE)$V1)
message("Proteins of interest: ", length(prot_bg))

# Comparison file for gene name lookup
comp_file <- read_tsv(COMPARISON_FILE)

# Select columns
cols <- if (WHAT_TO_COMPARE == 9) {
  CUSTOM_COLUMNS
} else {
  COLUMN_MAP[[as.character(WHAT_TO_COMPARE)]]
}
if (is.null(cols)) stop("Invalid WHAT_TO_COMPARE value: ", WHAT_TO_COMPARE)

data_sub <- ordered_expr[rownames(ordered_expr) %in% prot_bg, cols]
message("Matrix dimensions: ", nrow(data_sub), " proteins x ", ncol(data_sub), " samples")


# =============================================================================
# SCALE & CLUSTER
# =============================================================================

message("Scaling and clustering...")
mat         <- scale_rows(data_sub, center = TRUE, scale = TRUE)
q           <- quantile(mat, probs = c(0.01, 0.99), na.rm = TRUE)
col_scale   <- colorRamp2(c(q[1], 0, q[2]), c("blue", "white", "red"))

# Hierarchical clustering of proteins (pearson correlation distance)
dist_pearson <- as.dist(1 - cor(t(mat), use = "pairwise.complete.obs",
                                 method = "pearson"))
hclust_genes <- hclust(dist_pearson, method = "average")


# =============================================================================
# BUILD HEATMAP
# =============================================================================

message("Building heatmap...")
panel_name <- paste0("(", nrow(mat), " Proteins)")
ht <- Heatmap(
  mat,
  name            = panel_name,
  col             = col_scale,
  row_labels      = comp_file[rownames(mat), "Genes"],
  row_names_gp    = gpar(fontsize = max(2, min(8, 200 / nrow(mat)))),
  column_names_gp = gpar(fontsize = max(4, min(10, 200 / ncol(mat)))),
  cluster_rows    = hclust_genes,
  cluster_columns = FALSE,
  show_row_names  = TRUE,
  show_column_names = TRUE
)


# =============================================================================
# EXPORT
# =============================================================================

create_dir(FIGURES_DIR)
create_dir(RESULTS_DIR)

pdf_height <- pdf_height_from_nrow(nrow(mat))
pdf_path   <- file.path(FIGURES_DIR, paste0(panel_name, ".pdf"))
pdf(pdf_path, width = 10, height = pdf_height)
print(ht)
dev.off()
message("Heatmap saved: ", pdf_path)

write_tsv(as.data.frame(mat),
          file.path(RESULTS_DIR, paste0(panel_name, ".tsv")))

message("Done.")
