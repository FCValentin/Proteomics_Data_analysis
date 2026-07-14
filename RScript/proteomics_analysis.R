# =============================================================================
# proteomics_analysis.R
# -----------------------------------------------------------------------------
# Master analysis script for quantitative proteomics (Sperm vs Progenitor).
# Compiles all sub-analyses:
#   1. Differential expression heatmaps (filtered + all samples)
#   2. DE protein barplot summary
#   3. Custom candidate gene heatmaps (EPI/TF list)
#   4. Pathway enrichment analysis (pathfindR / GO-BP)
#   5. Transcriptomic / proteomic integration (volcano plots + FACS-like)
#
# Usage:
#   Rscript proteomics_analysis.R <data_dir>
#   Or set DATA_DIR below and source() in RStudio.
#
# Author  : Valentin FRANCOIS--CAMPION, PhD
# Contact : valentin.francoiscampion@gmail.com
# GitHub  : https://github.com/FCValentin/proteomics-quantitative-analysis
# Project : Quantitative proteomics — sperm vs progenitor cells
# Date    : 2023 (CR2TI, UMR 1064, Nantes Universite)
# =============================================================================


# =============================================================================
# I. DEPENDENCIES
# =============================================================================

source("proteomics_utils.R")

.load_pkg("circlize")
.load_pkg("ComplexHeatmap")
.load_pkg("Matrix")

# Optional — only needed for pathway enrichment (Section IV)
PATHFINDR_AVAILABLE <- requireNamespace("pathfindR", quietly = TRUE)
if (!PATHFINDR_AVAILABLE)
  message("Note: 'pathfindR' not installed — pathway enrichment will be skipped.")


# =============================================================================
# II. PARAMETERS — edit here
# =============================================================================

DATA_DIR      <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(DATA_DIR) || DATA_DIR == "") DATA_DIR <- "."

EXPR_FILE       <- file.path(DATA_DIR, "data", "dataMatrix_wash_groups_newSettings.tsv")
ANNOT_FILE      <- file.path(DATA_DIR, "data", "SampleAnnot.tsv")
ORTHOLOG_FILE   <- file.path(DATA_DIR, "data", "Orthologs.tsv")
COMPARISONS_DIR <- file.path(DATA_DIR, "Comparaison")
FIGURES_DIR     <- file.path(DATA_DIR, "Figures")
RESULTS_DIR     <- file.path(DATA_DIR, "results")
FIGS_DIR        <- file.path(DATA_DIR, "figs")

# Differential expression thresholds
PVAL       <- 0.05
FOLD_CHANGE <- 0.58     # absolute log2 FC threshold

# Pathway enrichment parameters
PEA_PVAL    <- 0.05
PEA_DATABASE <- "GO-BP"
PEA_START_IDX <- 28     # index of first comparison file to run PEA on

# Annotation colour mapping
ANNOT_COLOURS <- list(
  Batch       = c("1" = "green", "2" = "red",   "3" = "blue", "4" = "yellow"),
  Feader      = c("DE" = "pink", "ELB" = "brown"),
  Replication = c("WithGeminin" = "grey",  "NoGeminin" = "black"),
  CellType    = c("Sperm" = "lightgreen", "Progenitor" = "cyan")
)


# =============================================================================
# III. INTERNAL HELPER FUNCTIONS
# =============================================================================

#' Build a ComplexHeatmap column annotation bar
#' @param sample_ids Character vector. Sample identifiers (must be rownames of annot).
#' @param annot data.frame. Full sample annotation.
#' @param colours Named list. Colour mapping per annotation variable.
#' @return HeatmapAnnotation object.
build_col_annotation <- function(sample_ids, annot, colours = ANNOT_COLOURS) {
  HeatmapAnnotation(
    Batch       = annot[sample_ids, "Batch"],
    Feader      = annot[sample_ids, "Feader"],
    Replication = annot[sample_ids, "Geminin"],
    CellType    = annot[sample_ids, "Tissue"],
    annotation_legend_param = list(
      Batch       = list(title = "Batch",
                         at = c(1, 2, 3, 4), labels = c("1", "2", "3", "4"),
                         color_bar = "discrete"),
      Feader      = list(title = "Feeder",
                         at = c("DE", "ELB"), labels = c("DE", "ELB"),
                         color_bar = "discrete"),
      Replication = list(title = "Replication",
                         at = c("WithGeminin", "NoGeminin"),
                         labels = c("+Gem", "-Gem"),
                         color_bar = "discrete"),
      CellType    = list(title = "Cell Type",
                         at = c("Sperm", "Progenitor"),
                         labels = c("Sperm", "Progenitor"),
                         color_bar = "discrete")
    ),
    col = colours
  )
}


#' Build a scaled + clustered Heatmap object
#'
#' @param data_mat     Matrix (proteins x samples). Expression data.
#' @param annot        data.frame. Full sample annotation.
#' @param file_annot   data.frame. Comparison file with 'Genes' column.
#' @param panel_name   Character. Heatmap legend name.
#' @param cluster_cols Logical or hclust. Column clustering specification.
#' @return Heatmap object.
build_heatmap <- function(data_mat, annot, file_annot, panel_name,
                           cluster_cols = TRUE) {
  mat        <- scale_rows(data_mat)
  q          <- quantile(unlist(mat), probs = c(0.01, 0.99), na.rm = TRUE)
  col_scale  <- colorRamp2(c(q[1], 0, q[2]), c("blue", "white", "red"))

  # Protein clustering (pearson)
  dist_p     <- as.dist(1 - cor(t(mat), use = "pairwise.complete.obs",
                                 method = "pearson"))
  gene_clust <- hclust(dist_p, method = "average")

  # Sample clustering (euclidean) — only when cluster_cols is TRUE
  if (isTRUE(cluster_cols)) {
    dist_e      <- dist(t(mat), method = "euclidean")
    sample_clust <- hclust(dist_e, method = "complete")
    sample_order <- sample_clust$labels[sample_clust$order]
    col_clust    <- sample_clust
  } else {
    sample_order <- colnames(mat)
    col_clust    <- FALSE
  }

  col_annot <- build_col_annotation(sample_order, annot)

  row_labels <- if (!is.null(file_annot) && "Genes" %in% colnames(file_annot)) {
    file_annot[rownames(mat), "Genes"]
  } else {
    rownames(mat)
  }

  font_r <- max(2, min(8, 200 / nrow(mat)))
  font_c <- max(4, min(10, 200 / ncol(mat)))

  Heatmap(
    mat,
    name               = panel_name,
    col                = col_scale,
    row_labels         = row_labels,
    row_names_gp       = gpar(fontsize = font_r),
    column_names_gp    = gpar(fontsize = font_c),
    top_annotation     = col_annot,
    cluster_rows       = gene_clust,
    cluster_columns    = col_clust,
    show_row_names     = TRUE,
    show_column_names  = TRUE
  )
}


#' Export a Heatmap to PDF with adaptive height
#' @param ht     Heatmap or HeatmapList.
#' @param path   Character. PDF output path.
#' @param n_rows Integer. Number of protein rows (used for height scaling).
export_heatmap_pdf <- function(ht, path, n_rows) {
  h <- pdf_height_from_nrow(n_rows)
  pdf(path, width = 10, height = h)
  print(ht)
  dev.off()
  message("Heatmap saved: ", path)
}


#' Plot labelled scatter (volcano / FACS-like)
#' @param x, y     Numeric vectors. X and Y coordinates.
#' @param labels   Character vector. Point labels.
#' @param highlight Character vector. Point IDs to highlight in red.
#' @param xlab, ylab, main Character. Axis / title labels.
#' @param xlim, ylim Numeric. Axis limits (NULL = auto).
plot_labelled_scatter <- function(x, y, labels, highlight = NULL,
                                   xlab = "", ylab = "", main = "",
                                   xlim = NULL, ylim = NULL) {
  col_vec <- ifelse(labels %in% highlight, "red", "black")
  plot(x, y, type = "n", xlab = xlab, ylab = ylab, main = main,
       xlim = xlim, ylim = ylim)
  text(x, y, labels = labels, cex = 0.2, col = col_vec)
}


# =============================================================================
# IV. LOAD DATA
# =============================================================================

message("=== Quantitative Proteomics Analysis ===")
message("[1/6] Loading data...")

expr   <- read_tsv(EXPR_FILE)
annot  <- read_tsv(ANNOT_FILE)
check_sample_names(expr, annot)
ordered_expr <- expr[order(rownames(expr)), rownames(annot)]

# Comparison files
filenames <- list.files(COMPARISONS_DIR, pattern = "\\.tsv$",
                        full.names = TRUE, recursive = FALSE)
if (length(filenames) == 0) stop("No comparison .tsv files found in: ", COMPARISONS_DIR)
message("Comparison files: ", length(filenames))

# Sample subsets per comparison
sample_to_keep  <- read_tsv(file.path(DATA_DIR, "data", "SamplePerCond.tsv"))

create_dir(FIGURES_DIR)
create_dir(RESULTS_DIR)
create_dir(FIGS_DIR)


# =============================================================================
# V. SECTION 1 — DE HEATMAPS PER COMPARISON
# =============================================================================

message("[2/6] Building DE heatmaps for all comparisons...")

ht_all      <- list()   # all samples
ht_filtered <- list()   # filtered (SamplePerCond subset)
n_de_prot   <- c()

for (i in seq_along(filenames)) {
  message("  [", i, "/", length(filenames), "] ", basename(filenames[i]))
  comp        <- read_tsv(filenames[i])
  comp_label  <- parse_comparison_label(filenames[i])
  de_proteins <- rownames(comp)[comp$Qvalue < PVAL &
                                  comp$AbsoluteAVGLog2Ratio > FOLD_CHANGE]

  # All samples
  data_all <- ordered_expr[de_proteins, , drop = FALSE]
  if (nrow(data_all) > 1) {
    panel_all     <- paste0(comp_label, " (", nrow(data_all), " Proteins)")
    ht_all[[i]]  <- build_heatmap(data_all, annot, comp, panel_all,
                                   cluster_cols = TRUE)
    write_tsv(as.data.frame(scale_rows(data_all)),
              file.path(RESULTS_DIR,
                        paste0("MatrixHeatmap_", comp_label,
                               "_Pval_", PVAL, "_FC_", FOLD_CHANGE, ".tsv")))
  }

  # Filtered subset (SamplePerCond)
  keep_cols <- unlist(sample_to_keep[filenames[i], ])
  keep_cols <- keep_cols[keep_cols %in% colnames(ordered_expr)]
  data_filt <- ordered_expr[de_proteins, keep_cols, drop = FALSE]
  if (nrow(data_filt) > 1 && ncol(data_filt) > 0) {
    panel_filt        <- paste0(comp_label, "_filtered (", nrow(data_filt), " Proteins)")
    ht_filtered[[i]] <- build_heatmap(data_filt, annot, comp, panel_filt,
                                       cluster_cols = TRUE)
  }

  n_de_prot <- c(n_de_prot, nrow(data_all))
}
names(n_de_prot) <- sapply(filenames, parse_comparison_label)

# DE summary barplot
pdf(file.path(FIGURES_DIR,
              paste0("Distribution_Barplot_NbrProtDE_Pval_", PVAL,
                     "_FC_", FOLD_CHANGE, ".pdf")),
    width = 10, height = 10)
par(mar = c(15, 4, 4, 2))
barplot(n_de_prot, main = "# DE Proteins per Comparison",
        cex.names = 0.7, las = 2, ylim = c(0, max(n_de_prot) * 1.1),
        names.arg = n_de_prot)
abline(h = nrow(ordered_expr), lty = 2, col = "grey50")
dev.off()
write_tsv(as.data.frame(n_de_prot), file.path(RESULTS_DIR, "BarplotValues.tsv"))

# Export heatmaps
valid_all  <- Filter(Negate(is.null), ht_all)
valid_filt <- Filter(Negate(is.null), ht_filtered)

export_heatmap_pdf(
  valid_all,
  file.path(FIGURES_DIR,
            paste0("Heatmap_AllSamples_Pval_", PVAL, "_FC_", FOLD_CHANGE, ".pdf")),
  n_rows = max(n_de_prot)
)
export_heatmap_pdf(
  valid_filt,
  file.path(FIGURES_DIR,
            paste0("Heatmap_Filtered_Pval_", PVAL, "_FC_", FOLD_CHANGE, ".pdf")),
  n_rows = max(n_de_prot)
)


# =============================================================================
# VI. SECTION 2 — CUSTOM CANDIDATE GENE HEATMAPS
# =============================================================================

message("[3/6] Custom candidate gene heatmaps...")

PVAL2        <- 0.001
FOLD_CHANGE2 <- 1.0
sample_to_keep2 <- read_tsv(file.path(DATA_DIR, "data", "SamplePerCondHeatmap3.tsv"))

de_prot_union <- c()
for (i in seq_len(nrow(sample_to_keep2))) {
  f            <- read_tsv(rownames(sample_to_keep2)[i])
  de_prot_union <- c(de_prot_union,
                     rownames(f)[f$Qvalue < PVAL2 & f$AbsoluteAVGLog2Ratio > FOLD_CHANGE2])
}
prot_bg <- unique(de_prot_union)

data_cand   <- ordered_expr[rownames(ordered_expr) %in% prot_bg,
                              colnames(ordered_expr) %in% unlist(unique(sample_to_keep2)),
                              drop = FALSE]
# Column orderings (original exploration — uncomment desired layout)
# Layout 1 : DataOrder <- data_cand[, c(13,14,15,16,9,10,12,11,1,4,3,2,5,7,6,8)]
# Layout 2 : DataOrder <- data_cand[, c(13,14,15,16,9,10,12,11,1,2,5,6,3,4,7,8)]
# Layout 3 : DataOrder <- data_cand[, c(1,7,5,4,9,13,11,15,2,8,6,3,10,14,12,16)]
# Layout 4 (default):
DataOrder <- data_cand[, c(2, 3, 9, 10, 5, 7, 13, 15, 1, 4, 11, 12, 6, 8, 14, 16)]

file_ref <- read_tsv(filenames[min(42, length(filenames))])
ht_cand  <- build_heatmap(DataOrder, annot, file_ref,
                           paste0("(", nrow(DataOrder), " Proteins) — ELB Replication"),
                           cluster_cols = FALSE)

export_heatmap_pdf(
  ht_cand,
  file.path(FIGURES_DIR, paste0("HeatmapReplicationELB_Pval_",
                                 PVAL2, "_FC_", FOLD_CHANGE2, ".pdf")),
  n_rows = nrow(DataOrder)
)
write_tsv(as.data.frame(scale_rows(DataOrder)),
          file.path(RESULTS_DIR, paste0("MatrixHeatmapReplicationELB_Pval_",
                                         PVAL2, "_FC_", FOLD_CHANGE2, ".tsv")))

# EPI / TF candidate list heatmap
epi_list  <- unique(read.table(
  file.path(DATA_DIR,
            "ListOfDemandes/list_orthologs_EPI_and_TF_SEGM_DE/list_orthologs_EPI_and_TF_SEGM_DE.txt")
)$V1)
data_epi  <- ordered_expr[rownames(ordered_expr) %in% epi_list,
                            c(37,38,39,40,33,34,36,35,1,11,5,4,13,19,15,21,
                              2,12,6,3,14,20,16,22), drop = FALSE]
ht_epi    <- build_heatmap(data_epi, annot, file_ref,
                            paste0("(", nrow(data_epi), " Proteins) — EPI/TF"),
                            cluster_cols = FALSE)

export_heatmap_pdf(
  ht_epi,
  file.path(DATA_DIR,
            "ListOfDemandes/Heatmap-list_orthologs_EPI_and_TF_SEGM_DE.pdf"),
  n_rows = nrow(data_epi)
)
write_tsv(as.data.frame(scale_rows(data_epi)),
          file.path(DATA_DIR,
                    "ListOfDemandes/Matrix-list_orthologs_EPI_and_TF_SEGM_DE.tsv"))


# =============================================================================
# VII. SECTION 3 — PATHWAY ENRICHMENT ANALYSIS (pathfindR)
# =============================================================================

if (PATHFINDR_AVAILABLE) {
  message("[4/6] Running pathway enrichment analysis (", PEA_DATABASE, ")...")
  library(pathfindR)
  orthologs <- read_tsv(ORTHOLOG_FILE)

  for (i in PEA_START_IDX:length(filenames)) {
    message("  PEA [", i, "/", length(filenames), "] ", basename(filenames[i]))
    comp        <- read_tsv(filenames[i])
    comp_label  <- parse_comparison_label(filenames[i])
    input_df    <- comp[rownames(comp) %in% rownames(orthologs),
                        c("Genes", "AVGLog2Ratio", "Qvalue")]
    input_df    <- input_df[rownames(orthologs)[rownames(orthologs) %in% rownames(input_df)], ]
    input_df$Genes <- orthologs$HumanName

    run_pea <- function(df, label) {
      tryCatch(run_pathfindR(df, p_val_threshold = PEA_PVAL, gene_sets = PEA_DATABASE),
               error = function(e) { message("  PEA failed: ", e$message); NULL })
    }

    out_all  <- run_pea(input_df, "all")
    out_up   <- run_pea(input_df[input_df$AVGLog2Ratio > 0, ], "up")
    out_down <- run_pea(input_df[input_df$AVGLog2Ratio < 0, ], "down")

    pdf_path <- file.path(RESULTS_DIR,
                           paste0("GoEnrichment_", comp_label, "_Pval_",
                                  PEA_PVAL, "_", PEA_DATABASE, ".pdf"))
    pdf(pdf_path, width = 10, height = 10)
    for (res in list(out_all, out_up, out_down)) {
      if (!is.null(res) && nrow(res) > 0) {
        print(enrichment_chart(result_df = res, top_terms = 30))
        print(term_gene_heatmap(result_df = res, genes_df = input_df))
        print(term_gene_graph(result_df = res, use_description = TRUE))
        print(UpSet_plot(result_df = res, genes_df = input_df))
      }
    }
    dev.off()
    message("  Saved: ", pdf_path)
  }
} else {
  message("[4/6] Skipping pathway enrichment — pathfindR not installed.")
}


# =============================================================================
# VIII. SECTION 4 — TRANSCRIPTOMIC / PROTEOMIC INTEGRATION
# =============================================================================

message("[5/6] Transcriptomic / proteomic integration...")

# Load RNA-seq values (expected in RNA-Seq/results/MatrixVolcanoPlotPerso.tsv)
rnaseq_dir       <- file.path(DATA_DIR, "RNA-Seq")
fc_matrix_file   <- file.path(rnaseq_dir, "results", "MatrixVolcanoPlotPerso.tsv")
logfc_threshold  <- 0.58

if (!file.exists(fc_matrix_file)) {
  message("  FC matrix not found at: ", fc_matrix_file, " — skipping integration.")
} else {
  fc_matrix  <- read_tsv(fc_matrix_file)

  # Protein mean expression per condition
  file_main      <- read_tsv(filenames[min(37, length(filenames))])
  rownames(ordered_expr) <- file_main[rownames(ordered_expr), "Genes"]

  pr_prot <- rowMeans(ordered_expr[, c(33, 34, 36, 35), drop = FALSE], na.rm = TRUE)
  sp_prot <- rowMeans(ordered_expr[, c(37, 38, 39, 40), drop = FALSE], na.rm = TRUE)

  # RNA-seq fold-change
  rna_expr <- tryCatch(
    read_tsv(file.path(rnaseq_dir, "dataMatrix.tsv")),
    error = function(e) NULL
  )

  prot_de_file <- file.path(DATA_DIR,
    "ListOfDemandes/list_ortholog_up_down_sperm_vs_prog/list_ortholog_up_down_sperm_vs_prog.txt")
  prot_de <- if (file.exists(prot_de_file)) {
    rownames(read.table(prot_de_file, stringsAsFactors = FALSE))
  } else { c() }

  common_prots <- intersect(rownames(fc_matrix),
                             intersect(names(pr_prot), names(sp_prot)))
  fc_sub       <- fc_matrix[common_prots, , drop = FALSE]
  fc_sub       <- cbind(fc_sub,
                        PrProtValues = pr_prot[common_prots],
                        SpProtValues = sp_prot[common_prots])
  write_tsv(fc_sub, file.path(RESULTS_DIR, "MatrixVolcanoPlotPerso.tsv"))

  # ── Volcano plots ──────────────────────────────────────────────────────────
  pdf(file.path(FIGS_DIR, "VolcanoPlots.pdf"), width = 10, height = 10)

  plot_labelled_scatter(
    x = fc_sub[, 2], y = fc_sub[, 1],
    labels = file_main[rownames(fc_sub), "Genes"],
    highlight = prot_de,
    xlab = "log2FC(Sp/Pr) Protein", ylab = "Progenitor RNA-Seq",
    main = "Volcano plot — FC vs RNA-Seq"
  )
  abline(v = c(-logfc_threshold, logfc_threshold), lty = 2, col = "grey50")

  plot_labelled_scatter(
    x = file_main[rownames(fc_sub), "AVGLog2Ratio"],
    y = fc_sub[, 1],
    labels = file_main[rownames(fc_sub), "Genes"],
    highlight = prot_de,
    xlab = "AVG Log2 Ratio (Protein)", ylab = "Progenitor RNA-Seq",
    main = "Volcano plot — AVG Log2 Ratio vs RNA-Seq"
  )
  abline(v = c(-logfc_threshold, logfc_threshold), lty = 2, col = "grey50")

  plot_labelled_scatter(
    x = fc_sub$PrProtValues, y = fc_sub[, 1],
    labels = file_main[rownames(fc_sub), "Genes"],
    highlight = prot_de,
    xlab = "Prog Protein", ylab = "Progenitor RNA-Seq",
    main = "Progenitor Protein vs RNA-Seq"
  )
  dev.off()
  message("  Volcano plots saved.")

  # ── FACS-like quadrant plots ───────────────────────────────────────────────
  .facs_quadrants <- function(fc_sub, x_col, xlab, file_pdf, xlim_val) {
    pdf(file_pdf, width = 10, height = 10)
    quads <- list(
      list(
        subset = fc_sub[fc_sub$RNAValues < 10 & fc_sub$PrProtValues < 7, ],
        main   = "Double Negative (low RNA + low Pr Prot)"
      ),
      list(
        subset = fc_sub[fc_sub$RNAValues > 14 & fc_sub$PrProtValues > 12, ],
        main   = "Double Positive (high RNA + high Pr Prot)"
      ),
      list(
        subset = fc_sub[fc_sub$RNAValues > 14 & fc_sub$PrProtValues < 7, ],
        main   = "High RNA / Low Pr Prot"
      ),
      list(
        subset = fc_sub[fc_sub$RNAValues < 10 & fc_sub$PrProtValues > 12, ],
        main   = "Low RNA / High Pr Prot"
      )
    )
    for (q in quads) {
      plot_labelled_scatter(
        x         = q$subset[, x_col],
        y         = q$subset[, 1],
        labels    = file_main[rownames(q$subset), "Genes"],
        highlight = prot_de,
        xlab      = xlab,
        ylab      = "Progenitor RNA-Seq",
        main      = q$main,
        xlim      = xlim_val,
        ylim      = c(0, 20)
      )
    }
    dev.off()
    message("  FACS-like saved: ", file_pdf)
  }

  .facs_quadrants(fc_sub, "SpProtValues",
                  "Sperm Protein (mean expr)",
                  file.path(FIGS_DIR, "FACSlike_SpProt.pdf"),
                  xlim_val = c(0, 20))

  .facs_quadrants(fc_sub, 2,
                  "log2FC(Sp/Pr) Protein",
                  file.path(FIGS_DIR, "FACSlike_FC.pdf"),
                  xlim_val = c(-2, 2))
}


# =============================================================================
# IX. DONE
# =============================================================================

message("[6/6] Analysis complete.")
message("  Figures  : ", FIGURES_DIR)
message("  Results  : ", RESULTS_DIR)
message("  FACS/Volcano: ", FIGS_DIR)
