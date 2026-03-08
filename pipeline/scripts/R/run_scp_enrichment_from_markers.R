#!/usr/bin/env Rscript
################################################################################
# run_scp_enrichment_from_markers.R  (FULL, corrected)
#
# Inputs:
#   --pivot      scp_pivot.tsv (must contain Genes column; background = all genes detected anywhere)
#   --annotation scp_annotation.tsv (Run=stem, Condition)
#   --detection  detection_markers.tsv (Cluster, Genes, det_rate_in/out, qvalue, etc.)
#   --intensity  intensity_markers_detected_only.tsv (Cluster, Genes, log2FC_detected_only, qvalue, etc.)
#   --aucell     aucell_scores.tsv (Pathway + sample columns)
#   --outdir     output folder
#
# Outputs:
#   tables/aucell/aucell_pathway_markers.tsv
#   plots/aucell/aucell_pathways_<cluster>.pdf
#   tables/go/<cluster>_(DET_UP|INT_UP|INT_DOWN)_(ALL|FDR|P|SKIPPED|NONE).*
#   plots/go/go_<cluster>_(DET_UP|INT_UP|INT_DOWN).pdf
#
# Notes:
# - Intended for SCP where missingness is high.
# - Background is dataset-specific: all genes observed anywhere in pivot$Genes.
# - GO uses universe=background Entrez IDs (not "genome").
# - AUCell pathway markers are Wilcoxon cluster vs rest.
# - Robust option handling: uses opt[["flag-name"]] for hyphen flags.
# - Robust plotting: safe even if df empty or max_terms invalid.
################################################################################

suppressPackageStartupMessages({
  if (!require("optparse", quietly = TRUE)) stop("Package 'optparse' required.")
  library(optparse)
})

option_list <- list(
  make_option("--pivot", type = "character", help = "scp_pivot.tsv (Genes column present)"),
  make_option("--annotation", type = "character", help = "scp_annotation.tsv (Run, Condition)"),
  make_option("--detection", type = "character", help = "detection_markers.tsv"),
  make_option("--intensity", type = "character", help = "intensity_markers_detected_only.tsv"),
  make_option("--aucell", type = "character", help = "aucell_scores.tsv (Pathway x samples)"),
  make_option("--outdir", type = "character", help = "Output directory"),

  make_option("--species", type = "character", default = "auto",
              help = "auto|human|mouse [default %default]"),

  make_option("--go-ont", type = "character", default = "BP",
              help = "GO ontology [default %default]"),
  make_option("--fdr", type = "double", default = 0.05,
              help = "FDR threshold for reporting [default %default]"),
  make_option("--pval", type = "double", default = 0.05,
              help = "Nominal p threshold for extra tables [default %default]"),

  make_option("--det-delta", type = "double", default = 0.10,
              help = "Min det_rate_in - det_rate_out [default %default]"),
  make_option("--det-q", type = "double", default = 0.10,
              help = "Max qvalue for detection markers [default %default]"),

  make_option("--int-q", type = "double", default = 0.05,
              help = "Max qvalue for intensity markers [default %default]"),
  make_option("--int-fc", type = "double", default = 0.25,
              help = "Min abs(log2FC_detected_only) [default %default]"),

  make_option("--max-terms", type = "integer", default = 20,
              help = "Max terms/pathways to plot [default %default]"),
  make_option("--seed", type = "integer", default = 42,
              help = "Random seed [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

# ---- resolve hyphenated options safely ----
FDR_THR   <- as.numeric(opt[["fdr"]])
P_THR     <- as.numeric(opt[["pval"]])
DET_DELTA <- as.numeric(opt[["det-delta"]])
DET_Q     <- as.numeric(opt[["det-q"]])
INT_Q     <- as.numeric(opt[["int-q"]])
INT_FC    <- as.numeric(opt[["int-fc"]])
MAX_TERMS <- as.integer(opt[["max-terms"]])
if (!is.finite(FDR_THR)) FDR_THR <- 0.05
if (!is.finite(P_THR)) P_THR <- 0.05
if (!is.finite(DET_DELTA)) DET_DELTA <- 0.10
if (!is.finite(DET_Q)) DET_Q <- 0.10
if (!is.finite(INT_Q)) INT_Q <- 0.05
if (!is.finite(INT_FC)) INT_FC <- 0.25
if (is.na(MAX_TERMS) || MAX_TERMS < 1) MAX_TERMS <- 20L

# ---- soft-fail package checks (no auto-install) ----
has_pkg <- function(x) suppressWarnings(requireNamespace(x, quietly = TRUE))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "tables", "go"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "tables", "aucell"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "plots", "go"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "plots", "aucell"), recursive = TRUE, showWarnings = FALSE)

soft_fail <- function(msg) {
  writeLines(msg, con = file.path(opt$outdir, "SKIPPED_missing_packages.txt"))
  cat(msg, "\n")
  quit(status = 0)
}

req <- c("data.table", "ggplot2")
miss <- req[!vapply(req, has_pkg, logical(1))]
if (length(miss) > 0) soft_fail(paste0("SCP enrichment skipped: missing packages: ", paste(miss, collapse = ", ")))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

can_go <- has_pkg("clusterProfiler") && has_pkg("AnnotationDbi") &&
  ((opt$species %in% c("auto", "human") && has_pkg("org.Hs.eg.db")) ||
     (opt$species == "mouse" && has_pkg("org.Mm.eg.db")))

OrgDb <- NULL
if (can_go) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(AnnotationDbi)
  })
} else {
  writeLines(
    "GO skipped (missing clusterProfiler/AnnotationDbi/org.* packages).",
    con = file.path(opt$outdir, "tables", "go", "SKIPPED_missing_packages.txt")
  )
}

# ---- helpers ----
truncate_terms <- function(x, max_chars = 35) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  too <- nchar(x) > max_chars
  x[too] <- paste0(substr(x[too], 1, max_chars - 3), "...")
  x
}

split_genes <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  toks <- unlist(strsplit(x, ";", fixed = TRUE))
  toks <- trimws(toks)
  toks <- toks[toks != ""]
  unique(toks)
}

detect_species <- function(genes) {
  g <- head(unique(genes), 200)
  g <- g[!is.na(g) & g != ""]
  if (length(g) == 0) return("human")
  n_upper <- sum(grepl("^[A-Z][A-Z0-9-]+$", g))
  if (n_upper / length(g) > 0.6) "human" else "mouse"
}

sym2entrez <- function(sym, OrgDb) {
  sym <- unique(sym[!is.na(sym) & sym != ""])
  if (length(sym) == 0) return(character(0))
  m <- AnnotationDbi::select(OrgDb, keys = sym, columns = c("ENTREZID"), keytype = "SYMBOL")
  unique(as.character(m$ENTREZID[!is.na(m$ENTREZID)]))
}

stem_from_header <- function(h) {
  h <- as.character(h)
  h2 <- basename(gsub("\\\\", "/", h))
  h2 <- sub("\\.(raw|d|mzml|wiff)$", "", h2, ignore.case = TRUE)
  h2
}

plot_dot <- function(df, out_pdf, xcol, ycol, sizecol, colorcol, xlab, colorlab, max_terms = 20) {
  # Always create a note file if nothing to plot.
  if (is.null(df) || nrow(df) == 0) {
    writeLines("No results.", con = sub("\\.pdf$", ".txt", out_pdf))
    return(invisible(NULL))
  }

  # Guard max_terms
  if (is.null(max_terms) || !is.finite(max_terms) || max_terms < 1) {
    max_terms <- min(20L, nrow(df))
  } else {
    max_terms <- as.integer(max_terms)
  }

  # Ensure needed columns exist
  for (c in c(xcol, ycol, sizecol, colorcol)) {
    if (!(c %in% colnames(df))) stop(sprintf("plot_dot missing required column '%s'", c))
  }

  # Coerce numeric columns for plotting
  df[[xcol]] <- suppressWarnings(as.numeric(df[[xcol]]))
  df[[sizecol]] <- suppressWarnings(as.numeric(df[[sizecol]]))
  df[[colorcol]] <- suppressWarnings(as.numeric(df[[colorcol]]))

  df <- df[is.finite(df[[xcol]]) & is.finite(df[[sizecol]]) & is.finite(df[[colorcol]]), , drop = FALSE]
  if (nrow(df) == 0) {
    writeLines("No finite results to plot.", con = sub("\\.pdf$", ".txt", out_pdf))
    return(invisible(NULL))
  }

  # Select top terms by colorcol (typically -log10(FDR)) then limit
  df <- df[order(df[[colorcol]], decreasing = TRUE), , drop = FALSE]
  df <- head(df, max_terms)

  df[[ycol]] <- truncate_terms(df[[ycol]], 35)
  # ensure unique labels after truncation (prevents duplicated factor levels)
  df[[ycol]] <- make.unique(as.character(df[[ycol]]))
  df[[ycol]] <- factor(df[[ycol]], levels = rev(unique(df[[ycol]])))


  n <- nrow(df)
  w <- 6.2
  h <- max(4.8, 0.32 * n + 2.2)

  p <- ggplot(df, aes(x = .data[[xcol]], y = .data[[ycol]])) +
    geom_point(aes(size = .data[[sizecol]], color = .data[[colorcol]]), alpha = 0.9) +
    scale_color_gradient(low = "#2166AC", high = "#B2182B", name = colorlab) +
    scale_size_continuous(range = c(2.2, 7.5), name = "Size") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9),
      plot.margin = grid::unit(c(5.5, 18, 5.5, 5.5), "pt"),
      legend.position = "right"
    ) +
    labs(title = "", x = xlab, y = "")

  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(out_pdf, width = w, height = h, useDingbats = FALSE)
  print(p)
  dev.off()
  invisible(TRUE)
}

# ---- load inputs ----
pivot <- fread(opt$pivot, sep = "\t", header = TRUE, data.table = FALSE)
if (!("Genes" %in% colnames(pivot))) stop("Pivot must have Genes column.")
pivot$Genes <- as.character(pivot$Genes)

ann <- fread(opt$annotation, sep = "\t", header = TRUE, data.table = FALSE)
if (!all(c("Run", "Condition") %in% colnames(ann))) stop("Annotation must have Run and Condition.")
ann$Run <- as.character(ann$Run)
ann$Condition <- as.character(ann$Condition)

det <- fread(opt$detection, sep = "\t", header = TRUE, data.table = FALSE)
int <- fread(opt$intensity, sep = "\t", header = TRUE, data.table = FALSE)

auc <- fread(opt$aucell, sep = "\t", header = TRUE, data.table = FALSE)
if (!("Pathway" %in% colnames(auc))) stop("AUCell table must have Pathway column.")

# ---- background genes: all genes detected anywhere ----
bg_syms <- unique(unlist(lapply(pivot$Genes, split_genes)))
bg_syms <- bg_syms[!is.na(bg_syms) & bg_syms != ""]
cat(sprintf("Background genes (detected anywhere): %d\n", length(bg_syms)))

species <- opt$species
if (species == "auto") {
  species <- detect_species(bg_syms)
  cat(sprintf("Auto-detected species: %s\n", species))
}

if (can_go) {
  if (species == "human") {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    OrgDb <- org.Hs.eg.db
  } else if (species == "mouse") {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    OrgDb <- org.Mm.eg.db
  } else {
    can_go <- FALSE
    OrgDb <- NULL
  }
}

bg_entrez <- character(0)
if (can_go && !is.null(OrgDb)) {
  bg_entrez <- sym2entrez(bg_syms, OrgDb)
  cat(sprintf("Background Entrez IDs: %d\n", length(bg_entrez)))
}

# ---- clusters ----
clusters <- sort(unique(ann$Condition))

# ---- coerce numeric columns in marker tables ----
if ("det_rate_in" %in% colnames(det)) det$det_rate_in <- as.numeric(det$det_rate_in)
if ("det_rate_out" %in% colnames(det)) det$det_rate_out <- as.numeric(det$det_rate_out)
if ("qvalue" %in% colnames(det)) det$qvalue <- as.numeric(det$qvalue)
if ("pvalue" %in% colnames(det)) det$pvalue <- as.numeric(det$pvalue)

if ("qvalue" %in% colnames(int)) int$qvalue <- as.numeric(int$qvalue)
if ("pvalue" %in% colnames(int)) int$pvalue <- as.numeric(int$pvalue)
if ("log2FC_detected_only" %in% colnames(int)) int$log2FC_detected_only <- as.numeric(int$log2FC_detected_only)

# =============================================================================
# AUCell pathway markers (cluster vs rest, Wilcoxon)
# =============================================================================
auc_cols <- setdiff(colnames(auc), "Pathway")
auc_long <- melt(as.data.table(auc), id.vars = "Pathway",
                 variable.name = "SampleColumn", value.name = "AUC")
auc_long[, Run := vapply(SampleColumn, stem_from_header, character(1))]
auc_long <- merge(auc_long, as.data.table(ann), by = "Run", all.x = TRUE)

path_marks <- list()
for (cl in clusters) {
  sub_in <- auc_long[Condition == cl]
  sub_out <- auc_long[Condition != cl]
  if (nrow(sub_in) < 3 || nrow(sub_out) < 3) next

  pw <- auc_long[, .(
    mean_in = mean(AUC[Condition == cl], na.rm = TRUE),
    mean_out = mean(AUC[Condition != cl], na.rm = TRUE),
    delta = mean(AUC[Condition == cl], na.rm = TRUE) - mean(AUC[Condition != cl], na.rm = TRUE),
    pvalue = tryCatch(wilcox.test(AUC[Condition == cl], AUC[Condition != cl])$p.value, error = function(e) NA_real_)
  ), by = .(Pathway)]

  pw[, qvalue := p.adjust(pvalue, method = "BH")]
  pw[, Cluster := cl]
  path_marks[[cl]] <- pw
}
path_marks_df <- rbindlist(path_marks, fill = TRUE)
fwrite(path_marks_df, file.path(opt$outdir, "tables", "aucell", "aucell_pathway_markers.tsv"), sep = "\t")

# Plot AUCell pathway delta dotplot per cluster
for (cl in clusters) {
  sub <- path_marks_df[Cluster == cl & is.finite(qvalue), ]
  if (nrow(sub) == 0) next
  sub <- sub[order(qvalue), ]
  sub$logFDR <- -log10(sub$qvalue)
  sub$PathwayClean <- truncate_terms(sub$Pathway, 35)

  plot_dot(
    df = as.data.frame(sub),
    out_pdf = file.path(opt$outdir, "plots", "aucell", paste0("aucell_pathways_", cl, ".pdf")),
    xcol = "delta",
    ycol = "PathwayClean",
    sizecol = "logFDR",
    colorcol = "logFDR",
    xlab = "Mean AUC delta (in - out)",
    colorlab = "-log10(FDR)",
    max_terms = MAX_TERMS
  )
}

# =============================================================================
# GO enrichment (ORA) per cluster and marker type
# =============================================================================
if (can_go && !is.null(OrgDb) && length(bg_entrez) > 0) {
  for (cl in clusters) {

    # ---- detection UP ----
    det$det_delta <- det$det_rate_in - det$det_rate_out
    det_cl <- det[det$Cluster == cl &
                    is.finite(det$qvalue) & det$qvalue <= DET_Q &
                    is.finite(det$det_delta) & det$det_delta >= DET_DELTA, , drop = FALSE]
    det_syms <- unique(unlist(lapply(det_cl$Genes, split_genes)))
    det_syms <- det_syms[det_syms %in% bg_syms]
    det_ent <- sym2entrez(det_syms, OrgDb)

    # ---- intensity UP / DOWN ----
    int_cl <- int[int$Cluster == cl & is.finite(int$qvalue) & int$qvalue <= INT_Q, , drop = FALSE]
    up_syms <- unique(unlist(lapply(int_cl$Genes[is.finite(int_cl$log2FC_detected_only) & int_cl$log2FC_detected_only >= INT_FC], split_genes)))
    dn_syms <- unique(unlist(lapply(int_cl$Genes[is.finite(int_cl$log2FC_detected_only) & int_cl$log2FC_detected_only <= -INT_FC], split_genes)))
    up_syms <- up_syms[up_syms %in% bg_syms]
    dn_syms <- dn_syms[dn_syms %in% bg_syms]
    up_ent <- sym2entrez(up_syms, OrgDb)
    dn_ent <- sym2entrez(dn_syms, OrgDb)

    do_go <- function(ent, tag) {
      if (length(ent) < 5) {
        writeLines("Not enough genes.", con = file.path(opt$outdir, "tables", "go", paste0(cl, "_", tag, "_SKIPPED.txt")))
        return(invisible(NULL))
      }

      enr <- tryCatch(
        clusterProfiler::enrichGO(
          gene = ent,
          universe = bg_entrez,
          OrgDb = OrgDb,
          ont = opt[["go-ont"]],
          keyType = "ENTREZID",
          readable = TRUE
        ),
        error = function(e) NULL
      )

      if (is.null(enr) || nrow(as.data.frame(enr)) == 0) {
        writeLines("No enriched terms.", con = file.path(opt$outdir, "tables", "go", paste0(cl, "_", tag, "_NONE.txt")))
        return(invisible(NULL))
      }

      df <- as.data.frame(enr)
      fwrite(as.data.table(df), file.path(opt$outdir, "tables", "go", paste0(cl, "_", tag, "_ALL.tsv")), sep = "\t")

      df_fdr <- df[is.finite(df$p.adjust) & df$p.adjust < FDR_THR, , drop = FALSE]
      df_p <- df[is.finite(df$pvalue) & df$pvalue < P_THR, , drop = FALSE]
      fwrite(as.data.table(df_fdr), file.path(opt$outdir, "tables", "go", paste0(cl, "_", tag, "_FDR.tsv")), sep = "\t")
      fwrite(as.data.table(df_p), file.path(opt$outdir, "tables", "go", paste0(cl, "_", tag, "_P.tsv")), sep = "\t")

      # Plot set: if no FDR terms, fall back to top terms by p.adjust
      df_plot <- df_fdr
      if (nrow(df_plot) == 0) {
        df_plot <- df[order(df$p.adjust), , drop = FALSE]
        df_plot <- head(df_plot, MAX_TERMS)
      } else {
        df_plot <- df_plot[order(df_plot$p.adjust), , drop = FALSE]
      }

      # GeneRatio "x/y" -> numeric
      df_plot$GeneRatioNum <- as.numeric(sapply(strsplit(as.character(df_plot$GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])))
      df_plot$logFDR <- -log10(df_plot$p.adjust)
      df_plot$Description <- truncate_terms(df_plot$Description, 35)

      plot_dot(
        df = df_plot,
        out_pdf = file.path(opt$outdir, "plots", "go", paste0("go_", cl, "_", tag, ".pdf")),
        xcol = "GeneRatioNum",
        ycol = "Description",
        sizecol = "Count",
        colorcol = "logFDR",
        xlab = "Gene ratio",
        colorlab = "-log10(FDR)",
        max_terms = MAX_TERMS
      )

      invisible(enr)
    }

    do_go(det_ent, "DET_UP")
    do_go(up_ent, "INT_UP")
    do_go(dn_ent, "INT_DOWN")
  }
}

cat("\n✓ SCP enrichment complete\n")
