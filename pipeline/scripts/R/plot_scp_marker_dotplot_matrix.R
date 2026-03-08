#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  if (!require("optparse", quietly = TRUE)) stop("Package 'optparse' required.")
  library(optparse)
})

option_list <- list(
  make_option("--detection", type="character", help="detection_markers.tsv"),
  make_option("--intensity", type="character", help="intensity_markers_detected_only.tsv"),
  make_option("--out", type="character", help="Output PDF path"),
  make_option("--top-n", type="integer", default=30, help="Top N per cluster [default %default]"),
  make_option("--label-col", type="character", default="Genes",
              help="Preferred label column (fallback to Protein) [default %default]"),
  make_option("--max-label-chars", type="integer", default=35),
  make_option("--min-size", type="double", default=0.5),
  make_option("--max-size", type="double", default=5.0)
)
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

truncate_txt <- function(x, n=35) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  too <- nchar(x) > n
  x[too] <- paste0(substr(x[too], 1, n - 3), "...")
  x
}

pick_label <- function(df, preferred="Genes") {
  if (preferred %in% colnames(df) && any(df[[preferred]] != "" & !is.na(df[[preferred]]))) {
    lab <- df[[preferred]]
  } else if ("Gene" %in% colnames(df) && any(df[["Gene"]] != "" & !is.na(df[["Gene"]]))) {
    lab <- df[["Gene"]]
  } else if ("Protein" %in% colnames(df)) {
    lab <- df[["Protein"]]
  } else {
    lab <- rep("", nrow(df))
  }
  lab <- as.character(lab)
  lab[is.na(lab) | lab == ""] <- as.character(df$Protein[is.na(lab) | lab == ""])
  truncate_txt(lab, opt$`max-label-chars`)
}

# ---- load ----
det <- fread(opt$detection, sep="\t", header=TRUE, data.table=FALSE)
int <- fread(opt$intensity, sep="\t", header=TRUE, data.table=FALSE)

if (!all(c("Cluster","Protein") %in% colnames(det))) stop("Detection table must have Cluster, Protein")
if (!all(c("Cluster","Protein") %in% colnames(int))) stop("Intensity table must have Cluster, Protein")

# labels
det$Label <- pick_label(det, opt$`label-col`)
int$Label <- pick_label(int, opt$`label-col`)

# ---- define ranking per cluster ----
# Detection ranking: qvalue asc, |log2_odds_ratio| desc
if (!("log2_odds_ratio" %in% colnames(det))) det$log2_odds_ratio <- NA_real_
if (!("qvalue" %in% colnames(det))) det$qvalue <- det$pvalue
det$rank_key <- abs(det$log2_odds_ratio)

# Intensity ranking: qvalue asc, |log2FC_detected_only| desc
if (!("log2FC_detected_only" %in% colnames(int))) int$log2FC_detected_only <- NA_real_
if (!("qvalue" %in% colnames(int))) int$qvalue <- int$pvalue
int$rank_key <- abs(int$log2FC_detected_only)

top_markers <- function(df, effect_col) {
  clusters <- sort(unique(df$Cluster))
  out <- list()
  for (cl in clusters) {
    sub <- df[df$Cluster == cl, , drop=FALSE]
    sub <- sub[!is.na(sub[[effect_col]]) & is.finite(sub[[effect_col]]), , drop=FALSE]
    sub <- sub[order(sub$qvalue, -sub$rank_key), , drop=FALSE]
    out[[cl]] <- head(sub, opt$`top-n`)
  }
  do.call(rbind, out)
}

det_top <- top_markers(det, "log2_odds_ratio")
int_top <- top_markers(int, "log2FC_detected_only")

# marker sets (union across clusters)
det_labels <- unique(det_top$Label)
int_labels <- unique(int_top$Label)

# ---- build dotplot matrix frame: for each (Label, Cluster) fill size+color ----
make_plot_df <- function(df, labels, size_col, color_col, fill_missing_as=0) {
  clusters <- sort(unique(df$Cluster))
  grid <- expand.grid(Label=labels, Cluster=clusters, stringsAsFactors = FALSE)
  key <- paste(df$Label, df$Cluster, sep="__")
  idx <- match(paste(grid$Label, grid$Cluster, sep="__"), key)

  # size: det_rate_in is best; if missing, use 0
  size <- rep(fill_missing_as, nrow(grid))
  if (size_col %in% colnames(df)) {
    v <- df[[size_col]]
    v <- suppressWarnings(as.numeric(v))
    size[!is.na(idx)] <- v[idx[!is.na(idx)]]
  }

  # color: effect size; missing -> 0
  colv <- rep(0, nrow(grid))
  if (color_col %in% colnames(df)) {
    v <- df[[color_col]]
    v <- suppressWarnings(as.numeric(v))
    colv[!is.na(idx)] <- v[idx[!is.na(idx)]]
  }

  grid$size <- size
  grid$effect <- colv
  grid
}

# We use det_rate_in for size if present; otherwise fallback to det_in/n_in if available
if (!("det_rate_in" %in% colnames(det))) {
  if (all(c("det_in","n_in") %in% colnames(det))) det$det_rate_in <- det$det_in / det$n_in
}
if (!("det_rate_in" %in% colnames(int))) {
  if (all(c("n_detected_in","n_in") %in% colnames(int))) int$det_rate_in <- int$n_detected_in / int$n_in
  # if not available, leave size as 0
}

det_plot_df <- make_plot_df(det, det_labels, "det_rate_in", "log2_odds_ratio")
int_plot_df <- make_plot_df(int, int_labels, "det_rate_in", "log2FC_detected_only")

# Order labels: show strongest overall first (sum abs effect across clusters)
order_labels <- function(df) {
  agg <- aggregate(abs(effect) ~ Label, data=df, FUN=sum)
  agg <- agg[order(agg[,2], decreasing=TRUE), , drop=FALSE]
  agg$Label
}
det_order <- order_labels(det_plot_df)
int_order <- order_labels(int_plot_df)

det_plot_df$Label <- factor(det_plot_df$Label, levels=rev(det_order))
int_plot_df$Label <- factor(int_plot_df$Label, levels=rev(int_order))

# ---- plotting helper ----
dotplot <- function(df, title, effect_limits=NULL) {
  # dynamic range for color
  lim <- max(abs(df$effect), na.rm=TRUE)
  if (!is.finite(lim) || lim == 0) lim <- 1
  if (!is.null(effect_limits)) lim <- effect_limits

  n_rows <- length(unique(df$Label))
  n_cols <- length(unique(df$Cluster))

  # Adaptive portrait height: compact but not cramped
  h <- max(5, 0.22 * n_rows + 2.0)
  w <- max(6, 0.65 * n_cols + 3.0)

  p <- ggplot(df, aes(x=Cluster, y=Label)) +
    geom_point(aes(size=size, color=effect), alpha=0.95) +
    scale_color_gradient2(
      low="#2166AC", mid="white", high="#B2182B",
      limits=c(-lim, lim), oob=scales::squish, name="Effect"
    ) +
    scale_size_continuous(range=c(opt$`min-size`, opt$`max-size`), name="Detection") +
    theme_bw(base_size=10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size=8),
      axis.text.x = element_text(size=8),
      plot.margin = grid::unit(c(6, 18, 6, 6), "pt"),
      legend.position = "right"
    ) +
    labs(title="", x="", y="")
  list(p=p, w=w, h=h)
}

# ---- write PDF (2 pages) ----
dir.create(dirname(opt$out), recursive=TRUE, showWarnings=FALSE)
grDevices::pdf(opt$out, width=8, height=10, useDingbats=FALSE)

# Detection page
res1 <- dotplot(det_plot_df, "Detection markers (log2 OR)")
print(res1$p + ggtitle("") )

# Intensity page
res2 <- dotplot(int_plot_df, "Intensity markers (detected-only log2FC)")
print(res2$p + ggtitle("") )

dev.off()
cat(sprintf("Saved dotplot PDF: %s\n", opt$out))
