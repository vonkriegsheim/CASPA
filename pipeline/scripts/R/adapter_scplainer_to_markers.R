suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

getArg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

in_da   <- getArg("--scplainer-da")
out_dir <- getArg("--outdir")
prefix  <- getArg("--prefix", default = "scplainer_intensity_markers")
q_cut   <- as.numeric(getArg("--q", default = "0.05"))
min_abs <- as.numeric(getArg("--min-abs-log2fc", default = "0"))
top_n   <- as.integer(getArg("--top-n", default = "30"))

if (any(is.na(c(in_da, out_dir)))) {
  stop("Usage: Rscript adapter_scplainer_to_markers.R ",
       "--scplainer-da <scplainer_cluster_DA.tsv> ",
       "--outdir <markers_out_dir> ",
       "[--prefix scplainer_intensity_markers] ",
       "[--q 0.05] [--min-abs-log2fc 0] [--top-n 30]")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

da <- fread(in_da)

req <- c("protein_group","contrast","cluster_num","cluster_den","log2FC","pvalue","qvalue","ok")
miss <- setdiff(req, colnames(da))
if (length(miss) > 0) {
  stop("Input DA file missing required columns: ", paste(miss, collapse = ", "))
}

# Ensure ok is logical (fread may read TRUE/FALSE as logical; but be defensive)
if (!is.logical(da$ok)) {
  da[, ok := as.logical(ok)]
}

# Build a marker-like long table compatible with your ecosystem
# Mapping choices:
# - Cluster = cluster_num
# - Protein = protein_group (stable ID)
# - log2FC / qvalue kept
# - Add optional annotations if present
markers <- da[
  ,
  .(
    Cluster = cluster_num,
    Contrast = contrast,
    Baseline = cluster_den,
    Protein = protein_group,
    log2FC = log2FC,
    pvalue = pvalue,
    qvalue = qvalue,
    ok = ok,
    n_in = if ("n_in" %in% names(da)) n_in else NA_integer_,
    n_out = if ("n_out" %in% names(da)) n_out else NA_integer_,
    Genes = if ("genes" %in% names(da)) genes else NA_character_,
    Gene = if ("gene_primary" %in% names(da)) gene_primary else NA_character_,
    Description = if ("description" %in% names(da)) description else NA_character_
  )
]

# Write full table (unfiltered)
out_all <- file.path(out_dir, paste0(prefix, ".tsv"))
fwrite(markers, out_all, sep = "\t")

# Filtered significant set for convenience
markers_sig <- markers[
  ok == TRUE &
    !is.na(qvalue) & qvalue <= q_cut &
    !is.na(log2FC) & abs(log2FC) >= min_abs
]

out_sig <- file.path(out_dir, paste0(prefix, "_significant.tsv"))
fwrite(markers_sig, out_sig, sep = "\t")

# Top-N per cluster by q-value, then absolute effect size
markers_sig[, abs_log2FC := abs(log2FC)]
setorder(markers_sig, Cluster, qvalue, -abs_log2FC)

topN <- markers_sig[, head(.SD, top_n), by = Cluster]

# Clean up helper column if you don't want it in outputs
markers_sig[, abs_log2FC := NULL]
topN[, abs_log2FC := NULL]


out_top <- file.path(out_dir, paste0(prefix, "_topN.tsv"))
fwrite(topN, out_top, sep = "\t")

cat("Wrote:\n")
cat("  ", out_all, "\n", sep = "")
cat("  ", out_sig, "\n", sep = "")
cat("  ", out_top, "\n", sep = "")
