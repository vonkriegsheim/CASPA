suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
  library(scp)
})

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

rds_path   <- getArg("--rds")
model_name <- getArg("--name", default = "scplainer")
out_path   <- getArg("--out")

if (is.na(rds_path) || is.na(out_path)) {
  stop("Usage: Rscript scplainer_contrasts_from_rds.R --rds <sce_scplainer_fit.rds> --out <out.tsv> [--name scplainer]")
}

sce <- readRDS(rds_path)

write_empty <- function(msg) {
  message(msg)
  message("Writing empty output: ", out_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  fwrite(data.table(protein_group=character(), contrast=character(), cluster_num=character(),
                    cluster_den=character(), log2FC=numeric(), t=numeric(), df=numeric(),
                    pvalue=numeric(), qvalue=numeric(), ok=logical(),
                    n_in=integer(), n_out=integer()),
         out_path, sep = "\t")
  quit(status = 0)
}

if (!("cluster" %in% colnames(colData(sce)))) {
  write_empty("No 'cluster' in colData(sce) — skipping DA (likely single-group data).")
}

# Extract effect matrices + residuals from scp
eff <- tryCatch(
  scpModelEffects(sce, name = model_name),
  error = function(e) {
    write_empty(paste("scpModelEffects failed:", conditionMessage(e)))
  }
)
if (!("cluster" %in% names(eff))) {
  write_empty(paste("Model has no 'cluster' effect. Effects available:",
                    paste(names(eff), collapse = ", ")))
}

R <- tryCatch(
  scpModelResiduals(sce, name = model_name),
  error = function(e) {
    write_empty(paste("scpModelResiduals failed:", conditionMessage(e)))
  }
)
E <- eff[["cluster"]]  # features x cells

# Sanity: align cells
if (!identical(colnames(E), colnames(R))) stop("E and R columns do not align")
if (!identical(colnames(E), colnames(sce))) {
  # In some objects colnames(sce) can be NULL; enforce alignment by order
  # but warn if they exist and differ.
  if (!is.null(colnames(sce))) warning("colnames(E) != colnames(sce); proceeding using E column order.")
}
if (!identical(rownames(E), rownames(R))) stop("E and R rows do not align")

cl <- droplevels(sce$cluster)
tab <- table(cl)
baseline <- names(tab)[which.max(tab)]
others <- setdiff(levels(cl), baseline)

message("Baseline cluster: ", baseline)
message("Other clusters: ", paste(others, collapse = ", "))

# Helper: compute mean cluster effect for a group (per protein)
mean_effect <- function(g) {
  cols <- which(cl == g)
  matrixStats::rowMeans2(E[, cols, drop = FALSE], na.rm = TRUE)
}

# Helper: residual variance per protein (using available residuals)
# Use unbiased variance where possible:
resid_var <- function(cols) {
  # variance of residuals across selected cells, per protein
  x <- R[, cols, drop = FALSE]
  # rowVars uses unbiased estimator by default
  matrixStats::rowVars(x, na.rm = TRUE)
}

base_cols <- which(cl == baseline)
base_mean <- mean_effect(baseline)

out_list <- list()

for (g in others) {
  g_cols <- which(cl == g)

  # Effect-based logFC (already adjusted for other model terms)
  lfc <- mean_effect(g) - base_mean

  # Residual variance pooled from the two groups (Welch-style is possible; keep pooled for stability)
  v1 <- resid_var(g_cols)
  v0 <- resid_var(base_cols)

  n1 <- rowSums(!is.na(R[, g_cols, drop = FALSE]))
  n0 <- rowSums(!is.na(R[, base_cols, drop = FALSE]))

  # Avoid division by zero / invalid df
  ok <- (n1 >= 2) & (n0 >= 2) & is.finite(v1) & is.finite(v0)

  # Pooled variance (simple)
  sp2 <- (v1 + v0) / 2

  se <- rep(NA_real_, length(lfc))
  se[ok] <- sqrt(sp2[ok] * (1 / n1[ok] + 1 / n0[ok]))

  tval <- rep(NA_real_, length(lfc))
  tval[ok] <- lfc[ok] / se[ok]

  # Approx df: n1 + n0 - 2 (per protein)
  df <- rep(NA_real_, length(lfc))
  df[ok] <- (n1[ok] + n0[ok] - 2)

  pval <- rep(NA_real_, length(lfc))
  pval[ok] <- 2 * stats::pt(abs(tval[ok]), df = df[ok], lower.tail = FALSE)

  res <- data.frame(
    protein_group = rownames(E),
    contrast = paste0(g, "_vs_", baseline),
    cluster_num = g,
    cluster_den = baseline,
    log2FC = lfc,
    t = tval,
    df = df,
    pvalue = pval,
    ok = ok,
    n_in = n1,
    n_out = n0,
    stringsAsFactors = FALSE
  )


  # BH within this contrast
  res$qvalue <- NA_real_
  if (any(!is.na(res$pvalue))) {
    res$qvalue <- p.adjust(res$pvalue, method = "BH")
  }

  out_list[[g]] <- res
}

da <- rbindlist(out_list, use.names = TRUE, fill = TRUE)
da$ok <- !is.na(da$pvalue)

# Add annotations if present
if ("genes" %in% colnames(rowData(sce))) {
  da$genes <- rowData(sce)$genes[match(da$protein_group, rownames(sce))]
}
if ("gene_primary" %in% colnames(rowData(sce))) {
  da$gene_primary <- rowData(sce)$gene_primary[match(da$protein_group, rownames(sce))]
}
if ("description" %in% colnames(rowData(sce))) {
  da$description <- rowData(sce)$description[match(da$protein_group, rownames(sce))]
}

# Restore original protein group names (semicolons) from rowData
if ("protein_group" %in% colnames(rowData(sce))) {
  pg_map <- setNames(rowData(sce)$protein_group, rownames(sce))
  da$protein_group <- pg_map[da$protein_group]
}

fwrite(da, out_path, sep = "\t")
message("Wrote: ", out_path)
