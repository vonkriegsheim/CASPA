suppressPackageStartupMessages({
  library(data.table)
  library(SingleCellExperiment)
  library(matrixStats)
  library(scp)
  library(S4Vectors)
})

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag, default = NA_character_) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) return(default)
  args[i + 1]
}

pg_path        <- getArg("--pg-matrix")
manifest_path  <- getArg("--manifest", default = NA_character_)
ann_path       <- getArg("--annotation", default = NA_character_)
run_col        <- getArg("--run-col", default = "sample_id")
outdir         <- getArg("--outdir")
model_name     <- getArg("--name", default = "scplainer")
min_cells_feat <- as.integer(getArg("--min-cells-feature", default = "2"))
min_feat_cell  <- as.integer(getArg("--min-features-cell", default = "50"))

if (any(is.na(c(pg_path, outdir)))) {
  stop("Usage: Rscript run_scplainer_from_pg_matrix.R ",
       "--pg-matrix <pg.tsv> --outdir <outdir> ",
       "[--manifest <manifest.csv>] [--annotation <scp_annotation.tsv>] ",
       "[--run-col sample_id]")
}
if (is.na(manifest_path) && is.na(ann_path)) {
  stop("Provide either --manifest (sample_file/batch columns) or --annotation (Run/Batch columns).")
}

# Helper: strip path and extension to get a run stem (mirrors Python stem_from_header)
stem_from_header_r <- function(h) {
  h <- as.character(h)
  if (grepl("[/\\\\]", h)) h <- basename(h)
  sub("\\.(raw|d|mzml|wiff)$", "", h, ignore.case = TRUE)
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

diag_path <- file.path(outdir, "scplainer_model_diagnostics.txt")
sink(diag_path)
on.exit(sink(), add = TRUE)

cat("==== scplainer (scp) run ====\n")
cat("pg_matrix:   ", pg_path, "\n")
cat("manifest:    ", manifest_path, "\n")
cat("annotation:  ", ann_path, "\n")
cat("run_col:     ", run_col, "\n")
cat("outdir:      ", outdir, "\n")
cat("model_name:  ", model_name, "\n")
cat("min_feat_cell: ", min_feat_cell, "\n")
cat("min_cells_feat:", min_cells_feat, "\n\n")

# inputs
pg  <- fread(pg_path, sep = "\t", data.table = FALSE)

anno_cols_pg <- intersect(colnames(pg), c(
  "Protein.Group","Protein.Names","Genes","First.Protein.Description",
  "N.Sequences","N.Proteotypic.Sequences"
))
all_intensity_cols <- setdiff(colnames(pg), anno_cols_pg)

if (!is.na(manifest_path)) {
  # Classic manifest format: sample_file, sample_id (or run_col), batch
  man <- fread(manifest_path, data.table = FALSE)
  req_cols <- c("sample_file", run_col, "batch")
  missing_req <- setdiff(req_cols, colnames(man))
  if (length(missing_req) > 0) stop("Manifest missing columns: ", paste(missing_req, collapse=", "))
} else {
  # Derive manifest from annotation (Run/Condition/Batch format) via stem matching
  ann_tmp <- fread(ann_path, sep = "\t", data.table = FALSE)
  if (!"Run" %in% colnames(ann_tmp)) stop("Annotation must have a Run column.")
  stems <- sapply(all_intensity_cols, stem_from_header_r, USE.NAMES = TRUE)
  wanted <- ann_tmp$Run
  matched_mask <- stems %in% wanted
  matched_cols  <- all_intensity_cols[matched_mask]
  matched_stems <- stems[matched_mask]
  if (length(matched_cols) == 0) stop("No pg_matrix columns matched annotation Run values.")
  batch_vals <- if ("Batch" %in% colnames(ann_tmp)) {
    batch_map <- setNames(as.character(ann_tmp$Batch), ann_tmp$Run)
    sapply(matched_stems, function(s) if (s %in% names(batch_map)) batch_map[[s]] else "1")
  } else rep("1", length(matched_cols))
  man <- data.frame(sample_file = matched_cols,
                    sample_id   = matched_stems,
                    batch       = batch_vals,
                    stringsAsFactors = FALSE)
  run_col <- "sample_id"
}

ann <- NULL
if (!is.na(ann_path) && file.exists(ann_path)) {
  ann <- fread(ann_path, sep = "\t", data.table = FALSE)
  if (!all(c("Run","Condition") %in% colnames(ann))) stop("Annotation must contain Run, Condition")
}

# matrix
if (!("Protein.Group" %in% colnames(pg))) stop("pg_matrix must have Protein.Group")

# Drop rows with missing or empty Protein.Group (trailing blank lines from TSV writers)
bad_pg <- is.na(pg$Protein.Group) | pg$Protein.Group == ""
if (any(bad_pg)) {
  cat("Dropping", sum(bad_pg), "row(s) with empty/NA Protein.Group\n")
  pg <- pg[!bad_pg, , drop = FALSE]
}

miss_cols <- setdiff(man$sample_file, all_intensity_cols)
if (length(miss_cols) > 0) stop("Manifest sample_file not in matrix columns. Example: ", miss_cols[1])

# pivot_pack already log2 + median-normalised; NaN/NA marks missing values
expr_log <- as.matrix(pg[, man$sample_file, drop = FALSE])
# Sanitize protein group names: semicolons break SimpleList[character] subscripting
# in scp/S4Vectors — replace with "|" (fully reversible, "|" never appears in UniProt IDs).
safe_pg   <- gsub(";", "|", pg$Protein.Group)
rownames(expr_log) <- safe_pg
colnames(expr_log) <- man$sample_id   # use stems as cell IDs (no path separators)
mode(expr_log) <- "numeric"
# NaN from pivot → NA for R; no need to transform zeros/negatives (valid after median shift)
expr_log[is.nan(expr_log)] <- NA_real_

meta <- DataFrame(
  sample_file = man$sample_file,
  sample_id   = man[[run_col]],
  batch       = factor(man$batch)
)
rownames(meta) <- man$sample_id   # match colnames(expr_log)

sce <- SingleCellExperiment(
  assays  = list(logint = expr_log),
  colData = meta,
  rowData = DataFrame(
    protein_group      = pg[["Protein.Group"]],   # original (with ";") for output
    protein_group_safe = safe_pg,                  # safe names used as rownames(sce)
    genes = if ("Genes" %in% colnames(pg)) pg[["Genes"]] else NA_character_,
    description = if ("First.Protein.Description" %in% colnames(pg)) pg[["First.Protein.Description"]] else NA_character_
  )
)

rowData(sce)$gene_primary <- if (!all(is.na(rowData(sce)$genes))) sub(";.*$", "", rowData(sce)$genes) else NA_character_
rowData(sce)$uniprot_primary <- sub(";.*$", "", rowData(sce)$protein_group)

x <- assay(sce, "logint")
sce$n_features  <- colSums(!is.na(x))
sce$norm_median <- as.numeric(scale(matrixStats::colMedians(x, na.rm = TRUE), center = TRUE, scale = FALSE))

cat("Cells before filters:    ", ncol(sce), "\n")
cat("Features before filters: ", nrow(sce), "\n")

sce <- sce[, sce$n_features >= min_feat_cell]
sce <- sce[rowSums(!is.na(assay(sce,"logint"))) >= min_cells_feat, ]

cat("Cells after filters:     ", ncol(sce), "\n")
cat("Features after filters:  ", nrow(sce), "\n\n")

# clusters: match by sample_id (Run is stem)
if (!is.null(ann)) {
  idx <- match(sce$sample_id, ann$Run)
  if (anyNA(idx)) {
    bad <- sce$sample_id[is.na(idx)]
    stop("Some sample_id not found in annotation$Run. Examples:\n", paste(head(bad, 20), collapse="\n"))
  }
  sce$cluster <- factor(ann$Condition[idx])
  cat("Clusters present: ", paste(levels(sce$cluster), collapse=", "), "\n\n")
}

# Impute missing values with per-protein medians before model fitting.
# scpModelEffects internally joins per-protein effect vectors back to a full
# proteins x cells matrix using cell names. When a protein is detected in only
# a subset of cells, its per-protein model produces a shorter named vector, and
# the join fails with "subscript contains invalid names" because the full cell
# list contains names absent from that vector.
# Imputing with per-protein medians ensures every protein's model uses all cells,
# giving the join consistent dimensions. Filters above already exclude low-coverage
# proteins and cells, so medians here are well-estimated.
n_missing <- sum(is.na(assay(sce, "logint")))
cat("Missing values before imputation: ", n_missing,
    sprintf(" (%.1f%% of matrix)\n", 100 * n_missing / (nrow(sce) * ncol(sce))))
if (n_missing > 0) {
  cat("Imputing with per-protein medians...\n")
  logint_imp <- assay(sce, "logint")
  prot_med   <- matrixStats::rowMedians(logint_imp, na.rm = TRUE)
  for (i in seq_len(nrow(logint_imp))) {
    m <- prot_med[i]
    if (is.finite(m)) logint_imp[i, is.na(logint_imp[i, ])] <- m
  }
  assay(sce, "logint") <- logint_imp
  cat("Remaining NAs after imputation: ", sum(is.na(assay(sce, "logint"))), "\n\n")
}

# drop single-level factors
factor_terms <- c("batch", "cluster")
present <- factor_terms[factor_terms %in% colnames(colData(sce))]
nlev <- sapply(present, function(v) nlevels(as.factor(colData(sce)[[v]])))
drop_terms <- names(nlev)[nlev < 2]
use_batch   <- ("batch" %in% present)   && !("batch" %in% drop_terms)
use_cluster <- ("cluster" %in% present) && !("cluster" %in% drop_terms)

# Drop norm_median if near-constant (data already median-normalised in pivot_pack)
nm_sd <- sd(sce$norm_median, na.rm = TRUE)
use_norm_median <- is.finite(nm_sd) && nm_sd > 1e-6
if (!use_norm_median) cat("norm_median near-constant (sd=", nm_sd, "), dropping from model\n")



cat("Factor levels:\n")
for (v in present) cat("  ", v, ": ", nlevels(as.factor(colData(sce)[[v]])), "\n")
if (length(drop_terms) > 0) cat("Dropping single-level terms: ", paste(drop_terms, collapse=", "), "\n")
cat("\n")

# Build formula from active terms
terms <- c()
if (use_norm_median) terms <- c(terms, "norm_median")
if (use_batch)       terms <- c(terms, "batch")
if (use_cluster)     terms <- c(terms, "cluster")
if (length(terms) > 0) {
  model_formula <- as.formula(paste("~ 1 +", paste(terms, collapse = " + ")))
} else {
  model_formula <- ~ 1
}
cat("Model formula: ", deparse(model_formula), "\n\n")

cat("Running scpModelWorkflow...\n")
sce <- scpModelWorkflow(
  object  = sce,
  formula = model_formula,
  i       = "logint",
  name    = model_name,
  verbose = TRUE
)
cat("...done\n\n")

# Persist fitted object FIRST — downstream steps need this even if variance analysis fails
rds_path <- file.path(outdir, "sce_scplainer_fit.rds")
saveRDS(sce, rds_path)
cat("Saved fitted SCE: ", rds_path, "\n\n")

# variance — wrapped in tryCatch because scpVarianceAnalysis can fail when
# some protein model fits are rank-deficient (missing intercept term)
cat("Running scpVarianceAnalysis...\n")
var_ok <- tryCatch({
  var <- scpVarianceAnalysis(sce, name = model_name)
  fwrite(as.data.frame(var), file.path(outdir, "scplainer_variance_explained.tsv"), sep="\t")
  cat("Wrote scplainer_variance_explained.tsv\n\n")
  TRUE
}, error = function(e) {
  cat("WARNING: scpVarianceAnalysis failed: ", conditionMessage(e), "\n")
  cat("Writing empty variance table as fallback.\n\n")
  fwrite(data.frame(feature = character(0), Residuals = numeric(0)),
         file.path(outdir, "scplainer_variance_explained.tsv"), sep="\t")
  FALSE
})

# effects (often effect matrices)
tryCatch({
  cat("Extracting scpModelEffects...\n")
  eff <- scpModelEffects(sce, name = model_name)
  cat("Effects available: ", paste(names(eff), collapse=", "), "\n\n")

  if (use_cluster && ("cluster" %in% names(eff))) {
    cl_mat <- eff[["cluster"]]
    cat("cluster effect class: ", paste(class(cl_mat), collapse=", "), "\n")
    cat("cluster effect dim:   ", paste(dim(cl_mat), collapse=" x "), "\n\n")

    cl_df <- data.frame(protein_group = rownames(cl_mat), cl_mat, check.names = FALSE)
    cl_df$genes <- rowData(sce)$genes[match(cl_df$protein_group, rowData(sce)$protein_group)]
    cl_df$gene_primary <- rowData(sce)$gene_primary[match(cl_df$protein_group, rowData(sce)$protein_group)]
    cl_df$description <- rowData(sce)$description[match(cl_df$protein_group, rowData(sce)$protein_group)]
    fwrite(cl_df, file.path(outdir, "scplainer_cluster_effect_matrix.tsv"), sep="\t")
    cat("Wrote scplainer_cluster_effect_matrix.tsv\n\n")
  }

  # Dump keys to help find coefficients/covariance
  keys_path <- file.path(outdir, "scplainer_model_keys.txt")
  k <- capture.output({
    cat("scpModelNames:\n"); print(scpModelNames(sce))
    cat("\nmetadata names:\n"); print(names(metadata(sce)))
    cat("\nint_metadata names:\n"); print(names(int_metadata(sce)))
    cat("\ncolData names:\n"); print(colnames(colData(sce)))
    cat("\nassays:\n"); print(assayNames(sce))
  })
  writeLines(k, keys_path)
  cat("Wrote keys: ", keys_path, "\n\n")
}, error = function(e) {
  cat("WARNING: scpModelEffects/keys extraction failed: ", conditionMessage(e), "\n")
  cat("Primary outputs (RDS + variance TSV) already written — continuing.\n\n")
})
cat("DONE\n")
