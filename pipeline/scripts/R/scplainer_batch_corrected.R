suppressPackageStartupMessages({
  library(data.table)
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
  stop("Usage: Rscript scplainer_batch_corrected.R --rds <sce.rds> --out <corrected.tsv> [--name scplainer]")
}

sce <- readRDS(rds_path)
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

# Original log2 + median-normalised expression (from pivot_pack)
expr <- assay(sce, "logint")  # proteins x cells
corrected <- expr  # start with original; subtract effects where available

# Extract model effects — these cover only the subset of proteins that converged
eff <- tryCatch(
  scpModelEffects(sce, name = model_name),
  error = function(e) {
    message("WARNING: scpModelEffects failed (", conditionMessage(e),
            ") — writing uncorrected expression as fallback.")
    NULL
  }
)

# Helper: subtract an effect matrix from corrected, aligning by rownames
subtract_effect <- function(corrected, eff_mat, label) {
  shared <- intersect(rownames(corrected), rownames(eff_mat))
  if (length(shared) == 0) {
    message("WARNING: no shared proteins for ", label, " effect — skipping")
    return(corrected)
  }
  sub <- eff_mat[shared, , drop = FALSE]
  prev <- corrected[shared, , drop = FALSE]
  result <- prev - sub
  result[is.na(sub)] <- prev[is.na(sub)]  # keep original where model didn't converge
  corrected[shared, ] <- result
  message("Subtracted ", label, " effect for ", length(shared), "/",
          nrow(corrected), " proteins")
  corrected
}

if (!is.null(eff)) {
  if ("batch" %in% names(eff)) {
    corrected <- subtract_effect(corrected, eff[["batch"]], "batch")
  } else {
    message("No batch term in model — no batch correction applied")
  }
  if ("norm_median" %in% names(eff)) {
    corrected <- subtract_effect(corrected, eff[["norm_median"]], "norm_median")
  }
}

# Build output table matching pivot_pack format:
# Protein.Group | Genes | <sample columns...>
# Restore original protein group names (semicolons) from rowData if available
pg_labels <- rownames(corrected)
if ("protein_group" %in% colnames(rowData(sce))) {
  pg_map    <- setNames(rowData(sce)$protein_group, rownames(sce))
  pg_labels <- pg_map[pg_labels]
}
out_df <- data.frame(Protein.Group = pg_labels, check.names = FALSE)

if ("genes" %in% colnames(rowData(sce))) {
  out_df$Genes <- rowData(sce)$genes[match(rownames(corrected), rowData(sce)$protein_group)]
}

# Sample columns: use sample_id (run stems) as headers for consistency with pivot_pack
stem_map <- setNames(as.character(colData(sce)$sample_id), colnames(corrected))
col_headers <- stem_map[colnames(corrected)]

expr_df <- as.data.frame(corrected, check.names = FALSE)
colnames(expr_df) <- col_headers
out_df <- cbind(out_df, expr_df)

fwrite(out_df, out_path, sep = "\t")
message("Wrote batch-corrected expression: ", out_path,
        " (", nrow(out_df), " proteins x ", ncol(expr_df), " cells)")
