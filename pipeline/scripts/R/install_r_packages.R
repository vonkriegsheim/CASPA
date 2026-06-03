#!/usr/bin/env Rscript
# =============================================================================
# CASPA — R / Bioconductor dependency installer (native R, any OS)
#
#   Rscript pipeline/scripts/R/install_r_packages.R
#
# Use this on Windows (and anywhere you install R from CRAN rather than conda).
# Bioconductor ships Windows binaries via BiocManager even though bioconda does
# not — so native R + this script is the supported Windows route.
#
# The package list is derived from the library() calls in
# pipeline/scripts/R/*.R, plus a few transitive deps that the newest R/Bioc
# combinations fail to auto-install (notably IHW, which `scp` Imports, and its
# deps slam / fdrtool / lpsymphony).
#
# Strategy: install each missing package as a BINARY first (fast, no compiler).
# If no binary exists for your R/Bioc version, fall back to a SOURCE install —
# fine without a toolchain for pure-R / data packages (e.g. org.Hs.eg.db, IHW).
# Low-level compiled deps (slam, fdrtool, lpsymphony) are listed first so they
# are present before the source packages that need them.
# =============================================================================

options(timeout = 3600)
options(Ncpus = max(1L, parallel::detectCores() - 1L))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cat("R", as.character(getRversion()),
    "| Bioconductor", as.character(BiocManager::version()), "\n\n")

# Order matters: compiled CRAN deps and lpsymphony come before IHW (which is
# pure-R but needs them), and IHW comes before scp (which Imports IHW).
required <- c(
  # CRAN (binaries on Windows; some are deps of the Bioc source packages below)
  "slam", "fdrtool", "data.table", "matrixStats", "scales",
  "ggplot2", "igraph", "msigdbr", "optparse", "pheatmap", "uwot",
  # Bioconductor
  "lpsymphony", "S4Vectors", "AnnotationDbi", "SingleCellExperiment",
  "org.Hs.eg.db", "org.Mm.eg.db", "fgsea", "AUCell", "clusterProfiler",
  "IHW", "scp"
)

installed <- rownames(installed.packages())
missing   <- setdiff(required, installed)

if (length(missing) == 0L) {
  message("All ", length(required), " R packages already installed.")
} else {
  message("Installing ", length(missing), " missing package(s): ",
          paste(missing, collapse = ", "), "\n")
  for (pkg in missing) {
    # 1) try binary
    ok <- tryCatch({
      suppressWarnings(
        BiocManager::install(pkg, update = FALSE, ask = FALSE, type = "binary"))
      requireNamespace(pkg, quietly = TRUE)
    }, error = function(e) FALSE)
    # 2) fall back to source (no binary for this R/Bioc — usually pure-R/data)
    if (!isTRUE(ok)) {
      message("  '", pkg, "': no usable binary, trying source...")
      tryCatch(
        BiocManager::install(pkg, update = FALSE, ask = FALSE, type = "source"),
        error = function(e)
          message("  source install of '", pkg, "' failed: ", conditionMessage(e)))
    }
  }
}

# ---- Verify ----
installed     <- rownames(installed.packages())
still_missing <- setdiff(required, installed)

cat("\n=== CASPA R dependency check ===\n")
for (pkg in required) {
  cat(sprintf("%-22s %s\n", pkg,
              if (pkg %in% installed) as.character(packageVersion(pkg)) else "MISSING"))
}

if (length(still_missing)) {
  cat("\nMISSING:", paste(still_missing, collapse = ", "), "\n")
  cat("If a missing package needs compilation on Windows, install Rtools for your\n",
      "R version (https://cran.r-project.org/bin/windows/Rtools/) and re-run.\n", sep = "")
  quit(status = 1L)
}

cat("\nAll CASPA R dependencies present. scp loads:",
    requireNamespace("scp", quietly = TRUE), "\n")
