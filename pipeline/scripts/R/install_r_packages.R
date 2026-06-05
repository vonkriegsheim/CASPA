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
# RECOMMENDED: R 4.6.0 (Bioconductor 3.23) — the combination CASPA is tested
# end-to-end against. If a package has no prebuilt binary for your R/Bioc, this
# script falls back to a source install (fine for pure-R / data packages such as
# org.Hs.eg.db, GO.db; slower).
#
# The package list is derived from the library() calls in
# pipeline/scripts/R/*.R, plus the transitive deps that some R/Bioc combinations
# fail to auto-install. These are easy to miss because the dependent package
# still appears "installed" by name while failing to LOAD — e.g.:
#   getopt   -> optparse        (every R rule)
#   GSEABase -> AUCell          (scp_aucell)
#   GO.db    -> clusterProfiler (scp_enrichment)
#   IHW (+ slam/fdrtool/lpsymphony) -> scp
# so we list them explicitly AND verify by LOADING each package, not by name.
#
# Strategy: install each missing package as a BINARY first (fast, no compiler).
# If no binary exists for your R/Bioc version, fall back to a SOURCE install —
# fine without a toolchain for pure-R / data packages (e.g. org.Hs.eg.db, IHW).
# Low-level compiled deps (slam, fdrtool, lpsymphony) and the transitive deps
# above are listed before the packages that need them.
# =============================================================================

options(timeout = 3600)
options(Ncpus = max(1L, parallel::detectCores() - 1L))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

cat("R", as.character(getRversion()),
    "| Bioconductor", as.character(BiocManager::version()), "\n\n")

# Order matters: each transitive dep precedes the package that needs it
# (getopt < optparse, GSEABase < AUCell, GO.db < clusterProfiler, the compiled
# CRAN deps + lpsymphony < IHW, IHW < scp).
required <- c(
  # CRAN (binaries on Windows; some are deps of the Bioc source packages below)
  "slam", "fdrtool", "data.table", "matrixStats", "scales",
  "ggplot2", "igraph", "msigdbr", "getopt", "optparse", "pheatmap", "uwot",
  "irlba", "RSpectra", "gplots", "isoband", "plotly",
  # Bioconductor
  "lpsymphony", "S4Vectors", "AnnotationDbi", "SingleCellExperiment",
  "GO.db", "org.Hs.eg.db", "org.Mm.eg.db", "sparseMatrixStats", "GSEABase",
  "fgsea", "AUCell", "clusterProfiler", "IHW", "scp"
)

# Verify by LOADING, not by name: a package can be present in installed.packages()
# yet fail to load because a dependency is missing (the false-green that hides
# GSEABase/GO.db/getopt gaps). loads() reinstalls those too.
loads   <- function(p) isTRUE(suppressWarnings(requireNamespace(p, quietly = TRUE)))
missing <- required[!vapply(required, loads, logical(1L))]

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

# ---- Close the dependency graph -------------------------------------------
# A package can install yet fail to load because a transitive dep
# (Depends/Imports/LinkingTo) was skipped — e.g. no binary for this Bioc
# version, so the per-package binary install silently dropped it. Repeatedly
# install any declared dep of the installed set that isn't present. This is
# generic: it fixes deps we never enumerated, not just today's known nine.
# `attempted` guarantees termination — a dep that simply can't be installed is
# tried once, then left for the load-check below to report.
attempted <- character(0L)
repeat {
  ip   <- installed.packages()
  have <- rownames(ip)
  base <- rownames(installed.packages(priority = "base"))
  need <- unique(unlist(lapply(have, function(p)
    setdiff(tools::package_dependencies(
              p, db = ip, which = c("Depends", "Imports", "LinkingTo"))[[1]],
            c(have, base)))))
  need <- setdiff(need, attempted)
  if (!length(need)) break
  message("Closing dependency gap: ", paste(need, collapse = ", "))
  attempted <- union(attempted, need)
  tryCatch(
    suppressWarnings(BiocManager::install(need, update = FALSE, ask = FALSE)),
    error = function(e) message("  some deps could not be installed: ",
                                conditionMessage(e)))
}

# ---- Verify by LOADING each package (catches installed-but-won't-load) ----
cat("\n=== CASPA R dependency check ===\n")
broken <- character(0L)
for (pkg in required) {
  ok <- loads(pkg)
  if (!ok) broken <- c(broken, pkg)
  cat(sprintf("%-22s %s\n", pkg,
              if (ok) as.character(packageVersion(pkg)) else "FAILS TO LOAD"))
}

if (length(broken)) {
  cat("\nFAILS TO LOAD:", paste(broken, collapse = ", "), "\n")
  cat("These are installed by name but a dependency is missing, or they need\n",
      "compilation. On Windows install Rtools for your R version\n",
      "(https://cran.r-project.org/bin/windows/Rtools/) and re-run.\n", sep = "")
  quit(status = 1L)
}

cat("\nAll CASPA R dependencies present and loadable.\n")
