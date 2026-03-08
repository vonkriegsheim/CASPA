#!/usr/bin/env Rscript
################################################################################
# run_auc_from_pivot.R  (AUCell + plots)
#
# AUCell pathway activity scoring for SCP without intensity imputation.
# Input is pivot-like (log2, NA for missing). Missing is treated as 0 for AUCell.
#
# Outputs:
# - aucell_scores.tsv (pathways x cells)
# - heatmap (top variable pathways), PCA, optional UMAP
#
# Styling consistent with pipeline:
# - diverging palette: #2166AC -> white -> #B2182B
# - pathway name cleaning + truncation (35 chars)
# - heatmap clustering by 1 - correlation (average linkage)
################################################################################

suppressPackageStartupMessages({
  if (!require("optparse", quietly = TRUE)) stop("Package 'optparse' required.")
  library(optparse)
})

option_list <- list(
  make_option("--pivot", type="character", help="scp_pivot.tsv (log2, NA for missing)"),
  make_option("--annotation", type="character", help="scp_annotation.tsv (Run=stem, Condition, Batch)"),
  make_option("--outdir", type="character", help="Output directory"),
  make_option("--species", type="character", default="auto", help="auto|human|mouse [default %default]"),
  make_option("--collections", type="character", default="hallmark,reactome",
              help="Comma-separated: hallmark,reactome,go_bp,kegg [default %default]"),
  make_option("--top-n-heatmap", type="integer", default=50,
              help="Top N pathways by MAD to show in heatmap [default %default]"),
  make_option("--min_frac_geneset_present", type="double", default=0.20,
              help="Filter out gene sets with < this fraction present [default %default]"),
  make_option("--seed", type="integer", default=42),   # <-- COMMA HERE

  make_option("--cluster-cells", action="store_true", default=FALSE,
              help="Run secondary Leiden clustering in AUCell pathway space."),
  make_option("--leiden-resolution", type="double", default=0.8),
  make_option("--n-neighbors", type="integer", default=15),
  make_option("--top-n-cluster", type="integer", default=300)
)
opt <- parse_args(OptionParser(option_list=option_list))
set.seed(opt$seed)

# ---- soft-fail dependency checks ----
has_pkg <- function(x) suppressWarnings(requireNamespace(x, quietly=TRUE))

dir.create(opt$outdir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$outdir, "tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$outdir, "plots"), recursive=TRUE, showWarnings=FALSE)

soft_fail <- function(msg) {
  writeLines(msg, con=file.path(opt$outdir, "SKIPPED_missing_packages.txt"))
  cat(msg, "\n")
  quit(status=0)
}

req <- c("data.table","msigdbr","AUCell","ggplot2")
miss <- req[!vapply(req, has_pkg, logical(1))]
if (length(miss) > 0) soft_fail(paste0("AUCell skipped: missing packages: ", paste(miss, collapse=", ")))

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(AUCell)
  library(ggplot2)
})

can_pheatmap <- has_pkg("pheatmap")
can_umap <- has_pkg("uwot")
if (!can_pheatmap) warning("Package 'pheatmap' not installed: heatmap will be skipped.")
if (!can_umap) warning("Package 'uwot' not installed: UMAP will be skipped.")

can_igraph <- has_pkg("igraph")
if (isTRUE(opt$`cluster-cells`) && !can_igraph) {
  warning("igraph not available: AUCell-space Leiden clustering will be skipped.")
}

# ---- helpers ----
stem_from_header <- function(h) {
  h <- as.character(h)
  h2 <- basename(gsub("\\\\", "/", h))
  h2 <- sub("\\.(raw|d|mzml|wiff)$", "", h2, ignore.case=TRUE)
  h2
}

detect_species <- function(genes) {
  g <- head(unique(genes), 200)
  g <- g[!is.na(g) & g != ""]
  if (length(g) == 0) return("human")
  n_upper <- sum(grepl("^[A-Z][A-Z0-9-]+$", g))
  if (n_upper/length(g) > 0.6) "human" else "mouse"
}

clean_gene <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("[;,].*$", "", x)
  x <- trimws(x)
  x
}

clean_pathway_names <- function(names, max_chars = 35) {
  n <- gsub("^HALLMARK_|^REACTOME_|^GOBP_|^GOCC_|^GOMF_|^KEGG_", "", names)
  n <- gsub("_", " ", n)
  n <- tolower(n)
  n <- paste0(toupper(substring(n, 1, 1)), substring(n, 2))
  too_long <- nchar(n) > max_chars
  n[too_long] <- paste0(substr(n[too_long], 1, max_chars - 3), "...")
  n
}

# ---- load inputs ----
pivot <- fread(opt$pivot, sep="\t", header=TRUE, data.table=FALSE)
ann <- fread(opt$annotation, sep="\t", header=TRUE, data.table=FALSE)

if (!("Run" %in% colnames(ann)) || !("Condition" %in% colnames(ann))) {
  stop("annotation must have Run and Condition columns")
}
ann$Run <- as.character(ann$Run)
ann$Condition <- as.character(ann$Condition)
ann <- ann[order(ann$Condition, ann$Run), , drop=FALSE]

if (!("Genes" %in% colnames(pivot))) stop("pivot must contain Genes column")
pivot$Genes <- clean_gene(pivot$Genes)
pivot <- pivot[!is.na(pivot$Genes) & pivot$Genes != "", , drop=FALSE]

# map pivot sample columns -> stems, match to annotation
all_cols <- colnames(pivot)
sample_cols <- all_cols[grepl("(?i)\\.raw$|\\\\|/|\\.d$|\\.mzml$|\\.wiff$", all_cols, perl=TRUE)]
if (length(sample_cols) > 0) {
  stems <- vapply(sample_cols, stem_from_header, character(1))
  stem_to_col <- setNames(sample_cols, stems)
} else {
  # pivot already has bare stems (e.g. from scp_pivot) — match directly
  sample_cols <- all_cols[all_cols %in% ann$Run]
  stem_to_col <- setNames(sample_cols, sample_cols)
}

run_cols <- stem_to_col[ann$Run[ann$Run %in% names(stem_to_col)]]
run_cols <- unname(run_cols[!is.na(run_cols)])
cat(sprintf("Matched run columns: %d\n", length(run_cols)))
if (length(run_cols) < 3) stop("Need >=3 run columns matched between pivot and annotation (via stem mapping).")

# ---- expression matrix (genes x cells) for AUCell ----
expr <- as.matrix(pivot[, run_cols, drop=FALSE])
expr <- apply(expr, 2, as.numeric)
rownames(expr) <- pivot$Genes

# Treat missing as 0 = not detected (appropriate for AUCell in SCP context)
expr[!is.finite(expr)] <- 0

# Collapse duplicate genes by max
if (any(duplicated(rownames(expr)))) {
  dt <- as.data.table(expr, keep.rownames="Gene")
  expr <- as.data.frame(dt[, lapply(.SD, max, na.rm=TRUE), by=Gene])
  rownames(expr) <- expr$Gene
  expr$Gene <- NULL
  expr <- as.matrix(expr)
}

cat(sprintf("Expression matrix: %d genes x %d cells\n", nrow(expr), ncol(expr)))

# ---- species + gene sets ----
species <- opt$species
if (species == "auto") {
  species <- detect_species(rownames(expr))
  cat(sprintf("Auto-detected species: %s\n", species))
}
species_msigdb <- ifelse(species=="human", "Homo sapiens", "Mus musculus")

collections <- trimws(strsplit(opt$collections, ",")[[1]])
gene_sets <- list()
for (coll in collections) {
  gs <- switch(coll,
    "hallmark" =, "msigdb_hallmark" = { m <- msigdbr(species=species_msigdb, collection="H"); split(m$gene_symbol, m$gs_name) },
    "reactome" =, "msigdb_reactome" = { m <- msigdbr(species=species_msigdb, collection="C2", subcollection="CP:REACTOME"); split(m$gene_symbol, m$gs_name) },
    "go_bp"    =, "msigdb_go_bp"    = { m <- msigdbr(species=species_msigdb, collection="C5", subcollection="GO:BP"); split(m$gene_symbol, m$gs_name) },
    "kegg"     =, "msigdb_kegg"     = { m <- msigdbr(species=species_msigdb, collection="C2", subcollection="CP:KEGG"); split(m$gene_symbol, m$gs_name) },
    NULL
  )
  if (!is.null(gs)) gene_sets <- c(gene_sets, gs)
}
cat(sprintf("Loaded gene sets: %d\n", length(gene_sets)))
if (length(gene_sets) == 0) stop("No gene sets loaded")

# Filter gene sets by fraction present
present_genes <- rownames(expr)
keep_sets <- vapply(gene_sets, function(gs) {
  frac <- sum(gs %in% present_genes) / max(length(gs), 1)
  frac >= opt$min_frac_geneset_present
}, logical(1))
gene_sets2 <- gene_sets[keep_sets]
cat(sprintf("Gene sets after presence filter (>=%.2f): %d\n", opt$min_frac_geneset_present, length(gene_sets2)))

if (length(gene_sets2) == 0) {
  fwrite(data.table(Pathway=character(0)), file.path(opt$outdir, "tables", "aucell_scores.tsv"), sep="\t")
  writeLines("No gene sets passed presence filter.", con=file.path(opt$outdir, "RUN_OK.txt"))
  quit(status=0)
}

# ---- AUCell ----
rankings <- AUCell_buildRankings(expr, plotStats=FALSE, verbose=TRUE)

det_per_cell <- colSums(expr > 0)
aucMaxRank <- max(50, round(median(det_per_cell) * 0.2))
aucMaxRank <- min(aucMaxRank, nrow(expr))
cat(sprintf("AUCell aucMaxRank=%d\n", aucMaxRank))

cellsAUC <- AUCell_calcAUC(gene_sets2, rankings, aucMaxRank=aucMaxRank, verbose=TRUE)
auc_mat <- as.matrix(getAUC(cellsAUC))
cat(sprintf("AUC matrix: %d pathways x %d cells\n", nrow(auc_mat), ncol(auc_mat)))

# ---- write tables ----
auc_df <- as.data.frame(auc_mat)
auc_df$Pathway <- rownames(auc_df)
auc_df <- auc_df[, c("Pathway", colnames(auc_mat)), drop=FALSE]
fwrite(as.data.table(auc_df), file.path(opt$outdir, "tables", "aucell_scores.tsv"), sep="\t")

long <- data.table::melt(as.data.table(auc_df), id.vars="Pathway", variable.name="SampleColumn", value.name="AUC")
# Add stem + condition
long[, Run := vapply(SampleColumn, stem_from_header, character(1))]
long <- merge(long, as.data.table(ann[, c("Run","Condition")]), by="Run", all.x=TRUE)
fwrite(long, file.path(opt$outdir, "tables", "aucell_scores_long.tsv"), sep="\t")

writeLines("AUCell run completed; outputs written.", con=file.path(opt$outdir, "RUN_OK.txt"))

# ---- plots ----
# HEATMAP (top variable)
if (can_pheatmap && nrow(auc_mat) >= 2) {
  suppressPackageStartupMessages(library(pheatmap))

  madv <- apply(auc_mat, 1, mad, na.rm=TRUE)
  top_n <- min(opt$`top-n-heatmap`, length(madv))
  top_idx <- order(madv, decreasing=TRUE)[seq_len(top_n)]
  hm <- auc_mat[top_idx, , drop=FALSE]

  fwrite(
    data.table(Pathway=rownames(hm), MAD=madv[top_idx]),
    file.path(opt$outdir, "tables", "aucell_heatmap_top_variable_pathways.tsv"),
    sep="\t"
  )

  hm_z <- t(scale(t(hm)))
  hm_z[!is.finite(hm_z)] <- 0
  rownames(hm_z) <- clean_pathway_names(rownames(hm_z), 35)

  # column stems for annotation
  col_stems2 <- vapply(colnames(hm_z), stem_from_header, character(1))
  annot_col <- data.frame(
    Condition = ann$Condition[match(col_stems2, ann$Run)],
    row.names = colnames(hm_z)
  )

  # cluster by 1 - corr
  row_cor <- cor(t(hm_z), method="pearson", use="pairwise.complete.obs")
  col_cor <- cor(hm_z, method="pearson", use="pairwise.complete.obs")
  row_cor[!is.finite(row_cor)] <- 0; diag(row_cor) <- 1
  col_cor[!is.finite(col_cor)] <- 0; diag(col_cor) <- 1

  row_hc <- hclust(as.dist(1 - row_cor), method="average")
  col_hc <- hclust(as.dist(1 - col_cor), method="average")

  pal <- colorRampPalette(c("#2166AC","white","#B2182B"))(100)
  breaks <- seq(-2, 2, length.out=101)

  hm_width <- min(10, ncol(hm_z)*0.18 + 3)   # many cells -> keep narrow
  hm_height <- max(8, nrow(hm_z)*0.25 + 2)

  out_hm <- file.path(opt$outdir, "plots", paste0("aucell_heatmap_top", top_n, ".pdf"))
  grDevices::pdf(out_hm, width=hm_width, height=hm_height, useDingbats=FALSE)
  pheatmap::pheatmap(
    hm_z,
    color=pal,
    breaks=breaks,
    cluster_rows=row_hc,
    cluster_cols=col_hc,
    annotation_col=annot_col,
    fontsize=9,
    fontsize_row=8,
    fontsize_col=6,
    angle_col=45,
    border_color=NA,
    main=""
  )
  dev.off()
}


# -------------------------
# PCA (robust): drop constant cells and constant pathways
# -------------------------

# Drop constant/NA-variance CELLS (columns of auc_mat)
cell_sd <- apply(auc_mat, 2, sd, na.rm = TRUE)
keep_cells <- is.finite(cell_sd) & cell_sd > 0

if (sum(!keep_cells) > 0) {
  cat(sprintf("PCA: dropping %d/%d cells with zero/NA variance in AUCell scores\n",
              sum(!keep_cells), length(keep_cells)))
}

auc_mat_pca <- auc_mat[, keep_cells, drop = FALSE]

# Drop constant/NA-variance PATHWAYS (rows of auc_mat_pca)
path_sd <- apply(auc_mat_pca, 1, sd, na.rm = TRUE)
keep_path <- is.finite(path_sd) & path_sd > 0

if (sum(!keep_path) > 0) {
  cat(sprintf("PCA: dropping %d/%d pathways with zero/NA variance after cell filtering\n",
              sum(!keep_path), length(keep_path)))
}

auc_mat_pca2 <- auc_mat_pca[keep_path, , drop = FALSE]

if (ncol(auc_mat_pca2) < 3 || nrow(auc_mat_pca2) < 3) {
  writeLines(
    sprintf("PCA skipped: matrix too small after filtering (pathways=%d, cells=%d).",
            nrow(auc_mat_pca2), ncol(auc_mat_pca2)),
    con = file.path(opt$outdir, "tables", "PCA_SKIPPED.txt")
  )
} else {
  pca <- prcomp(t(auc_mat_pca2), center = TRUE, scale. = TRUE)

  var_exp <- (pca$sdev^2) / sum(pca$sdev^2)
  pc1 <- var_exp[1] * 100
  pc2 <- var_exp[2] * 100

  pca_df <- data.frame(
    SampleColumn = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  pca_df$Run <- vapply(pca_df$SampleColumn, stem_from_header, character(1))
  pca_df$Condition <- ann$Condition[match(pca_df$Run, ann$Run)]

  fwrite(as.data.table(pca_df),
         file.path(opt$outdir, "tables", "aucell_pca_coords.tsv"),
         sep = "\t")

  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 2.4, alpha = 0.85) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = "right") +
    labs(title = "", x = sprintf("PC1 (%.1f%%)", pc1), y = sprintf("PC2 (%.1f%%)", pc2)) +
    coord_equal()

  ggsave(file.path(opt$outdir, "plots", "aucell_pca.pdf"), p_pca, width = 6, height = 6)
}


# UMAP (square; optional) — use filtered matrix
if (can_umap) {
  if (exists("auc_mat_pca2") && ncol(auc_mat_pca2) >= 3 && nrow(auc_mat_pca2) >= 3) {
    suppressPackageStartupMessages(library(uwot))
    set.seed(opt$seed)

    um <- uwot::umap(
      t(auc_mat_pca2),
      n_neighbors = min(15, ncol(auc_mat_pca2) - 1),
      min_dist = 0.1,
      metric = "euclidean"
    )

    um_df <- data.frame(
      SampleColumn = colnames(auc_mat_pca2),
      UMAP1 = um[, 1],
      UMAP2 = um[, 2]
    )
    um_df$Run <- vapply(um_df$SampleColumn, stem_from_header, character(1))
    um_df$Condition <- ann$Condition[match(um_df$Run, ann$Run)]

    fwrite(as.data.table(um_df), file.path(opt$outdir, "tables", "aucell_umap_coords.tsv"), sep="\t")

    p_um <- ggplot(um_df, aes(x = UMAP1, y = UMAP2, color = Condition)) +
      geom_point(size = 2.4, alpha = 0.85) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank(), legend.position = "right") +
      labs(title = "", x = "UMAP1", y = "UMAP2") +
      coord_equal()

    ggsave(file.path(opt$outdir, "plots", "aucell_umap.pdf"), p_um, width = 6, height = 6)
  } else {
    writeLines("UMAP skipped: matrix too small or PCA filtering object missing.",
               con = file.path(opt$outdir, "tables", "UMAP_SKIPPED.txt"))
  }
}

# =============================================================================
# Secondary clustering: AUCell-space PCA/UMAP/Leiden (optional)
# =============================================================================
if (isTRUE(opt$`cluster-cells`) && can_umap && can_igraph) {
  suppressPackageStartupMessages({
    library(uwot)
    library(igraph)
  })

  # Build a clustering matrix: use top variable pathways for stability
  madv_all <- apply(auc_mat, 1, mad, na.rm=TRUE)
  top_n <- min(opt$`top-n-cluster`, length(madv_all))
  top_idx <- order(madv_all, decreasing=TRUE)[seq_len(top_n)]
  A <- auc_mat[top_idx, , drop=FALSE]   # pathways x cells

  # Row z-score (pathway-wise) to focus on patterns not absolute AUC scale
  A_z <- t(scale(t(A)))
  A_z[!is.finite(A_z)] <- 0

  # PCA on cells
  pca2 <- prcomp(t(A_z), center=TRUE, scale.=TRUE)
  pcs <- pca2$x[, 1:20, drop=FALSE]

  # UMAP on PCs
  set.seed(opt$seed)
  nn <- min(opt$`n-neighbors`, nrow(pcs) - 1)
  um2 <- uwot::umap(pcs, n_neighbors=nn, min_dist=0.1, metric="euclidean")

  # KNN graph on PC space (cells x PCs)
  D <- as.matrix(dist(pcs))
  nn <- min(opt$`n-neighbors`, nrow(pcs) - 1)

  # For each cell i, get its nn nearest neighbors (excluding itself)
  knn_list <- lapply(seq_len(nrow(D)), function(i) {
    order(D[i, ])[2:(nn + 1)]
  })

  # Build edge list as a 2-column matrix (i, j)
  edge_mat <- do.call(
    rbind,
    lapply(seq_along(knn_list), function(i) {
      cbind(rep(i, length(knn_list[[i]])), knn_list[[i]])
    })
  )

  # Ensure integer matrix with exactly two columns
  edge_mat <- matrix(as.integer(edge_mat), ncol = 2)

  g <- igraph::graph_from_edgelist(edge_mat, directed = FALSE)
  g <- igraph::simplify(g)


  cl <- igraph::cluster_leiden(g, resolution = opt$`leiden-resolution`)
  leiden_ids <- as.integer(cl$membership) - 1L
  cond_auc <- paste0("AUC", leiden_ids)

  # Write assignments (map SampleColumn -> Run stem)
  sample_cols_used <- colnames(A_z)
  run_stems_used <- vapply(sample_cols_used, stem_from_header, character(1))

  out_df <- data.frame(
    SampleColumn = sample_cols_used,
    Run = run_stems_used,
    Condition_AUCell = cond_auc,
    umap_1 = um2[, 1],
    umap_2 = um2[, 2],
    PC1 = pcs[, 1],
    PC2 = pcs[, 2],
    stringsAsFactors = FALSE
  )

  fwrite(as.data.table(out_df),
         file.path(opt$outdir, "tables", "aucell_cluster_assignments.tsv"),
         sep="\t")

  # Cluster sizes
  sizes <- as.data.table(table(out_df$Condition_AUCell))
  setnames(sizes, c("Condition_AUCell", "n_cells"))
  fwrite(sizes, file.path(opt$outdir, "tables", "aucell_cluster_sizes.tsv"), sep="\t")

  # Contingency vs primary clustering
  # primary Condition comes from annotation (Run stem)
  prim <- ann[, c("Run", "Condition")]
  comb <- merge(out_df[, c("Run", "Condition_AUCell")], prim, by="Run", all.x=TRUE)
  tab <- as.data.table(table(comb$Condition_AUCell, comb$Condition))
  setnames(tab, c("Condition_AUCell", "Condition_primary", "n"))
  fwrite(tab, file.path(opt$outdir, "tables", "aucell_cluster_vs_primary.tsv"), sep="\t")

  # Plot UMAP coloured by AUCell cluster
  p_auc_um <- ggplot(out_df, aes(x=umap_1, y=umap_2, color=Condition_AUCell)) +
    geom_point(size=2.2, alpha=0.85) +
    theme_bw(base_size=11) +
    theme(panel.grid.minor=element_blank(), legend.position="right") +
    labs(title="", x="UMAP1", y="UMAP2") +
    coord_equal()

  ggsave(file.path(opt$outdir, "plots", "aucell_umap_by_cluster.pdf"),
         p_auc_um, width=6, height=6)
}

cat("\n✓ AUCell complete (tables + plots)\n")
