## Seurat label transfer: Baron et al. 2016 mouse pancreas scRNA-seq (GSE84133, 2 mice)
## as reference, CASPA caerulein-injured pancreas single-cell proteome (234 cells) as query.
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
set.seed(0)
SP  <- "C:/Users/Alex/AppData/Local/Temp/claude/a4dfa14bdcb4ac46/01001967-73e9-4235-b4b9-6d4c84b21701/scratchpad/seurat"
REF <- file.path(SP, "ref_download")
QRY <- "C:/Users/Alex/OneDrive - University of Edinburgh/Documents/Mark1LLMpaper/D2/pivot_pack.tsv"
OUT <- file.path(SP, "out"); dir.create(OUT, showWarnings = FALSE)

## ---- reference -----------------------------------------------------------
read_baron <- function(f, mouse) {
  dt <- fread(f)
  meta <- data.frame(cell = dt[[1]], celltype = dt$assigned_cluster, mouse = mouse,
                     row.names = dt[[1]], stringsAsFactors = FALSE)
  m <- as.matrix(dt[, -(1:3)]); rownames(m) <- dt[[1]]
  list(counts = Matrix(t(m), sparse = TRUE), meta = meta)
}
b1 <- read_baron(file.path(REF, "GSM2230761_mouse1_umifm_counts.csv.gz"), "mouse1")
b2 <- read_baron(file.path(REF, "GSM2230762_mouse2_umifm_counts.csv.gz"), "mouse2")
stopifnot(identical(rownames(b1$counts), rownames(b2$counts)))
ref <- CreateSeuratObject(cbind(b1$counts, b2$counts), meta.data = rbind(b1$meta, b2$meta),
                          project = "Baron2016_mouse")
cat("reference:", ncol(ref), "cells,", nrow(ref), "genes\n")
print(table(ref$celltype))
ref <- NormalizeData(ref, verbose = FALSE)
ref <- FindVariableFeatures(ref, nfeatures = 3000, verbose = FALSE)

## ---- query (proteome) -------------------------------------------------------
pv <- fread(QRY)
genes <- sub(";.*$", "", pv$Genes)
keep  <- !is.na(genes) & genes != "" & !duplicated(genes)
cat("query proteins:", nrow(pv), "; unique leading gene symbols:", sum(keep), "\n")
X <- as.matrix(pv[keep, -(1:6)]); rownames(X) <- genes[keep]
cells <- colnames(X)
## log2 normalised intensities with NA = not detected -> 0 (Seurat 'data' convention)
det_rate <- rowMeans(!is.na(X))
X[is.na(X)] <- 0
q <- CreateSeuratObject(Matrix(X, sparse = TRUE), project = "CASPA_pancreas_D2")
q <- SetAssayData(q, layer = "data", new.data = Matrix(X, sparse = TRUE))
q$detected_proteins <- colSums(X != 0)
q$library <- sub("^SC_[0-9]+_(s[0-9]+)-.*$", "\\1", cells)
cat("query:", ncol(q), "cells; libraries:", paste(names(table(q$library)), table(q$library), collapse = ", "), "\n")

## ---- shared feature space -----------------------------------------------------
shared <- intersect(rownames(ref), rownames(q))
cat("shared gene symbols:", length(shared), "\n")
feats <- shared[det_rate[shared] >= 0.10]              # detected in >=10 % of query cells
cat("transfer features (shared, detected in >=10% of proteome cells):", length(feats), "\n")
ref <- ScaleData(ref, features = feats, verbose = FALSE)
ref <- RunPCA(ref, features = feats, npcs = 30, verbose = FALSE)
q   <- ScaleData(q, features = feats, verbose = FALSE)

run_transfer <- function(reduction, tag) {
  anchors <- FindTransferAnchors(reference = ref, query = q, features = feats,
                                 normalization.method = "LogNormalize",
                                 reduction = reduction, dims = 1:30,
                                 k.filter = NA, k.anchor = 5, k.score = 30, verbose = FALSE)
  na <- nrow(anchors@anchors)
  cat(sprintf("[%s] anchors: %d\n", tag, na))
  pred <- TransferData(anchorset = anchors, refdata = ref$celltype,
                       weight.reduction = if (reduction == "cca") "cca" else "pcaproject",
                       dims = 1:30, k.weight = min(20, na), verbose = FALSE)
  pred$cell <- rownames(pred)
  pred$method <- tag
  pred
}
p_cca <- run_transfer("cca", "cca")
p_pca <- run_transfer("pcaproject", "pcaproject")

summ <- function(p, tag) {
  cat(sprintf("\n== %s: predicted labels (n cells) ==\n", tag))
  print(sort(table(p$predicted.id), decreasing = TRUE))
  cat(sprintf("median max prediction score: %.2f; cells with score < 0.5: %d of %d\n",
              median(p$prediction.score.max), sum(p$prediction.score.max < 0.5), nrow(p)))
}
summ(p_cca, "CCA"); summ(p_pca, "PCA projection")
cat("\nagreement between CCA and PCA-projection calls:",
    sprintf("%.1f%%", 100 * mean(p_cca$predicted.id == p_pca$predicted.id)), "\n")

out <- data.frame(cell = p_cca$cell, library = q$library[p_cca$cell],
                  detected_proteins = q$detected_proteins[p_cca$cell],
                  seurat_cca_label = p_cca$predicted.id, seurat_cca_score = round(p_cca$prediction.score.max, 3),
                  seurat_pcaproject_label = p_pca$predicted.id, seurat_pcaproject_score = round(p_pca$prediction.score.max, 3))
sc <- p_cca[, grep("^prediction.score\\.", colnames(p_cca))]
sc <- sc[, colnames(sc) != "prediction.score.max", drop = FALSE]
colnames(sc) <- paste0("cca_score_", sub("^prediction.score.", "", colnames(sc)))
out <- cbind(out, round(sc[out$cell, ], 3))
fwrite(out, file.path(OUT, "seurat_label_transfer_per_cell.tsv"), sep = "\t")
fwrite(data.frame(feature = feats, query_detection_rate = round(det_rate[feats], 3)),
       file.path(OUT, "transfer_features.tsv"), sep = "\t")
fwrite(as.data.frame(table(ref$celltype)), file.path(OUT, "reference_celltype_counts.tsv"), sep = "\t")
cat("\nwritten:", file.path(OUT, "seurat_label_transfer_per_cell.tsv"), "\n")
