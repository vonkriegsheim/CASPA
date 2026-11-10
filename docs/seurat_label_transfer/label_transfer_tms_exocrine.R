## Seurat label transfer: Tabula Muris Senis droplet pancreas (official cell_ontology_class
## annotations, all ages) as reference; CASPA caerulein-injured pancreas proteome (234 cells) as query.
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(data.table)})
set.seed(0)
SP  <- "C:/Users/Alex/AppData/Local/Temp/claude/a4dfa14bdcb4ac46/01001967-73e9-4235-b4b9-6d4c84b21701/scratchpad/seurat"
REF <- file.path(SP, "tms_ref")
QRY <- "C:/Users/Alex/OneDrive - University of Edinburgh/Documents/Mark1LLMpaper/D2/pivot_pack.tsv"
OUT <- file.path(SP, "out"); dir.create(OUT, showWarnings = FALSE)

## ---- reference -----------------------------------------------------------
cnt  <- as(readMM(file.path(REF, "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(REF, "genes.tsv"))
meta <- fread(file.path(REF, "meta.tsv"), data.table = FALSE); rownames(meta) <- meta$cell
rownames(cnt) <- genes; colnames(cnt) <- meta$cell
meta$celltype <- as.character(meta$cell_ontology_class)
## drop very rare classes (< 10 cells) so k.score / anchor weighting is stable
tab <- table(meta$celltype); keepct <- c("pancreatic acinar cell","pancreatic ductal cel","leukocyte","endothelial cell","pancreatic stellate cell")
meta <- meta[meta$celltype %in% keepct, ]; cnt <- cnt[, rownames(meta)]
ref <- CreateSeuratObject(cnt, meta.data = meta, project = "TMS_droplet_pancreas")
cat("reference:", ncol(ref), "cells,", nrow(ref), "genes\n"); print(sort(table(ref$celltype), decreasing = TRUE))
ref <- NormalizeData(ref, verbose = FALSE)

## ---- query (proteome) -------------------------------------------------------
pv <- fread(QRY)
g <- sub(";.*$", "", pv$Genes); keep <- !is.na(g) & g != "" & !duplicated(g)
X <- as.matrix(pv[keep, -(1:6)]); rownames(X) <- g[keep]
det_rate <- rowMeans(!is.na(X)); X[is.na(X)] <- 0
q <- CreateSeuratObject(Matrix(X, sparse = TRUE), project = "CASPA_pancreas_D2")
q <- SetAssayData(q, layer = "data", new.data = Matrix(X, sparse = TRUE))
q$detected_proteins <- colSums(X != 0)
q$library <- sub("^SC_[0-9]+_(s[0-9]+)-.*$", "\\1", colnames(X))

shared <- intersect(rownames(ref), rownames(q))
feats  <- shared[det_rate[shared] >= 0.10]
cat("shared symbols:", length(shared), "; transfer features:", length(feats), "\n")
ref <- ScaleData(ref, features = feats, verbose = FALSE)
ref <- RunPCA(ref, features = feats, npcs = 30, verbose = FALSE)
q   <- ScaleData(q, features = feats, verbose = FALSE)

run_transfer <- function(reduction, tag) {
  anchors <- FindTransferAnchors(reference = ref, query = q, features = feats,
                                 normalization.method = "LogNormalize", reduction = reduction,
                                 dims = 1:30, k.filter = NA, k.anchor = 5, k.score = 30, verbose = FALSE)
  na <- nrow(anchors@anchors); cat(sprintf("[%s] anchors: %d\n", tag, na))
  pred <- TransferData(anchorset = anchors, refdata = ref$celltype,
                       weight.reduction = if (reduction == "cca") "cca" else "pcaproject",
                       dims = 1:30, k.weight = min(20, na), verbose = FALSE)
  pred$cell <- rownames(pred); pred
}
p_cca <- run_transfer("cca", "cca")
p_pca <- run_transfer("pcaproject", "pcaproject")
for (x in list(list(p_cca, "CCA"), list(p_pca, "PCA projection"))) {
  cat(sprintf("\n== %s ==\n", x[[2]])); print(sort(table(x[[1]]$predicted.id), decreasing = TRUE))
  cat(sprintf("median max score %.2f; score < 0.5: %d of %d\n", median(x[[1]]$prediction.score.max),
              sum(x[[1]]$prediction.score.max < 0.5), nrow(x[[1]])))
}
cat("\nCCA vs PCA-projection agreement:", sprintf("%.1f%%", 100 * mean(p_cca$predicted.id == p_pca$predicted.id)), "\n")

out <- data.frame(cell = p_cca$cell, library = q$library[p_cca$cell], detected_proteins = q$detected_proteins[p_cca$cell],
                  tms_cca_label = p_cca$predicted.id, tms_cca_score = round(p_cca$prediction.score.max, 3),
                  tms_pcaproject_label = p_pca$predicted.id, tms_pcaproject_score = round(p_pca$prediction.score.max, 3))
sc <- p_cca[, setdiff(grep("^prediction.score\\.", colnames(p_cca), value = TRUE), "prediction.score.max"), drop = FALSE]
colnames(sc) <- paste0("cca_score_", sub("^prediction.score.", "", colnames(sc)))
out <- cbind(out, round(sc[out$cell, ], 3))
fwrite(out, file.path(OUT, "seurat_TMSexocrine_label_transfer_per_cell.tsv"), sep = "\t")
fwrite(data.frame(feature = feats, query_detection_rate = round(det_rate[feats], 3)), file.path(OUT, "transfer_features_TMSexocrine.tsv"), sep = "\t")
fwrite(as.data.frame(sort(table(ref$celltype), decreasing = TRUE)), file.path(OUT, "reference_TMSexocrine_celltype_counts.tsv"), sep = "\t")
cat("written:", file.path(OUT, "seurat_TMSexocrine_label_transfer_per_cell.tsv"), "\n")
