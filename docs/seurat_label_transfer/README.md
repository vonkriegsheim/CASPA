# Seurat label transfer onto the CASPA pancreas single-cell proteome (10 Nov 2026)

Supporting analysis for the Discussion paragraph on reference-based annotation tools and Supplementary Table 3.
The R scripts carry the absolute paths of the analysis PC (`SP`, `QRY`): edit those two lines to rerun. `export_tms.py`
converts the Tabula Muris Senis h5ad to MTX; `compare_to_clusters.py` builds the cluster x label table from
`scp/clustering/scp_cluster_assignments.tsv` of a pipeline run.


Query: the 234 QC-passing cells of the caerulein-injured mouse pancreas run (`scp/pivot_pack.tsv` of the published pancreas run; the input matrix is `testset/Pancreas-mouse-Spectronaut`;
log2 median-normalised protein intensities; 3,370 unique leading gene symbols; missing = not detected, set to 0).

References (scRNA-seq, mouse, official author annotations):
1. Baron et al. 2016 (GEO GSE84133, mouse1 + mouse2, 1,886 islet cells, 13 classes; NO acinar cells) - `label_transfer.R`
2. Tabula Muris Senis droplet pancreas (figshare 12654728, file 23872580; 6,201 cells, 9 classes incl. 2,141 acinar) - `label_transfer_tms.R`
3. As 2 but reference restricted to the five non-islet classes (acinar, ductal, leukocyte, endothelial, stellate) - `label_transfer_tms_exocrine.R`

Method: Seurat 5 (R 4.5.2). Reference LogNormalize; features = genes shared between reference and proteome and detected in >=10% of
proteome cells (1,666-1,694); ScaleData/RunPCA(30); FindTransferAnchors(dims 1:30, k.filter = NA, k.anchor 5, k.score 30) with
reduction = "cca" and, separately, "pcaproject"; TransferData(k.weight = 20). Per-cell labels and class scores in `results/`.

Headline numbers (cells of 234):
| reference | CCA top calls | PCA-projection top calls | CCA/PCA agreement |
|---|---|---|---|
| Baron islets | beta 131, ductal 32, macrophage 19 | macrophage 116, schwann 59, ductal 33 | 12.8% |
| TMS pancreas (all classes) | B (beta) cell 125, acinar 57, leukocyte 31 | leukocyte 145, ductal 50, stellate 35, acinar 0 | 16.2% |
| TMS non-islet classes only | acinar 153, leukocyte 42, ductal 21, endothelial 16, stellate 2 | leukocyte 140, ductal 70, stellate 20, acinar 0 | 26.5% |

Insulin (Ins1) is detected in 2.6% and glucagon (Gcg) in 6.0% of proteome cells, so the beta-cell majority under references 1-2 is a
reference-abundance artefact. Of the 77 proteome cells in the top tertile of acinar-enzyme intensity (Cpa1, Ctrb1, Cela2a, Prss2, Pnlip, Cel),
CCA/TMS-all labels 37 acinar, 21 beta, 14 leukocyte; PCA-projection labels 57 leukocyte, 17 ductal, 0 acinar.

Still to do once `scp/clustering/scp_cluster_assignments.tsv` of the published D2 run is available: cluster x Seurat-label contingency
table against the CASPA/LLM annotations of C0-C7 (`compare_to_clusters.py`, to be written).
