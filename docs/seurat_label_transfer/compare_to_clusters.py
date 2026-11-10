"""Cluster x Seurat-label contingency for the CASPA pancreas run.

Usage: python compare_to_clusters.py <scp_cluster_assignments.tsv of the published D2 run>
Writes results/cluster_by_seurat_label_*.tsv and prints the tables.
"""
import sys, os
sys.stdout.reconfigure(encoding="utf-8")
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
CASPA_LABELS = {  # Fig. 6d, Claude Sonnet 4.6 condition B (abbreviated as in regenerate_fig6_pancreas.py)
    "C0": "Injured acinar", "C1": "Macrophage (Mrc1+)", "C2": "Stellate / myeloid (low depth)",
    "C3": "Lymphocyte / low-complexity immune", "C4": "Stellate (activated)", "C5": "Myofibroblast",
    "C6": "Myeloid", "C7": "Ductal"}

ca = pd.read_csv(sys.argv[1], sep="\t", dtype=str)
cell_col = next(c for c in ca.columns if c.lower() in ("cell", "sample", "run", "cell_id", "sample_id", "column"))
clu_col = next(c for c in ca.columns if "cluster" in c.lower())
ca = ca[[cell_col, clu_col]].rename(columns={cell_col: "cell", clu_col: "cluster"})
ca["cluster"] = ca["cluster"].map(lambda x: x if str(x).startswith("C") else f"C{x}")
print("clusters:", ca.cluster.value_counts().sort_index().to_dict())

for tag in ("seurat_TMS", "seurat_TMSexocrine", "seurat"):
    pred = pd.read_csv(os.path.join(RES, f"{tag}_label_transfer_per_cell.tsv"), sep="\t")
    m = pred.merge(ca, on="cell", how="left")
    print(f"\n[{tag}] cells matched to a cluster: {m.cluster.notna().sum()} of {len(m)}")
    for method in ("cca", "pcaproject"):
        col = next(c for c in m.columns if c.endswith(f"_{method}_label"))
        ct = pd.crosstab(m.cluster, m[col])
        ct.insert(0, "CASPA annotation", [CASPA_LABELS.get(c, "") for c in ct.index])
        ct.to_csv(os.path.join(RES, f"cluster_by_{tag}_{method}.tsv"), sep="\t")
        print(f"\n-- {tag} / {method} --"); print(ct.to_string())
