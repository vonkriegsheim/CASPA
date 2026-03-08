#!/usr/bin/env python3
"""
plot_cluster_centroid_similarity.py

Compute pairwise Pearson correlation between cluster-mean expression profiles
and render as a heatmap with dendrogram. Helps assess whether clustering
resolution is appropriate and which clusters are biologically related.

Inputs
------
--pivot             scp_pivot.tsv (log2, NA for missing)
--annotation        scp_annotation.tsv (Run, Condition)
--outdir            output directory
--top-n-variable    proteins to include in centroid (default 500, by MAD)
--format            pdf or png
"""

import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

KNOWN_META_COLS = {
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences",
    "Protein.Q.Value", "Global.Q.Value", "Protein.Descriptions",
}


def stem_from_header(h: str) -> str:
    s = str(h)
    p = Path(s)
    if "\\" in s or "/" in s:
        return p.stem
    if "." in s and s.rsplit(".", 1)[1].lower() in ("raw", "d", "mzml", "wiff"):
        return p.stem
    return s


def is_sample_column(col_name: str) -> bool:
    if col_name in KNOWN_META_COLS:
        return False
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pivot", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-n-variable", type=int, default=500)
    ap.add_argument("--format", default="pdf", choices=["pdf", "png"])
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
    from scipy.spatial.distance import squareform

    outdir = Path(args.outdir) / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    ann = pd.read_csv(args.annotation, sep="\t", dtype=str)
    run_to_cond = dict(zip(ann["Run"].astype(str), ann["Condition"].astype(str)))
    clusters = sorted(set(run_to_cond.values()))

    pivot = pd.read_csv(args.pivot, sep="\t", dtype=str)
    sample_cols = [c for c in pivot.columns if is_sample_column(c)]
    header_map = {c: stem_from_header(c) for c in sample_cols}
    run_set = set(ann["Run"].astype(str))
    matched = [c for c in sample_cols if header_map.get(c, "") in run_set]
    run_to_col = {header_map[c]: c for c in matched}

    runs = [r for r in ann["Run"].astype(str).tolist() if r in run_to_col]
    cols = [run_to_col[r] for r in runs]
    cond_arr = np.array([run_to_cond[r] for r in runs])

    X = pivot[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    mads = np.nanmedian(np.abs(X - np.nanmedian(X, axis=1, keepdims=True)), axis=1)
    mads = np.where(np.isfinite(mads), mads, 0.0)
    n_use = min(args.top_n_variable, X.shape[0])
    top_prot_idx = np.argsort(mads)[::-1][:n_use]
    X_top = X[top_prot_idx, :]

    n_cl = len(clusters)
    centroids = np.full((n_use, n_cl), np.nan)
    for j, cl in enumerate(clusters):
        mask = (cond_arr == cl)
        if mask.sum() == 0:
            continue
        centroids[:, j] = np.nanmean(X_top[:, mask], axis=1)

    keep = np.isfinite(centroids).any(axis=1)
    centroids = centroids[keep, :]

    # Pairwise Pearson correlation
    R = np.full((n_cl, n_cl), 0.0)
    for i in range(n_cl):
        for j in range(n_cl):
            ci = centroids[:, i]
            cj = centroids[:, j]
            both = np.isfinite(ci) & np.isfinite(cj)
            if both.sum() >= 10:
                R[i, j] = np.corrcoef(ci[both], cj[both])[0, 1]

    # Force exact symmetry — floating-point can cause tiny asymmetries that
    # break scipy squareform's symmetry validation
    R = (R + R.T) / 2.0
    np.fill_diagonal(R, 1.0)

    tables_dir = Path(args.outdir) / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(R, index=clusters, columns=clusters).to_csv(
        tables_dir / "cluster_centroid_correlation.tsv", sep="\t")

    dist_mat = 1.0 - R
    np.fill_diagonal(dist_mat, 0.0)
    dist_mat = np.clip(dist_mat, 0.0, None)
    dist_mat = (dist_mat + dist_mat.T) / 2.0  # second symmetrisation after clip
    np.fill_diagonal(dist_mat, 0.0)

    if n_cl >= 3:
        Z = linkage(squareform(dist_mat), method="average")
        order = leaves_list(Z)
    else:
        order = list(range(n_cl))

    R_ord = R[np.ix_(order, order)]
    labels_ord = [clusters[i] for i in order]

    figsize = max(5, n_cl * 0.7 + 2)
    fig, ax = plt.subplots(figsize=(figsize, figsize * 0.85))
    im = ax.imshow(R_ord, cmap="RdBu_r", vmin=-1, vmax=1, aspect="equal", interpolation="none")

    ax.set_xticks(range(n_cl))
    ax.set_xticklabels(labels_ord, fontsize=10, rotation=45, ha="right")
    ax.set_yticks(range(n_cl))
    ax.set_yticklabels(labels_ord, fontsize=10)

    for i in range(n_cl):
        for j in range(n_cl):
            v = R_ord[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=max(7, min(11, 80 // n_cl)),
                    color="white" if abs(v) > 0.6 else "black")

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson r (cluster centroids)", fontsize=9)
    ax.set_title(f"Cluster centroid similarity\n(top {n_use} variable proteins)", fontsize=10)

    fig.tight_layout()
    out_path = outdir / f"cluster_centroid_similarity.{args.format}"
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")

    if n_cl >= 3:
        fig2, ax2 = plt.subplots(figsize=(max(4, n_cl * 0.6), 4))
        dendrogram(Z, labels=[clusters[i] for i in range(n_cl)],
                   ax=ax2, color_threshold=0.7 * np.max(Z[:, 2]))
        ax2.set_xlabel("Cluster", fontsize=10)
        ax2.set_ylabel("1 - Pearson r", fontsize=10)
        ax2.set_title("Cluster similarity dendrogram", fontsize=10)
        fig2.tight_layout()
        fig2.savefig(outdir / f"cluster_similarity_dendrogram.{args.format}", bbox_inches="tight")
        plt.close(fig2)
        print(f"Wrote cluster_similarity_dendrogram.{args.format}")

    print(f"Cluster centroid similarity complete -> {outdir}")


if __name__ == "__main__":
    main()
