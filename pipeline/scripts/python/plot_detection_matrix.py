#!/usr/bin/env python3
"""
plot_detection_matrix.py

Binary detection (presence/absence) heatmap for SCP data.

Proteins ordered by informativeness (cluster-specificity of detection).
Cells ordered by cluster. Uses a white/navy binary colormap.

Inputs
------
--pivot              scp_pivot.tsv (log2, NA for missing)
--annotation         scp_annotation.tsv (Run, Condition)
--detection-markers  detection_markers.tsv (optional, for protein ordering)
--outdir             output directory
--top-n-proteins     max proteins to display (default: 200)
--format             pdf or png
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


def cluster_specificity_score(det: np.ndarray, cond_arr: np.ndarray,
                               clusters: list[str]) -> np.ndarray:
    """
    For each protein, compute max(det_rate_in - det_rate_out) across clusters.
    Higher = more cluster-specific.
    """
    n_prot, n_cells = det.shape
    scores = np.zeros(n_prot, dtype=float)
    for cl in clusters:
        mask_in = (cond_arr == cl)
        mask_out = ~mask_in
        if mask_in.sum() == 0 or mask_out.sum() == 0:
            continue
        rate_in = det[:, mask_in].mean(axis=1)
        rate_out = det[:, mask_out].mean(axis=1)
        diff = rate_in - rate_out
        scores = np.maximum(scores, diff)
    return scores


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pivot", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--detection-markers", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--top-n-proteins", type=int, default=200)
    ap.add_argument("--min-det-fraction", type=float, default=0.05,
                    help="Min fraction of cells where protein must be detected to be included.")
    ap.add_argument("--format", default="pdf", choices=["pdf", "png"])
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.patches as mpatches

    outdir = Path(args.outdir) / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Load annotation ----
    ann = pd.read_csv(args.annotation, sep="\t", dtype=str)
    run_to_cond = dict(zip(ann["Run"].astype(str), ann["Condition"].astype(str)))
    clusters = sorted(set(run_to_cond.values()))

    # ---- Load pivot ----
    pivot = pd.read_csv(args.pivot, sep="\t", dtype=str)
    sample_cols = [c for c in pivot.columns if is_sample_column(c)]
    header_map = {c: stem_from_header(c) for c in sample_cols}
    run_set = set(ann["Run"].astype(str))
    matched = [c for c in sample_cols if header_map.get(c, "") in run_set]
    run_to_col = {header_map[c]: c for c in matched}

    # Order cells by cluster
    runs_by_cluster = []
    for cl in clusters:
        for r in ann["Run"].astype(str).tolist():
            if run_to_cond.get(r) == cl and r in run_to_col:
                runs_by_cluster.append(r)

    cols_ordered = [run_to_col[r] for r in runs_by_cluster]
    cond_arr = np.array([run_to_cond[r] for r in runs_by_cluster])

    X = pivot[cols_ordered].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    det = np.isfinite(X).astype(np.float32)  # 0/1 binary

    n_prot, n_cells = det.shape

    # Filter: minimum detection fraction
    det_frac = det.mean(axis=1)
    keep_mask = det_frac >= args.min_det_fraction
    det = det[keep_mask, :]
    prot_col = "Protein.Group" if "Protein.Group" in pivot.columns else pivot.columns[0]
    gene_col = "Genes" if "Genes" in pivot.columns else None
    labels_all = (pivot[gene_col].fillna("").astype(str) + "|" + pivot[prot_col].astype(str)
                  if gene_col else pivot[prot_col].astype(str)).tolist()
    labels = [labels_all[i] for i in np.where(keep_mask)[0]]

    print(f"Proteins passing min det fraction ({args.min_det_fraction:.2f}): {det.shape[0]} / {n_prot}")

    if det.shape[0] == 0 or n_cells == 0:
        print("No proteins or cells pass filters — writing empty sentinel and exiting.")
        outdir.mkdir(parents=True, exist_ok=True)
        return

    # Score and select top N proteins
    if args.detection_markers and Path(args.detection_markers).exists():
        dm = pd.read_csv(args.detection_markers, sep="\t", dtype=str)
        dm["log2_odds_ratio"] = pd.to_numeric(dm.get("log2_odds_ratio", pd.Series()), errors="coerce")
        dm["qvalue"] = pd.to_numeric(dm.get("qvalue", pd.Series()), errors="coerce")
        id_col = "Protein" if "Protein" in dm.columns else dm.columns[0]
        # Score: log2OR * significance per cluster, take max across clusters
        dm["score"] = dm["log2_odds_ratio"].abs() * (-np.log10(dm["qvalue"].clip(1e-300)))
        best_per_prot = dm.groupby(id_col)["score"].max()
        scores = np.array([float(best_per_prot.get(l.split("|")[-1], 0.0)) for l in labels])
    else:
        scores = cluster_specificity_score(det, cond_arr, clusters)

    top_n = min(args.top_n_proteins, det.shape[0])
    top_idx = np.argsort(scores)[::-1][:top_n]
    det_top = det[top_idx, :]
    labels_top = [labels[i] for i in top_idx]

    # Short labels (gene name only)
    def short_label(l):
        gene = l.split("|")[0].strip()
        if gene and gene not in ("", "nan"):
            return gene[:20]
        return l.split("|")[-1][:15]

    labels_short = [short_label(l) for l in labels_top]

    # ---- Figure ----
    cmap_tab = plt.get_cmap("tab10") if len(clusters) <= 10 else plt.get_cmap("tab20")
    cluster_colors = {c: cmap_tab(i % cmap_tab.N) for i, c in enumerate(clusters)}

    # Adaptive figure size
    cell_width = max(0.04, min(0.12, 8.0 / n_cells))
    prot_height = max(0.06, min(0.18, 10.0 / top_n))
    fig_w = max(8, n_cells * cell_width + 2.0)
    fig_h = max(6, top_n * prot_height + 1.5)

    fig = plt.figure(figsize=(fig_w, fig_h))
    # Axes layout: annotation bar on top, heatmap below
    ax_bar = fig.add_axes([0.1, 0.92, 0.8, 0.04])
    ax_hm = fig.add_axes([0.1, 0.08, 0.8, 0.82])

    # Cluster colour bar
    bar_arr = np.array([[cluster_colors[c] for c in cond_arr]])  # 1 x n_cells x 4 RGBA
    ax_bar.imshow(bar_arr, aspect="auto", interpolation="none")
    ax_bar.set_yticks([])
    ax_bar.set_xticks([])
    ax_bar.set_title("Detection matrix (presence/absence)", fontsize=10, pad=4)

    # Add cluster labels to bar
    boundaries = [0]
    prev = cond_arr[0]
    for i, c in enumerate(cond_arr[1:], 1):
        if c != prev:
            boundaries.append(i)
            prev = c
    boundaries.append(n_cells)
    for k in range(len(boundaries) - 1):
        mid = (boundaries[k] + boundaries[k + 1]) / 2.0
        cl = cond_arr[boundaries[k]]
        ax_bar.text(mid, 0, cl, ha="center", va="center", fontsize=6,
                    color="white", fontweight="bold")

    # Heatmap
    binary_cmap = mcolors.LinearSegmentedColormap.from_list("det", ["#FFFFFF", "#1A3A6C"])
    ax_hm.imshow(det_top, aspect="auto", interpolation="none",
                 cmap=binary_cmap, vmin=0, vmax=1, origin="upper")

    # Y-axis labels: show if few enough proteins
    if top_n <= 80:
        ax_hm.set_yticks(range(top_n))
        ax_hm.set_yticklabels(labels_short, fontsize=max(4, min(8, 600 // top_n)))
    else:
        ax_hm.set_yticks([])

    ax_hm.set_xticks([])
    ax_hm.set_xlabel(f"Cells (n={n_cells}, ordered by cluster)", fontsize=9)
    ax_hm.set_ylabel(f"Proteins (top {top_n} by cluster-specificity)", fontsize=9)

    # Add vertical separators between clusters
    for b in boundaries[1:-1]:
        ax_hm.axvline(b - 0.5, color="#888888", linewidth=0.8, alpha=0.6)
        ax_bar.axvline(b - 0.5, color="#ffffff", linewidth=0.8, alpha=0.8)

    # Legend
    patches = [mpatches.Patch(color=cluster_colors[c], label=c) for c in clusters]
    fig.legend(handles=patches, loc="lower right", fontsize=7, framealpha=0.8,
               bbox_to_anchor=(0.98, 0.01), ncol=max(1, len(clusters) // 5))

    out_path = outdir / f"detection_matrix.{args.format}"
    fig.savefig(out_path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
