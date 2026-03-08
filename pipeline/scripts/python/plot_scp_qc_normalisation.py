#!/usr/bin/env python3
"""
plot_scp_qc_normalisation.py

Before/after median normalisation QC plots for SCP pivot.

Inputs
------
--pivot         scp_pivot.tsv (log2, median-normalised)
--shifts        scp_pivot_shifts.tsv (RunStem, SampleColumn, median_shift_applied)
--annotation    scp_annotation.tsv (Run, Condition)
--outdir        output directory

Outputs
-------
plots/
  normalisation_violin.pdf   -- per-cell median: pre vs post, grouped by cluster
  normalisation_shifts.pdf   -- shift magnitude per cell, coloured by cluster
  normalisation_scatter.pdf  -- pre-norm median vs post-norm median scatter
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
    ap.add_argument("--pivot", required=True, help="scp_pivot.tsv (log2, normalised)")
    ap.add_argument("--shifts", required=True, help="scp_pivot_shifts.tsv")
    ap.add_argument("--annotation", required=True, help="scp_annotation.tsv")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--format", default="pdf", choices=["pdf", "png"])
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    outdir = Path(args.outdir) / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Load data ----
    ann = pd.read_csv(args.annotation, sep="\t", dtype=str)
    run_to_cond = dict(zip(ann["Run"].astype(str), ann["Condition"].astype(str)))

    shifts_df = pd.read_csv(args.shifts, sep="\t", dtype=str)
    shifts_df["median_shift_applied"] = pd.to_numeric(shifts_df["median_shift_applied"], errors="coerce")
    run_to_shift = dict(zip(shifts_df["RunStem"].astype(str),
                            shifts_df["median_shift_applied"].tolist()))

    pivot = pd.read_csv(args.pivot, sep="\t", dtype=str)
    all_sample_cols = [c for c in pivot.columns if is_sample_column(c)]
    header_map = {c: stem_from_header(c) for c in all_sample_cols}

    # Align columns to annotation order
    wanted = set(ann["Run"].astype(str))
    matched = [c for c in all_sample_cols if header_map.get(c, "") in wanted]
    col_to_stem = {c: header_map[c] for c in matched}
    run_to_col = {v: k for k, v in col_to_stem.items()}

    runs = [r for r in ann["Run"].astype(str).tolist() if r in run_to_col]
    cols = [run_to_col[r] for r in runs]
    conds = [run_to_cond.get(r, "unknown") for r in runs]
    shifts = np.array([run_to_shift.get(r, 0.0) for r in runs], dtype=float)

    X_norm = pivot[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    # Recover pre-norm: X_pre = X_norm - shift (shift was added to bring to global median)
    X_pre = X_norm - shifts[None, :]

    # Per-cell medians
    med_post = np.nanmedian(X_norm, axis=0)
    med_pre = np.nanmedian(X_pre, axis=0)

    clusters = sorted(set(conds))
    n_clusters = len(clusters)
    cmap = plt.get_cmap("tab10") if n_clusters <= 10 else plt.get_cmap("tab20")
    cluster_colors = {c: cmap(i % cmap.N) for i, c in enumerate(clusters)}
    cell_colors = [cluster_colors[c] for c in conds]

    # ---- Plot 1: violin of per-cell medians before vs after ----
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=False)

    for ax, (vals, title) in zip(axes, [(med_pre, "Pre-normalisation"), (med_post, "Post-normalisation")]):
        by_cluster = {c: vals[[i for i, cond in enumerate(conds) if cond == c]] for c in clusters}
        # Filter out clusters with no finite values (avoids zero-size array error)
        plot_clusters = [c for c in clusters if len(by_cluster[c]) > 0 and np.any(np.isfinite(by_cluster[c]))]
        if not plot_clusters:
            ax.set_title(f"{title}\n(no data)")
            continue
        positions = list(range(len(plot_clusters)))
        vp = ax.violinplot([by_cluster[c][np.isfinite(by_cluster[c])] for c in plot_clusters],
                           positions=positions, showmedians=True, showextrema=True)
        for i, (body, cl) in enumerate(zip(vp["bodies"], plot_clusters)):
            body.set_facecolor(cluster_colors[cl])
            body.set_alpha(0.7)
        vp["cmedians"].set_color("black")
        vp["cmedians"].set_linewidth(1.5)
        for j, cl in enumerate(plot_clusters):
            ys = by_cluster[cl]
            xs = np.random.normal(j, 0.07, size=len(ys))
            ax.scatter(xs, ys, s=10, alpha=0.6, c=[cluster_colors[cl]], zorder=3, linewidths=0)
        ax.set_xticks(positions)
        ax.set_xticklabels(plot_clusters, fontsize=9, rotation=45, ha="right")
        ax.set_xlabel("Cluster", fontsize=10)
        ax.set_ylabel("Per-cell log2 median intensity", fontsize=10)
        ax.set_title(title, fontsize=11)
        ax.grid(axis="y", alpha=0.3, linewidth=0.5)

    fig.tight_layout(pad=1.5)
    fig.savefig(outdir / f"normalisation_violin.{args.format}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote normalisation_violin.{args.format}")

    # ---- Plot 2: shift per cell, sorted and coloured by cluster ----
    fig, ax = plt.subplots(figsize=(max(8, len(runs) * 0.18), 4))

    order = np.argsort(shifts)
    xs = np.arange(len(order))
    bar_colors = [cell_colors[i] for i in order]
    ax.bar(xs, shifts[order], color=bar_colors, width=0.8, alpha=0.85, linewidth=0)
    ax.axhline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Cell (sorted by shift)", fontsize=10)
    ax.set_ylabel("Shift applied (log2)", fontsize=10)
    ax.set_title("Median normalisation shifts per cell", fontsize=11)
    ax.tick_params(bottom=False, labelbottom=False)
    legend_patches = [mpatches.Patch(color=cluster_colors[c], label=c) for c in clusters]
    ax.legend(handles=legend_patches, fontsize=8, loc="upper left", framealpha=0.8)
    ax.grid(axis="y", alpha=0.3, linewidth=0.5)

    fig.tight_layout()
    fig.savefig(outdir / f"normalisation_shifts.{args.format}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote normalisation_shifts.{args.format}")

    # ---- Plot 3: pre vs post scatter ----
    fig, ax = plt.subplots(figsize=(6, 6))
    finite_pre = med_pre[np.isfinite(med_pre)]
    finite_post = med_post[np.isfinite(med_post)]
    if len(finite_pre) == 0 or len(finite_post) == 0:
        ax.set_title("Effect of median normalisation\n(no data)")
        fig.tight_layout()
        fig.savefig(outdir / f"normalisation_scatter.{args.format}", bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote normalisation_scatter.{args.format} (empty)")
        print(f"Normalisation QC complete -> {outdir}")
        return
    sc = ax.scatter(med_pre, med_post, c=cell_colors, s=25, alpha=0.85, linewidths=0)
    lims = [min(np.nanmin(finite_pre), np.nanmin(finite_post)) - 0.2,
            max(np.nanmax(finite_pre), np.nanmax(finite_post)) + 0.2]
    ax.plot(lims, lims, "k--", linewidth=0.8, alpha=0.5, label="y = x")
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Pre-norm per-cell median (log2)", fontsize=10)
    ax.set_ylabel("Post-norm per-cell median (log2)", fontsize=10)
    ax.set_title("Effect of median normalisation", fontsize=11)
    ax.set_aspect("equal")
    legend_patches = [mpatches.Patch(color=cluster_colors[c], label=c) for c in clusters]
    ax.legend(handles=legend_patches, fontsize=8, framealpha=0.8)

    fig.tight_layout()
    fig.savefig(outdir / f"normalisation_scatter.{args.format}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote normalisation_scatter.{args.format}")

    print(f"Normalisation QC complete -> {outdir}")


if __name__ == "__main__":
    main()
