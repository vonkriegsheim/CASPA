#!/usr/bin/env python3
"""
plot_harmony_diagnostics.py

Before/after Harmony batch correction diagnostic plots.

Requires scp_cluster_assignments.tsv to contain:
  umap_1, umap_2                      -- post-correction UMAP (always present)
  umap_1_uncorrected, umap_2_uncorrected  -- pre-correction UMAP (present when batch_applied=True)
  Library                             -- batch/library label
  cluster_id                          -- Leiden cluster

Also reads scp_clustering_report.json to annotate whether Harmony was applied.

Inputs
------
--assignments       scp_cluster_assignments.tsv
--report            scp_clustering_report.json  (optional, for metadata)
--outdir            output directory
--format            pdf or png

Outputs
-------
plots/
  harmony_before_after.pdf   -- 2x2 grid: pre/post UMAP × coloured by Library/Cluster
  harmony_batch_summary.pdf  -- bar chart of cell counts per library × cluster
"""

import argparse
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


def scatter_panel(ax, x, y, labels, title, cmap_name="tab10", s=20.0, show_legend=True):
    """Categorical scatter on ax."""
    import matplotlib.pyplot as plt
    cats = sorted(set(str(l) for l in labels if pd.notna(l)))
    if not cats:
        ax.set_title(title, fontsize=9)
        return
    cmap = plt.get_cmap(cmap_name) if len(cats) <= 10 else plt.get_cmap("tab20")
    colors = {c: cmap(i % cmap.N) for i, c in enumerate(cats)}
    for cat in cats:
        mask = np.array([str(l) == cat for l in labels])
        ax.scatter(x[mask], y[mask], c=[colors[cat]], s=s, alpha=0.85,
                   linewidths=0, label=cat, zorder=2)
    if show_legend:
        ax.legend(fontsize=7, markerscale=1.0, framealpha=0.7,
                  loc="best", ncol=max(1, len(cats) // 8))
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("UMAP1", fontsize=8)
    ax.set_ylabel("UMAP2", fontsize=8)
    ax.set_aspect("equal", adjustable="datalim")
    ax.tick_params(labelsize=7)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignments", required=True)
    ap.add_argument("--report",      default=None,  help="scp_clustering_report.json")
    ap.add_argument("--outdir",      required=True)
    ap.add_argument("--format",      default="pdf", choices=["pdf", "png"])
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = Path(args.outdir) / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Load assignments ----
    df = pd.read_csv(args.assignments, sep="\t", dtype=str)
    for col in ["umap_1", "umap_2", "umap_1_uncorrected", "umap_2_uncorrected"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    x_post = df["umap_1"].to_numpy()
    y_post = df["umap_2"].to_numpy()
    has_pre = "umap_1_uncorrected" in df.columns and "umap_2_uncorrected" in df.columns
    x_pre = df["umap_1_uncorrected"].to_numpy() if has_pre else x_post
    y_pre = df["umap_2_uncorrected"].to_numpy() if has_pre else y_post

    library = df.get("Library", pd.Series(["unknown"] * len(df))).fillna("unknown").tolist()
    cluster = df.get("cluster_id", pd.Series(["unknown"] * len(df))).fillna("unknown").tolist()

    n_cells = len(df)
    s = max(6.0, min(40.0, 40.0 * np.sqrt(12.0 / n_cells)))

    # ---- Report metadata ----
    batch_applied = False
    n_libraries   = len(set(library))
    if args.report and Path(args.report).exists():
        try:
            rpt = json.loads(Path(args.report).read_text(encoding="utf-8"))
            batch_applied = bool(rpt.get("batch_applied", False))
            n_libraries   = int(rpt.get("n_libraries", n_libraries))
        except Exception:
            pass

    # ---- Plot 1: 2×2 before/after grid ----
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))

    if has_pre and batch_applied:
        scatter_panel(axes[0][0], x_pre,  y_pre,  library, "Before Harmony — by Library",   s=s)
        scatter_panel(axes[0][1], x_pre,  y_pre,  cluster, "Before Harmony — by Cluster",    s=s)
        scatter_panel(axes[1][0], x_post, y_post, library, "After Harmony — by Library",     s=s)
        scatter_panel(axes[1][1], x_post, y_post, cluster, "After Harmony — by Cluster",     s=s)
        fig.suptitle(f"Harmony batch correction diagnostics  "
                     f"(n_libraries={n_libraries}, n_cells={n_cells})",
                     fontsize=11, y=1.01)
    else:
        label = "No Harmony (single library)" if n_libraries <= 1 else "Harmony not applied"
        scatter_panel(axes[0][0], x_post, y_post, library, f"{label} — by Library", s=s)
        scatter_panel(axes[0][1], x_post, y_post, cluster, f"{label} — by Cluster", s=s)
        axes[1][0].set_visible(False)
        axes[1][1].set_visible(False)
        fig.suptitle(f"UMAP — batch diagnostics  "
                     f"(n_libraries={n_libraries}, n_cells={n_cells})",
                     fontsize=11, y=1.01)

    fig.tight_layout(pad=1.5)
    out1 = outdir / f"harmony_before_after.{args.format}"
    fig.savefig(out1, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out1}")

    # ---- Plot 2: Cell counts per Library × Cluster (stacked bar) ----
    df2 = pd.DataFrame({"Library": library, "Cluster": cluster})
    if df2["Library"].nunique() > 1:
        cross = (df2.groupby(["Library", "Cluster"])
                    .size()
                    .unstack(fill_value=0))
        n_lib = len(cross)
        n_cl  = cross.shape[1]
        cmap  = plt.get_cmap("tab10") if n_cl <= 10 else plt.get_cmap("tab20")
        colors = [cmap(i % cmap.N) for i in range(n_cl)]

        fig2, ax2 = plt.subplots(figsize=(max(5, n_lib * 0.8), 5))
        bottom = np.zeros(n_lib)
        for k, cl in enumerate(cross.columns):
            vals = cross[cl].to_numpy(dtype=float)
            ax2.bar(range(n_lib), vals, bottom=bottom, color=colors[k],
                    label=cl, width=0.7)
            bottom += vals

        ax2.set_xticks(range(n_lib))
        ax2.set_xticklabels(cross.index.tolist(), rotation=35, ha="right", fontsize=9)
        ax2.set_ylabel("Cell count", fontsize=10)
        ax2.set_title("Cell composition per library / batch", fontsize=10)
        ax2.legend(title="Cluster", fontsize=8, loc="upper right", framealpha=0.8,
                   ncol=max(1, n_cl // 6))
        ax2.grid(axis="y", alpha=0.3, linewidth=0.5)
        fig2.tight_layout()
        out2 = outdir / f"harmony_batch_summary.{args.format}"
        fig2.savefig(out2, bbox_inches="tight")
        plt.close(fig2)
        print(f"Wrote {out2}")
    else:
        print("Only one library detected — batch summary plot skipped.")

    print(f"Harmony diagnostics complete -> {outdir}")


if __name__ == "__main__":
    main()
