#!/usr/bin/env python3
"""
plot_scp_marker_volcanos.py

Volcano plots for all three SCP marker modalities.

  Detection markers:  x = log2_odds_ratio,     y = -log10(qvalue)
  Intensity markers:  x = log2FC_detected_only, y = -log10(qvalue)
  Scplainer markers:  x = log2FC,               y = -log10(qvalue)

Each modality × contrast gets one panel.  Panels are collected into a
multi-page PDF (one page per modality).

Inputs
------
--detection-markers   detection_markers.tsv   (optional)
--intensity-markers   intensity_markers_detected_only.tsv (optional)
--scplainer-markers   scplainer_intensity_markers.tsv (optional)
--outdir              output directory
--q-threshold         FDR threshold (default 0.05)
--fc-threshold        |FC| threshold for detection (default 1.0)
--top-n-labels        gene labels to show (default 8)
--format              pdf or png (default pdf)
"""

import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


def _neg_log10(arr):
    arr = np.asarray(arr, dtype=float)
    return -np.log10(np.clip(arr, 1e-300, 1.0))


def volcano_panel(ax, fc: np.ndarray, pval: np.ndarray, labels: list[str],
                  cluster: str, title: str,
                  q_thresh: float, fc_thresh: float, top_n: int,
                  color_up="#B2182B", color_down="#2166AC", color_ns="#999999"):
    """Draw a single volcano panel on ax."""
    y = _neg_log10(pval)
    finite = np.isfinite(fc) & np.isfinite(y)

    sig_up   = finite & (pval <= q_thresh) & (fc >  fc_thresh)
    sig_down = finite & (pval <= q_thresh) & (fc < -fc_thresh)
    ns       = finite & ~(sig_up | sig_down)

    ax.scatter(fc[ns], y[ns],   s=6,  alpha=0.4, c=color_ns,   linewidths=0, zorder=1)
    ax.scatter(fc[sig_down], y[sig_down], s=10, alpha=0.75, c=color_down, linewidths=0, zorder=2)
    ax.scatter(fc[sig_up],   y[sig_up],   s=10, alpha=0.75, c=color_up,   linewidths=0, zorder=2)

    # Threshold lines
    y_thresh = _neg_log10(np.array([q_thresh]))[0]
    xlim_ext = max(1.0, np.nanmax(np.abs(fc[finite])) * 1.05) if finite.any() else 2.0
    ax.axhline(y_thresh, color="black", linewidth=0.7, linestyle="--", alpha=0.6)
    ax.axvline( fc_thresh, color="black", linewidth=0.7, linestyle="--", alpha=0.4)
    ax.axvline(-fc_thresh, color="black", linewidth=0.7, linestyle="--", alpha=0.4)
    ax.axvline(0, color="black", linewidth=0.5, alpha=0.3)

    # Labels for top hits
    if top_n > 0 and (sig_up.any() or sig_down.any()):
        sig_mask = sig_up | sig_down
        sig_score = y.copy()
        sig_score[~sig_mask] = -np.inf
        top_idx = np.argsort(sig_score)[::-1][:top_n]
        try:
            from adjustText import adjust_text
            texts = []
            for i in top_idx:
                if sig_mask[i]:
                    texts.append(ax.text(fc[i], y[i], labels[i], fontsize=6,
                                         fontstyle="italic", ha="center"))
            adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))
        except ImportError:
            for i in top_idx:
                if sig_mask[i]:
                    ax.text(fc[i], y[i], labels[i], fontsize=6,
                            fontstyle="italic", ha="center", va="bottom")

    n_up   = int(sig_up.sum())
    n_down = int(sig_down.sum())
    ax.text(0.97, 0.97, f"↑{n_up}  ↓{n_down}", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, color="black")

    ax.set_title(f"{title}\n{cluster}", fontsize=8)
    ax.set_xlabel("log₂FC (or log₂ odds ratio)", fontsize=7)
    ax.set_ylabel("−log₁₀(q-value)", fontsize=7)
    ax.tick_params(labelsize=6)
    ax.set_xlim(-xlim_ext, xlim_ext)
    y_max = np.nanmax(y[finite]) * 1.1 if finite.any() else 5.0
    ax.set_ylim(-0.1, max(y_max, y_thresh * 1.5))


def plot_modality(df: pd.DataFrame, fc_col: str, pval_col: str, label_col: str,
                  cluster_col: str, modality_title: str,
                  q_thresh: float, fc_thresh: float, top_n: int,
                  out_path: Path, fmt: str):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    clusters = sorted(df[cluster_col].unique())
    n_cl = len(clusters)
    if n_cl == 0:
        print(f"  No clusters in data for {out_path.name} — skipping.")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.touch()
        return
    ncols = min(3, n_cl)
    nrows = int(np.ceil(n_cl / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows), squeeze=False)

    for k, cl in enumerate(clusters):
        r, c = divmod(k, ncols)
        sub = df[df[cluster_col] == cl].copy()
        fc   = pd.to_numeric(sub[fc_col], errors="coerce").to_numpy()
        pval = pd.to_numeric(sub[pval_col], errors="coerce").to_numpy()
        lbls = sub[label_col].astype(str).tolist()
        volcano_panel(axes[r][c], fc, pval, lbls, cl, modality_title,
                      q_thresh, fc_thresh, top_n)

    for k in range(n_cl, nrows * ncols):
        r, c = divmod(k, ncols)
        axes[r][c].set_visible(False)

    fig.tight_layout(pad=0.8)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--detection-markers",  default=None)
    ap.add_argument("--intensity-markers",  default=None)
    ap.add_argument("--scplainer-markers",  default=None)
    ap.add_argument("--outdir",             required=True)
    ap.add_argument("--q-threshold",        type=float, default=0.05)
    ap.add_argument("--fc-threshold",       type=float, default=0.5)
    ap.add_argument("--top-n-labels",       type=int,   default=8)
    ap.add_argument("--format",             default="pdf", choices=["pdf", "png"])
    args = ap.parse_args()

    outdir = Path(args.outdir) / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    def gene_label(df: pd.DataFrame) -> list[str]:
        if "Genes" in df.columns:
            return df["Genes"].fillna("").astype(str).tolist()
        if "Gene" in df.columns:
            return df["Gene"].fillna("").astype(str).tolist()
        return df.iloc[:, 0].astype(str).tolist()

    # ---- Detection markers ----
    if args.detection_markers and Path(args.detection_markers).exists():
        dm = pd.read_csv(args.detection_markers, sep="\t", dtype=str)
        dm["label"] = gene_label(dm)
        if "Cluster" in dm.columns and "log2_odds_ratio" in dm.columns and "qvalue" in dm.columns:
            plot_modality(dm, "log2_odds_ratio", "qvalue", "label", "Cluster",
                          "Detection markers\n(log₂ odds ratio)",
                          args.q_threshold, args.fc_threshold, args.top_n_labels,
                          outdir / f"volcano_detection_markers.{args.format}", args.format)

    # ---- Intensity markers ----
    if args.intensity_markers and Path(args.intensity_markers).exists():
        im = pd.read_csv(args.intensity_markers, sep="\t", dtype=str)
        im["label"] = gene_label(im)
        if "Cluster" in im.columns and "log2FC_detected_only" in im.columns and "qvalue" in im.columns:
            plot_modality(im, "log2FC_detected_only", "qvalue", "label", "Cluster",
                          "Intensity markers\n(detected-only log₂FC)",
                          args.q_threshold, args.fc_threshold, args.top_n_labels,
                          outdir / f"volcano_intensity_markers.{args.format}", args.format)

    # ---- Scplainer markers ----
    if args.scplainer_markers and Path(args.scplainer_markers).exists():
        sm = pd.read_csv(args.scplainer_markers, sep="\t", dtype=str)
        sm["label"] = gene_label(sm)
        # scplainer markers use Cluster column; filter ok==TRUE if present
        if "ok" in sm.columns:
            sm = sm[sm["ok"].astype(str).str.upper().isin(["TRUE", "1"])].copy()
        cluster_col = "Cluster" if "Cluster" in sm.columns else "contrast"
        if cluster_col in sm.columns and "log2FC" in sm.columns and "qvalue" in sm.columns:
            plot_modality(sm, "log2FC", "qvalue", "label", cluster_col,
                          "scplainer markers\n(model-based log₂FC)",
                          args.q_threshold, args.fc_threshold, args.top_n_labels,
                          outdir / f"volcano_scplainer_markers.{args.format}", args.format)

    print(f"Marker volcanos complete -> {outdir}")


if __name__ == "__main__":
    main()
