#!/usr/bin/env python
"""
plot_scplainer_volcano.py

Volcano plots directly from scplainer_cluster_DA.tsv.

Input schema (minimum):
  protein_group, contrast, log2FC, pvalue, qvalue
Optional:
  gene_primary, genes, ok

Outputs:
  volcano_scplainer_{contrast}.pdf (or png/svg) in outdir

Notes:
  - Uses qvalue as adjusted p-value on y-axis: -log10(qvalue)
  - Points with missing p/q or log2FC are dropped
  - If ok column exists, only ok==True are kept (unless --keep-not-ok)
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="Volcano plots from scplainer_cluster_DA.tsv")
    ap.add_argument("--input", required=True, help="scplainer_cluster_DA.tsv")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--format", default="pdf", choices=["pdf", "png", "svg"])
    ap.add_argument("--dpi", type=int, default=300)

    ap.add_argument("--q-threshold", type=float, default=0.05, help="q-value threshold")
    ap.add_argument("--fc-threshold", type=float, default=0.25, help="abs(log2FC) threshold")
    ap.add_argument("--top-n-labels", type=int, default=10, help="Top N labels per direction")
    ap.add_argument("--keep-not-ok", action="store_true", help="Do not filter ok==False rows")
    ap.add_argument("--label-col", default="auto",
                    help="Label column: auto|gene_primary|genes|protein_group")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t")
    print(f"Loaded {len(df)} rows from {args.input}")

    required = {"contrast", "log2FC", "pvalue", "qvalue", "protein_group"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {sorted(missing)}")

    # Filter to "ok" if available
    if "ok" in df.columns and not args.keep_not_ok:
        # ok may be TRUE/FALSE strings or booleans
        ok = df["ok"]
        if ok.dtype == object:
            ok = ok.astype(str).str.lower().isin(["true", "t", "1", "yes"])
        df = df[ok].copy()
        print(f"Filtered to ok rows: {len(df)} remaining")

    # Pick label column
    if args.label_col != "auto":
        label_col = args.label_col
    else:
        if "gene_primary" in df.columns and df["gene_primary"].notna().any():
            label_col = "gene_primary"
        elif "genes" in df.columns and df["genes"].notna().any():
            label_col = "genes"
        else:
            label_col = "protein_group"
    print(f"Using label column: {label_col}")

    # Clean numeric cols
    df["log2FC"] = pd.to_numeric(df["log2FC"], errors="coerce")
    df["pvalue"] = pd.to_numeric(df["pvalue"], errors="coerce")
    df["qvalue"] = pd.to_numeric(df["qvalue"], errors="coerce")

    # Drop unusable rows
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["log2FC", "pvalue", "qvalue"]).copy()
    print(f"After dropping NA/Inf rows: {len(df)}")

    # y-axis: -log10(qvalue)
    df["neg_log10_q"] = -np.log10(np.clip(df["qvalue"].values, 1e-300, 1.0))

    # significance categories
    sig = (df["qvalue"] <= args.q_threshold) & (df["log2FC"].abs() >= args.fc_threshold)
    up = sig & (df["log2FC"] > 0)
    down = sig & (df["log2FC"] < 0)

    df["significance"] = "NS"
    df.loc[up, "significance"] = "Up"
    df.loc[down, "significance"] = "Down"

    colors = {"Up": "#E41A1C", "Down": "#377EB8", "NS": "#BDBDBD"}

    # Plot per contrast
    for ctr, sub in df.groupby("contrast", sort=False):
        print(f"Plotting {ctr} ({len(sub)} rows)")

        fig, ax = plt.subplots(figsize=(5, 5))

        # background first
        for cat in ["NS", "Down", "Up"]:
            s = sub[sub["significance"] == cat]
            if s.empty:
                continue
            ax.scatter(
                s["log2FC"],
                s["neg_log10_q"],
                c=colors[cat],
                s=12 if cat == "NS" else 16,
                alpha=0.5 if cat == "NS" else 0.85,
                edgecolors="none",
                rasterized=True,
            )

        # threshold lines
        ax.axhline(-np.log10(args.q_threshold), color="#666666", linewidth=0.6, linestyle="--")
        ax.axvline(-args.fc_threshold, color="#666666", linewidth=0.6, linestyle="--")
        ax.axvline(args.fc_threshold, color="#666666", linewidth=0.6, linestyle="--")

        ax.set_xlabel(r"$\log_2$ fold change", fontsize=10)
        ax.set_ylabel(r"$-\log_{10}$ adjusted p-value (q)", fontsize=10)
        ax.tick_params(axis="both", labelsize=8)

        # symmetric x-axis
        x_max = max(abs(sub["log2FC"].min()), abs(sub["log2FC"].max()))
        x_max = max(x_max * 1.1, args.fc_threshold * 1.5)
        ax.set_xlim(-x_max, x_max)

        # labels: top N per direction by qvalue then abs(fc)
        lab = sub[sub["significance"] != "NS"].copy()
        if not lab.empty and args.top_n_labels > 0:
            # rank within up/down separately
            lab["abs_fc"] = lab["log2FC"].abs()
            up_hits = lab[lab["log2FC"] > 0].sort_values(["qvalue", "abs_fc"], ascending=[True, False]).head(args.top_n_labels)
            dn_hits = lab[lab["log2FC"] < 0].sort_values(["qvalue", "abs_fc"], ascending=[True, False]).head(args.top_n_labels)
            hits = pd.concat([up_hits, dn_hits], axis=0)

            texts = []
            for _, row in hits.iterrows():
                lbl = str(row.get(label_col, "")).strip()
                if not lbl or lbl.lower() == "nan":
                    continue
                t = ax.annotate(
                    lbl,
                    xy=(row["log2FC"], row["neg_log10_q"]),
                    fontsize=7,
                    fontstyle="italic",
                    ha="center",
                    va="bottom",
                    zorder=5,
                )
                texts.append(t)

            # optional: adjustText
            try:
                from adjustText import adjust_text
                if texts:
                    adjust_text(
                        texts,
                        ax=ax,
                        arrowprops=dict(arrowstyle="-", color="#444444", linewidth=0.4),
                        expand_points=(1.4, 1.4),
                        force_text=(0.5, 0.7),
                    )
            except ImportError:
                pass

        # spines
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(0.6)
        ax.spines["bottom"].set_linewidth(0.6)

        safe = str(ctr).replace("/", "_").replace("\\", "_").replace(" ", "_")
        out_path = outdir / f"volcano_scplainer_{safe}.{args.format}"
        fig.savefig(out_path, dpi=args.dpi, bbox_inches="tight", pad_inches=0.15, facecolor="white")
        plt.close(fig)
        print(f"  Saved: {out_path}")

    print(f"Done. Volcano plots saved to {outdir}")


if __name__ == "__main__":
    main()
