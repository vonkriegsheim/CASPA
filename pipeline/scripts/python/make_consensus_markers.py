#!/usr/bin/env python3
"""
make_consensus_markers.py

Aggregate detection, intensity, and scplainer markers into a single
consensus ranking per cluster using a score-based Borda count.

Score per modality per protein = abs(FC) * -log10(qvalue), normalised
to [0, 1] within cluster.  Proteins not present in a modality receive 0
for that modality.  Consensus score = mean of available modality scores.

Inputs
------
--detection-markers   detection_markers.tsv
--intensity-markers   intensity_markers_detected_only.tsv
--scplainer-markers   scplainer_intensity_markers.tsv
--out                 output TSV path
--q-threshold         pre-filter: only include rows with qvalue <= threshold
                      before ranking (default 1.0 = no pre-filter)
--top-n               rows to keep per cluster in output (default 50, 0 = all)

Output columns
--------------
Cluster, Protein, Genes, consensus_score,
score_detection, score_intensity, score_scplainer,
n_modalities,
det_rate_in, log2OR, qvalue_detection,
log2FC_intensity, qvalue_intensity,
log2FC_scplainer, qvalue_scplainer
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def score_modality(df: pd.DataFrame, fc_col: str, q_col: str) -> pd.Series:
    """Compute abs(FC) * -log10(qvalue) score, return Series indexed by (Cluster, Protein).
    Deduplicates by keeping the max score per (Cluster, Protein)."""
    df = df.copy()
    df["_fc"] = pd.to_numeric(df[fc_col], errors="coerce").abs()
    df["_q"]  = pd.to_numeric(df[q_col],  errors="coerce")
    df["_score"] = df["_fc"] * (-np.log10(df["_q"].clip(1e-300, 1.0)))
    df["_score"] = df["_score"].fillna(0.0)
    return df.groupby(["Cluster", "Protein"])["_score"].max()


def normalise_within_cluster(score_series: pd.Series) -> pd.Series:
    """Normalise scores to [0, 1] within each cluster."""
    def _norm(grp):
        mn, mx = grp.min(), grp.max()
        if mx == mn:
            return grp * 0.0
        return (grp - mn) / (mx - mn)
    return score_series.groupby(level="Cluster", group_keys=False).apply(_norm)


def best_row_per_prot(df: pd.DataFrame, cluster_col: str, prot_col: str,
                      score_col: str) -> pd.DataFrame:
    """Keep the row with highest score per (cluster, protein)."""
    df = df.copy()
    df["_s"] = pd.to_numeric(df.get(score_col, pd.Series(0.0, index=df.index)), errors="coerce").fillna(0.0)
    return (df.sort_values("_s", ascending=False)
              .drop_duplicates(subset=[cluster_col, prot_col])
              .drop(columns=["_s"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--detection-markers",  default=None)
    ap.add_argument("--intensity-markers",  default=None)
    ap.add_argument("--scplainer-markers",  default=None)
    ap.add_argument("--out",                required=True)
    ap.add_argument("--q-threshold",        type=float, default=1.0)
    ap.add_argument("--top-n",              type=int,   default=50)
    args = ap.parse_args()

    # ---- Load each modality ----
    loaded = {}

    if args.detection_markers and Path(args.detection_markers).exists():
        dm = pd.read_csv(args.detection_markers, sep="\t", dtype=str)
        if "Cluster" in dm.columns and "Protein" in dm.columns and \
           "log2_odds_ratio" in dm.columns and "qvalue" in dm.columns:
            dm["qvalue"] = pd.to_numeric(dm["qvalue"], errors="coerce")
            if args.q_threshold < 1.0:
                dm = dm[dm["qvalue"] <= args.q_threshold]
            loaded["detection"] = dm
            print(f"Detection markers: {len(dm)} rows, "
                  f"{dm['Cluster'].nunique()} clusters")

    if args.intensity_markers and Path(args.intensity_markers).exists():
        im = pd.read_csv(args.intensity_markers, sep="\t", dtype=str)
        if "Cluster" in im.columns and "Protein" in im.columns and \
           "log2FC_detected_only" in im.columns and "qvalue" in im.columns:
            im["qvalue"] = pd.to_numeric(im["qvalue"], errors="coerce")
            if args.q_threshold < 1.0:
                im = im[im["qvalue"] <= args.q_threshold]
            loaded["intensity"] = im
            print(f"Intensity markers: {len(im)} rows, "
                  f"{im['Cluster'].nunique()} clusters")

    if args.scplainer_markers and Path(args.scplainer_markers).exists():
        sm = pd.read_csv(args.scplainer_markers, sep="\t", dtype=str)
        # scplainer markers use Cluster = cluster_num
        if "Cluster" not in sm.columns and "contrast" in sm.columns:
            sm["Cluster"] = sm["contrast"]
        if "Cluster" in sm.columns and "Protein" in sm.columns and \
           "log2FC" in sm.columns and "qvalue" in sm.columns:
            sm["qvalue"] = pd.to_numeric(sm["qvalue"], errors="coerce")
            if "ok" in sm.columns:
                sm = sm[sm["ok"].astype(str).str.upper().isin(["TRUE", "1"])]
            if args.q_threshold < 1.0:
                sm = sm[sm["qvalue"] <= args.q_threshold]
            loaded["scplainer"] = sm
            print(f"Scplainer markers: {len(sm)} rows, "
                  f"{sm['Cluster'].nunique()} clusters")

    if not loaded:
        print("No valid marker inputs found. Exiting.")
        return

    # ---- Collect all (Cluster, Protein) combinations ----
    all_keys = set()
    for name, df in loaded.items():
        for _, row in df[["Cluster", "Protein"]].drop_duplicates().iterrows():
            all_keys.add((str(row["Cluster"]), str(row["Protein"])))

    all_keys_df = pd.DataFrame(list(all_keys), columns=["Cluster", "Protein"])

    # ---- Score each modality ----
    scores = {}
    if "detection" in loaded:
        raw = score_modality(loaded["detection"], "log2_odds_ratio", "qvalue")
        scores["score_detection"] = normalise_within_cluster(raw)
    if "intensity" in loaded:
        raw = score_modality(loaded["intensity"], "log2FC_detected_only", "qvalue")
        scores["score_intensity"] = normalise_within_cluster(raw)
    if "scplainer" in loaded:
        raw = score_modality(loaded["scplainer"], "log2FC", "qvalue")
        scores["score_scplainer"] = normalise_within_cluster(raw)

    # ---- Build consensus DataFrame ----
    result = all_keys_df.copy()
    result = result.set_index(["Cluster", "Protein"])

    for col, ser in scores.items():
        result[col] = ser
    result = result.fillna(0.0).reset_index()

    score_cols = [c for c in ["score_detection", "score_intensity", "score_scplainer"]
                  if c in result.columns]
    result["n_modalities"] = (result[score_cols] > 0).sum(axis=1)
    result["consensus_score"] = result[score_cols].mean(axis=1)

    # ---- Merge annotation columns from each modality ----
    def pick_gene(df):
        if "Genes" in df.columns:
            return df[["Cluster", "Protein", "Genes"]].drop_duplicates(subset=["Cluster", "Protein"])
        return None

    gene_dfs = [pick_gene(df) for df in loaded.values() if pick_gene(df) is not None]
    if gene_dfs:
        genes = pd.concat(gene_dfs, ignore_index=True).drop_duplicates(subset=["Cluster", "Protein"])
        result = result.merge(genes, on=["Cluster", "Protein"], how="left")

    # Merge per-modality stats
    if "detection" in loaded:
        dm_sel = best_row_per_prot(loaded["detection"], "Cluster", "Protein", "log2_odds_ratio")
        det_cols = {c: f"{c}_det" for c in ["det_rate_in", "det_rate_out", "log2_odds_ratio", "qvalue"]
                    if c in dm_sel.columns}
        dm_sel = dm_sel.rename(columns=det_cols)[["Cluster", "Protein"] + list(det_cols.values())]
        result = result.merge(dm_sel, on=["Cluster", "Protein"], how="left")

    if "intensity" in loaded:
        im_sel = best_row_per_prot(loaded["intensity"], "Cluster", "Protein", "log2FC_detected_only")
        int_cols = {c: f"{c}_int" for c in ["log2FC_detected_only", "qvalue"]
                    if c in im_sel.columns}
        im_sel = im_sel.rename(columns=int_cols)[["Cluster", "Protein"] + list(int_cols.values())]
        result = result.merge(im_sel, on=["Cluster", "Protein"], how="left")

    if "scplainer" in loaded:
        sm_sel = best_row_per_prot(loaded["scplainer"], "Cluster", "Protein", "log2FC")
        scp_cols = {c: f"{c}_scp" for c in ["log2FC", "qvalue"]
                    if c in sm_sel.columns}
        sm_sel = sm_sel.rename(columns=scp_cols)[["Cluster", "Protein"] + list(scp_cols.values())]
        result = result.merge(sm_sel, on=["Cluster", "Protein"], how="left")

    # ---- Sort and trim ----
    result = result.sort_values(["Cluster", "consensus_score"], ascending=[True, False])

    if args.top_n > 0:
        result = result.groupby("Cluster", group_keys=False).head(args.top_n)

    # ---- Write ----
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out}: {len(result)} rows across {result['Cluster'].nunique()} clusters")


if __name__ == "__main__":
    main()
