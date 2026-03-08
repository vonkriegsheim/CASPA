#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

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
    cl = col_name.lower()
    if any(ext in cl for ext in (".raw", ".d", ".mzml", ".wiff")):
        return True
    if "\\" in col_name or "/" in col_name:
        return True
    if "." not in col_name:
        return True
    return False

def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    m = np.isfinite(p)
    if not m.any():
        return out
    pv = p[m]
    order = np.argsort(pv)
    ranked = pv[order]
    n = ranked.size
    q = ranked * (n / (np.arange(1, n + 1)))
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty_like(q)
    tmp[order] = q
    out[m] = tmp
    return out

def bonferroni(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan)
    m = np.isfinite(p)
    if not m.any():
        return out
    pv = p[m]
    out[m] = np.clip(pv * pv.size, 0, 1)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, help="SCP pivot TSV (log2, NA for missing)")
    ap.add_argument("--annotation", required=True, help="scp_annotation.tsv (Run=stem, Condition)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--protein_col", default=None)
    ap.add_argument("--min_cells_detected", type=int, default=5)
    ap.add_argument("--p_adjust", choices=["BH", "bonferroni"], default="BH")
    args = ap.parse_args()

    ann = pd.read_csv(args.annotation, sep="\t", dtype=str).fillna("")
    if "Run" not in ann.columns or "Condition" not in ann.columns:
        raise SystemExit("annotation must contain Run and Condition")

    run_to_cond = dict(zip(ann["Run"].astype(str), ann["Condition"].astype(str)))
    runs_order = ann["Run"].astype(str).tolist()

    mat = pd.read_csv(args.matrix, sep="\t", dtype=str)
    header_map = {c: stem_from_header(c) for c in mat.columns}
    sample_cols_all = [c for c in mat.columns if is_sample_column(c)]
    sample_cols = [c for c in sample_cols_all if header_map.get(c, "") in run_to_cond]

    if not sample_cols:
        raise SystemExit("No sample columns matched between matrix and annotation.")

    # reorder
    col_to_run = {c: header_map[c] for c in sample_cols}
    run_to_col = {}
    for c in sample_cols:
        run_to_col.setdefault(col_to_run[c], c)
    runs = [r for r in runs_order if r in run_to_col]
    cols = [run_to_col[r] for r in runs]

    prot_col = args.protein_col
    if prot_col is None:
        prot_col = "Protein.Group" if "Protein.Group" in mat.columns else ("Protein" if "Protein" in mat.columns else mat.columns[0])

    meta_cols = [c for c in ["Genes", "Protein.Names", "First.Protein.Description"] if c in mat.columns]
    meta = mat[[prot_col] + meta_cols].copy().rename(columns={prot_col: "Protein"}).drop_duplicates(subset=["Protein"])

    X = mat[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    det = np.isfinite(X).astype(np.int8)

    det_counts = det.sum(axis=1)
    keep = det_counts >= args.min_cells_detected
    det = det[keep, :]
    prot_ids = mat.loc[keep, prot_col].astype(str).tolist()

    conds = np.array([run_to_cond[r] for r in runs], dtype=object)
    clusters = sorted(pd.unique(conds).tolist())

    rows = []
    for cl in clusters:
        in_mask = (conds == cl)
        out_mask = ~in_mask
        n_in = int(in_mask.sum())
        n_out = int(out_mask.sum())
        if n_in < 2 or n_out < 2:
            continue

        a = det[:, in_mask].sum(axis=1)
        b = n_in - a
        c = det[:, out_mask].sum(axis=1)
        d = n_out - c

        pvals = np.full(det.shape[0], np.nan)
        ors = np.full(det.shape[0], np.nan)

        for i in range(det.shape[0]):
            table = np.array([[a[i], b[i]], [c[i], d[i]]], dtype=int)
            ors[i] = ((a[i] + 0.5) * (d[i] + 0.5)) / ((b[i] + 0.5) * (c[i] + 0.5))
            pvals[i] = fisher_exact(table, alternative="two-sided")[1]

        adj = bh_fdr(pvals) if args.p_adjust == "BH" else bonferroni(pvals)
        det_in = a / max(n_in, 1)
        det_out = c / max(n_out, 1)

        for i, pid in enumerate(prot_ids):
            rows.append({
                "Cluster": cl,
                "Protein": pid,
                "det_rate_in": float(det_in[i]),
                "det_rate_out": float(det_out[i]),
                "odds_ratio": float(ors[i]),
                "log2_odds_ratio": float(np.log2(ors[i])) if ors[i] > 0 else np.nan,
                "pvalue": float(pvals[i]),
                "qvalue": float(adj[i]),
                "n_in": n_in,
                "n_out": n_out,
                "det_in": int(a[i]),
                "det_out": int(c[i]),
            })

    out_df = pd.DataFrame(rows)
    out_df = out_df.merge(meta, on="Protein", how="left")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)

if __name__ == "__main__":
    main()
