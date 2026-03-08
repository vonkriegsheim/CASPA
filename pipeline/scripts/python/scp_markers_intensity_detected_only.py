#!/usr/bin/env python3
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ttest_ind

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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True, help="SCP pivot TSV (log2 + median-normalised, NaN for missing)")
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--protein_col", default=None)
    ap.add_argument("--test", choices=["wilcoxon", "welch"], default="wilcoxon")
    ap.add_argument("--min_detected_in", type=int, default=5)
    ap.add_argument("--min_detected_out", type=int, default=5)
    args = ap.parse_args()

    ann = pd.read_csv(args.annotation, sep="\t", dtype=str).fillna("")
    if "Run" not in ann.columns or "Condition" not in ann.columns:
        raise SystemExit("annotation must contain Run and Condition columns.")
    run_to_cond = dict(zip(ann["Run"].astype(str), ann["Condition"].astype(str)))
    runs_order = ann["Run"].astype(str).tolist()

    mat = pd.read_csv(args.matrix, sep="\t", dtype=str)
    header_map = {c: stem_from_header(c) for c in mat.columns}
    sample_cols_all = [c for c in mat.columns if is_sample_column(c)]
    sample_cols = [c for c in sample_cols_all if header_map.get(c, "") in run_to_cond]
    if not sample_cols:
        raise SystemExit("No sample columns matched between pivot and annotation.")

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

    # pivot_pack.tsv is already log2 + median-normalised; NaN = missing
    X = mat[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    conds = np.array([run_to_cond[r] for r in runs], dtype=object)
    clusters = sorted(pd.unique(conds).tolist())

    rows = []
    for cl in clusters:
        in_mask = (conds == cl)
        out_mask = ~in_mask
        if in_mask.sum() < 2 or out_mask.sum() < 2:
            continue

        pvals = []
        tmp_rows = []
        for i in range(X.shape[0]):
            xin = X[i, in_mask]
            xout = X[i, out_mask]
            xin = xin[np.isfinite(xin)]
            xout = xout[np.isfinite(xout)]

            if xin.size < args.min_detected_in or xout.size < args.min_detected_out:
                continue

            med_in = float(np.median(xin))
            med_out = float(np.median(xout))
            log2fc = med_in - med_out

            if args.test == "wilcoxon":
                p = mannwhitneyu(xin, xout, alternative="two-sided").pvalue
            else:
                p = ttest_ind(xin, xout, equal_var=False).pvalue

            tmp_rows.append({
                "Cluster": cl,
                "Protein": str(mat.iloc[i][prot_col]),
                "n_detected_in": int(xin.size),
                "n_detected_out": int(xout.size),
                "median_in": med_in,
                "median_out": med_out,
                "log2FC_detected_only": float(log2fc),
                "pvalue": float(p),
            })
            pvals.append(p)

        if tmp_rows:
            q = bh_fdr(np.array(pvals, dtype=float))
            for rr, qq in zip(tmp_rows, q):
                rr["qvalue"] = float(qq)
                rows.append(rr)

    out_df = pd.DataFrame(rows)
    if out_df.empty:
        out_df = pd.DataFrame(columns=[
            "Cluster","Protein","n_detected_in","n_detected_out","median_in","median_out",
            "log2FC_detected_only","pvalue","qvalue"
        ])
    out_df = out_df.merge(meta, on="Protein", how="left")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, sep="\t", index=False)

if __name__ == "__main__":
    main()
