#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


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


def per_run_median_normalise_log(X_log: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    run_meds = np.nanmedian(X_log, axis=0)
    global_med = float(np.nanmedian(run_meds))
    shifts = global_med - run_meds
    Xn = X_log + shifts[None, :]
    return Xn, shifts, global_med


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pg-matrix", required=True, help="Reduced DIA-NN-like report.pg_matrix.tsv (linear)")
    ap.add_argument("--manifest", required=True, help="Reduced manifest (CSV/TSV)")
    ap.add_argument("--run-col", default="sample_id", help="Manifest column holding stem run ids")
    ap.add_argument("--out-pivot", required=True, help="Output pivot TSV (log2+median-normalised)")
    ap.add_argument("--out-shifts", required=True, help="Output per-run shift TSV")
    ap.add_argument("--out-report", required=True, help="Output JSON report")
    ap.add_argument("--zero-is-missing", action="store_true", default=True)
    ap.add_argument("--log-transform", action="store_true", default=True, help="Apply log2(x+1) to linear pg_matrix")
    ap.add_argument("--median-normalise", action="store_true", default=True)
    args = ap.parse_args()

    # manifest: autodetect delimiter
    man = pd.read_csv(args.manifest, sep=None, engine="python", dtype=str).fillna("")
    if args.run_col not in man.columns:
        raise SystemExit(f"Manifest missing run-col '{args.run_col}'. Columns: {man.columns.tolist()}")

    man["RunStem"] = man[args.run_col].astype(str).str.strip()
    man = man[man["RunStem"] != ""].copy()
    wanted = set(man["RunStem"].tolist())

    mat = pd.read_csv(args.pg_matrix, sep="\t", dtype=str)
    header_map = {c: stem_from_header(c) for c in mat.columns}

    all_sample_cols = [c for c in mat.columns if is_sample_column(c)]
    sample_cols = [c for c in all_sample_cols if header_map.get(c, "") in wanted]

    if not sample_cols:
        raise SystemExit("No sample columns matched between pg_matrix and manifest.")

    # reorder to manifest order
    col_to_stem = {c: header_map[c] for c in sample_cols}
    stem_to_col = {}
    for c in sample_cols:
        stem_to_col.setdefault(col_to_stem[c], c)

    stems = [s for s in man["RunStem"].tolist() if s in stem_to_col]
    cols = [stem_to_col[s] for s in stems]

    # meta columns = everything not a sample column structurally
    meta_cols = [c for c in mat.columns if c not in all_sample_cols]

    X = mat[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float, copy=True)
    if args.zero_is_missing:
        X[X == 0] = np.nan

    # Drop rows with no detected values in the selected sample columns.
    # These arise when a protein passed global FDR in the full matrix but was
    # never quantified in the cells included in this manifest subset.
    n_det_per_row = np.sum(np.isfinite(X), axis=1)
    empty_mask = n_det_per_row == 0
    if empty_mask.any():
        n_dropped = int(empty_mask.sum())
        print(f"Dropping {n_dropped} protein row(s) with no detected values in selected samples.")
        keep = ~empty_mask
        X = X[keep]
        mat = mat.iloc[keep].reset_index(drop=True)

    # Auto-detect intensity scale and apply log2(x+1) if needed
    finite_vals = X[np.isfinite(X)]
    p95 = float(np.percentile(finite_vals, 95)) if finite_vals.size > 0 else 0.0
    if args.log_transform:
        if finite_vals.size > 0 and p95 > 25:
            print(f"Detected linear-scale intensities (p95={p95:.1f}), applying log2(x+1).")
            X_log = np.where(np.isnan(X), np.nan, np.log2(X + 1.0))
            autodetect_decision = "linear_to_log2"
        else:
            print(f"Intensities appear log-transformed (p95={p95:.1f}), using as-is.")
            X_log = X.copy()
            autodetect_decision = "already_log"
    else:
        X_log = X
        autodetect_decision = "log_transform_disabled"

    shifts = np.zeros(X_log.shape[1], dtype=float)
    global_med = float("nan")
    if args.median_normalise:
        X_norm, shifts, global_med = per_run_median_normalise_log(X_log)
    else:
        X_norm = X_log

    out = mat[meta_cols].copy()
    out = pd.concat([out, pd.DataFrame(X_norm, columns=cols)], axis=1)

    # write pivot
    out_path = Path(args.out_pivot)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)

    # write shifts
    shifts_df = pd.DataFrame({
        "RunStem": stems,
        "SampleColumn": cols,
        "median_shift_applied": shifts,
    })
    shifts_df.to_csv(Path(args.out_shifts), sep="\t", index=False)

    report = {
        "tool": "scp_make_pivot_from_pg_matrix",
        "pg_matrix": args.pg_matrix,
        "manifest": args.manifest,
        "n_proteins_raw": int(n_det_per_row.size),
        "n_proteins_dropped_empty": int(empty_mask.sum()),
        "n_proteins": int(X.shape[0]),
        "n_runs": int(X.shape[1]),
        "log_transform": bool(args.log_transform),
        "autodetect_decision": autodetect_decision,
        "p95_raw": p95,
        "median_normalise": bool(args.median_normalise),
        "global_median_after_log": None if not np.isfinite(global_med) else float(global_med),
    }
    Path(args.out_report).write_text(json.dumps(report, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
