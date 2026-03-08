#!/usr/bin/env python3
"""
scp_qc_filter_cells_pg_matrix.py

SCP-oriented QC filtering + QC plots from DIA-NN report.pg_matrix.tsv (wide matrix).

Key features:
- Forced exclusion of library/control runs by Run-name regex (e.g. contains 'cells' or starts with 'library_').
- Bottom-N median * multiplier cutoff on number of detected proteins per run/cell.
- Writes backups and reduced inputs for downstream:
    - report.complete.pg_matrix.tsv (backup)
    - report.pg_matrix.tsv          (reduced: meta cols + kept sample cols only)
    - sample_sheet.complete.csv     (backup)
    - sample_sheet.csv              (reduced manifest, kept runs only; optional overwrite)
- QC plots (6):
    01 detected proteins per run (scatter)
    02 missingness per run (scatter)
    03 violin of intensity distributions + individual points (strip) per "Class"
    04 correlation heatmap (subset)
    05 PCA (subset)
    06 mean vs missingness scatter
"""

import argparse
import json
import re
import shutil
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
    if "\\" in s or "/" in s:
        return Path(s).stem
    if "." in s and s.rsplit(".", 1)[1].lower() in ("raw", "d", "mzml", "wiff"):
        return Path(s).stem
    return s


def is_sample_column(col_name: str) -> bool:
    if col_name in KNOWN_META_COLS:
        return False
    col_lower = col_name.lower()
    if any(ext in col_lower for ext in (".raw", ".d", ".mzml", ".wiff")):
        return True
    if "\\" in col_name or "/" in col_name:
        return True
    if "." not in col_name:
        return True
    return False


def compute_detected_cutoff(n_detected: pd.Series, bottom_n: int, multiplier: float, min_detected: int):
    x = pd.to_numeric(n_detected, errors="coerce").dropna().values
    if x.size == 0:
        return float(min_detected), 0.0, float(min_detected)

    xs = np.sort(x)
    b = xs[: min(bottom_n, xs.size)]
    bottom_med = float(np.median(b)) if b.size else 0.0
    cutoff_dynamic = float(multiplier) * bottom_med
    cutoff_final = max(float(min_detected), cutoff_dynamic)
    return cutoff_final, bottom_med, cutoff_dynamic


def savefig(fig, outpath: Path, fmt: str, dpi: int):
    outpath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpath.with_suffix(f".{fmt}"), dpi=dpi, bbox_inches="tight", facecolor="white")


def backup_copy(src: Path, dst: Path):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def choose_class_from_run(run: str) -> str:
    """
    Used for plot 3 grouping after forced exclusion.
    Returns things like '1cell', '2cell', ... if present; else 'other'.
    """
    s = str(run).strip()
    m = re.match(r"(?i)^(?P<n>\d+)cell\b", s)
    if m:
        return f"{m.group('n')}cell"
    # If someone has '1cells' etc (rare), normalize:
    m = re.match(r"(?i)^(?P<n>\d+)cells\b", s)
    if m:
        return f"{m.group('n')}cells"
    if re.match(r"(?i)^library\b", s):
        return "library"
    return "other"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pg-matrix", required=True, help="report.pg_matrix.tsv")
    ap.add_argument("--manifest", required=True, help="sample_sheet.csv (or SCP manifest)")
    ap.add_argument("--outdir", required=True, help="work/<dataset>/scp/qc")
    ap.add_argument("--run-col", default="sample_id", help="Manifest column used as Run ID (default sample_id)")

    # QC cutoff parameters
    ap.add_argument("--bottom-n", type=int, default=5)
    ap.add_argument("--multiplier", type=float, default=1.7)
    ap.add_argument("--min-detected", type=int, default=400)

    # Forced library exclusion (by run name)
    ap.add_argument(
        "--exclude-run-regex",
        default=r"(^library[_-]|\bcells\b)",
        help="Regex applied to Run name to force-exclude runs before protein-count QC. Case-insensitive.",
    )

    # Reduced-input writing/backups
    ap.add_argument("--write-reduced-inputs", action="store_true", help="Write backups + reduced pg_matrix + reduced manifest.")
    ap.add_argument("--backup-tag", default="complete", help="Backup tag; writes report.<tag>.pg_matrix.tsv and sample_sheet.<tag>.csv")
    ap.add_argument("--overwrite-pg-matrix", action="store_true", help="Overwrite input pg_matrix with reduced matrix.")
    ap.add_argument("--overwrite-manifest", action="store_true", help="Overwrite input manifest with reduced manifest.")
    ap.add_argument("--reduced-pg-name", default="report.filtered.pg_matrix.tsv", help="Reduced pg_matrix filename (if not overwriting).")
    ap.add_argument("--reduced-manifest-name", default="sample_sheet.filtered.csv", help="Reduced manifest filename (if not overwriting).")

    # Plot params
    ap.add_argument("--format", choices=["pdf", "png"], default="pdf")
    ap.add_argument("--dpi", type=int, default=250)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--max-proteins-for-corr", type=int, default=15000)
    ap.add_argument("--max-cells-heatmap", type=int, default=120)
    ap.add_argument("--max-points-per-run", type=int, default=2000, help="For plot 3 strip points: cap points per run.")
    ap.add_argument("--max-total-points", type=int, default=200000, help="For plot 3: cap total points across runs.")

    # Explicit output paths (Snakemake-friendly)
    ap.add_argument("--out-filtered-pg-matrix", default=None,
                    help="Explicit output path for filtered pg_matrix (Snakemake-friendly)")
    ap.add_argument("--out-filtered-manifest", default=None,
                    help="Explicit output path for filtered manifest (Snakemake-friendly)")

    # Debug
    ap.add_argument("--debug", action="store_true")

    args = ap.parse_args()
    np.random.seed(args.seed)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap

    outdir = Path(args.outdir)
    plots_dir = outdir / "plots"
    tables_dir = outdir / "tables"
    plots_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    pg_path = Path(args.pg_matrix)
    man_path = Path(args.manifest)

    # ---- Load manifest ----
    sep = "\t" if man_path.suffix.lower() in (".tsv", ".txt") else ","
    man_raw = pd.read_csv(man_path, sep=sep, dtype=str).fillna("")
    if args.run_col not in man_raw.columns:
        raise SystemExit(f"Manifest missing --run-col '{args.run_col}'. Columns: {man_raw.columns.tolist()}")

    man = man_raw.copy()
    man["Run"] = man[args.run_col].astype(str).str.strip()
    man = man[man["Run"] != ""].copy()

    if args.debug:
        print("[DEBUG] running:", __file__)
        print("[DEBUG] manifest rows:", man.shape[0])
        print("[DEBUG] run-col:", args.run_col)
        print("[DEBUG] example runs:", man["Run"].head(5).tolist())

    # ---- Forced exclusion by Run regex ----
    forced_excl = pd.DataFrame()
    if args.exclude_run_regex:
        m = man["Run"].astype(str).str.contains(args.exclude_run_regex, case=False, regex=True, na=False)
        forced_excl = man[m].copy()
        man = man[~m].copy()
        forced_excl["reason"] = f"excluded_by_run_regex: {args.exclude_run_regex}"
        forced_excl.to_csv(tables_dir / "scp_exclusions_forced_library.tsv", sep="\t", index=False)

    if man.empty:
        raise SystemExit("All runs were excluded by --exclude-run-regex. Adjust the regex.")

    wanted_runs = set(man["Run"].tolist())

    # ---- Load pg_matrix ----
    mat = pd.read_csv(pg_path, sep="\t", dtype=str)
    header_map = {c: stem_from_header(c) for c in mat.columns}

    # Identify ALL sample columns in matrix (structural)
    all_sample_cols = [c for c in mat.columns if is_sample_column(c)]
    meta_cols = [c for c in mat.columns if c not in all_sample_cols]

    # Identify sample columns that match the manifest runs (post forced exclusion)
    run_cols = [c for c in all_sample_cols if header_map.get(c, "") in wanted_runs]
    if not run_cols:
        raise SystemExit(
            "No matching run columns found in pg_matrix for the provided manifest AFTER forced exclusion.\n"
            f"Manifest runs (first 10): {man['Run'].head(10).tolist()}\n"
            f"Matrix sample cols (first 20): {all_sample_cols[:20]}"
        )

    # Reorder to manifest order
    col_to_run = {c: header_map[c] for c in run_cols}
    run_to_col = {}
    for c in run_cols:
        run_to_col.setdefault(col_to_run[c], c)

    ordered_runs = [r for r in man["Run"].tolist() if r in run_to_col]
    ordered_cols = [run_to_col[r] for r in ordered_runs]

    # Numeric matrix
    X = mat[ordered_cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float, copy=True)
    X[X == 0] = np.nan  # treat zeros as missing

    n_prot, n_cells = X.shape
    print(f"SCP QC: matrix {n_prot} proteins × {n_cells} runs (after forced exclusion)")

    present = np.isfinite(X)
    n_detected = present.sum(axis=0)
    missing_frac = 1.0 - (n_detected / max(n_prot, 1))

    X_log = np.where(np.isnan(X), np.nan, np.log2(X + 1.0))
    mean_log = np.nanmean(X_log, axis=0)
    median_log = np.nanmedian(X_log, axis=0)

    qc = pd.DataFrame({
        "Run": ordered_runs,
        "n_detected_proteins": n_detected,
        "missing_fraction": missing_frac,
        "mean_log2_intensity": mean_log,
        "median_log2_intensity": median_log,
    })

    # For plotting groups (post forced exclusion)
    qc["Class"] = qc["Run"].map(choose_class_from_run)

    # ---- Apply cutoff on remaining runs ----
    cutoff_final, bottom_med, cutoff_dynamic = compute_detected_cutoff(
        qc["n_detected_proteins"],
        bottom_n=args.bottom_n,
        multiplier=args.multiplier,
        min_detected=args.min_detected,
    )
    qc["keep"] = qc["n_detected_proteins"] >= cutoff_final

    qc.to_csv(tables_dir / "scp_qc_metrics.tsv", sep="\t", index=False)

    excl_cutoff = qc[~qc["keep"]].copy()
    excl_cutoff["reason"] = (
        f"n_detected < cutoff; cutoff_final={cutoff_final:.1f} "
        f"(max(min_detected={args.min_detected}, multiplier={args.multiplier} * bottom{args.bottom_n}_median={bottom_med:.1f}))"
    )
    excl_cutoff.to_csv(tables_dir / "scp_exclusions_cutoff.tsv", sep="\t", index=False)

    kept = qc[qc["keep"]].copy()
    kept_runs = set(kept["Run"].tolist())

    # ---- Reduced manifest contents (kept runs only) ----
    man_kept = man[man["Run"].isin(kept_runs)].copy()
    man_kept.to_csv(tables_dir / "scp_manifest.filtered.tsv", sep="\t", index=False)

    # ---- Reduced pg_matrix contents (meta cols + kept sample cols only) ----
    kept_cols = [run_to_col[r] for r in ordered_runs if r in kept_runs]
    reduced_mat = mat[meta_cols + kept_cols].copy()

    # ---- Write backups + reduced inputs if requested ----
    pg_backup = None
    man_backup = None
    pg_reduced_path = None
    man_reduced_path = None

    if args.write_reduced_inputs:
        # backups
        pg_backup = pg_path.with_name(f"report.{args.backup_tag}.pg_matrix.tsv")
        man_backup = man_path.with_name(f"sample_sheet.{args.backup_tag}.csv")
        backup_copy(pg_path, pg_backup)
        backup_copy(man_path, man_backup)

        # choose output names
        if args.overwrite_pg_matrix:
            pg_reduced_path = pg_path
        else:
            pg_reduced_path = pg_path.with_name(args.reduced_pg_name)

        if args.overwrite_manifest:
            man_reduced_path = man_path
        else:
            man_reduced_path = man_path.with_name(args.reduced_manifest_name)

        # write reduced pg_matrix
        reduced_mat.to_csv(pg_reduced_path, sep="\t", index=False)

        # write reduced manifest in the SAME column schema as original manifest (drop added 'Run')
        out_cols = list(man_raw.columns)
        man_kept[out_cols].to_csv(man_reduced_path, index=False)

    # ---- Explicit output paths (Snakemake) ----
    if args.out_filtered_pg_matrix:
        Path(args.out_filtered_pg_matrix).parent.mkdir(parents=True, exist_ok=True)
        reduced_mat.to_csv(args.out_filtered_pg_matrix, sep="\t", index=False)

    if args.out_filtered_manifest:
        Path(args.out_filtered_manifest).parent.mkdir(parents=True, exist_ok=True)
        out_cols = list(man_raw.columns)
        man_kept[out_cols].to_csv(args.out_filtered_manifest, sep="\t", index=False)

    # ---- JSON report ----
    report = {
        "tool": "scp_qc_filter_cells_pg_matrix",
        "parameters": {
            "bottom_n": args.bottom_n,
            "multiplier": args.multiplier,
            "min_detected": args.min_detected,
            "exclude_run_regex": args.exclude_run_regex,
            "write_reduced_inputs": bool(args.write_reduced_inputs),
            "overwrite_pg_matrix": bool(args.overwrite_pg_matrix),
            "overwrite_manifest": bool(args.overwrite_manifest),
        },
        "paths": {
            "pg_matrix_in": str(pg_path),
            "manifest_in": str(man_path),
            "pg_matrix_backup": str(pg_backup) if pg_backup else None,
            "manifest_backup": str(man_backup) if man_backup else None,
            "pg_matrix_reduced": str(pg_reduced_path) if pg_reduced_path else None,
            "manifest_reduced": str(man_reduced_path) if man_reduced_path else None,
        },
        "matrix_shape_after_forced_exclusion": {"n_proteins": int(n_prot), "n_runs": int(n_cells)},
        "cutoff": {
            "bottom_median": float(bottom_med),
            "cutoff_dynamic": float(cutoff_dynamic),
            "cutoff_final": float(cutoff_final),
        },
        "kept_after_cutoff": int(kept.shape[0]),
        "excluded_by_cutoff": int(excl_cutoff.shape[0]),
        "excluded_by_run_regex": int(forced_excl.shape[0]),
    }
    (outdir / "scp_qc_report.json").write_text(json.dumps(report, indent=2))

    # ---- Plotting palettes ----
    class_order = sorted(qc["Class"].unique().tolist())
    class_pal = dict(zip(class_order, sns.color_palette("Set2", n_colors=max(3, len(class_order)))))

    qc_plot = qc.sort_values(["Class", "n_detected_proteins"]).reset_index(drop=True)

    # ---- Plot 1: proteins detected per run ----
    fig, ax = plt.subplots(figsize=(max(10, len(qc_plot) * 0.03), 4.0))
    ax.scatter(
        np.arange(len(qc_plot)), qc_plot["n_detected_proteins"],
        c=[class_pal.get(x, (0.5, 0.5, 0.5)) for x in qc_plot["Class"]],
        s=12, alpha=0.85, linewidths=0
    )
    ax.axhline(cutoff_final, color="#444444", linestyle="--", linewidth=0.8)
    ax.set_title("SCP QC: proteins detected per run (post forced exclusion)")
    ax.set_xlabel("Runs (ordered by Class, n_detected)")
    ax.set_ylabel("# detected proteins")
    savefig(fig, plots_dir / "01_detected_proteins_per_run", args.format, args.dpi)
    plt.close(fig)

    # ---- Plot 2: missingness per run ----
    fig, ax = plt.subplots(figsize=(max(10, len(qc_plot) * 0.03), 4.0))
    ax.scatter(
        np.arange(len(qc_plot)), qc_plot["missing_fraction"] * 100.0,
        c=[class_pal.get(x, (0.5, 0.5, 0.5)) for x in qc_plot["Class"]],
        s=12, alpha=0.85, linewidths=0
    )
    ax.set_title("SCP QC: missingness per run (post forced exclusion)")
    ax.set_xlabel("Runs (ordered by Class, n_detected)")
    ax.set_ylabel("% missing proteins")
    savefig(fig, plots_dir / "02_missingness_per_run", args.format, args.dpi)
    plt.close(fig)

    # ---- Plot 3: violin of identified protein groups per run (n_detected) ----
    # Use kept runs to reflect post-QC dataset; switch to `qc` if you prefer all runs post-forced-exclusion.
    plot3_df = kept.copy()
    if plot3_df.empty:
        # Fallback: if nothing kept, plot from qc to avoid crash
        plot3_df = qc.copy()

    fig, ax = plt.subplots(figsize=(max(6, len(class_order) * 1.2), 4.4))
    sns.violinplot(
        data=plot3_df, x="Class", y="n_detected_proteins",
        palette=class_pal, cut=0, inner="quartile", linewidth=0.7, ax=ax
    )
    # Overlay individual runs as points
    sns.stripplot(
        data=plot3_df, x="Class", y="n_detected_proteins",
        color="black", alpha=0.35, size=3.0, jitter=0.25, ax=ax
    )
    ax.set_title("Protein groups identified per run (post-QC kept runs)")
    ax.set_xlabel("")
    ax.set_ylabel("# detected protein groups")
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    savefig(fig, plots_dir / "03_detected_proteins_violin_by_class_with_points", args.format, args.dpi)
    plt.close(fig)


    # ---- Plot 4: correlation heatmap on kept runs (subset) ----
    kept_list = kept["Run"].tolist()
    max_cells_hm = args.max_cells_heatmap
    if len(kept_list) > max_cells_hm:
        kept_list = sorted(
            kept_list,
            key=lambda r: kept.loc[kept["Run"] == r, "n_detected_proteins"].iloc[0],
            reverse=True
        )[:max_cells_hm]
        print(f"Heatmap: using top {max_cells_hm} runs by n_detected for readability.")

    col_idx = [ordered_runs.index(r) for r in kept_list if r in ordered_runs]
    Xc = X_log[:, col_idx].copy()

    # Impute per-column median for correlation
    col_meds = np.nanmedian(Xc, axis=0)
    col_meds = np.where(np.isnan(col_meds), 0.0, col_meds)
    miss = ~np.isfinite(Xc)
    if miss.any():
        Xc[miss] = np.take(col_meds, np.where(miss)[1])

    # Subsample proteins for speed
    if Xc.shape[0] > args.max_proteins_for_corr:
        ridx = np.random.choice(np.arange(Xc.shape[0]), size=args.max_proteins_for_corr, replace=False)
        Xc = Xc[ridx, :]

    corr = np.corrcoef(Xc.T)
    corr_df = pd.DataFrame(corr, index=kept_list, columns=kept_list)

    corr_vals = corr_df.to_numpy().copy()
    np.fill_diagonal(corr_vals, np.nan)
    vmin = float(np.nanquantile(corr_vals, 0.02))
    vmax = float(np.nanquantile(corr_vals, 0.98))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
        vmin, vmax = 0.0, 1.0
    elif (vmax - vmin) < 0.05:
        mid = 0.5 * (vmin + vmax)
        vmin = max(0.0, mid - 0.05)
        vmax = min(1.0, mid + 0.05)

    cmap = LinearSegmentedColormap.from_list("white_red", ["#FFFFFF", "#B2182B"])

    fig, ax = plt.subplots(figsize=(9.0, 8.0))
    sns.heatmap(corr_df, cmap=cmap, vmin=vmin, vmax=vmax, square=True, cbar_kws={"label": "Pearson r"}, ax=ax)
    ax.set_title(f"Run correlation (subset), scale [{vmin:.2f}, {vmax:.2f}]")
    ax.tick_params(axis="x", rotation=90, labelsize=5)
    ax.tick_params(axis="y", rotation=0, labelsize=5)
    savefig(fig, plots_dir / "04_run_correlation_heatmap_subset", args.format, args.dpi)
    plt.close(fig)

    # ---- Plot 5: PCA of kept runs (subset) ----
    Xp = X_log[:, col_idx].copy()
    prot_meds = np.nanmedian(Xp, axis=1)
    prot_meds = np.where(np.isnan(prot_meds), 0.0, prot_meds)
    miss = ~np.isfinite(Xp)
    if miss.any():
        Xp[miss] = prot_meds[np.where(miss)[0]]

    Xp = Xp - Xp.mean(axis=1, keepdims=True)
    sd = Xp.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    Xp = Xp / sd

    U, S, Vt = np.linalg.svd(Xp.T, full_matrices=False)
    pcs = U[:, :2] * S[:2]

    pca_df = pd.DataFrame({
        "Run": kept_list,
        "PC1": pcs[:, 0],
        "PC2": pcs[:, 1],
        "Class": [choose_class_from_run(r) for r in kept_list],
    })

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    for cls in sorted(pca_df["Class"].unique().tolist()):
        sub = pca_df[pca_df["Class"] == cls]
        ax.scatter(sub["PC1"], sub["PC2"], s=35, alpha=0.85, label=cls, color=class_pal.get(cls, "#777777"),
                   edgecolor="white", linewidth=0.5)
    ax.set_title("PCA of runs (subset; PCA-only median impute)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.legend(fontsize=7, title="Class", title_fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, plots_dir / "05_pca_runs_subset", args.format, args.dpi)
    plt.close(fig)

    # ---- Plot 6: mean vs missingness ----
    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    for cls in sorted(qc["Class"].unique().tolist()):
        sub = qc[qc["Class"] == cls]
        ax.scatter(sub["missing_fraction"] * 100.0, sub["mean_log2_intensity"], s=28, alpha=0.85,
                   label=cls, color=class_pal.get(cls, "#777777"), edgecolor="white", linewidth=0.4)
    ax.set_title("Mean intensity vs missingness (all post-forced-exclusion runs)")
    ax.set_xlabel("% missing proteins")
    ax.set_ylabel("Mean log2 intensity")
    ax.legend(fontsize=7, title="Class", title_fontsize=8, bbox_to_anchor=(1.02, 1), loc="upper left")
    savefig(fig, plots_dir / "06_mean_vs_missingness", args.format, args.dpi)
    plt.close(fig)

    print(
        f"SCP QC complete. "
        f"Forced-excluded(library)={report['excluded_by_run_regex']} "
        f"Cutoff-excluded={report['excluded_by_cutoff']} "
        f"Kept={report['kept_after_cutoff']} "
        f"Cutoff={cutoff_final:.1f}"
    )
    print(f"Outputs: {outdir}")


if __name__ == "__main__":
    main()
