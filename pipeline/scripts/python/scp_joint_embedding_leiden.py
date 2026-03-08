#!/usr/bin/env python3
"""
scp_joint_embedding_leiden.py

SCP joint embedding + Leiden clustering from a DIA-NN-like wide pg_matrix.

Inputs
------
--pg-matrix     work/<dataset>/search/diann/report.pg_matrix.tsv
--manifest      staging/<dataset>/sample_sheet.csv (or filtered SCP manifest TSV)
--run-col       column in manifest that matches run/cell IDs (default: sample_id)
--outdir        work/<dataset>/scp/clustering

Key methods
-----------
- Builds X_int (log intensities, missing=NaN) and X_det (0/1)
- Median normalisation per run (on log intensities, detected-only)
- PCA on intensity (embedding-only impute by protein median)
- PCA on detection (SVD on centred detection)
- Joint embedding = [PC_int, PC_det] with block-wise standardisation
- Optional Harmony batch correction in embedding space (Library)
- Neighbors/UMAP/Leiden
- Writes annotation + contrasts compatible with bulk pipeline style

Outputs
-------
outdir/
  scp_cluster_assignments.tsv
  scp_annotation.tsv
  scp_contrasts.tsv
  scp_clustering_report.json
  qc/
    umap_by_cluster.png
    umap_by_library.png
    umap_by_n_detected.png
    pca_intensity.png
    pca_detection.png
"""

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

def adaptive_dot_size(n: int, s_max_at_12: float = 36.0, s_min: float = 6.0) -> float:
    """
    Return a scatter marker size that decreases as n increases.
    - s_max_at_12 is the dot size around n=12.
      (We set this to ~half the previous 'large' look.)
    - s_min is the floor for very large n.
    """
    n = max(int(n), 1)
    # Smooth decay: size ~ 1/sqrt(n)
    s = s_max_at_12 * np.sqrt(12.0 / n)
    return float(np.clip(s, s_min, s_max_at_12))


def detect_library_column(man: pd.DataFrame) -> tuple[str | None, pd.Series]:
    """
    Detect a manifest column containing values like 'library-*' (case-insensitive).
    Returns (source_col, Library_series). If not found, returns library-1 for all.
    """
    for c in man.columns:
        s = man[c].astype(str)
        if s.str.contains(r"(?i)\blibrary[-_]", regex=True, na=False).any():
            lib = s.fillna("").astype(str).str.strip()
            lib[lib == ""] = "library-1"
            return c, lib
    # fallback by column name
    for c in man.columns:
        if c.strip().lower() in ("library", "batch"):
            lib = man[c].astype(str).fillna("").str.strip()
            lib[lib == ""] = "library-1"
            return c, lib
    return None, pd.Series(["library-1"] * len(man), index=man.index)


def batch_mixing_entropy(cluster_labels: list[str], batch_labels: list[str]) -> tuple[float, dict]:
    """
    Compute mean normalised Shannon entropy of batch distribution per cluster.

    For each cluster, measure how evenly batches are represented.
    Returns (mean_entropy, per_cluster_dict) where entropy is in [0, 1].
    1.0 = perfectly mixed (all batches equally represented in every cluster).
    0.0 = completely segregated (each cluster is one batch).
    """
    df = pd.DataFrame({"cluster": cluster_labels, "batch": batch_labels})
    n_batches_global = df["batch"].nunique()
    if n_batches_global <= 1:
        return 1.0, {}

    max_ent = np.log2(n_batches_global)
    per_cluster = {}
    for cl, grp in df.groupby("cluster"):
        counts = grp["batch"].value_counts().to_numpy(dtype=float)
        props = counts / counts.sum()
        ent = -np.sum(props * np.log2(props + 1e-300))
        per_cluster[cl] = float(ent / max_ent)

    # Weight by cluster size
    sizes = df["cluster"].value_counts()
    weighted = sum(per_cluster[cl] * sizes[cl] for cl in per_cluster) / len(df)
    return float(weighted), per_cluster


def pca_svd(samples_x_features: np.ndarray, n_components: int) -> tuple[np.ndarray, np.ndarray]:
    """
    Simple PCA via SVD on samples x features.
    Returns (scores, singular_values) for first n_components.
    """
    X = samples_x_features
    X = X - X.mean(axis=0, keepdims=True)
    U, S, Vt = np.linalg.svd(X, full_matrices=False)
    k = min(n_components, U.shape[1])
    scores = U[:, :k] * S[:k]
    return scores, S[:k]


def standardise_block(Z: np.ndarray) -> np.ndarray:
    """
    Standardise each column to unit variance (and zero mean) to balance blocks.
    """
    Z = Z - Z.mean(axis=0, keepdims=True)
    sd = Z.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    return Z / sd


def make_contrasts(conditions: list[str]) -> pd.DataFrame:
    """
    Baseline = largest condition (should be passed in separately if desired),
    but here we just create empty; actual baseline chosen later from sizes.
    """
    return pd.DataFrame(columns=["Contrast", "Numerator", "Denominator"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pg-matrix", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--run-col", default="sample_id")
    ap.add_argument("--outdir", required=True)

    ap.add_argument("--n_pcs_int", type=int, default=20)
    ap.add_argument("--n_pcs_det", type=int, default=10)
    ap.add_argument("--n_neighbors", type=int, default=15)
    ap.add_argument("--leiden_resolution", type=float, default=0.8)
    ap.add_argument("--min_cluster_size", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)

    ap.add_argument(
        "--batch_correct",
        choices=["auto", "on", "off"],
        default="auto",
        help="Harmony correction on joint PCs using Library (if available).",
    )
    ap.add_argument("--harmony_theta", type=float, default=2.0)
    ap.add_argument("--harmony_theta_max", type=float, default=5.0,
                    help="Maximum theta for adaptive Harmony. 0 = disable adaptive loop.")
    ap.add_argument("--harmony_entropy_target", type=float, default=0.6,
                    help="Target mean normalised batch entropy (0-1). "
                         "Below this, theta is increased iteratively.")

    # Plot controls (optional override; otherwise adaptive sizing is used)
    ap.add_argument(
        "--plot_umap_size",
        type=float,
        default=None,
        help="Override dot size for UMAP/PCA (else auto).",
    )

    args = ap.parse_args()


    np.random.seed(args.seed)

    # Heavy imports
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import scanpy as sc

    outdir = Path(args.outdir)
    qc_dir = outdir / "qc"
    outdir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    # ---- load manifest ----
    man = pd.read_csv(args.manifest, sep=None, engine="python", dtype=str).fillna("")
    if args.run_col not in man.columns:
        raise SystemExit(f"Manifest missing run-col '{args.run_col}'. Columns: {man.columns.tolist()}")

    man["Run"] = man[args.run_col].astype(str).str.strip()
    man = man[man["Run"] != ""].copy()

    lib_src, lib = detect_library_column(man)
    man["Library"] = lib.astype(str)

    # ---- load pg_matrix ----
    mat = pd.read_csv(args.pg_matrix, sep="\t", dtype=str)
    header_map = {c: stem_from_header(c) for c in mat.columns}

    wanted_runs = set(man["Run"].tolist())
    all_sample_cols = [c for c in mat.columns if is_sample_column(c)]
    run_cols = [c for c in all_sample_cols if header_map.get(c, "") in wanted_runs]

    if not run_cols:
        raise SystemExit("No matching run columns found between pg_matrix and manifest.")

    # Reorder columns to manifest order (critical for consistent output)
    col_to_run = {c: header_map[c] for c in run_cols}
    run_to_col = {}
    for c in run_cols:
        run_to_col.setdefault(col_to_run[c], c)

    runs = [r for r in man["Run"].tolist() if r in run_to_col]
    cols = [run_to_col[r] for r in runs]

    # pivot_pack.tsv is already log2 + median-normalised; NaN = missing
    X_norm = mat[cols].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float, copy=True)

    n_prot, n_cells = X_norm.shape
    print(f"SCP clustering: {n_prot} proteins × {n_cells} cells (pre-normalised pivot)")

    # QC metrics to carry into obs
    present = np.isfinite(X_norm)
    n_detected = present.sum(axis=0)
    missing_frac = 1.0 - (n_detected / max(n_prot, 1))

    # ---- Build PCA blocks ----
    # Intensity PCA: impute protein-wise medians for PCA only
    X_int = X_norm.copy()  # proteins x cells
    prot_meds = np.nanmedian(X_int, axis=1)
    prot_meds = np.where(np.isnan(prot_meds), 0.0, prot_meds)
    miss = ~np.isfinite(X_int)
    if miss.any():
        X_int[miss] = prot_meds[np.where(miss)[0]]

    # z-score proteins (features) across cells for PCA
    X_int = X_int - X_int.mean(axis=1, keepdims=True)
    sd = X_int.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    X_int = X_int / sd

    # Detection matrix: proteins x cells binary
    X_det = present.astype(float)
    # centre proteins
    X_det = X_det - X_det.mean(axis=1, keepdims=True)

    # PCA in sample space (cells x proteins)
    pc_int, sv_int = pca_svd(X_int.T, n_components=min(args.n_pcs_int, n_cells - 1))
    pc_det, sv_det = pca_svd(X_det.T, n_components=min(args.n_pcs_det, n_cells - 1))

    # balance blocks
    pc_int = standardise_block(pc_int)
    pc_det = standardise_block(pc_det)

    Z = np.concatenate([pc_int, pc_det], axis=1)
    Z = standardise_block(Z)

    # ---- optional Harmony batch correction (adaptive theta) ----
    batch_applied = False
    Z_use = Z
    libs = man.set_index("Run").loc[runs, "Library"].astype(str).tolist()
    n_libs = len(set(libs))

    # Save pre-Harmony embedding for diagnostics
    Z_precorrect = Z.copy()

    harmony_log = []  # track theta iterations

    def _orient_Z(Z_corr, Z_ref):
        """Ensure Z_corr is (cells x dims) matching Z_ref shape."""
        if Z_corr.shape == (Z_ref.shape[1], Z_ref.shape[0]):
            return Z_corr.T
        if Z_corr.shape == Z_ref.shape:
            return Z_corr
        raise ValueError(
            f"Unexpected Harmony Z_corr shape {Z_corr.shape}; "
            f"expected {Z_ref.shape} or {(Z_ref.shape[1], Z_ref.shape[0])}"
        )

    def _run_harmony_and_cluster(Z_in, theta):
        """Run Harmony → neighbors → UMAP → Leiden. Returns (Z_corrected, adata, clusters)."""
        import harmonypy as hm
        ho = hm.run_harmony(
            Z_in, pd.DataFrame({"Library": libs}), "Library",
            theta=theta, max_iter_harmony=50,
        )
        Z_out = _orient_Z(ho.Z_corr, Z_in)
        if Z_out.shape[0] != len(runs) and Z_out.shape[1] == len(runs):
            Z_out = Z_out.T

        ad = sc.AnnData(X=Z_out)
        ad.obs_names = pd.Index(runs)
        ad.obs["Library"] = pd.Categorical(libs)
        n_neigh = min(args.n_neighbors, max(3, n_cells - 1))
        sc.pp.neighbors(ad, n_neighbors=n_neigh, use_rep="X", random_state=args.seed)
        sc.tl.umap(ad, random_state=args.seed)
        sc.tl.leiden(ad, resolution=args.leiden_resolution, random_state=args.seed,
                     flavor="igraph", n_iterations=2, directed=False)
        cl = ["C" + c for c in ad.obs["leiden"].astype(str).tolist()]
        return Z_out, ad, cl

    if args.batch_correct in ("auto", "on") and n_libs > 1:
        try:
            theta = args.harmony_theta
            theta_max = args.harmony_theta_max
            entropy_target = args.harmony_entropy_target
            adaptive = theta_max > theta and entropy_target > 0

            Z_use, adata_tmp, clusters_tmp = _run_harmony_and_cluster(Z, theta)
            mean_ent, per_cl_ent = batch_mixing_entropy(clusters_tmp, libs)
            harmony_log.append({"theta": theta, "entropy": mean_ent, "n_clusters": len(set(clusters_tmp))})
            print(f"Harmony theta={theta:.1f}: batch mixing entropy={mean_ent:.3f} "
                  f"(target={entropy_target:.2f}, clusters={len(set(clusters_tmp))})")

            # Adaptive loop: increase theta if mixing is insufficient
            while adaptive and mean_ent < entropy_target and theta < theta_max:
                theta = min(theta + 1.0, theta_max)
                Z_use, adata_tmp, clusters_tmp = _run_harmony_and_cluster(Z, theta)
                mean_ent, per_cl_ent = batch_mixing_entropy(clusters_tmp, libs)
                harmony_log.append({"theta": theta, "entropy": mean_ent, "n_clusters": len(set(clusters_tmp))})
                print(f"Harmony theta={theta:.1f}: batch mixing entropy={mean_ent:.3f} "
                      f"(clusters={len(set(clusters_tmp))})")

            batch_applied = True

            # Flag batch-dominated clusters (per-cluster entropy < 0.3)
            batch_dominated = {cl: ent for cl, ent in per_cl_ent.items() if ent < 0.3}
            if batch_dominated:
                print(f"\nWARNING: {len(batch_dominated)} cluster(s) remain batch-dominated "
                      f"after theta={theta:.1f}:")
                cl_batch_df = pd.DataFrame({"cluster": clusters_tmp, "batch": libs})
                for cl, ent in sorted(batch_dominated.items()):
                    top_batch = (cl_batch_df[cl_batch_df["cluster"] == cl]["batch"]
                                 .value_counts().head(1))
                    n_in_cl = (cl_batch_df["cluster"] == cl).sum()
                    print(f"  {cl}: entropy={ent:.3f}, n={n_in_cl}, "
                          f"top batch={top_batch.index[0]} ({top_batch.iloc[0]}/{n_in_cl} = "
                          f"{top_batch.iloc[0]/n_in_cl:.0%})")
                print("  These may be real batch-specific biology or technical artifacts.")
                print("  Review harmony_before_after.pdf to decide.\n")
            else:
                print(f"All clusters well-mixed (per-cluster entropy >= 0.3).")

            # Store in log for report JSON
            harmony_log.append({
                "final_theta": theta,
                "batch_dominated_clusters": {
                    cl: {"entropy": ent} for cl, ent in batch_dominated.items()
                } if batch_dominated else {},
            })

            print(f"Harmony applied (final theta={theta:.1f}, entropy={mean_ent:.3f}, "
                  f"n_libraries={n_libs}).")

        except Exception as e:
            if args.batch_correct == "on":
                raise
            print(f"WARNING: Harmony not applied (harmonypy missing or failed): {e}")
            adata_tmp = None
            clusters_tmp = None
    else:
        print(f"Batch correction skipped (mode={args.batch_correct}, n_libraries={n_libs}).")
        adata_tmp = None
        clusters_tmp = None

    # Final orientation guard
    if Z_use.shape[0] != len(runs) and Z_use.shape[1] == len(runs):
        Z_use = Z_use.T
    if Z_use.shape[0] != len(runs):
        raise ValueError(f"Z_use rows ({Z_use.shape[0]}) must equal n_cells ({len(runs)}).")

    # ---- Build AnnData (reuse from adaptive loop or create fresh) ----
    if adata_tmp is not None:
        adata = adata_tmp
    else:
        adata = sc.AnnData(X=Z_use)
        adata.obs_names = pd.Index(runs)
        n_neighbors = min(args.n_neighbors, max(3, n_cells - 1))
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X", random_state=args.seed)
        sc.tl.umap(adata, random_state=args.seed)
        sc.tl.leiden(adata, resolution=args.leiden_resolution, random_state=args.seed,
                     flavor="igraph", n_iterations=2, directed=False)
        clusters_tmp = ["C" + c for c in adata.obs["leiden"].astype(str).tolist()]

    adata.obs["Run"] = runs
    if "Library" not in adata.obs.columns:
        adata.obs["Library"] = pd.Categorical(libs)
    adata.obs["n_detected_proteins"] = n_detected
    adata.obs["missing_fraction"] = missing_frac

    raw_clusters = adata.obs["leiden"].astype(str).tolist()
    clusters = ["C" + c for c in raw_clusters]
    adata.obs["cluster_id"] = clusters

    # Merge tiny clusters
    sizes = pd.Series(clusters).value_counts()
    small = set(sizes[sizes < args.min_cluster_size].index.tolist())
    if small:
        adata.obs.loc[adata.obs["cluster_id"].isin(small), "cluster_id"] = "C_other"

    # Baseline = largest cluster
    final_sizes = adata.obs["cluster_id"].value_counts()
    baseline = final_sizes.index[0]

    # Contrasts = each cluster vs baseline
    conds = sorted(adata.obs["cluster_id"].unique().tolist())
    contrast_rows = []
    for c in conds:
        if c == baseline:
            continue
        contrast_rows.append({
            "Contrast": f"inferred_{c}_vs_{baseline}",
            "Numerator": c,
            "Denominator": baseline,
        })
    contrasts_df = pd.DataFrame(contrast_rows, columns=["Contrast", "Numerator", "Denominator"])

    # Annotation table (bulk-compatible)
    # Normalise batch: use Library values but cap batch levels relative to cell count.
    # This prevents unique-per-cell batch values (e.g. run IDs) from making the
    # scplainer design matrix underdetermined, while allowing genuine acquisition batches.
    _libs_series = pd.Series(libs)
    _n_libs = _libs_series.nunique()
    # Allow up to n_cells/10 batch levels (minimum 20), ensuring enough cells per batch
    max_batch_levels = max(20, n_cells // 10)
    if _n_libs > max_batch_levels:
        batch_vals = ["1"] * len(libs)
    else:
        batch_vals = libs
    ann_out = pd.DataFrame({
        "Run": runs,
        "Condition": adata.obs["cluster_id"].tolist(),
        "BioRep": (adata.obs.groupby("cluster_id").cumcount() + 1).tolist(),
        "Batch": batch_vals,
    })

    # Compute UMAP on pre-Harmony (uncorrected) embedding for diagnostics
    if batch_applied:
        import scanpy as _sc2
        _adata_pre = _sc2.AnnData(X=Z_precorrect)
        _adata_pre.obs_names = pd.Index(runs)
        _sc2.pp.neighbors(_adata_pre, n_neighbors=min(args.n_neighbors, max(3, n_cells - 1)),
                          use_rep="X", random_state=args.seed)
        _sc2.tl.umap(_adata_pre, random_state=args.seed)
        um_pre = _adata_pre.obsm["X_umap"]
    else:
        um_pre = adata.obsm["X_umap"]  # same as corrected when no batch correction

    # Assignments — include pre-Harmony UMAP coords for diagnostic comparison
    um = adata.obsm["X_umap"]
    assign = pd.DataFrame({
        "Run": runs,
        "Library": libs,
        "cluster_id": adata.obs["cluster_id"].tolist(),
        "cluster_id_raw": raw_clusters,
        "umap_1": um[:, 0],
        "umap_2": um[:, 1],
        "umap_1_uncorrected": um_pre[:, 0],
        "umap_2_uncorrected": um_pre[:, 1],
        "n_detected_proteins": n_detected,
        "missing_fraction": missing_frac,
    })

    # ---- write outputs ----
    assign.to_csv(outdir / "scp_cluster_assignments.tsv", sep="\t", index=False)
    ann_out.to_csv(outdir / "scp_annotation.tsv", sep="\t", index=False)
    contrasts_df.to_csv(outdir / "scp_contrasts.tsv", sep="\t", index=False)

    report = {
        "tool": "scp_joint_embedding_leiden",
        "status": "OK",
        "inputs": {"pg_matrix": args.pg_matrix, "manifest": args.manifest},
        "n_proteins": int(n_prot),
        "n_cells": int(n_cells),
        "params": {
            "n_pcs_int": args.n_pcs_int,
            "n_pcs_det": args.n_pcs_det,
            "n_neighbors": min(args.n_neighbors, max(3, n_cells - 1)),
            "leiden_resolution": args.leiden_resolution,
            "min_cluster_size": args.min_cluster_size,
            "seed": args.seed,
            "batch_correct": args.batch_correct,
            "harmony_theta": args.harmony_theta,
        },
        "library_source_column": lib_src,
        "n_libraries": int(n_libs),
        "batch_applied": bool(batch_applied),
        "harmony_adaptive_log": harmony_log,
        "cluster_sizes": final_sizes.to_dict(),
        "baseline_cluster": baseline,
        "n_contrasts": int(len(contrasts_df)),
        "small_clusters_merged": sorted(list(small)),
    }
    (outdir / "scp_clustering_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")

    # ---- QC plots (PDF, square, adaptive dot size) ----
    import matplotlib.pyplot as plt

    # PDF-only
    plot_ext = "pdf"

    # Square figure size (inches)
    figsize_sq = (6.0, 6.0)

    # Adaptive dot size (half-ish of previous)
    s_auto = adaptive_dot_size(n_cells, s_max_at_12=50.0, s_min=12.0)
    s_use = float(args.plot_umap_size) if args.plot_umap_size is not None else s_auto

    # Use a larger size for scanpy UMAP than matplotlib scatter
    umap_size = s_use * 4.0  # tune 2.5–4.0 if needed

    # UMAP by cluster
    plt.figure(figsize=figsize_sq)
    sc.pl.umap(adata, color="cluster_id", show=False, size=umap_size)
    plt.savefig(qc_dir / f"umap_by_cluster.{plot_ext}", bbox_inches="tight")
    plt.close()

    # UMAP by library
    plt.figure(figsize=figsize_sq)
    sc.pl.umap(adata, color="Library", show=False, size=umap_size)
    plt.savefig(qc_dir / f"umap_by_library.{plot_ext}", bbox_inches="tight")
    plt.close()

    # UMAP by detected proteins
    plt.figure(figsize=figsize_sq)
    sc.pl.umap(adata, color="n_detected_proteins", show=False, size=umap_size)
    plt.savefig(qc_dir / f"umap_by_n_detected.{plot_ext}", bbox_inches="tight")
    plt.close()

    # PCA intensity block (PC1/PC2) - square, adaptive dots
    fig, ax = plt.subplots(figsize=figsize_sq)
    ax.scatter(pc_int[:, 0], pc_int[:, 1], s=s_use, alpha=0.80, linewidths=0)
    ax.set_title("Intensity PCA block (PC1/PC2)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_aspect("equal", adjustable="datalim")
    plt.savefig(qc_dir / f"pca_intensity.{plot_ext}", bbox_inches="tight")
    plt.close(fig)

    # PCA detection block
    fig, ax = plt.subplots(figsize=figsize_sq)
    ax.scatter(pc_det[:, 0], pc_det[:, 1], s=s_use, alpha=0.80, linewidths=0)
    ax.set_title("Detection PCA block (PC1/PC2)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_aspect("equal", adjustable="datalim")
    plt.savefig(qc_dir / f"pca_detection.{plot_ext}", bbox_inches="tight")
    plt.close(fig)


    print(f"Done. Clusters={len(final_sizes)} baseline={baseline}; outputs -> {outdir}")


if __name__ == "__main__":
    main()
