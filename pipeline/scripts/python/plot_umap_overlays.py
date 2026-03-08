#!/usr/bin/env python3
"""
plot_umap_overlays.py

Generate a comprehensive set of UMAP overlay plots for SCP data.

Inputs
------
--assignments   scp_cluster_assignments.tsv
--pivot         scp_pivot.tsv (log2, NA for missing)
--annotation    scp_annotation.tsv (Run, Condition)
--aucell        aucell_scores.tsv  -- optional
--detection-markers   detection_markers.tsv -- optional, guides auto protein selection
--proteins      comma-separated gene names to overlay, e.g. INS,GCG,SST,KRT19
                  OR path to a plain-text file with one gene per line
--outdir        output directory
--top-n-proteins  number of auto-selected proteins (default 9). Ignored if
                  --proteins is supplied and enough matches are found.
--top-n-pathways  number of AUCell pathways (default 9)

Outputs
-------
plots/
  umap_qc_grid.pdf
  umap_batch.pdf
  umap_top_proteins.pdf        -- auto-selected (MAD or detection markers)
  umap_custom_proteins.pdf     -- only when --proteins supplied
  umap_aucell_pathways.pdf     -- only when --aucell supplied
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
    "Protein.Q.Value", "Global.Q.Value", "Protein.Descriptions", "Pathway",
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
    # Anything not in the known metadata set is a sample column — covers
    # file-like headers (with extensions/paths) AND clean stem names from
    # SCP pivots (e.g. "SC_01_s1-a1_1_7515").
    return True


def build_cell_df(assign_df: pd.DataFrame, ann_df: pd.DataFrame) -> pd.DataFrame:
    df = assign_df.copy()
    if "Condition" not in df.columns:
        df = df.merge(ann_df[["Run", "Condition"]].drop_duplicates(), on="Run", how="left")
    return df


def load_pivot_matrix(pivot_path: str, cell_df: pd.DataFrame):
    pivot = pd.read_csv(pivot_path, sep="\t", dtype=str)
    header_map = {c: stem_from_header(c) for c in pivot.columns}
    sample_cols = [c for c in pivot.columns if is_sample_column(c)]

    run_set = set(cell_df["Run"].astype(str))
    matched_cols = [c for c in sample_cols if header_map.get(c, "") in run_set]
    col_to_run = {c: header_map[c] for c in matched_cols}
    run_to_col = {v: k for k, v in col_to_run.items()}

    runs_ordered = cell_df["Run"].astype(str).tolist()
    cols_ordered = [run_to_col.get(r) for r in runs_ordered]
    valid = [(i, c) for i, c in enumerate(cols_ordered) if c is not None]
    if not valid:
        return None, [], []

    _, valid_cols = zip(*valid)
    X = pivot[list(valid_cols)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    prot_col = "Protein.Group" if "Protein.Group" in pivot.columns else pivot.columns[0]
    gene_col = "Genes" if "Genes" in pivot.columns else None
    if gene_col:
        protein_ids = (pivot["Genes"].fillna("").astype(str) + "|" +
                       pivot[prot_col].astype(str)).tolist()
    else:
        protein_ids = pivot[prot_col].astype(str).tolist()

    return X, protein_ids, runs_ordered


def parse_protein_list(proteins_arg: str) -> list[str]:
    """
    Accept either a comma-separated string of gene names
    or a path to a file with one gene per line.
    Returns a list of stripped, non-empty strings.
    """
    p = Path(proteins_arg)
    if p.exists() and p.is_file():
        genes = [l.strip() for l in p.read_text(encoding="utf-8").splitlines()]
    else:
        genes = [g.strip() for g in proteins_arg.split(",")]
    return [g for g in genes if g]


def find_protein_indices(wanted_genes: list[str], protein_ids: list[str]) -> tuple[list[int], list[str]]:
    """
    Match gene names (case-insensitive) against the gene|protein_id labels.
    Returns (indices, display_labels) in the order of wanted_genes.
    Missing genes are reported but skipped.
    """
    indices = []
    labels_out = []
    pid_lower = [p.split("|")[0].strip().lower() for p in protein_ids]

    for gene in wanted_genes:
        gl = gene.lower()
        # Exact match first
        hits = [i for i, p in enumerate(pid_lower) if p == gl]
        if not hits:
            # Prefix match (handles isoforms like INS1, INS2)
            hits = [i for i, p in enumerate(pid_lower) if p.startswith(gl)]
        if hits:
            indices.append(hits[0])
            labels_out.append(gene)
        else:
            print(f"  WARNING: gene '{gene}' not found in pivot — skipping")

    return indices, labels_out


def pick_top_proteins_auto(X: np.ndarray, protein_ids: list[str],
                           n: int, marker_path: str = None) -> list[int]:
    """Auto-select top N by detection markers or MAD."""
    if marker_path and Path(marker_path).exists():
        dm = pd.read_csv(marker_path, sep="\t", dtype=str)
        gene_col = "Genes" if "Genes" in dm.columns else None
        id_col = "Protein" if "Protein" in dm.columns else dm.columns[0]
        if gene_col and "det_rate_in" in dm.columns:
            dm["det_rate_in"] = pd.to_numeric(dm["det_rate_in"], errors="coerce")
            top_genes = (dm.sort_values("det_rate_in", ascending=False)
                         .drop_duplicates(subset=[id_col])
                         .head(n)[gene_col].astype(str).tolist())
            indices = []
            for g in top_genes:
                for i, pid in enumerate(protein_ids):
                    if pid.startswith(g + "|") or pid == g:
                        indices.append(i)
                        break
            if len(indices) >= n // 2:
                return indices[:n]

    mads = np.nanmedian(np.abs(X - np.nanmedian(X, axis=1, keepdims=True)), axis=1)
    mads = np.where(np.isfinite(mads), mads, 0.0)
    return np.argsort(mads)[::-1][:n].tolist()


def clean_pathway_name(name: str, max_len: int = 35) -> str:
    import re
    n = re.sub(r"^(HALLMARK_|REACTOME_|GOBP_|GOCC_|GOMF_|KEGG_)", "", name)
    n = n.replace("_", " ").lower()
    n = n[0].upper() + n[1:] if n else n
    return n[:max_len - 3] + "..." if len(n) > max_len else n


def scatter_continuous(ax, x, y, values, label: str, cmap: str = "viridis",
                       s: float = 20, absent_marker: bool = False):
    """
    Plot continuous values on a UMAP.
    absent_marker=True: missing/NaN cells drawn as grey × markers so they
    are unambiguously distinct from low-but-detected values.
    absent_marker=False: missing cells drawn as small filled grey dots
    (appropriate for QC metrics where absence just means low/zero).
    """
    import matplotlib.pyplot as plt
    finite = np.isfinite(values)
    n_absent = int((~finite).sum())

    if absent_marker and n_absent > 0:
        # Cross/x markers — clearly "not detected", not confused with any colormap colour
        ax.scatter(x[~finite], y[~finite],
                   marker="x", c="#aaaaaa", s=max(s * 0.6, 8),
                   alpha=0.5, linewidths=0.8, zorder=1,
                   label=f"not detected (n={n_absent})")
    elif n_absent > 0:
        ax.scatter(x[~finite], y[~finite],
                   c="#cccccc", s=s * 0.5, alpha=0.5, linewidths=0, zorder=1)

    sc = ax.scatter(x[finite], y[finite], c=values[finite],
                    cmap=cmap, s=s, alpha=0.85, linewidths=0, zorder=2)
    plt.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)

    if absent_marker and n_absent > 0:
        ax.legend(fontsize=6, loc="lower right", framealpha=0.7,
                  markerscale=1.0, handlelength=1.0)

    ax.set_title(label, fontsize=9)
    ax.set_xlabel("UMAP1", fontsize=8)
    ax.set_ylabel("UMAP2", fontsize=8)
    ax.set_aspect("equal", adjustable="datalim")
    ax.tick_params(labelsize=7)


def scatter_categorical(ax, x, y, labels, title: str, s: float = 20,
                        centroid_labels: bool = False):
    import matplotlib.pyplot as plt
    cats = sorted(set(str(l) for l in labels))
    cmap = plt.get_cmap("tab10") if len(cats) <= 10 else plt.get_cmap("tab20")
    colors = {c: cmap(i % cmap.N) for i, c in enumerate(cats)}
    for cat in cats:
        mask = np.array([str(l) == cat for l in labels])
        ax.scatter(x[mask], y[mask], c=[colors[cat]], s=s, alpha=0.85,
                   linewidths=0, label=cat, zorder=2)
        if centroid_labels:
            cx, cy = np.median(x[mask]), np.median(y[mask])
            ax.text(cx, cy, cat, fontsize=6, fontweight="bold",
                    ha="center", va="center", zorder=3,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec="none", alpha=0.7))
    ax.legend(fontsize=7, markerscale=1.2, framealpha=0.7,
              loc="best", ncol=max(1, len(cats) // 8))
    ax.set_title(title, fontsize=9)
    ax.set_xlabel("UMAP1", fontsize=8)
    ax.set_ylabel("UMAP2", fontsize=8)
    ax.set_aspect("equal", adjustable="datalim")
    ax.tick_params(labelsize=7)


def plot_protein_grid(X, protein_indices, display_labels, umap1, umap2, s, fmt, out_path):
    import matplotlib.pyplot as plt
    n_top = len(protein_indices)
    if n_top == 0:
        print(f"  No proteins to plot for {out_path.name}")
        return
    ncols = min(3, n_top)
    nrows = int(np.ceil(n_top / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 5.0 * nrows), squeeze=False)

    for k, (idx, lbl) in enumerate(zip(protein_indices, display_labels)):
        r, c = divmod(k, ncols)
        scatter_continuous(axes[r][c], umap1, umap2, X[idx, :],
                           lbl[:30], cmap="YlOrRd", s=s, absent_marker=True)
    for k in range(n_top, nrows * ncols):
        r, c = divmod(k, ncols)
        axes[r][c].set_visible(False)

    fig.tight_layout(pad=1.0)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--assignments",       required=True)
    ap.add_argument("--pivot",             required=True)
    ap.add_argument("--annotation",        required=True)
    ap.add_argument("--shifts",            default=None,
                    help="pivot_shifts.tsv with per-run median normalisation shifts")
    ap.add_argument("--corrected-pivot",   default=None, dest="corrected_pivot",
                    help="Batch-corrected expression TSV from scplainer. "
                         "When provided, protein overlays use these values instead of pivot_pack.")
    ap.add_argument("--aucell",            default=None)
    ap.add_argument("--detection-markers", default=None)
    ap.add_argument("--proteins",          default=None,
                    help="Comma-separated gene names (e.g. INS,GCG,SST) "
                         "or path to a text file with one gene per line. "
                         "Produces umap_custom_proteins.pdf in addition to "
                         "the auto-selected overlay.")
    ap.add_argument("--custom-proteins",   default=None,
                    help="Additional comma-separated gene names or path to "
                         "a text file. Merged with --proteins (deduplicated) "
                         "into umap_custom_proteins.pdf.")
    ap.add_argument("--outdir",            required=True)
    ap.add_argument("--top-n-proteins",    type=int, default=9)
    ap.add_argument("--top-n-pathways",    type=int, default=9)
    ap.add_argument("--cell-type-annotations", default=None,
                    help="Path to cluster_cell_type_annotations.tsv from LLM. "
                         "Produces umap_cell_types.pdf with centroid labels.")
    ap.add_argument("--dot-size",          type=float, default=None)
    ap.add_argument("--format",            default="pdf", choices=["pdf", "png"])
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = Path(args.outdir) / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    # ---- Load core data ----
    assign   = pd.read_csv(args.assignments, sep="\t", dtype=str)
    ann      = pd.read_csv(args.annotation,  sep="\t", dtype=str)
    cell_df  = build_cell_df(assign, ann)

    umap1   = pd.to_numeric(cell_df["umap_1"], errors="coerce").to_numpy()
    umap2   = pd.to_numeric(cell_df["umap_2"], errors="coerce").to_numpy()
    n_cells = len(cell_df)
    s = args.dot_size if args.dot_size else max(6.0, min(36.0, 36.0 * np.sqrt(12.0 / n_cells)))

    # ---- Load pivot ----
    X, protein_ids, _ = load_pivot_matrix(args.pivot, cell_df)
    total_intensity = np.nansum(X, axis=0) if X is not None else None

    # Batch-corrected expression for protein overlays (if available)
    X_overlay = X  # default: use raw pivot
    overlay_ids = protein_ids
    if args.corrected_pivot and Path(args.corrected_pivot).exists():
        X_corr, corr_ids, _ = load_pivot_matrix(args.corrected_pivot, cell_df)
        if X_corr is not None and len(corr_ids) > 0:
            X_overlay = X_corr
            overlay_ids = corr_ids
            print(f"Using batch-corrected expression for protein overlays "
                  f"({X_corr.shape[0]} proteins)")
        else:
            print("WARNING: corrected pivot loaded but empty, falling back to pivot_pack")

    # ---- Plot 1: QC overlays ----
    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes = axes.flatten()

    n_det     = pd.to_numeric(cell_df.get("n_detected_proteins",
                              pd.Series(np.nan, index=cell_df.index)), errors="coerce").to_numpy()
    miss      = pd.to_numeric(cell_df.get("missing_fraction",
                              pd.Series(np.nan, index=cell_df.index)), errors="coerce").to_numpy()

    # Load median normalisation shifts from pivot_shifts.tsv (produced by scp_pivot)
    med_shift = np.full(n_cells, np.nan)
    if args.shifts and Path(args.shifts).exists():
        shifts_df = pd.read_csv(args.shifts, sep="\t", dtype=str)
        shift_map = dict(zip(shifts_df["RunStem"].astype(str),
                             pd.to_numeric(shifts_df["median_shift_applied"], errors="coerce")))
        med_shift = np.array([shift_map.get(str(r), np.nan) for r in cell_df["Run"]])

    scatter_continuous(axes[0], umap1, umap2, n_det,     "n detected proteins",           cmap="viridis",  s=s)
    scatter_continuous(axes[1], umap1, umap2, miss,      "missing fraction",              cmap="Reds",     s=s)
    scatter_continuous(axes[2], umap1, umap2, med_shift, "median normalisation shift",    cmap="coolwarm", s=s)
    if total_intensity is not None and len(total_intensity) == n_cells:
        scatter_continuous(axes[3], umap1, umap2, total_intensity, "total log2 intensity", cmap="plasma", s=s)
    else:
        axes[3].set_visible(False)

    fig.tight_layout(pad=1.5)
    fig.savefig(outdir / f"umap_qc_grid.{args.format}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote umap_qc_grid.{args.format}")

    # ---- Plot 2: Batch/library ----
    fig, ax = plt.subplots(figsize=(6, 6))
    lib_labels = cell_df.get("Library", cell_df.get("Batch", pd.Series(["unknown"] * n_cells)))
    scatter_categorical(ax, umap1, umap2, lib_labels.tolist(), "Batch / Library", s=s)
    fig.tight_layout()
    fig.savefig(outdir / f"umap_batch.{args.format}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote umap_batch.{args.format}")

    # ---- Plot 2b: Cluster UMAP (same format as overlays) ----
    fig, ax = plt.subplots(figsize=(6, 6))
    cluster_labels = cell_df.get("cluster_id",
                     cell_df.get("Condition", pd.Series(["unknown"] * n_cells)))
    scatter_categorical(ax, umap1, umap2, cluster_labels.tolist(), "Clusters", s=s)
    fig.tight_layout()
    fig.savefig(outdir / f"umap_clusters.{args.format}", bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote umap_clusters.{args.format}")

    # ---- Plot 3a: Auto-selected protein overlays ----
    if X_overlay is not None and len(overlay_ids) > 0:
        auto_idx = pick_top_proteins_auto(X_overlay, overlay_ids, args.top_n_proteins,
                                          args.detection_markers)
        auto_labels = []
        for idx in auto_idx:
            lbl = overlay_ids[idx].split("|")[0] if "|" in overlay_ids[idx] else overlay_ids[idx]
            auto_labels.append(lbl[:30])

        plot_protein_grid(X_overlay, auto_idx, auto_labels, umap1, umap2, s,
                          args.format, outdir / f"umap_top_proteins.{args.format}")

        # ---- Plot 3b: Custom protein overlays ----
        combined_proteins = []
        if args.proteins:
            combined_proteins += parse_protein_list(args.proteins)
        if args.custom_proteins:
            combined_proteins += parse_protein_list(args.custom_proteins)
        # Deduplicate while preserving order
        seen = set()
        wanted = []
        for g in combined_proteins:
            gl = g.lower()
            if gl not in seen:
                seen.add(gl)
                wanted.append(g)
        if wanted:
            print(f"Custom protein list: {wanted}")
            custom_idx, custom_labels = find_protein_indices(wanted, overlay_ids)
            print(f"  Matched {len(custom_idx)}/{len(wanted)} genes in pivot")
            plot_protein_grid(X_overlay, custom_idx, custom_labels, umap1, umap2, s,
                              args.format, outdir / f"umap_custom_proteins.{args.format}")

    # ---- Plot 4: AUCell pathway overlays ----
    if args.aucell and Path(args.aucell).exists():
        auc = pd.read_csv(args.aucell, sep="\t", dtype=str)
        pathway_col = "Pathway" if "Pathway" in auc.columns else auc.columns[0]
        auc_sample_cols = [c for c in auc.columns if c != pathway_col]
        auc_stem_map = {c: stem_from_header(c) for c in auc_sample_cols}
        run_set = set(cell_df["Run"].astype(str))
        matched_auc = [c for c in auc_sample_cols if auc_stem_map.get(c, "") in run_set]

        if matched_auc:
            run_to_auc = {auc_stem_map[c]: c for c in matched_auc}
            runs_ordered = cell_df["Run"].astype(str).tolist()
            auc_cols_ord = [run_to_auc.get(r) for r in runs_ordered]
            valid_auc = [(i, c) for i, c in enumerate(auc_cols_ord) if c is not None]

            if valid_auc:
                _, valid_auc_cols = zip(*valid_auc)
                A = auc[list(valid_auc_cols)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
                mads = np.nanmedian(np.abs(A - np.nanmedian(A, axis=1, keepdims=True)), axis=1)
                mads = np.where(np.isfinite(mads), mads, 0.0)
                top_path_idx = np.argsort(mads)[::-1][:args.top_n_pathways].tolist()

                n_top = len(top_path_idx)
                if n_top == 0:
                    print("  top-n-pathways=0, skipping AUCell overlay")
                else:
                    ncols = min(3, n_top)
                    nrows = int(np.ceil(n_top / ncols))
                    fig, axes = plt.subplots(nrows, ncols,
                                             figsize=(5.5 * ncols, 5.0 * nrows), squeeze=False)

                    for k, idx in enumerate(top_path_idx):
                        r, c = divmod(k, ncols)
                        pname = clean_pathway_name(str(auc.iloc[idx][pathway_col]))
                        scatter_continuous(axes[r][c], umap1, umap2, A[idx, :],
                                           pname, cmap="OrRd", s=s)
                    for k in range(n_top, nrows * ncols):
                        r, c = divmod(k, ncols)
                        axes[r][c].set_visible(False)

                    fig.tight_layout(pad=1.0)
                    fig.savefig(outdir / f"umap_aucell_pathways.{args.format}", bbox_inches="tight")
                    plt.close(fig)
                    print(f"Wrote umap_aucell_pathways.{args.format}")

    # ---- Plot 5: Cell-type annotation UMAP ----
    if args.cell_type_annotations and Path(args.cell_type_annotations).exists():
        ct = pd.read_csv(args.cell_type_annotations, sep="\t", dtype=str)
        # Expect columns: cluster_id (or Cluster), cell_type
        ct_cluster_col = "cluster_id" if "cluster_id" in ct.columns else "Cluster"
        ct_type_col = "cell_type" if "cell_type" in ct.columns else ct.columns[1]
        ct_map = dict(zip(ct[ct_cluster_col].astype(str), ct[ct_type_col].astype(str)))

        cl_ids = cell_df.get("cluster_id",
                             cell_df.get("Condition", pd.Series(["unknown"] * n_cells)))
        ct_labels = [f"{c}: {ct_map.get(str(c), '?')}" for c in cl_ids]
        ct_labels = [l[:40] + "..." if len(l) > 40 else l for l in ct_labels]

        fig, ax = plt.subplots(figsize=(8, 8))
        scatter_categorical(ax, umap1, umap2, ct_labels, "Cell types (LLM)",
                            s=s, centroid_labels=True)
        fig.tight_layout()
        fig.savefig(outdir / f"umap_cell_types.{args.format}", bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote umap_cell_types.{args.format}")

    print(f"UMAP overlays complete -> {outdir}")


if __name__ == "__main__":
    main()
