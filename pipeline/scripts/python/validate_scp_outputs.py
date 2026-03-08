#!/usr/bin/env python3
"""
validate_scp_outputs.py

Automated validation of single-cell proteomics pipeline outputs.
Reads a workdir (containing scp/ subdirectory) and checks every SCP output
for file existence, numerical sanity, cross-table consistency, statistical
validity, batch correction quality, and LLM annotation quality.

Returns pass/fail with details.

Usage
-----
    python validate_scp_outputs.py --workdir /path/to/dataset
    python validate_scp_outputs.py --workdir /path/to/dataset --known-markers markers.json
    python validate_scp_outputs.py --workdir /path/to/dataset --strict
"""

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class CheckResult:
    """Single validation check result."""
    def __init__(self, name: str, passed: bool, detail: str = "",
                 severity: str = "FAIL"):
        self.name = name
        self.passed = passed
        self.detail = detail
        self.severity = severity  # FAIL or WARN

    def __repr__(self):
        tag = "PASS" if self.passed else self.severity
        s = f"[{tag}] {self.name}"
        if self.detail:
            s += f"  -- {self.detail}"
        return s


def load_tsv(path: Path) -> pd.DataFrame | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    return pd.read_csv(path, sep="\t", dtype=str)


def to_float_col(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col], errors="coerce")


# ---------------------------------------------------------------------------
# 3a. File existence checks
# ---------------------------------------------------------------------------

# Files that MUST exist and be non-empty
REQUIRED_FILES = [
    "scp/clustering/scp_clustering_report.json",
    "scp/clustering/scp_cluster_assignments.tsv",
    "scp/clustering/scp_annotation.tsv",
    "scp/clustering/scp_contrasts.tsv",
    "scp/markers/detection_markers.tsv",
    "scp/markers/intensity_markers_detected_only.tsv",
    "scp/markers/consensus_markers.tsv",
    "scp/llm/cluster_llm_prompt.md",
    "scp/llm/cluster_cell_type_annotations.tsv",
    "scp/pivot_pack.tsv",
]

# Files that SHOULD exist (warn if missing)
OPTIONAL_FILES = [
    "scp/qc/scp_qc_report.json",
    "scp/markers/scplainer_intensity_markers.tsv",
    "scp/scplainer/scplainer_variance_explained.tsv",
    "scp/aucell/tables/aucell_scores.tsv",
    "scp/clustering/qc/umap_by_cluster.pdf",
    "scp/clustering/qc/umap_by_library.pdf",
    "scp/clustering/qc/pca_intensity.pdf",
    "scp/llm/plots/umap_cell_types.pdf",
    "scp/llm/cluster_summary.tsv",
]


def check_file_existence(workdir: Path) -> list[CheckResult]:
    results = []
    for rel in REQUIRED_FILES:
        p = workdir / rel
        if not p.exists():
            results.append(CheckResult(f"file_exists:{rel}", False,
                                       "MISSING (required)"))
        elif p.stat().st_size == 0:
            results.append(CheckResult(f"file_exists:{rel}", False,
                                       "EXISTS but EMPTY"))
        else:
            results.append(CheckResult(f"file_exists:{rel}", True))

    for rel in OPTIONAL_FILES:
        p = workdir / rel
        if not p.exists():
            results.append(CheckResult(f"file_exists:{rel}", False,
                                       "MISSING (optional)", severity="WARN"))
        elif p.stat().st_size == 0:
            results.append(CheckResult(f"file_exists:{rel}", False,
                                       "EXISTS but EMPTY", severity="WARN"))
        else:
            results.append(CheckResult(f"file_exists:{rel}", True))
    return results


# ---------------------------------------------------------------------------
# 3b. Numerical sanity
# ---------------------------------------------------------------------------

def check_numerical_sanity(workdir: Path) -> list[CheckResult]:
    results = []
    scp = workdir / "scp"

    # --- Intensity markers: log2FC range ---
    im_path = scp / "markers/intensity_markers_detected_only.tsv"
    im = load_tsv(im_path)
    if im is not None and "log2FC_detected_only" in im.columns:
        fc = to_float_col(im, "log2FC_detected_only").dropna()
        if fc.size > 0:
            mn, mx = float(fc.min()), float(fc.max())
            ok = mn >= -15 and mx <= 15
            results.append(CheckResult(
                "intensity_log2FC_range", ok,
                f"range=[{mn:.2f}, {mx:.2f}]; expect [-15,+15]"
            ))
        # p-values in [0,1]
        if "pvalue" in im.columns:
            pv = to_float_col(im, "pvalue").dropna()
            ok = pv.size > 0 and float(pv.min()) >= 0 and float(pv.max()) <= 1
            results.append(CheckResult("intensity_pvalue_range", ok,
                                       f"[{float(pv.min()):.2e}, {float(pv.max()):.2e}]"))
        if "qvalue" in im.columns:
            qv = to_float_col(im, "qvalue").dropna()
            ok = qv.size > 0 and float(qv.min()) >= 0 and float(qv.max()) <= 1
            results.append(CheckResult("intensity_qvalue_range", ok,
                                       f"[{float(qv.min()):.2e}, {float(qv.max()):.2e}]"))

    # --- Detection markers: det_rate in [0,1] ---
    dm_path = scp / "markers/detection_markers.tsv"
    dm = load_tsv(dm_path)
    if dm is not None:
        for col in ["det_rate_in", "det_rate_out"]:
            if col in dm.columns:
                vals = to_float_col(dm, col).dropna()
                ok = vals.size > 0 and float(vals.min()) >= 0 and float(vals.max()) <= 1
                results.append(CheckResult(
                    f"detection_{col}_range", ok,
                    f"[{float(vals.min()):.3f}, {float(vals.max()):.3f}]"
                ))
        if "pvalue" in dm.columns:
            pv = to_float_col(dm, "pvalue").dropna()
            ok = pv.size > 0 and float(pv.min()) >= 0 and float(pv.max()) <= 1
            results.append(CheckResult("detection_pvalue_range", ok))
        if "qvalue" in dm.columns:
            qv = to_float_col(dm, "qvalue").dropna()
            ok = qv.size > 0 and float(qv.min()) >= 0 and float(qv.max()) <= 1
            results.append(CheckResult("detection_qvalue_range", ok))

    # --- AUCell scores in [0,1] ---
    ac_path = scp / "aucell/tables/aucell_scores.tsv"
    ac = load_tsv(ac_path)
    if ac is not None:
        first_col = ac.columns[0]
        num_cols = [c for c in ac.columns if c != first_col]
        if num_cols:
            vals = ac[num_cols].apply(pd.to_numeric, errors="coerce")
            flat = vals.values.flatten()
            flat = flat[np.isfinite(flat)]
            if flat.size > 0:
                ok = float(flat.min()) >= -0.01 and float(flat.max()) <= 1.01
                results.append(CheckResult(
                    "aucell_score_range", ok,
                    f"[{float(flat.min()):.4f}, {float(flat.max()):.4f}]"
                ))

    # --- Clustering report: n_cells > 0, n_proteins > 0 ---
    cr_path = scp / "clustering/scp_clustering_report.json"
    if cr_path.exists():
        cr = json.loads(cr_path.read_text(encoding="utf-8"))
        n_cells = cr.get("n_cells", 0)
        n_prot = cr.get("n_proteins", 0)
        results.append(CheckResult("clustering_n_cells", n_cells > 0,
                                   f"n_cells={n_cells}"))
        results.append(CheckResult("clustering_n_proteins", n_prot > 0,
                                   f"n_proteins={n_prot}"))

    # --- Consensus markers row count ---
    cm_path = scp / "markers/consensus_markers.tsv"
    cm = load_tsv(cm_path)
    if cm is not None and "Cluster" in cm.columns:
        n_clusters = cm["Cluster"].nunique()
        max_per = cm.groupby("Cluster").size().max()
        results.append(CheckResult(
            "consensus_top_n", max_per <= 50,
            f"{n_clusters} clusters, max {max_per} per cluster (expect <=50)"
        ))

    return results


# ---------------------------------------------------------------------------
# 3c. Cross-table consistency
# ---------------------------------------------------------------------------

def check_cross_table_consistency(workdir: Path) -> list[CheckResult]:
    results = []
    scp = workdir / "scp"

    # Load clustering report for reference cluster IDs
    cr_path = scp / "clustering/scp_clustering_report.json"
    if not cr_path.exists():
        return results
    cr = json.loads(cr_path.read_text(encoding="utf-8"))
    ref_clusters = set(cr.get("cluster_sizes", {}).keys())
    if not ref_clusters:
        return results

    # Cluster IDs from annotation
    ann = load_tsv(scp / "clustering/scp_annotation.tsv")
    if ann is not None and "Condition" in ann.columns:
        ann_clusters = set(ann["Condition"].unique())
        ok = ann_clusters == ref_clusters
        results.append(CheckResult(
            "cluster_ids_annotation_vs_report", ok,
            f"annotation={sorted(ann_clusters)}, report={sorted(ref_clusters)}"
        ))

    # Cluster IDs from marker files
    for name, rel in [
        ("detection", "markers/detection_markers.tsv"),
        ("intensity", "markers/intensity_markers_detected_only.tsv"),
        ("consensus", "markers/consensus_markers.tsv"),
    ]:
        df = load_tsv(scp / rel)
        if df is not None and "Cluster" in df.columns:
            mk_clusters = set(df["Cluster"].unique())
            # markers may not have C_other if it's tiny, so check subset
            missing = ref_clusters - mk_clusters
            extra = mk_clusters - ref_clusters
            ok = len(missing) <= 1 and len(extra) == 0  # allow baseline missing
            results.append(CheckResult(
                f"cluster_ids_{name}_vs_report", ok,
                f"missing={sorted(missing)}, extra={sorted(extra)}"
            ))

    # Cluster IDs from LLM annotations
    llm = load_tsv(scp / "llm/cluster_cell_type_annotations.tsv")
    if llm is not None and "Cluster" in llm.columns:
        llm_clusters = set(llm["Cluster"].unique())
        missing = ref_clusters - llm_clusters
        extra = llm_clusters - ref_clusters
        ok = len(missing) == 0  # LLM should annotate every cluster
        results.append(CheckResult(
            "cluster_ids_llm_vs_report", ok,
            f"missing={sorted(missing)}, extra={sorted(extra)}"
        ))

    # Number of cells: assignments vs clustering report
    assign = load_tsv(scp / "clustering/scp_cluster_assignments.tsv")
    n_cells_report = cr.get("n_cells", 0)
    if assign is not None:
        n_cells_assign = len(assign)
        ok = n_cells_assign == n_cells_report
        results.append(CheckResult(
            "n_cells_assignments_vs_report", ok,
            f"assignments={n_cells_assign}, report={n_cells_report}"
        ))

    # AUCell columns should match cells
    ac = load_tsv(scp / "aucell/tables/aucell_scores.tsv")
    if ac is not None:
        first_col = ac.columns[0]
        n_aucell_cells = len([c for c in ac.columns if c != first_col])
        # Allow some tolerance (some cells may be filtered)
        ok = abs(n_aucell_cells - n_cells_report) <= max(5, n_cells_report * 0.1)
        results.append(CheckResult(
            "n_cells_aucell_vs_report", ok,
            f"aucell={n_aucell_cells}, report={n_cells_report}",
            severity="WARN" if not ok else "FAIL"
        ))

    # Protein IDs in markers exist in pivot
    pivot = load_tsv(scp / "pivot_pack.tsv")
    if pivot is not None:
        prot_col = "Protein.Group" if "Protein.Group" in pivot.columns else pivot.columns[0]
        pivot_prots = set(pivot[prot_col].unique())

        for name, rel in [
            ("detection", "markers/detection_markers.tsv"),
            ("intensity", "markers/intensity_markers_detected_only.tsv"),
        ]:
            df = load_tsv(scp / rel)
            if df is not None and "Protein" in df.columns:
                mk_prots = set(df["Protein"].unique())
                missing = mk_prots - pivot_prots
                frac_missing = len(missing) / max(len(mk_prots), 1)
                ok = frac_missing < 0.05  # <5% missing is OK
                results.append(CheckResult(
                    f"proteins_{name}_in_pivot", ok,
                    f"{len(missing)}/{len(mk_prots)} not in pivot ({frac_missing:.1%})"
                ))

    return results


# ---------------------------------------------------------------------------
# 3d. Statistical validity
# ---------------------------------------------------------------------------

def check_statistical_validity(workdir: Path) -> list[CheckResult]:
    results = []
    scp = workdir / "scp"

    im = load_tsv(scp / "markers/intensity_markers_detected_only.tsv")
    if im is None or "pvalue" not in im.columns or "Cluster" not in im.columns:
        return results

    im["_pvalue"] = to_float_col(im, "pvalue")
    im["_qvalue"] = to_float_col(im, "qvalue") if "qvalue" in im.columns else np.nan

    from scipy import stats

    # p-value uniformity under null (KS test per cluster)
    ks_warnings = []
    for cl, grp in im.groupby("Cluster"):
        pv = grp["_pvalue"].dropna().values
        if pv.size < 20:
            continue
        ks_stat, ks_p = stats.kstest(pv, "uniform")
        if ks_p < 0.01:
            ks_warnings.append(f"{cl}(KS_p={ks_p:.2e})")

    if ks_warnings:
        results.append(CheckResult(
            "pvalue_uniformity", False,
            f"Non-uniform p-values (potential bias): {', '.join(ks_warnings)}",
            severity="WARN"
        ))
    else:
        results.append(CheckResult("pvalue_uniformity", True,
                                   "All clusters pass KS test (p>=0.01)"))

    # BH FDR monotonicity: q-values sorted by p-value should be non-decreasing
    fdr_issues = []
    if "qvalue" in im.columns:
        for cl, grp in im.groupby("Cluster"):
            sub = grp[["_pvalue", "_qvalue"]].dropna().sort_values("_pvalue")
            if sub.shape[0] < 5:
                continue
            qvals = sub["_qvalue"].values
            # Check non-decreasing (allow tiny floating point violations)
            diffs = np.diff(qvals)
            violations = (diffs < -1e-10).sum()
            if violations > 0:
                fdr_issues.append(f"{cl}({violations} violations)")

    if fdr_issues:
        results.append(CheckResult(
            "fdr_monotonicity", False,
            f"BH q-values not monotone: {', '.join(fdr_issues)}"
        ))
    else:
        results.append(CheckResult("fdr_monotonicity", True))

    # log2FC distribution: symmetric around 0
    if "log2FC_detected_only" in im.columns:
        fc = to_float_col(im, "log2FC_detected_only").dropna().values
        if fc.size > 10:
            mean_abs = float(np.mean(np.abs(fc)))
            ok = mean_abs < 2.0
            results.append(CheckResult(
                "log2FC_symmetry", ok,
                f"mean|log2FC|={mean_abs:.2f} (warn if >2.0)",
                severity="WARN" if not ok else "FAIL"
            ))

    return results


# ---------------------------------------------------------------------------
# 3e. Batch correction quality
# ---------------------------------------------------------------------------

def check_batch_correction(workdir: Path) -> list[CheckResult]:
    results = []
    scp = workdir / "scp"

    cr_path = scp / "clustering/scp_clustering_report.json"
    if not cr_path.exists():
        return results
    cr = json.loads(cr_path.read_text(encoding="utf-8"))

    batch_applied = cr.get("batch_applied", False)
    n_libs = cr.get("n_libraries", 1)

    if n_libs <= 1:
        results.append(CheckResult("batch_correction", True,
                                   f"Single library (n={n_libs}), no batch correction needed"))
        return results

    if not batch_applied:
        results.append(CheckResult(
            "batch_correction", False,
            f"n_libraries={n_libs} but batch correction not applied",
            severity="WARN"
        ))
        return results

    # Extract entropy from harmony log
    harmony_log = cr.get("harmony_adaptive_log", [])
    final_entropy = None
    for entry in harmony_log:
        if "entropy" in entry:
            final_entropy = entry["entropy"]

    if final_entropy is not None:
        ok = final_entropy >= 0.7
        results.append(CheckResult(
            "harmony_entropy", ok,
            f"entropy={final_entropy:.3f} (target >=0.7)",
            severity="WARN" if not ok else "FAIL"
        ))
    else:
        results.append(CheckResult("harmony_entropy", True,
                                   "No entropy recorded (batch applied=True)", severity="WARN"))

    # Check for batch-dominated clusters
    for entry in harmony_log:
        dominated = entry.get("batch_dominated_clusters", {})
        if dominated:
            for cl, info in dominated.items():
                ent = info.get("entropy", 0)
                results.append(CheckResult(
                    f"batch_dominated_{cl}", False,
                    f"entropy={ent:.3f} (<0.3 = batch-dominated)",
                    severity="WARN"
                ))

    return results


# ---------------------------------------------------------------------------
# 3f. LLM annotation quality
# ---------------------------------------------------------------------------

def check_llm_annotations(workdir: Path) -> list[CheckResult]:
    results = []
    scp = workdir / "scp"

    llm = load_tsv(scp / "llm/cluster_cell_type_annotations.tsv")
    if llm is None:
        return results

    if "Cluster" not in llm.columns or "cell_type" not in llm.columns:
        results.append(CheckResult("llm_columns", False,
                                   f"Missing Cluster or cell_type column. Cols: {llm.columns.tolist()}"))
        return results

    # All clusters annotated
    n_annotated = llm["Cluster"].nunique()
    cr_path = scp / "clustering/scp_clustering_report.json"
    if cr_path.exists():
        cr = json.loads(cr_path.read_text(encoding="utf-8"))
        n_expected = len(cr.get("cluster_sizes", {}))
        ok = n_annotated >= n_expected
        results.append(CheckResult(
            "llm_all_clusters_annotated", ok,
            f"annotated={n_annotated}, expected={n_expected}"
        ))

    # No empty annotations
    empty = llm["cell_type"].isna() | (llm["cell_type"].str.strip() == "")
    n_empty = int(empty.sum())
    results.append(CheckResult("llm_no_empty", n_empty == 0,
                               f"{n_empty} empty annotations"))

    # Unknown annotations (warn only)
    unknown = llm["cell_type"].str.lower().str.contains("unknown", na=False)
    n_unknown = int(unknown.sum())
    if n_unknown > 0:
        results.append(CheckResult(
            "llm_no_unknown", False,
            f"{n_unknown} clusters annotated as 'Unknown'",
            severity="WARN"
        ))

    # Confidence present
    if "confidence" in llm.columns:
        conf = llm["confidence"].fillna("").str.strip()
        n_missing = int((conf == "").sum())
        results.append(CheckResult(
            "llm_confidence_present", n_missing == 0,
            f"{n_missing} missing confidence values",
            severity="WARN" if n_missing > 0 else "FAIL"
        ))

    return results


# ---------------------------------------------------------------------------
# Phase 4. Known-marker recovery benchmark
# ---------------------------------------------------------------------------

KNOWN_MARKERS = {
    "BrainSCPpx": {
        "species": "human",
        "expected_types": {
            "neurons": ["RBFOX3", "SYN1", "MAP2", "NEFL", "SNAP25"],
            "astrocytes": ["GFAP", "AQP4", "SLC1A3", "ALDH1L1", "S100B"],
            "oligodendrocytes": ["MBP", "PLP1", "MOG", "CNP", "MAG"],
            "microglia": ["AIF1", "CX3CR1", "CSF1R", "ITGAM", "CD68"],
            "endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "FLT1"],
        },
    },
    "NeutrophilGBMpaper": {
        "species": "human",
        "expected_types": {
            "mature_neutrophils": ["ELANE", "MPO", "LCN2", "CTSG", "PRTN3"],
            "neutrophil_activation": ["S100A8", "S100A9", "MMP9", "CAMP", "LTF"],
            "tumor": ["EGFR", "VIM", "FN1", "COL1A1", "SPARC"],
        },
    },
    "Day7Caerulein": {
        "species": "mouse",
        "expected_types": {
            "acinar": ["Amy2a", "Cpa1", "Cela1", "Prss1", "Ctrb1"],
            "ductal": ["Krt19", "Krt8", "Sox9", "Epcam", "Spp1"],
            "immune": ["Ptprc", "Cd68", "Lyz2", "Lcp1", "Itgam"],
            "stellate": ["Vim", "Col1a1", "Des", "Acta2", "Sparc"],
        },
    },
}


def check_known_marker_recovery(workdir: Path, dataset_key: str | None = None,
                                custom_markers: dict | None = None) -> list[CheckResult]:
    results = []
    scp = workdir / "scp"

    # Determine which marker set to use
    markers_spec = None
    if custom_markers:
        markers_spec = custom_markers
    elif dataset_key:
        markers_spec = KNOWN_MARKERS.get(dataset_key)
    else:
        # Try to auto-detect from workdir name
        wname = workdir.name.lower()
        for key in KNOWN_MARKERS:
            if key.lower() in wname:
                markers_spec = KNOWN_MARKERS[key]
                dataset_key = key
                break

    if markers_spec is None:
        results.append(CheckResult("known_markers", True,
                                   "No known-marker set for this dataset (skipped)",
                                   severity="WARN"))
        return results

    # Load consensus markers
    cm = load_tsv(scp / "markers/consensus_markers.tsv")
    # Load intensity markers
    im = load_tsv(scp / "markers/intensity_markers_detected_only.tsv")
    # Load detection markers
    dm = load_tsv(scp / "markers/detection_markers.tsv")

    if cm is None and im is None and dm is None:
        results.append(CheckResult("known_markers_files", False,
                                   "No marker files found"))
        return results

    # Get all gene names from consensus top-50
    consensus_genes = set()
    if cm is not None and "Genes" in cm.columns:
        consensus_genes = set(
            g.strip()
            for genes_str in cm["Genes"].dropna()
            for g in str(genes_str).split(";")
            if g.strip()
        )

    # Get significant genes from intensity markers (q < 0.05)
    sig_intensity_genes = set()
    if im is not None and "Genes" in im.columns and "qvalue" in im.columns:
        im_filt = im[to_float_col(im, "qvalue") < 0.05]
        sig_intensity_genes = set(
            g.strip()
            for genes_str in im_filt["Genes"].dropna()
            for g in str(genes_str).split(";")
            if g.strip()
        )

    # Get significant genes from detection markers (q < 0.05)
    sig_detection_genes = set()
    if dm is not None and "Genes" in dm.columns and "qvalue" in dm.columns:
        dm_filt = dm[to_float_col(dm, "qvalue") < 0.05]
        sig_detection_genes = set(
            g.strip()
            for genes_str in dm_filt["Genes"].dropna()
            for g in str(genes_str).split(";")
            if g.strip()
        )

    all_recovered = consensus_genes | sig_intensity_genes | sig_detection_genes

    # Case-insensitive matching
    all_recovered_lower = {g.lower() for g in all_recovered}

    expected = markers_spec.get("expected_types", {})
    total_expected = 0
    total_found = 0

    for cell_type, gene_list in expected.items():
        found = []
        missing = []
        for g in gene_list:
            if g.lower() in all_recovered_lower:
                found.append(g)
            else:
                missing.append(g)
        total_expected += len(gene_list)
        total_found += len(found)

        frac = len(found) / max(len(gene_list), 1)
        ok = frac >= 0.4  # At least 40% of markers recovered per cell type
        results.append(CheckResult(
            f"known_markers_{dataset_key}_{cell_type}", ok,
            f"{len(found)}/{len(gene_list)} recovered ({frac:.0%}): "
            f"found={found}, missing={missing}",
            severity="WARN" if not ok else "FAIL"
        ))

    overall_frac = total_found / max(total_expected, 1)
    ok = overall_frac >= 0.5  # At least 50% overall
    results.append(CheckResult(
        f"known_markers_{dataset_key}_overall", ok,
        f"{total_found}/{total_expected} total markers recovered ({overall_frac:.0%})"
    ))

    return results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Validate SCP pipeline outputs")
    ap.add_argument("--workdir", required=True, help="Dataset working directory containing scp/")
    ap.add_argument("--dataset-key", default=None,
                    help="Known-marker dataset key (BrainSCPpx, NeutrophilGBMpaper, Day7Caerulein)")
    ap.add_argument("--known-markers", default=None,
                    help="JSON file with custom known markers (overrides built-in)")
    ap.add_argument("--strict", action="store_true",
                    help="Treat WARNings as failures for exit code")
    ap.add_argument("--json-out", default=None,
                    help="Write structured results to JSON file")
    args = ap.parse_args()

    workdir = Path(args.workdir)
    if not workdir.exists():
        print(f"ERROR: workdir does not exist: {workdir}")
        sys.exit(2)
    if not (workdir / "scp").exists():
        print(f"ERROR: no scp/ directory in {workdir}")
        sys.exit(2)

    custom_markers = None
    if args.known_markers:
        custom_markers = json.loads(Path(args.known_markers).read_text(encoding="utf-8"))

    print(f"Validating SCP outputs: {workdir}")
    print(f"{'='*70}")

    all_results: list[CheckResult] = []

    # Phase 3a: File existence
    print("\n## File Existence")
    res = check_file_existence(workdir)
    all_results.extend(res)
    for r in res:
        if not r.passed:
            print(f"  {r}")

    # Phase 3b: Numerical sanity
    print("\n## Numerical Sanity")
    res = check_numerical_sanity(workdir)
    all_results.extend(res)
    for r in res:
        if not r.passed:
            print(f"  {r}")

    # Phase 3c: Cross-table consistency
    print("\n## Cross-table Consistency")
    res = check_cross_table_consistency(workdir)
    all_results.extend(res)
    for r in res:
        if not r.passed:
            print(f"  {r}")

    # Phase 3d: Statistical validity
    print("\n## Statistical Validity")
    res = check_statistical_validity(workdir)
    all_results.extend(res)
    for r in res:
        if not r.passed:
            print(f"  {r}")

    # Phase 3e: Batch correction
    print("\n## Batch Correction Quality")
    res = check_batch_correction(workdir)
    all_results.extend(res)
    for r in res:
        if not r.passed:
            print(f"  {r}")

    # Phase 3f: LLM annotations
    print("\n## LLM Annotation Quality")
    res = check_llm_annotations(workdir)
    all_results.extend(res)
    for r in res:
        if not r.passed:
            print(f"  {r}")

    # Phase 4: Known-marker recovery
    print("\n## Known-marker Recovery")
    res = check_known_marker_recovery(workdir, args.dataset_key, custom_markers)
    all_results.extend(res)
    for r in res:
        print(f"  {r}")

    # Summary
    n_pass = sum(1 for r in all_results if r.passed)
    n_fail = sum(1 for r in all_results if not r.passed and r.severity == "FAIL")
    n_warn = sum(1 for r in all_results if not r.passed and r.severity == "WARN")
    n_total = len(all_results)

    print(f"\n{'='*70}")
    print(f"SUMMARY: {n_pass}/{n_total} passed, {n_fail} FAIL, {n_warn} WARN")

    if n_fail > 0:
        print("\nFAILURES:")
        for r in all_results:
            if not r.passed and r.severity == "FAIL":
                print(f"  {r}")

    if n_warn > 0:
        print("\nWARNINGS:")
        for r in all_results:
            if not r.passed and r.severity == "WARN":
                print(f"  {r}")

    # JSON output
    if args.json_out:
        out = {
            "workdir": str(workdir),
            "n_pass": n_pass,
            "n_fail": n_fail,
            "n_warn": n_warn,
            "n_total": n_total,
            "overall": "PASS" if n_fail == 0 and (not args.strict or n_warn == 0) else "FAIL",
            "checks": [
                {
                    "name": r.name,
                    "passed": r.passed,
                    "severity": r.severity,
                    "detail": r.detail,
                }
                for r in all_results
            ],
        }
        Path(args.json_out).write_text(json.dumps(out, indent=2), encoding="utf-8")
        print(f"\nJSON report -> {args.json_out}")

    if n_fail > 0:
        print("\nRESULT: FAIL")
        sys.exit(1)
    elif args.strict and n_warn > 0:
        print("\nRESULT: FAIL (strict mode, warnings present)")
        sys.exit(1)
    else:
        print("\nRESULT: PASS")
        sys.exit(0)


if __name__ == "__main__":
    main()
