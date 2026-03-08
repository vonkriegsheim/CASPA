#!/usr/bin/env python3
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
"""
make_cluster_summary.py

Generate a per-cluster summary table combining cell counts, QC metrics,
top markers from all three modalities, and top AUCell pathways.

Inputs
------
--annotation          scp_annotation.tsv (Run, Condition)
--assignments         scp_cluster_assignments.tsv (n_detected_proteins etc.)
--detection-markers   detection_markers.tsv       (optional)
--intensity-markers   intensity_markers.tsv        (optional)
--scplainer-markers   scplainer_intensity_markers.tsv (optional)
--aucell-scores       aucell_scores.tsv            (optional)
--consensus-markers   consensus_markers.tsv        (optional)
--out                 output TSV path
--top-n-markers       how many top markers to include per cluster (default 5)
--q-threshold         significance cutoff for marker inclusion (default 0.05)
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def stem_from_header(h: str) -> str:
    s = str(h)
    p = Path(s)
    if "\\" in s or "/" in s:
        return p.stem
    if "." in s and s.rsplit(".", 1)[1].lower() in ("raw", "d", "mzml", "wiff"):
        return p.stem
    return s


def gene_name(df: pd.DataFrame, row_idx) -> str:
    for col in ["Genes", "Gene", "gene_primary"]:
        if col in df.columns:
            v = str(df.iloc[row_idx].get(col, ""))
            if v and v not in ("", "nan"):
                return v.split(";")[0].split("|")[0].strip()
    prot_col = "Protein" if "Protein" in df.columns else df.columns[0]
    return str(df.iloc[row_idx][prot_col])[:15]


def top_markers_str(df: pd.DataFrame, cluster: str, cluster_col: str,
                    score_col: str, q_col: str, q_thresh: float, top_n: int) -> str:
    sub = df[df[cluster_col].astype(str) == str(cluster)].copy()
    sub[score_col] = pd.to_numeric(sub.get(score_col, pd.Series(dtype=float)), errors="coerce")
    sub[q_col]     = pd.to_numeric(sub.get(q_col, pd.Series(dtype=float)), errors="coerce")
    sub = sub[sub[q_col] <= q_thresh].sort_values(score_col, ascending=False).head(top_n)
    parts = []
    for i in range(len(sub)):
        name = gene_name(sub, i)
        score = sub.iloc[i][score_col]
        q = sub.iloc[i][q_col]
        parts.append(f"{name} ({score_col}={score:.2f}, q={q:.1e})")
    return "; ".join(parts) if parts else ""


def top_pathways_str(auc_mat: np.ndarray, pathway_names: list[str],
                     cell_cols: list[str], cluster_cell_mask: np.ndarray,
                     top_n: int) -> str:
    if not cluster_cell_mask.any():
        return ""
    in_scores  = auc_mat[:, cluster_cell_mask].mean(axis=1)
    out_scores = auc_mat[:, ~cluster_cell_mask].mean(axis=1) if (~cluster_cell_mask).any() else np.zeros(auc_mat.shape[0])
    delta = in_scores - out_scores
    top_idx = np.argsort(delta)[::-1][:top_n]
    parts = []
    for i in top_idx:
        name = pathway_names[i][:35]
        parts.append(f"{name} (AUC={in_scores[i]:.3f}, delta={delta[i]:+.3f})")
    return "; ".join(parts)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--annotation",        required=True)
    ap.add_argument("--assignments",       required=True)
    ap.add_argument("--detection-markers", default=None)
    ap.add_argument("--intensity-markers", default=None)
    ap.add_argument("--scplainer-markers", default=None)
    ap.add_argument("--aucell-scores",     default=None)
    ap.add_argument("--consensus-markers", default=None)
    ap.add_argument("--out",               required=True)
    ap.add_argument("--top-n-markers",     type=int,   default=5)
    ap.add_argument("--q-threshold",       type=float, default=0.05)
    ap.add_argument("--llm-export",        default=None,
                    help="Path for LLM-ready markdown export (default: cluster_llm_prompt.md beside --out)")
    ap.add_argument("--custom-proteins",   default=None,
                    help="Comma-separated gene names (or path to txt file, one gene per line) "
                         "that are biologically meaningful to the researcher. "
                         "Included verbatim in the LLM prompt as researcher-provided markers.")
    ap.add_argument("--config", action="append", default=[],
                    help="dataset.json (or multiple configs). Experiment context "
                         "(project.name, project.description, project.species_label) "
                         "is included in the LLM prompt.")
    args = ap.parse_args()

    # ---- Merge config(s) for experiment context ----
    exp_cfg = {}
    for cfg_path in args.config:
        p = Path(cfg_path)
        if p.exists():
            with open(p, encoding="utf-8") as f:
                c = json.load(f)
            # Deep-merge top-level dicts
            for k, v in c.items():
                if isinstance(v, dict) and isinstance(exp_cfg.get(k), dict):
                    exp_cfg[k].update(v)
                else:
                    exp_cfg[k] = v

    # ---- Core data ----
    ann    = pd.read_csv(args.annotation,  sep="\t", dtype=str)
    assign = pd.read_csv(args.assignments, sep="\t", dtype=str)

    run_to_cond = dict(zip(ann["Run"].astype(str), ann["Condition"].astype(str)))
    clusters = sorted(set(run_to_cond.values()))

    # Merge assignment QC metrics
    assign["Run"] = assign["Run"].astype(str)
    assign["n_detected_proteins"] = pd.to_numeric(assign.get("n_detected_proteins"), errors="coerce")
    assign["missing_fraction"]    = pd.to_numeric(assign.get("missing_fraction"),    errors="coerce")
    assign["Condition"] = assign["Run"].map(run_to_cond)

    # Optional inputs
    dm = pd.read_csv(args.detection_markers,  sep="\t", dtype=str) if (
        args.detection_markers and Path(args.detection_markers).exists()) else None
    im = pd.read_csv(args.intensity_markers,  sep="\t", dtype=str) if (
        args.intensity_markers and Path(args.intensity_markers).exists()) else None
    sm = pd.read_csv(args.scplainer_markers,  sep="\t", dtype=str) if (
        args.scplainer_markers and Path(args.scplainer_markers).exists()) else None
    cm = pd.read_csv(args.consensus_markers,  sep="\t", dtype=str) if (
        args.consensus_markers and Path(args.consensus_markers).exists()) else None

    # AUCell: pathways x cells wide table
    auc_mat = None
    auc_path_names = []
    auc_run_order  = []
    if args.aucell_scores and Path(args.aucell_scores).exists():
        auc_df = pd.read_csv(args.aucell_scores, sep="\t", dtype=str)
        path_col = "Pathway" if "Pathway" in auc_df.columns else auc_df.columns[0]
        auc_path_names = auc_df[path_col].astype(str).tolist()
        auc_sample_cols = [c for c in auc_df.columns if c != path_col]
        auc_stem_map = {c: stem_from_header(c) for c in auc_sample_cols}
        # Align to annotation
        run_set = set(ann["Run"].astype(str))
        matched = [c for c in auc_sample_cols if auc_stem_map.get(c, "") in run_set]
        if matched:
            auc_run_order = [auc_stem_map[c] for c in matched]
            auc_mat = auc_df[matched].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    # ---- Build summary ----
    rows = []
    for cl in clusters:
        cl_cells = assign[assign["Condition"] == cl]
        n_cells = len(cl_cells)
        med_det = float(cl_cells["n_detected_proteins"].median()) if not cl_cells["n_detected_proteins"].isna().all() else float("nan")
        med_miss = float(cl_cells["missing_fraction"].median()) if not cl_cells["missing_fraction"].isna().all() else float("nan")

        row = {
            "Cluster":              cl,
            "n_cells":              n_cells,
            "median_n_detected":    round(med_det, 1) if np.isfinite(med_det) else "",
            "median_missing_frac":  round(med_miss, 3) if np.isfinite(med_miss) else "",
        }

        # Detection markers: top by det_rate_in (how specific to this cluster)
        if dm is not None and "Cluster" in dm.columns and "det_rate_in" in dm.columns:
            row["top_detection_markers"] = top_markers_str(
                dm, cl, "Cluster", "det_rate_in", "qvalue", args.q_threshold, args.top_n_markers)
        else:
            row["top_detection_markers"] = ""

        # Intensity markers: top by abs(log2FC)
        if im is not None and "Cluster" in im.columns and "log2FC_detected_only" in im.columns:
            im_sub = im[im["Cluster"].astype(str) == str(cl)].copy()
            im_sub["_abs"] = pd.to_numeric(im_sub["log2FC_detected_only"], errors="coerce").abs()
            im_sub["log2FC_detected_only"] = pd.to_numeric(im_sub["log2FC_detected_only"], errors="coerce")
            im_sub["qvalue"] = pd.to_numeric(im_sub["qvalue"], errors="coerce")
            im_top = im_sub[im_sub["qvalue"] <= args.q_threshold].sort_values("_abs", ascending=False).head(args.top_n_markers)
            parts = []
            for i in range(len(im_top)):
                name = gene_name(im_top, i)
                fc = im_top.iloc[i]["log2FC_detected_only"]
                q = im_top.iloc[i]["qvalue"]
                parts.append(f"{name} (log2FC={fc:+.2f}, q={q:.1e})")
            row["top_intensity_markers"] = "; ".join(parts)
        else:
            row["top_intensity_markers"] = ""

        # Scplainer markers
        if sm is not None:
            sm_cl_col = "Cluster" if "Cluster" in sm.columns else "contrast"
            if sm_cl_col in sm.columns and "log2FC" in sm.columns:
                sm_sub = sm[sm[sm_cl_col].astype(str) == str(cl)].copy()
                sm_sub["_abs"] = pd.to_numeric(sm_sub["log2FC"], errors="coerce").abs()
                sm_sub["log2FC"] = pd.to_numeric(sm_sub["log2FC"], errors="coerce")
                sm_sub["qvalue"] = pd.to_numeric(sm_sub["qvalue"], errors="coerce")
                if "ok" in sm_sub.columns:
                    sm_sub = sm_sub[sm_sub["ok"].astype(str).str.upper().isin(["TRUE", "1"])]
                sm_top = sm_sub[sm_sub["qvalue"] <= args.q_threshold].sort_values("_abs", ascending=False).head(args.top_n_markers)
                parts = []
                for i in range(len(sm_top)):
                    name = gene_name(sm_top, i)
                    fc = sm_top.iloc[i]["log2FC"]
                    q = sm_top.iloc[i]["qvalue"]
                    parts.append(f"{name} (log2FC={fc:+.2f}, q={q:.1e})")
                row["top_scplainer_markers"] = "; ".join(parts)
            else:
                row["top_scplainer_markers"] = ""
        else:
            row["top_scplainer_markers"] = ""

        # Consensus top markers
        if cm is not None and "Cluster" in cm.columns and "consensus_score" in cm.columns:
            cm_sub = cm[cm["Cluster"].astype(str) == str(cl)].copy()
            cm_sub["consensus_score"] = pd.to_numeric(cm_sub["consensus_score"], errors="coerce")
            cm_sub["n_modalities"] = pd.to_numeric(cm_sub.get("n_modalities", 0), errors="coerce")
            cm_sub = cm_sub.sort_values("consensus_score", ascending=False).head(args.top_n_markers)
            parts = []
            for i in range(len(cm_sub)):
                name = gene_name(cm_sub, i)
                score = cm_sub.iloc[i]["consensus_score"]
                n_mod = int(cm_sub.iloc[i]["n_modalities"]) if pd.notna(cm_sub.iloc[i]["n_modalities"]) else 0
                parts.append(f"{name} (score={score:.2f}, {n_mod}/3 modalities)")
            row["top_consensus_markers"] = "; ".join(parts)
        else:
            row["top_consensus_markers"] = ""

        # AUCell top pathways
        if auc_mat is not None:
            cl_run_mask = np.array([run_to_cond.get(r, "") == cl for r in auc_run_order])
            row["top_aucell_pathways"] = top_pathways_str(
                auc_mat, auc_path_names, auc_run_order, cl_run_mask, args.top_n_markers)
        else:
            row["top_aucell_pathways"] = ""

        rows.append(row)

    summary = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out}")
    print(summary[["Cluster", "n_cells", "median_n_detected"]].to_string(index=False))

    # ---- Parse custom protein list ----
    custom_proteins = []
    if args.custom_proteins:
        cp = args.custom_proteins.strip()
        from pathlib import Path as _Path
        if _Path(cp).exists():
            custom_proteins = [l.strip() for l in _Path(cp).read_text(encoding="utf-8").splitlines() if l.strip()]
        else:
            custom_proteins = [g.strip() for g in cp.split(",") if g.strip()]

    # ---- LLM-ready markdown export ----
    llm_path = Path(args.llm_export) if args.llm_export else Path(args.out).parent / "cluster_llm_prompt.md"
    llm_path.parent.mkdir(parents=True, exist_ok=True)

    lines = []
    lines.append("# Cluster marker summary for cell type annotation")
    lines.append("")

    # Experiment context from config
    proj = exp_cfg.get("project", {})
    exp_name = proj.get("name", "")
    exp_desc = proj.get("description", "")
    exp_species = proj.get("species_label", "")
    if exp_name or exp_desc or exp_species:
        lines.append("## Experiment context")
        lines.append("")
        if exp_name:
            lines.append(f"- **Experiment**: {exp_name}")
        if exp_species:
            lines.append(f"- **Species**: {exp_species}")
        if exp_desc:
            lines.append(f"- **Description**: {exp_desc}")
        lines.append("")

    # -- Change 1: Run metadata block (placeholders for calling script) --
    lines.append("## Run metadata (auto-populated at runtime)")
    lines.append("")
    lines.append("- **Prompt version**: v2")
    lines.append("- **Model**: {model_name}")
    lines.append("- **Temperature**: 0")
    lines.append("- **Run timestamp**: {iso_timestamp}")
    lines.append("- **Prompt hash (SHA256)**: {sha256_of_prompt_text}")
    lines.append("")
    lines.append("*These fields are populated by the calling script and stored alongside the")
    lines.append("raw LLM output for audit trail purposes.*")
    lines.append("")

    # -- Change 2: Expanded annotation instructions --
    lines.append(
        "You are a single-cell proteomics expert. Below is a per-cluster summary from a "
        "label-free single-cell proteomics (SCP) experiment. For each cluster you are "
        "given statistical evidence across three complementary modalities: detection "
        "specificity (Fisher exact), intensity fold-change (Wilcoxon, detected-only), "
        "and scplainer (linear mixed-effects, correcting for technical variation). "
        "Consensus markers are ranked across all three modalities."
    )
    lines.append("")
    lines.append("**Your task for each cluster:**")
    lines.append("")
    lines.append("1. Propose the most likely **cell type label**, drawing on the marker evidence "
                 "and the experiment context above.")
    lines.append("")
    lines.append("2. List the **key supporting markers** -- prioritise proteins with concordant "
                 "evidence across >=2 modalities and/or strong detection specificity "
                 "(det_rate_in >> det_rate_out, low q-value).")
    lines.append("")
    lines.append("3. List any **contradictory or confounding signals** that weaken the call "
                 "(e.g. immune markers in a nominally stromal cluster, low n_detected, "
                 "absence of expected canonical markers).")
    lines.append("")
    lines.append("4. Assign a **structured confidence score** using the rubric below.")
    lines.append("")
    lines.append("5. If confidence is Low or if the cluster is ambiguous, explicitly state what "
                 "additional marker evidence would resolve the annotation, citing specific "
                 "proteins by name.")
    lines.append("")
    lines.append("**Confidence rubric (use exactly these labels):**")
    lines.append("")
    lines.append("| Label | Criteria |")
    lines.append("|---|---|")
    lines.append("| High | >=3 canonical lineage markers concordant across >=2 modalities; no major contradictions |")
    lines.append("| Medium-high | 2-3 canonical markers concordant; minor contradictions explainable by biology |")
    lines.append("| Medium | 1-2 canonical markers; contradictions present but cell type still most parsimonious |")
    lines.append("| Low | Canonical markers absent or contradicted; cluster is ambiguous or likely a doublet/artefact |")
    lines.append("")
    lines.append("**Output format -- return a Markdown table with exactly these columns:**")
    lines.append("")
    lines.append("| Cluster | Suggested cell type | Key supporting markers | Contradictions / caveats | Confidence | Resolving markers (if Low/ambiguous) |")
    lines.append("")
    lines.append("Then, after the table, provide a **per-cluster narrative** (one paragraph each) "
                 "with quantitative citations (e.g. \"Cd68 det_rate_in=1.00 vs 0.25, q=0.0000\"). "
                 "Do not assert a cell type without citing at least one specific statistic.")
    lines.append("")
    lines.append("**Important caveats for this dataset:**")
    lines.append("")
    lines.append("- Proteins detected across >=90% of all clusters should be treated as "
                 "constitutively expressed or background and are not lineage-informative on "
                 "their own.")
    lines.append("- Depleted detection markers (proteins specifically *absent* in a cluster) are "
                 "biologically informative and should be used to exclude lineages, not only "
                 "positive markers to assign them.")
    lines.append("")
    lines.append("**Contamination vs ambient carryover -- be conservative:**")
    lines.append("")
    lines.append("- In single-cell proteomics, keratins (KRT5, KRT14, KRT10, TGM1), serum "
                 "proteins (ALB, TF, HP, A2M), and hemoglobin (HBB, HBA1) are commonly "
                 "detected as **ambient carryover** from sample preparation rather than true "
                 "contaminant cell populations.")
    lines.append("- **Do NOT label a cluster as 'Contaminant' solely because it expresses "
                 "keratins, ALB, or HBB.** These proteins are frequently detected in genuine "
                 "cell populations alongside their true lineage markers.")
    lines.append("- Only call contamination when ALL of: (a) the cluster lacks any coherent "
                 "lineage-specific markers across all modalities, AND (b) contaminant markers "
                 "are overwhelmingly dominant (>80% of top markers), AND (c) no alternative "
                 "biological explanation is plausible given the experiment context.")
    lines.append("- When contaminant-associated proteins are present alongside genuine lineage "
                 "markers, annotate the cluster by its **biological identity** and note "
                 "ambient contamination as a caveat in the Contradictions column.")
    lines.append("- For ambiguous clusters with mixed contaminant and lineage signals, prefer "
                 "'Low confidence' biological annotation over a 'Contaminant' label. The "
                 "researcher can decide post-hoc whether to filter these clusters.")
    lines.append("")

    # -- Change 3: Marker coverage score instruction --
    lines.append("## Marker coverage score")
    lines.append("")
    lines.append("For each cluster, the calling script will pre-compute a **marker coverage "
                 "score** defined as:")
    lines.append("")
    lines.append("    coverage = (number of recommended panel markers detected at det_rate >= 0.20")
    lines.append("                in this cluster) / (total recommended panel markers for the")
    lines.append("                assigned cell type)")
    lines.append("")
    lines.append("This score (0-1) will be appended to the annotation table post-hoc by the "
                 "pipeline. You do not need to compute it, but you **must** assign a cell type "
                 "label that maps unambiguously to one of the panels defined in "
                 "the recommended markers output so the script can look it up. If the best label "
                 "is ambiguous between two panels, name both separated by \" / \" (e.g. "
                 "\"Macrophage (Mrc1+) / Inflammatory myeloid\").")
    lines.append("")

    # -- Change 4: Pass structure --
    lines.append("## Pass structure")
    lines.append("")
    lines.append("This is **Pass 1**. Your output will be reviewed and a supplemental marker "
                 "query will be generated for any cluster where:")
    lines.append("- Confidence is Low or Medium, or")
    lines.append("- A contradictory signal was flagged, or")
    lines.append("- The suggested cell type requires confirmation of a specific marker not yet "
                 "in the data.")
    lines.append("")
    lines.append("Pass 2 will re-present the same cluster data with additional researcher-"
                 "requested marker statistics appended. In Pass 2, update only clusters where "
                 "new evidence changes the call or confidence; carry forward unchanged calls "
                 "explicitly stating \"No change from Pass 1.\"")
    lines.append("")

    # -- Column definitions (retained from v1 but reformatted) --
    lines.append("## Data columns provided per cluster")
    lines.append("")
    lines.append("- **n_cells**: number of single cells in the cluster")
    lines.append("- **median_n_detected**: median number of proteins detected per cell")
    lines.append("- **detection markers**: proteins most specifically *detected* in this cluster "
                 "(Fisher exact test on presence/absence). Columns: det_rate_in, det_rate_out, qvalue")
    lines.append("- **depleted detection markers**: proteins specifically *absent/under-detected* "
                 "in this cluster vs others")
    lines.append("- **intensity markers (up-regulated)**: proteins with highest *intensity* in detected cells "
                 "(Wilcoxon, detected-only). Columns: log2FC (positive), qvalue")
    lines.append("- **depleted intensity markers (down-regulated)**: proteins with *lowest* intensity "
                 "in this cluster vs others (negative log2FC)")
    lines.append("- **scplainer markers (up/down-regulated)**: linear mixed-effects model, "
                 "accounting for technical variation. Columns: log2FC, qvalue")
    lines.append("- **consensus markers**: proteins ranked consistently across all three modalities. "
                 "Columns: score, n_modalities (how many of 3 agreed)")
    lines.append("- **AUCell pathways**: most active MSigDB Hallmark pathways. "
                 "Columns: AUC_in, delta_vs_rest")
    lines.append("- **researcher markers of interest**: biologically relevant markers provided by "
                 "the researcher with per-cluster stats (may not reach significance everywhere)")
    lines.append("")
    lines.append("---")
    lines.append("")

    for _, row in summary.iterrows():
        lines.append(f"## Cluster {row['Cluster']}")
        lines.append("")
        lines.append(f"- **n_cells**: {row['n_cells']}")
        lines.append(f"- **median_n_detected**: {row.get('median_n_detected', 'N/A')}")

        cl_str = str(row["Cluster"])

        # -- Consensus markers table --
        if cm is not None and "Cluster" in cm.columns and "consensus_score" in cm.columns:
            cm_sub = cm[cm["Cluster"].astype(str) == cl_str].copy()
            cm_sub["consensus_score"] = pd.to_numeric(cm_sub["consensus_score"], errors="coerce")
            cm_sub["n_modalities"] = pd.to_numeric(cm_sub.get("n_modalities", 0), errors="coerce")
            cm_top = cm_sub.sort_values("consensus_score", ascending=False).head(args.top_n_markers)
            if len(cm_top):
                lines.append("- **consensus markers**:")
                lines.append("")
                lines.append("  | Gene | score | n_modalities |")
                lines.append("  |------|-------|--------------|")
                for i in range(len(cm_top)):
                    g = gene_name(cm_top, i)
                    sc = cm_top.iloc[i]["consensus_score"]
                    nm = int(cm_top.iloc[i]["n_modalities"]) if pd.notna(cm_top.iloc[i]["n_modalities"]) else 0
                    lines.append(f"  | {g} | {sc:.2f} | {nm}/3 |")
                lines.append("")

        # -- Detection markers table (enriched) --
        if dm is not None and "Cluster" in dm.columns and "det_rate_in" in dm.columns:
            dm_sub = dm[dm["Cluster"].astype(str) == cl_str].copy()
            dm_sub["det_rate_in"] = pd.to_numeric(dm_sub["det_rate_in"], errors="coerce")
            dm_sub["det_rate_out"] = pd.to_numeric(dm_sub.get("det_rate_out", pd.Series(dtype=float)), errors="coerce")
            dm_sub["qvalue"] = pd.to_numeric(dm_sub["qvalue"], errors="coerce")
            dm_sig = dm_sub[dm_sub["qvalue"] <= args.q_threshold]
            dm_top = dm_sig.sort_values("det_rate_in", ascending=False).head(args.top_n_markers)
            if len(dm_top):
                lines.append("- **detection markers**:")
                lines.append("")
                lines.append("  | Gene | det_rate_in | det_rate_out | qvalue |")
                lines.append("  |------|-------------|--------------|--------|")
                for i in range(len(dm_top)):
                    g = gene_name(dm_top, i)
                    dr_in = dm_top.iloc[i]["det_rate_in"]
                    dr_out = dm_top.iloc[i]["det_rate_out"]
                    q = dm_top.iloc[i]["qvalue"]
                    dr_out_s = f"{dr_out:.2f}" if pd.notna(dr_out) else "—"
                    lines.append(f"  | {g} | {dr_in:.2f} | {dr_out_s} | {q:.1e} |")
                lines.append("")

            # Depleted: proteins detected much less in this cluster than elsewhere
            dm_dep = dm_sig[dm_sig["det_rate_out"] > dm_sig["det_rate_in"]].copy()
            dm_dep["_depletion"] = dm_dep["det_rate_out"] - dm_dep["det_rate_in"]
            dm_dep = dm_dep.sort_values("_depletion", ascending=False).head(args.top_n_markers)
            if len(dm_dep):
                lines.append("- **depleted detection markers** (under-detected vs other clusters):")
                lines.append("")
                lines.append("  | Gene | det_rate_in | det_rate_out | qvalue |")
                lines.append("  |------|-------------|--------------|--------|")
                for i in range(len(dm_dep)):
                    g = gene_name(dm_dep, i)
                    dr_in = dm_dep.iloc[i]["det_rate_in"]
                    dr_out = dm_dep.iloc[i]["det_rate_out"]
                    q = dm_dep.iloc[i]["qvalue"]
                    dr_out_s = f"{dr_out:.2f}" if pd.notna(dr_out) else "—"
                    lines.append(f"  | {g} | {dr_in:.2f} | {dr_out_s} | {q:.1e} |")
                lines.append("")

        # -- Intensity markers table (up-regulated) --
        if im is not None and "Cluster" in im.columns and "log2FC_detected_only" in im.columns:
            im_sub = im[im["Cluster"].astype(str) == cl_str].copy()
            im_sub["log2FC_detected_only"] = pd.to_numeric(im_sub["log2FC_detected_only"], errors="coerce")
            im_sub["qvalue"] = pd.to_numeric(im_sub["qvalue"], errors="coerce")
            im_sig = im_sub[im_sub["qvalue"] <= args.q_threshold]

            im_up = im_sig[im_sig["log2FC_detected_only"] > 0].sort_values("log2FC_detected_only", ascending=False).head(args.top_n_markers)
            if len(im_up):
                lines.append("- **intensity markers** (up-regulated):")
                lines.append("")
                lines.append("  | Gene | log2FC | qvalue |")
                lines.append("  |------|--------|--------|")
                for i in range(len(im_up)):
                    g = gene_name(im_up, i)
                    fc = im_up.iloc[i]["log2FC_detected_only"]
                    q = im_up.iloc[i]["qvalue"]
                    lines.append(f"  | {g} | {fc:+.2f} | {q:.1e} |")
                lines.append("")

            im_down = im_sig[im_sig["log2FC_detected_only"] < 0].sort_values("log2FC_detected_only", ascending=True).head(args.top_n_markers)
            if len(im_down):
                lines.append("- **depleted intensity markers** (down-regulated):")
                lines.append("")
                lines.append("  | Gene | log2FC | qvalue |")
                lines.append("  |------|--------|--------|")
                for i in range(len(im_down)):
                    g = gene_name(im_down, i)
                    fc = im_down.iloc[i]["log2FC_detected_only"]
                    q = im_down.iloc[i]["qvalue"]
                    lines.append(f"  | {g} | {fc:+.2f} | {q:.1e} |")
                lines.append("")

        # -- Scplainer markers table (up / down) --
        if sm is not None:
            sm_cl_col = "Cluster" if "Cluster" in sm.columns else "contrast"
            if sm_cl_col in sm.columns and "log2FC" in sm.columns:
                sm_sub = sm[sm[sm_cl_col].astype(str) == cl_str].copy()
                sm_sub["log2FC"] = pd.to_numeric(sm_sub["log2FC"], errors="coerce")
                sm_sub["qvalue"] = pd.to_numeric(sm_sub["qvalue"], errors="coerce")
                if "ok" in sm_sub.columns:
                    sm_sub = sm_sub[sm_sub["ok"].astype(str).str.upper().isin(["TRUE", "1"])]
                sm_sig = sm_sub[sm_sub["qvalue"] <= args.q_threshold]

                sm_up = sm_sig[sm_sig["log2FC"] > 0].sort_values("log2FC", ascending=False).head(args.top_n_markers)
                if len(sm_up):
                    lines.append("- **scplainer markers** (up-regulated):")
                    lines.append("")
                    lines.append("  | Gene | log2FC | qvalue |")
                    lines.append("  |------|--------|--------|")
                    for i in range(len(sm_up)):
                        g = gene_name(sm_up, i)
                        fc = sm_up.iloc[i]["log2FC"]
                        q = sm_up.iloc[i]["qvalue"]
                        lines.append(f"  | {g} | {fc:+.2f} | {q:.1e} |")
                    lines.append("")

                sm_down = sm_sig[sm_sig["log2FC"] < 0].sort_values("log2FC", ascending=True).head(args.top_n_markers)
                if len(sm_down):
                    lines.append("- **depleted scplainer markers** (down-regulated):")
                    lines.append("")
                    lines.append("  | Gene | log2FC | qvalue |")
                    lines.append("  |------|--------|--------|")
                    for i in range(len(sm_down)):
                        g = gene_name(sm_down, i)
                        fc = sm_down.iloc[i]["log2FC"]
                        q = sm_down.iloc[i]["qvalue"]
                        lines.append(f"  | {g} | {fc:+.2f} | {q:.1e} |")
                    lines.append("")

        # -- AUCell pathways table --
        if auc_mat is not None:
            cl_run_mask = np.array([run_to_cond.get(r, "") == cl_str for r in auc_run_order])
            if cl_run_mask.any():
                in_scores  = auc_mat[:, cl_run_mask].mean(axis=1)
                out_scores = auc_mat[:, ~cl_run_mask].mean(axis=1) if (~cl_run_mask).any() else np.zeros(auc_mat.shape[0])
                delta = in_scores - out_scores
                top_idx = np.argsort(delta)[::-1][:args.top_n_markers]
                lines.append("- **AUCell pathways**:")
                lines.append("")
                lines.append("  | Pathway | AUC_in | delta_vs_rest |")
                lines.append("  |---------|--------|---------------|")
                for i in top_idx:
                    name = auc_path_names[i][:40]
                    lines.append(f"  | {name} | {in_scores[i]:.3f} | {delta[i]:+.3f} |")
                lines.append("")

        # Researcher-provided markers: look up actual stats from marker files
        if custom_proteins:
            lines.append("- **researcher markers of interest** (per-marker stats):")
            lines.append("")
            lines.append("  | Gene | det_rate_in | det_rate_out | log2_OR | median_log2int_in | log2FC_int | det_qval | int_qval |")
            lines.append("  |------|-------------|--------------|---------|-------------------|------------|----------|----------|")

            for gene in custom_proteins:
                gl = gene.lower()

                # Look up detection markers
                det_row = None
                if dm is not None and "Cluster" in dm.columns:
                    sub = dm[dm["Cluster"].astype(str) == cl_str]
                    for gcol in ["Genes", "Gene"]:
                        if gcol in sub.columns:
                            match = sub[sub[gcol].astype(str).str.lower().str.split(";").str[0].str.strip() == gl]
                            if not match.empty:
                                det_row = match.iloc[0]
                                break

                # Look up intensity markers
                int_row = None
                if im is not None and "Cluster" in im.columns:
                    sub = im[im["Cluster"].astype(str) == cl_str]
                    for gcol in ["Genes", "Gene"]:
                        if gcol in sub.columns:
                            match = sub[sub[gcol].astype(str).str.lower().str.split(";").str[0].str.strip() == gl]
                            if not match.empty:
                                int_row = match.iloc[0]
                                break

                def fmt(val, col, decimals=3):
                    if val is None:
                        return "—"
                    v = val.get(col, None)
                    if v is None or (isinstance(v, float) and v != v):
                        return "—"
                    try:
                        return f"{float(v):.{decimals}f}"
                    except (ValueError, TypeError):
                        return "—"

                det_rate_in  = fmt(det_row, "det_rate_in", 2)
                det_rate_out = fmt(det_row, "det_rate_out", 2)
                log2_or      = fmt(det_row, "log2_odds_ratio", 2)
                det_q        = fmt(det_row, "qvalue", 4)
                med_int_in   = fmt(int_row, "median_in", 2)
                log2fc_int   = fmt(int_row, "log2FC_detected_only", 2)
                int_q        = fmt(int_row, "qvalue", 4)

                lines.append(f"  | {gene} | {det_rate_in} | {det_rate_out} | {log2_or} | {med_int_in} | {log2fc_int} | {det_q} | {int_q} |")

            lines.append("")
            lines.append(
                "  *det_rate: fraction of cells with protein detected (0–1); "
                "log2_OR: log2 odds ratio for detection specificity; "
                "median_log2int_in: median log2(intensity+1) in detected cells of this cluster; "
                "log2FC_int: log2 fold change vs other clusters (detected-only, positive=up in cluster); "
                "— = not tested (below min_cells threshold or gene not matched)*"
            )

        lines.append("")

    llm_path.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote LLM export -> {llm_path}")


if __name__ == "__main__":
    main()
