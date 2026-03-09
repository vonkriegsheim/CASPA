#!/usr/bin/env python3
import sys
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
"""
annotate_clusters_llm.py

Two-pass LLM cell type annotation for single-cell proteomics clusters.

Pass 1: Annotate clusters from cluster_llm_prompt.md. Ask the LLM which
        additional markers it needs that were missing (all "—") from the data.

Pivot lookup: For markers the LLM requests, compute raw per-cluster detection
        rates directly from scp_pivot.tsv (no statistical threshold applied).
        This recovers signal from proteins that were too sparse for the marker
        tests but are still biologically informative.

Pass 2: Re-annotate with the supplemental detection rates. Ask the LLM to:
        (a) confirm/revise cell type assignments
        (b) output a recommended marker panel for the NEXT SCP experiment,
            accounting for single-cell LFQ abundance constraints.

Outputs
-------
cluster_cell_type_annotations.tsv   -- final parsed annotation table
cluster_annotation_raw_pass1.txt    -- raw pass 1 response
cluster_annotation_raw_pass2.txt    -- raw pass 2 response
scp_recommended_markers.txt         -- LLM-suggested panel for next experiment

Requires:  pip install openai
Env vars:  ELM_API_KEY or OPENAI_API_KEY, ELM_BASE_URL
"""

import argparse
import csv
import datetime
import hashlib
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ---- Condition B: generic amendments + Round 0 ------------------------------

_AMENDMENT_1_DEVELOPMENTAL = """
## Developmental context

When the experiment context specifies a developmental stage (prenatal, embryonic,
neonatal, postnatal), constrain cell type labels to populations that are known to
exist at that stage. Do not assign mature/terminally differentiated cell type labels
(e.g., "astrocyte", "oligodendrocyte", "mature neuron") unless the cluster expresses
>=2 canonical markers of the mature state AND lacks progenitor markers (e.g., SOX2,
NES, PCNA, TOP2A). When mature markers are absent, prefer progenitor/precursor
terminology (e.g., "astroglial progenitor" rather than "astrocyte", "OPC" rather
than "oligodendrocyte", "neuroblast" rather than "neuron").

This applies regardless of what marker databases (e.g., PanglaoDB) associate with
a given protein — database entries typically reflect adult tissue, not developmental
stages.
""".strip()

_AMENDMENT_2_MECHANISTIC = """
## Mechanistic annotation rules

Do not assign a mechanistic or pathway-level qualifier (e.g., "IFN-activated",
"hypoxic", "senescent", "apoptotic") to a cell type label based on fewer than
2 concordant markers from that pathway, each detected in >=10% of the cluster's
cells. A single pathway-associated marker at low detection is insufficient
evidence for a mechanistic label — report it in the Contradictions/caveats
column instead.

When a cluster shows coordinated depletion of biosynthetic machinery (ribosomal
proteins, mitochondrial complex subunits, replication factors) without a coherent
activation signature, the most parsimonious explanation is sample handling or
dissociation damage. Describe the phenotype factually (e.g., "biosynthetically
depleted", "low-complexity proteome") rather than inferring a specific activation
mechanism.
""".strip()

_AMENDMENT_3_NONSELF = """
## Interpreting non-self proteins in mixed clusters

Before labelling any cluster as "contaminant" or "non-[expected cell type]",
consider whether the source cell type has the biological capacity to acquire
non-self proteins through any of the following mechanisms:

1. **Phagocytosis / efferocytosis**: Innate immune cells (neutrophils, macrophages,
   dendritic cells, microglia) routinely internalise dead cells, debris, and
   environmental material. A neutrophil or macrophage cluster showing keratin,
   desmosome, immunoglobulin, or tissue-specific markers may have phagocytosed
   this material rather than being a contaminant. Supporting evidence for
   phagocytic uptake includes: enrichment of lysosomal proteins (CTSD, LAMP1,
   LAMP2, SCARB2), proteasome subunits (PSMA/PSMB family), or
   immunosuppressive markers (ARG1).

2. **Lytic cell death / NETosis**: Neutrophils undergoing lytic NETosis expel
   their nuclear and granular contents, producing a protein-depleted remnant
   that may acquire complement and coagulation proteins from the vascular
   environment. A cluster showing coordinated depletion of histones, granule
   proteins, and metabolic enzymes — especially if combined with enrichment of
   complement or coagulation factors — is more consistent with a terminal
   functional state than with debris or contamination.

3. **Trogocytosis / membrane transfer**: Immune cells in vascular or tissue
   contact can acquire membrane fragments from neighbouring cells. An immune
   cell cluster showing tissue-specific markers (e.g., brain proteins in a
   tumour-infiltrating immune cell) at low-to-moderate levels may reflect
   tissue contact rather than tissue identity.

4. **Ambient carryover**: Highly abundant tissue proteins (digestive enzymes
   in pancreas, hemoglobin in vascularised tissue, keratins from sample
   handling) are detected ubiquitously at similar levels across all clusters.
   These are not informative for cell identity when detection rates are
   >90% across all clusters — use intensity fold-change to distinguish
   genuine enrichment from ambient background.

When non-self proteins co-occur with evidence of the mechanisms above, annotate
the cluster by its **functional identity** (e.g., "phagocytic neutrophil with
ingested epithelial material") and note the non-self proteins as acquired
cargo in the Contradictions/caveats column. Reserve the "contaminant" label
for clusters where no functional identity can be established AND no
acquisition mechanism is plausible.
""".strip()

_ROUND0_SYSTEM = """\
You are a single-cell proteomics expert preparing to annotate clusters from the \
following experiment. Before seeing any data, generate dataset-specific analytical \
constraints based on the experiment description alone.

Based on this experimental context, generate exactly four sections. Be specific and \
cite biological reasoning, but do NOT reference specific cluster IDs or per-cluster \
statistics (you have not seen the data yet). Be concise — maximum 600 words total.

### 1. Developmental / tissue vocabulary constraints
What mature cell types should NOT be used as labels given this tissue and developmental \
stage? What progenitor/precursor terminology should be used instead? What cell types \
ARE expected?

### 2. Expected ambient signals
Given the tissue source and sample preparation, what proteins are likely to be detected \
as ambient carryover in ALL clusters regardless of cell identity? Which are \
non-discriminative for cell type?

### 3. Non-self protein acquisition
Given the cell types expected, which populations can acquire non-self proteins? \
Through what mechanisms (phagocytosis, NETosis, trogocytosis)? What would this look \
like in the data?

### 4. Expected artefacts
What technical artefacts are common in this experiment type (dissociation stress, \
batch effects, doublets)? How would these manifest in the clustering?"""


def _extract_experiment_context(prompt_text: str) -> str:
    """Extract experiment context block (everything before cluster data / task instructions)."""
    lines = prompt_text.splitlines()
    context_lines = []
    for line in lines:
        if line.strip().startswith("You are a single-cell proteomics expert"):
            break
        context_lines.append(line)
    return "\n".join(context_lines).strip()


def _run_round0(client, model: str, context: str, max_tokens: int,
                thinking_budget: int = 0) -> str:
    """Call the LLM with only experiment context (no cluster data) to get Round 0 constraints."""
    print("[Round 0] Generating blind context constraints...")
    text, finish = _call_llm(
        client, model, _ROUND0_SYSTEM, context,
        max_tokens=max_tokens, temperature=0,
        thinking_budget=thinking_budget,
    )
    print(f"[Round 0] Complete ({len(text):,} chars, finish={finish})")
    return text


def _build_condition_b_prompt(prompt_text: str, round0_output: str) -> str:
    """
    Inject Amendments 1-3 and the Round 0 constraints into the base v2 prompt.

    Injection points (same as mark2/mark3 experiments):
      - Amendment 1 (developmental context): after run metadata sentinel line
      - Round 0 block: immediately after Amendment 1
      - Amendments 2+3 (mechanistic + non-self): after contamination block
    """
    # Amendment 1 + Round 0 after run metadata sentinel
    meta_sentinel = (
        "*These fields are populated by the calling script and stored alongside the\n"
        "raw LLM output for audit trail purposes.*"
    )
    r0_block = (
        "## Dataset-specific constraints (auto-generated — Round 0 blind reasoning)\n\n"
        "The following constraints were generated from the experiment description "
        "before any cluster data was examined. They capture expected biology, "
        "ambient signals, and artefacts specific to this experiment type.\n\n"
        + round0_output
    )
    if meta_sentinel in prompt_text:
        prompt_text = prompt_text.replace(
            meta_sentinel,
            meta_sentinel + "\n\n" + _AMENDMENT_1_DEVELOPMENTAL + "\n\n" + r0_block,
            1,
        )
    else:
        # Fallback: prepend to task instructions
        task_start = "You are a single-cell proteomics expert."
        prompt_text = prompt_text.replace(
            task_start,
            _AMENDMENT_1_DEVELOPMENTAL + "\n\n" + r0_block + "\n\n" + task_start,
            1,
        )

    # Amendments 2+3 after contamination block
    contam_end = "The researcher can decide post-hoc whether to filter these clusters."
    if contam_end in prompt_text:
        prompt_text = prompt_text.replace(
            contam_end,
            contam_end + "\n\n" + _AMENDMENT_2_MECHANISTIC + "\n\n" + _AMENDMENT_3_NONSELF,
            1,
        )

    return prompt_text


# ---- shared utilities -------------------------------------------------------

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
    return True


# ---- pivot lookup -----------------------------------------------------------

def load_pivot_gene_index(pivot_path: str) -> tuple[pd.DataFrame, list[str], list[str]]:
    """
    Load pivot and build a gene->row index.
    Returns (pivot_df, sample_cols, gene_labels_per_row)
    gene_labels_per_row: lowercase first gene name for each protein row.
    """
    pivot = pd.read_csv(pivot_path, sep="\t", dtype=str)
    sample_cols = [c for c in pivot.columns if is_sample_column(c)]
    gene_col = "Genes" if "Genes" in pivot.columns else None

    if gene_col:
        gene_labels = (pivot[gene_col]
                       .fillna("")
                       .astype(str)
                       .str.split(";")
                       .str[0]
                       .str.strip()
                       .str.lower()
                       .tolist())
    else:
        pg_col = "Protein.Group" if "Protein.Group" in pivot.columns else pivot.columns[0]
        gene_labels = pivot[pg_col].str.lower().tolist()

    return pivot, sample_cols, gene_labels


def cluster_detection_rates_for_genes(
    genes: list[str],
    pivot_path: str,
    annotation_path: str,
) -> dict[str, dict[str, str]]:
    """
    For each requested gene, compute per-cluster detection rate directly from
    the pivot (no statistical threshold). Returns:
        {gene: {cluster: "det_rate (n_detected/n_cells)"}}
    """
    pivot, sample_cols, gene_labels = load_pivot_gene_index(pivot_path)

    ann = pd.read_csv(annotation_path, sep="\t", dtype=str)
    run_to_cond = dict(zip(ann["Run"].astype(str), ann["Condition"].astype(str)))
    clusters = sorted(set(run_to_cond.values()))

    header_map = {c: stem_from_header(c) for c in sample_cols}
    run_set = set(ann["Run"].astype(str))
    matched = [c for c in sample_cols if header_map.get(c, "") in run_set]
    run_to_col = {header_map[c]: c for c in matched}
    runs = ann["Run"].astype(str).tolist()
    cols_ordered = [run_to_col.get(r) for r in runs]
    cond_arr = [run_to_cond.get(r, "") for r in runs]

    valid_pairs = [(i, c) for i, c in enumerate(cols_ordered) if c is not None]
    if not valid_pairs:
        return {}
    valid_idx, valid_cols = zip(*valid_pairs)
    valid_cond = [cond_arr[i] for i in valid_idx]

    X = pivot[list(valid_cols)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    result = {}
    for gene in genes:
        gl = gene.lower().strip()
        # Find matching rows (exact first, then prefix)
        row_indices = [i for i, g in enumerate(gene_labels) if g == gl]
        if not row_indices:
            row_indices = [i for i, g in enumerate(gene_labels) if g.startswith(gl)]
        if not row_indices:
            result[gene] = {cl: "not_found" for cl in clusters}
            continue

        # Use first matching row
        row_vals = X[row_indices[0], :]
        is_detected = np.isfinite(row_vals) & (row_vals > 0)

        cl_rates = {}
        for cl in clusters:
            cl_mask = np.array([c == cl for c in valid_cond])
            n_cells = cl_mask.sum()
            n_det = is_detected[cl_mask].sum() if n_cells > 0 else 0
            rate = n_det / n_cells if n_cells > 0 else 0.0
            cl_rates[cl] = f"{rate:.2f} ({int(n_det)}/{int(n_cells)})"

        result[gene] = cl_rates

    return result


# ---- LLM response parsing ---------------------------------------------------

def parse_annotation_table(text: str) -> list[dict]:
    """
    Parse the LLM annotation table.  Handles both v1 (4-col) and v2 (6-col) formats:
      v1: | Cluster | Suggested cell type | Key evidence | Confidence |
      v2: | Cluster | Suggested cell type | Key supporting markers | Contradictions / caveats | Confidence | Resolving markers (if Low/ambiguous) |
    """
    rows = []
    lines = text.splitlines()
    in_table = False
    header_passed = False
    n_cols = 0

    for line in lines:
        line = line.strip()
        if not line.startswith("|"):
            if in_table and rows:
                break
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]

        if not in_table:
            low = [c.lower() for c in cells]
            if any("cluster" in c for c in low) and any("cell" in c or "type" in c for c in low):
                in_table = True
                n_cols = len(cells)
                continue

        if in_table and not header_passed:
            if all(re.match(r"^[-: ]+$", c) for c in cells if c):
                header_passed = True
                continue

        if in_table and header_passed and len(cells) >= 3:
            row = {
                "Cluster":      cells[0].strip("* ").replace("**", ""),
                "cell_type":    cells[1].strip("* ").replace("**", ""),
            }
            if n_cols >= 6:
                # v2 format: 6 columns
                row["key_evidence"]      = cells[2] if len(cells) > 2 else ""
                row["contradictions"]    = cells[3] if len(cells) > 3 else ""
                row["confidence"]        = cells[4] if len(cells) > 4 else ""
                row["resolving_markers"] = cells[5] if len(cells) > 5 else ""
            else:
                # v1 format: 4 columns
                row["key_evidence"]      = cells[2] if len(cells) > 2 else ""
                row["contradictions"]    = ""
                row["confidence"]        = cells[3] if len(cells) > 3 else ""
                row["resolving_markers"] = ""
            rows.append(row)

    return rows


def parse_requested_markers(text: str) -> list[str]:
    """
    Extract the list of additional markers the LLM says it needs.
    Looks for a section tagged REQUESTED_MARKERS: gene1, gene2, ...
    """
    match = re.search(
        r"REQUESTED_MARKERS:\s*([^\n]+)", text, re.IGNORECASE
    )
    if not match:
        return []
    raw = match.group(1)
    genes = [g.strip().strip("*").strip() for g in re.split(r"[,;]+", raw)]
    return [g for g in genes if g and re.match(r"^[A-Za-z][A-Za-z0-9_\-]*$", g)]


def parse_recommended_panel(text: str) -> str:
    """
    Extract the RECOMMENDED_PANEL section from the pass 2 response.
    """
    match = re.search(
        r"RECOMMENDED_PANEL:(.*?)(?:\n[A-Z_]+:|$)", text, re.IGNORECASE | re.DOTALL
    )
    if match:
        return match.group(1).strip()
    return ""


# ---- main -------------------------------------------------------------------

def _call_llm(client, model: str, system: str, user: str,
              max_tokens: int, temperature: float,
              thinking_budget: int = 0) -> tuple[str, str]:
    """
    Call the LLM using whichever endpoint the model supports.
    Returns (response_text, finish_reason).
    GPT-5 series uses the Responses API (v1/responses).
    Claude models use the Anthropic Messages API.
    All other models use Chat Completions (v1/chat/completions).
    thinking_budget > 0 enables extended thinking for Claude models.
    """
    # Detect Anthropic client by type
    try:
        import anthropic as _anthropic_mod
        _is_anthropic = isinstance(client, _anthropic_mod.Anthropic)
    except ImportError:
        _is_anthropic = False

    if _is_anthropic:
        kwargs = dict(
            model=model,
            max_tokens=max_tokens,
            system=system,
            messages=[{"role": "user", "content": user}],
        )
        if thinking_budget > 0:
            kwargs["thinking"] = {"type": "enabled", "budget_tokens": thinking_budget}
            # Extended thinking requires temperature=1
        else:
            kwargs["temperature"] = temperature
        r = client.messages.create(**kwargs)
        # Extract only text blocks (skip thinking blocks)
        text = "".join(
            block.text for block in r.content
            if hasattr(block, "type") and block.type == "text"
        )
        return text, r.stop_reason

    use_responses_api = model.startswith("gpt-5") or model.startswith("o")

    if use_responses_api:
        r = client.responses.create(
            model=model,
            instructions=system,
            input=user,
            max_output_tokens=max_tokens,
        )
        return r.output_text, r.status
    else:
        r = client.chat.completions.create(
            model=model,
            max_tokens=max_tokens,
            temperature=temperature,
            messages=[
                {"role": "system", "content": system},
                {"role": "user",   "content": user},
            ],
        )
        return r.choices[0].message.content, r.choices[0].finish_reason


def parse_recommended_panel_genes(text: str) -> list[str]:
    """
    Extract a flat deduplicated gene list from the RECOMMENDED_PANEL table.
    Table rows: Cell type | Recommended markers | Rationale
    """
    panel_text = parse_recommended_panel(text)
    if not panel_text:
        return []
    genes = []
    for line in panel_text.splitlines():
        line = line.strip()
        if not line.startswith("|") or "---" in line:
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) >= 2:
            marker_cell = cells[1]
            if "recommended" in marker_cell.lower() or "marker" in marker_cell.lower():
                continue  # skip header row
            for g in re.split(r"[,;]+", marker_cell):
                g = g.strip().strip("*").strip()
                if g and re.match(r"^[A-Za-z][A-Za-z0-9_\-]*$", g) and "note" not in g.lower():
                    genes.append(g)
    seen = set()
    return [g for g in genes if not (g in seen or seen.add(g))]


def load_panglaodb(path: str, species: str = "mouse") -> dict[str, list[str]]:
    """
    Load PanglaoDB markers TSV(.gz) and return {cell_type_lower: [Gene1, Gene2, ...]}.
    Filters to specified species (Mm for mouse, Hs for human).
    """
    import gzip
    p = Path(path)
    if not p.exists():
        return {}

    opener = gzip.open if p.suffix == ".gz" else open
    cell_markers: dict[str, list[str]] = {}
    sp_code = "Mm" if species == "mouse" else "Hs"

    with opener(p, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sp = row.get("species", "")
            if sp_code not in sp:
                continue
            ct = row.get("cell type", "").strip().lower()
            gene = row.get("official gene symbol", "").strip()
            if ct and gene:
                # Convert to species-appropriate nomenclature
                if species == "mouse":
                    gene = gene[0].upper() + gene[1:].lower() if len(gene) > 1 else gene.upper()
                cell_markers.setdefault(ct, []).append(gene)

    return cell_markers


def compute_marker_coverage(
    annotations: list[dict],
    panglaodb_markers: dict[str, list[str]],
    pivot_path: str,
    annotation_path: str,
    det_rate_threshold: float = 0.20,
) -> list[dict]:
    """
    For each annotated cluster, compute marker coverage score:
      coverage = (markers detected at >= det_rate_threshold) / (total panel markers)
    Returns the annotations list with 'marker_coverage' and 'panel_matched' added.
    """
    if not panglaodb_markers or not Path(pivot_path).exists():
        for row in annotations:
            row["marker_coverage"] = ""
            row["panel_matched"] = ""
        return annotations

    pivot, sample_cols, gene_labels = load_pivot_gene_index(pivot_path)
    ann = pd.read_csv(annotation_path, sep="\t", dtype=str)
    run_to_cond = dict(zip(ann["Run"].astype(str), ann["Condition"].astype(str)))

    header_map = {c: stem_from_header(c) for c in sample_cols}
    run_set = set(ann["Run"].astype(str))
    matched = [c for c in sample_cols if header_map.get(c, "") in run_set]
    run_to_col = {header_map[c]: c for c in matched}
    runs = ann["Run"].astype(str).tolist()
    cols_ordered = [run_to_col.get(r) for r in runs]
    cond_arr = [run_to_cond.get(r, "") for r in runs]

    valid_pairs = [(i, c) for i, c in enumerate(cols_ordered) if c is not None]
    if not valid_pairs:
        for row in annotations:
            row["marker_coverage"] = ""
            row["panel_matched"] = ""
        return annotations
    valid_idx, valid_cols = zip(*valid_pairs)
    valid_cond = [cond_arr[i] for i in valid_idx]
    X = pivot[list(valid_cols)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)

    # Build gene name -> row index map
    gene_to_row = {}
    for i, g in enumerate(gene_labels):
        if g and g not in gene_to_row:
            gene_to_row[g] = i

    for row in annotations:
        ct = row.get("cell_type", "").strip().lower()
        cluster = row.get("Cluster", "")

        # Fuzzy match: find best PanglaoDB cell type for this annotation
        best_ct = None
        best_n = 0
        for pdb_ct, genes in panglaodb_markers.items():
            # Check if annotation label contains the panglaodb cell type or vice versa
            if pdb_ct in ct or ct in pdb_ct or any(
                w in ct for w in pdb_ct.split() if len(w) > 3
            ):
                if len(genes) > best_n:
                    best_ct = pdb_ct
                    best_n = len(genes)

        if not best_ct:
            row["marker_coverage"] = ""
            row["panel_matched"] = ""
            continue

        panel_genes = panglaodb_markers[best_ct]
        cl_mask = np.array([c == cluster for c in valid_cond])
        n_cells = cl_mask.sum()

        if n_cells == 0:
            row["marker_coverage"] = ""
            row["panel_matched"] = best_ct
            continue

        n_detected = 0
        n_total = len(panel_genes)
        for gene in panel_genes:
            gl = gene.lower()
            ridx = gene_to_row.get(gl)
            if ridx is None:
                continue
            vals = X[ridx, :][cl_mask]
            det = np.sum(np.isfinite(vals) & (vals > 0))
            rate = det / n_cells
            if rate >= det_rate_threshold:
                n_detected += 1

        coverage = n_detected / n_total if n_total > 0 else 0.0
        row["marker_coverage"] = f"{coverage:.2f} ({n_detected}/{n_total})"
        row["panel_matched"] = best_ct

    return annotations


def _nomenclature_instruction(species: str) -> str:
    if species == "mouse":
        return ("Use mouse gene nomenclature (title case: Chga, Vim, Acta2, Krt19). "
                "Do NOT use human all-caps nomenclature (CHGA, VIM, ACTA2, KRT19).")
    return ("Use human gene nomenclature (all-caps: CHGA, VIM, ACTA2, KRT19). "
            "Do NOT use mouse title-case nomenclature (Chga, Vim, Acta2, Krt19).")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prompt-md",    required=True)
    ap.add_argument("--pivot",        required=True,
                    help="scp_pivot.tsv — used for sub-threshold marker lookup")
    ap.add_argument("--annotation",   required=True,
                    help="scp_annotation.tsv — Run/Condition mapping")
    ap.add_argument("--model",        default="gpt-4-turbo",
                    help="Model name as exposed by ELM (default: gpt-4-turbo)")
    ap.add_argument("--base-url",     default="",
                    help="ELM API base URL. Can also be set via ELM_BASE_URL env var. "
                         "Try https://elm.edina.ac.uk/ or https://elm.edina.ac.uk/v1")
    ap.add_argument("--max-tokens",   type=int, default=8000)
    ap.add_argument("--temperature",  type=float, default=0.0)
    ap.add_argument("--out-tsv",      required=True)
    ap.add_argument("--species", default="human", choices=["human", "mouse"],
                    help="Species for gene nomenclature (human=all-caps, mouse=title-case)")
    ap.add_argument("--experiment-context", default=None,
                    help="Free-text description of the experiment to inform marker selection. E.g. 'Murine pancreas, 2-day caerulein pancreatitis, single-cell LFQ proteomics (Bruker timsTOF)'. Injected into both LLM system prompts.")
    ap.add_argument("--out-markers-txt", default=None,
                    help="Path for flat gene list from recommended panel "
                         "(one gene per line, for --proteins in plot_umap_overlays.py). "
                         "Default: scp_recommended_markers_flat.txt beside --out-tsv.")
    ap.add_argument("--out-dir",      default=None,
                    help="Directory for auxiliary outputs. Defaults to --out-tsv directory.")
    ap.add_argument("--panglaodb",   default=None,
                    help="Path to PanglaoDB_markers.tsv.gz for cross-validation coverage scoring. "
                         "Default: pipeline/assets/panglaodb_markers.tsv.gz beside script.")
    ap.add_argument("--provider",   default="openai", choices=["openai", "claude"],
                    help="LLM provider: 'openai' (default) or 'claude' (Anthropic).")
    ap.add_argument("--thinking-budget", type=int, default=0,
                    help="Extended thinking token budget for Claude models (0 = disabled). "
                         "Recommended: 8000-16000 for annotation tasks.")
    ap.add_argument("--condition-b", action="store_true",
                    help="Enable Condition B prompting: inject 3 generic amendments "
                         "(developmental context, mechanistic threshold, non-self proteins) "
                         "plus a Round 0 blind context reasoning call into the prompt "
                         "before annotation. Round 0 uses the same model and API key.")
    args = ap.parse_args()

    if args.provider == "claude":
        api_key = os.environ.get("ANTHROPIC_API_KEY", "").strip()
        if not api_key:
            # Also try machine.json scp.anthropic_key
            machine_json = Path(__file__).parent.parent.parent.parent / "machine.json"
            if machine_json.exists():
                import json as _json
                m = _json.load(open(machine_json))
                api_key = m.get("scp", {}).get("anthropic_key", "").strip()
        if not api_key:
            print("WARNING: No ANTHROPIC_API_KEY found. Skipping LLM annotation.",
                  file=sys.stderr)
            Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
            Path(args.out_tsv).write_text(
                "Cluster\tcell_type\tkey_evidence\tcontradictions\tconfidence\tresolving_markers\n",
                encoding="utf-8")
            sys.exit(0)
    else:
        api_key = os.environ.get("ELM_API_KEY", "").strip()
        if not api_key:
            api_key = os.environ.get("OPENAI_API_KEY", "").strip()
        if not api_key:
            print("WARNING: No API key found (ELM_API_KEY or OPENAI_API_KEY). "
                  "Skipping LLM annotation.", file=sys.stderr)
            # Write empty output so Snakemake is satisfied
            Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
            Path(args.out_tsv).write_text(
                "Cluster\tcell_type\tkey_evidence\tcontradictions\tconfidence\tresolving_markers\n",
                encoding="utf-8")
            sys.exit(0)

    # ---- Build experiment context block ----
    nomenclature = _nomenclature_instruction(args.species)
    context_block = ""
    if args.experiment_context:
        context_block = (
            f"\n\nEXPERIMENT CONTEXT:\n{args.experiment_context}\n\n"
            "Use this context to:\n"
            "- Anticipate which cell types are expected and in what proportions\n"
            "- Account for disease/treatment-induced proteomic shifts\n"
            "- Prioritise markers known to be detectable by single-cell LFQ "
              "in this specific tissue and condition\n"
            "- Flag markers that are canonical scRNA-seq markers but are "
              "unlikely to be detected by LFQ proteomics at single-cell depth "
              "(e.g. low-abundance hormones, transcription factors)\n"
            f"NOMENCLATURE: {nomenclature}\n"
        )

    if args.provider == "claude":
        try:
            import anthropic as _anthropic_mod
        except ImportError:
            print("ERROR: pip install anthropic", file=sys.stderr)
            sys.exit(1)
        print(f"Using Anthropic Claude endpoint (model={args.model}, "
              f"thinking_budget={args.thinking_budget})")
        client = _anthropic_mod.Anthropic(api_key=api_key)
    else:
        try:
            from openai import OpenAI
        except ImportError:
            print("ERROR: pip install openai", file=sys.stderr)
            sys.exit(1)

        base_url = (os.environ.get("ELM_BASE_URL", "").strip()
                    or args.base_url.strip())

        if base_url:
            print(f"Using custom base URL: {base_url}")
            client = OpenAI(api_key=api_key, base_url=base_url)
        else:
            print("Using default OpenAI endpoint (api.openai.com)")
            client = OpenAI(api_key=api_key)

    # Guard: when thinking is enabled max_tokens must exceed the budget;
    # bump silently so the Anthropic API doesn't reject the call.
    if args.thinking_budget > 0 and args.max_tokens <= args.thinking_budget:
        args.max_tokens = args.thinking_budget + 4096
        print(f"INFO: max_tokens bumped to {args.max_tokens} "
              f"(thinking_budget={args.thinking_budget} + 4096 output headroom)")

    out_dir = Path(args.out_dir) if args.out_dir else Path(args.out_tsv).parent
    out_dir.mkdir(parents=True, exist_ok=True)

    prompt_text = Path(args.prompt_md).read_text(encoding="utf-8")
    n_clusters = prompt_text.count("## Cluster")
    print(f"Loaded prompt: {n_clusters} clusters, {len(prompt_text):,} chars")

    round0_text = ""  # populated below if condition_b; used in Pass 2 preamble

    # ---- Condition B: Round 0 + generic amendments ----
    if args.condition_b:
        context = _extract_experiment_context(prompt_text)
        round0_text = _run_round0(
            client, args.model, context,
            max_tokens=args.max_tokens,
            thinking_budget=args.thinking_budget,
        )
        round0_path = out_dir / "round0_constraints.md"
        round0_path.write_text(
            f"# Round 0 constraints (auto-generated)\n\n"
            f"*Generated: {datetime.datetime.now(datetime.timezone.utc).isoformat()}*\n"
            f"*Model: {args.model}*\n\n"
            + round0_text,
            encoding="utf-8",
        )
        print(f"Wrote Round 0 constraints -> {round0_path}")
        prompt_text = _build_condition_b_prompt(prompt_text, round0_text)
        print(f"Condition B prompt built ({len(prompt_text):,} chars)")

    # ---- v2: Compute prompt hash and fill metadata placeholders ----
    run_timestamp = datetime.datetime.now(datetime.timezone.utc).isoformat()
    prompt_hash = hashlib.sha256(prompt_text.encode("utf-8")).hexdigest()

    # Fill metadata placeholders in the prompt template
    prompt_text = (prompt_text
                   .replace("{model_name}", args.model)
                   .replace("{iso_timestamp}", run_timestamp)
                   .replace("{sha256_of_prompt_text}", prompt_hash))

    # Write metadata sidecar JSON for audit trail
    run_metadata = {
        "prompt_version": "v2_condB" if args.condition_b else "v2",
        "model": args.model,
        "temperature": args.temperature,
        "condition_b": args.condition_b,
        "run_timestamp": run_timestamp,
        "prompt_sha256": prompt_hash,
        "prompt_file": str(Path(args.prompt_md).resolve()),
        "prompt_chars": len(prompt_text),
        "n_clusters": n_clusters,
    }
    metadata_path = out_dir / "annotation_run_metadata.json"
    metadata_path.write_text(json.dumps(run_metadata, indent=2), encoding="utf-8")
    print(f"Wrote run metadata -> {metadata_path}")

    # =========================================================================
    # PASS 1 — annotate + ask for missing markers
    # =========================================================================
    system_p1 = (
        "You are a single-cell proteomics expert. "
        "Your input is quantitative protein detection and intensity data from label-free "
        "single-cell proteomics. Proteins shown as '—' were not detected in enough cells "
        "to run statistics — this is a fundamental abundance limitation of single-cell LFQ "
        "proteomics, NOT necessarily evidence of absence.\n\n"
        "CONTAMINATION GUIDANCE — BE CONSERVATIVE:\n"
        "In single-cell proteomics, keratins (KRT5, KRT14, KRT10), serum proteins (ALB, "
        "TF, HP), and hemoglobin (HBB, HBA1) are commonly detected as ambient carryover "
        "from sample preparation, NOT as true contaminant cell populations. Do NOT label "
        "a cluster as 'Contaminant' solely because it expresses these proteins. Only call "
        "contamination when the cluster lacks ANY coherent lineage markers across ALL "
        "modalities AND contaminant markers dominate >80% of top markers AND no biological "
        "explanation is plausible. When contaminant-associated proteins co-occur with "
        "genuine lineage markers, annotate the cluster by its biological identity and note "
        "ambient contamination as a caveat. Prefer a low-confidence biological annotation "
        "over a contaminant label.\n\n"
        "Respond with THREE sections in this exact order:\n\n"
        "1. An annotation table with EXACTLY this header and one row per cluster:\n"
        "   | Cluster | Suggested cell type | Key supporting markers | "
        "Contradictions / caveats | Confidence | Resolving markers (if Low/ambiguous) |\n"
        "   Confidence: High / Medium-high / Medium / Low\n"
        "   (use the rubric in the prompt)\n\n"
        "2. A per-cluster narrative (one paragraph each) with quantitative citations "
        "   (e.g. \"Cd68 det_rate_in=1.00 vs 0.25, q=0.0000\"). "
        "   Do not assert a cell type without citing at least one specific statistic.\n\n"
        "3. A single line listing markers you would need to improve low-confidence "
        "   assignments. Format EXACTLY as:\n"
        "   REQUESTED_MARKERS: Gene1, Gene2, Gene3\n"
        f"   CRITICAL: {nomenclature} "
        "   Only list genes not already shown in the data or shown as '—'. "
        "   Prioritise proteins known to be detectable by LFQ proteomics at single-cell "
        "   level (i.e. high-abundance structural or metabolic proteins, not "
        "   low-abundance signalling peptides or transcription factors).\n\n"
        "Do not include any preamble before the table."
    )

    # Save prompts for reproducibility / debugging
    (out_dir / "prompt_pass1_system.txt").write_text(system_p1, encoding="utf-8")
    (out_dir / "prompt_pass1_user.txt").write_text(context_block + prompt_text, encoding="utf-8")

    print(f"\nPass 1: calling {args.model}...")
    text_p1, finish_p1 = _call_llm(
        client, args.model, system_p1,
        context_block + prompt_text,
        args.max_tokens, args.temperature,
        thinking_budget=args.thinking_budget,
    )
    (out_dir / "cluster_annotation_raw_pass1.txt").write_text(text_p1, encoding="utf-8")
    print(f"Pass 1 complete ({len(text_p1):,} chars, finish={finish_p1})")

    rows_p1 = parse_annotation_table(text_p1)
    print(f"Parsed {len(rows_p1)} annotations from pass 1")

    requested = parse_requested_markers(text_p1)
    print(f"LLM requested {len(requested)} additional markers: {requested}")

    # =========================================================================
    # PIVOT LOOKUP — compute raw detection rates for requested markers
    # =========================================================================
    supplemental_text = ""
    if requested and Path(args.pivot).exists() and Path(args.annotation).exists():
        print("Looking up sub-threshold detection rates in pivot...")
        rates = cluster_detection_rates_for_genes(requested, args.pivot, args.annotation)

        lines = ["\n\n## Supplemental: raw per-cluster detection rates for requested markers"]
        lines.append(
            "*These proteins were below the min_cells_detected threshold for statistical "
            "testing. Rates computed directly from the pivot matrix (detected = log2 value "
            "present and > 0). No statistical test applied.*\n"
        )

        # Get cluster list from annotation
        ann = pd.read_csv(args.annotation, sep="\t", dtype=str)
        clusters = sorted(set(ann["Condition"].astype(str)))

        header = "| Gene | " + " | ".join(clusters) + " |"
        sep    = "|------|" + "|".join(["------"] * len(clusters)) + "|"
        lines.append(header)
        lines.append(sep)

        for gene in requested:
            cl_vals = rates.get(gene, {})
            row = f"| {gene} | " + " | ".join(cl_vals.get(cl, "—") for cl in clusters) + " |"
            lines.append(row)

        lines.append(
            "\n*Format: detection_rate (n_detected/n_cells_in_cluster). "
            "not_found = gene name not matched in pivot.*"
        )
        supplemental_text = "\n".join(lines)
        print(supplemental_text[:500])

    # =========================================================================
    # PASS 2 — refine with supplemental data + generate recommended panel
    # =========================================================================
    system_p2 = (
        "You are a single-cell proteomics expert. "
        "You previously annotated clusters from label-free single-cell proteomics data. "
        "You are now provided with supplemental raw detection rates for the markers you "
        "requested. These rates bypass the statistical detection threshold — a detection "
        "rate of even 10-30% in a cluster can be biologically informative for low-abundance "
        "cell-type markers in SCP data.\n\n"
        "CONTAMINATION GUIDANCE — BE CONSERVATIVE:\n"
        "Keratins, ALB, and HBB are common ambient carryover in SCP, not proof of "
        "contaminant cells. Only call contamination when no coherent lineage markers exist "
        "and contaminant markers overwhelmingly dominate. If a cluster you previously "
        "labelled as contaminant now shows any lineage-specific signal in the supplemental "
        "data, reclassify it with a biological label (even at Low confidence).\n\n"
        "Respond with THREE sections in this exact order:\n\n"
        "1. A REVISED annotation table with EXACTLY this header:\n"
        "   | Cluster | Suggested cell type | Key supporting markers | "
        "Contradictions / caveats | Confidence | Resolving markers (if Low/ambiguous) |\n"
        "   For unchanged clusters, write \"No change from Pass 1\" in the Key supporting "
        "   markers column and carry forward the original values.\n\n"
        "2. Brief rationale per cluster, noting where supplemental marker data "
        "   changed or confirmed your original assessment. Cite specific detection rates.\n\n"
        "3. A recommended marker panel for the NEXT single-cell proteomics experiment "
        "   on this tissue type. This panel should:\n"
        "   - Prioritise proteins confirmed as detectable in this dataset\n"
        "   - Include the most discriminative markers per cell type identified\n"
        "   - Note where canonical markers (e.g. low-abundance hormones) are unlikely "
        "     to be detected by LFQ and suggest higher-abundance surrogates\n"
        "   - Flag any cell types that remain ambiguous and need additional markers\n"
        "   Format the panel as a markdown table with EXACTLY this structure:\n"
        "   RECOMMENDED_PANEL:\n"
        "   | Cell type | Recommended markers | Rationale |\n"
        "   | --- | --- | --- |\n"
        "   | Stellate activated | Vim, Des, S100a6 | Confirmed detectable at 60-90%... |\n"
        "   (one data row per cell type identified)\n"
        f"   CRITICAL: {nomenclature}\n\n"
        "Do not include any preamble before the table."
    )

    # Condition B: re-surface R0 constraints at the front of Pass 2 so they are
    # not lost in the middle of the ~50 KB re-sent cluster prompt (primacy effect).
    r0_preamble = ""
    if args.condition_b and round0_text:
        r0_preamble = (
            "## Key constraints for revision (Round 0 — read before interpreting supplemental data)\n\n"
            "These constraints were generated from the experiment context alone, before any "
            "cluster data was seen. Apply them when deciding whether supplemental detection "
            "rates change any annotation:\n\n"
            + round0_text
            + "\n\n---\n\n"
        )

    pass2_user = (
        f"{r0_preamble}"
        f"Here is my original cluster data:\n\n{context_block}{prompt_text}"
        f"{supplemental_text}\n\n"
        f"Your pass 1 annotations were:\n{text_p1}\n\n"
        "Please now provide your refined annotation and recommended marker panel."
    ) if supplemental_text else (
        f"{r0_preamble}"
        f"Here is my original cluster data:\n\n{context_block}{prompt_text}\n\n"
        f"Your pass 1 annotations were:\n{text_p1}\n\n"
        "No additional markers could be retrieved from the pivot. "
        "Please finalise your annotation and provide the recommended panel."
    )

    (out_dir / "prompt_pass2_system.txt").write_text(system_p2, encoding="utf-8")
    (out_dir / "prompt_pass2_user.txt").write_text(pass2_user, encoding="utf-8")

    print(f"\nPass 2: calling {args.model}...")
    text_p2, finish_p2 = _call_llm(
        client, args.model, system_p2, pass2_user,
        args.max_tokens, args.temperature,
        thinking_budget=args.thinking_budget,
    )
    (out_dir / "cluster_annotation_raw_pass2.txt").write_text(text_p2, encoding="utf-8")
    print(f"Pass 2 complete ({len(text_p2):,} chars, finish={finish_p2})")

    rows_p2 = parse_annotation_table(text_p2)
    final_rows = rows_p2 if rows_p2 else rows_p1
    print(f"Parsed {len(final_rows)} annotations from pass 2")

    # =========================================================================
    # WRITE OUTPUTS
    # =========================================================================

    # Annotation TSV
    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["Cluster", "cell_type", "key_evidence", "contradictions",
                  "confidence", "resolving_markers"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in final_rows:
            writer.writerow(row)
    print(f"Wrote annotations -> {out_path}")

    # ---- PanglaoDB cross-validation: marker coverage score ----
    panglaodb_path = args.panglaodb
    if not panglaodb_path:
        # Default: look beside script in pipeline/assets/
        script_dir = Path(__file__).resolve().parent
        default_path = script_dir.parent.parent / "assets" / "panglaodb_markers.tsv.gz"
        if default_path.exists():
            panglaodb_path = str(default_path)

    if panglaodb_path and Path(panglaodb_path).exists():
        print(f"Loading PanglaoDB markers for cross-validation: {panglaodb_path}")
        pdb_markers = load_panglaodb(panglaodb_path, species=args.species)
        print(f"  Loaded {sum(len(v) for v in pdb_markers.values())} markers across "
              f"{len(pdb_markers)} cell types")
        final_rows = compute_marker_coverage(
            final_rows, pdb_markers, args.pivot, args.annotation,
        )
        # Re-write TSV with coverage columns
        fieldnames_ext = fieldnames + ["marker_coverage", "panel_matched"]
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames_ext, delimiter="\t")
            writer.writeheader()
            for row in final_rows:
                writer.writerow(row)
        print(f"Updated annotations with marker coverage -> {out_path}")
    else:
        print("NOTE: PanglaoDB markers not found; skipping coverage scoring")

    # Update metadata sidecar with final results
    run_metadata["n_annotations_parsed"] = len(final_rows)
    run_metadata["pass2_used"] = bool(rows_p2)
    run_metadata["panglaodb_used"] = bool(panglaodb_path and Path(panglaodb_path).exists())
    metadata_path.write_text(json.dumps(run_metadata, indent=2), encoding="utf-8")

    # Recommended panel
    panel = parse_recommended_panel(text_p2)
    if panel:
        panel_path = out_dir / "scp_recommended_markers.txt"
        panel_path.write_text(panel, encoding="utf-8")
        print(f"Wrote recommended panel -> {panel_path}")
        # Also write flat gene list for direct use by plot_umap_overlays.py --proteins
        flat_genes = parse_recommended_panel_genes(text_p2)
        if flat_genes:
            markers_txt_path = (
                Path(args.out_markers_txt)
                if args.out_markers_txt
                else Path(args.out_tsv).parent / "scp_recommended_markers_flat.txt"
            )
            markers_txt_path.write_text("\n".join(flat_genes), encoding="utf-8")
            print(f"Wrote flat gene list ({len(flat_genes)} genes) -> {markers_txt_path}")
        else:
            print("NOTE: could not extract gene list from recommended panel")
    else:
        print("NOTE: could not parse RECOMMENDED_PANEL section — check raw pass 2 response")

    # Summary to stdout
    print("\nFinal annotation summary:")
    for row in final_rows:
        contra = row.get("contradictions", "")[:30]
        resolv = row.get("resolving_markers", "")[:30]
        print(f"  {row['Cluster']:8s} | {row['cell_type'][:40]:40s} | {row['confidence']:12s} | {contra}")


if __name__ == "__main__":
    main()
