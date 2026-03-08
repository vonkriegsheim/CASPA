# ==============================================================================
# CASPA — Context-Aware Single Cell Proteomic Analysis
# Snakemake >= 7.0
#
# INVOCATION (from your workdir):
#   python /path/to/caspa/run.py --workdir /path/to/MyExperiment --cores 30
#
#   OR directly:
#   snakemake \
#     --snakefile /path/to/CASPA/Snakefile \
#     --directory  /path/to/MyExperiment \
#     --cores      30 \
#     --keep-going
#
# CONFIG LOAD ORDER (last file wins):
#   1. caspa_defaults.json   — bundled lab defaults (scp, plots, enrichment)
#   2. config/caspa.json     — per-experiment overrides
#
# INPUT OPTIONS:
#   A) Pre-exported pg_matrix (DIA-NN, FragPipe):
#      config/caspa.json -> input.pg_matrix = "/path/to/report.pg_matrix.tsv"
#   B) Spectronaut long-format TSV:
#      config/caspa.json -> input.spectronaut_tsv = "/path/to/spectronaut_export.tsv"
#      CASPA will convert it to pg_matrix format automatically.
#
# REQUIRED INPUTS:
#   config/caspa.json        — per-experiment config (see config/caspa.json.template)
#   config/ms_inputs.tsv     — sample sheet: sample_id | sample_file | batch
#
# OUTPUT:
#   scp/                     — all analysis outputs
# ==============================================================================

import os

# ── Config loading ──────────────────────────────────────────────────────────────
configfile: os.path.join(workflow.basedir, "caspa_defaults.json")
configfile: "config/caspa.json"

# ── Script roots — relative to Snakefile via workflow.basedir ──────────────────
PY  = os.path.join(workflow.basedir, "pipeline/scripts/python")
RSC = os.path.join(workflow.basedir, "pipeline/scripts/R")

# ── Config convenience helpers ──────────────────────────────────────────────────
def _c(*keys, default=None):
    v = config
    for k in keys:
        if not isinstance(v, dict):
            return default
        v = v.get(k)
        if v is None:
            return default
    return v if v is not None else default

def _fmt():
    """First entry of plots.formats list, e.g. 'pdf'."""
    fmts = _c("plots", "formats", default=["pdf"])
    return fmts[0] if isinstance(fmts, list) else str(fmts)

def _gsea_collections():
    """Comma-separated string of GSEA collections for CLI flags."""
    gs = _c("enrichment", "gsea", "genesets", default=["hallmark", "reactome"])
    return ",".join(gs) if isinstance(gs, list) else str(gs)

# ── Sentinel helper ─────────────────────────────────────────────────────────────
def sentinel(path):
    return os.path.join(path, ".done")

# ── Input routing: direct pg_matrix OR Spectronaut → conversion ─────────────────
_INPUT_IS_SPECTRONAUT = bool(_c("input", "spectronaut_tsv"))
INPUT_PG_MATRIX = (
    "caspa_input/pg_matrix.tsv" if _INPUT_IS_SPECTRONAUT
    else _c("input", "pg_matrix")
)
INPUT_MANIFEST = _c("input", "sample_sheet", default="config/ms_inputs.tsv")

# ── All-targets collector ───────────────────────────────────────────────────────
def all_targets():
    t = []
    t.append("scp/scplainer/scplainer_cluster_DA.tsv")
    t.append(sentinel("scp/scplainer/plots"))
    t.append("scp/markers/scplainer_intensity_markers.tsv")
    t.append("scp/llm/cluster_cell_type_annotations.tsv")
    t.append("scp/llm/plots/umap_cell_types.pdf")
    t.append("scp/viz/plots/harmony_before_after.pdf")
    t.append(sentinel("scp/viz/plots/normalisation"))
    t.append(sentinel("scp/viz/plots/marker_volcanos"))
    t.append("scp/viz/tables/cluster_centroid_correlation.tsv")
    t.append(sentinel("scp/viz/plots/detection_matrix"))
    t.append("scp/markers/marker_dotplot_matrix.pdf")
    t.append(sentinel("scp/enrichment"))
    return t

# ── Master target ───────────────────────────────────────────────────────────────
rule all:
    input: all_targets()


# ==============================================================================
# OPTIONAL: Spectronaut → pg_matrix conversion
# ==============================================================================

if _INPUT_IS_SPECTRONAUT:
    rule caspa_spectronaut_convert:
        input:
            spectronaut_tsv = _c("input", "spectronaut_tsv"),
        output:
            pg_matrix = "caspa_input/pg_matrix.tsv",
        shell:
            "python {PY}/spectronaut_to_diann_pg_matrix.py"
            " --input  {input.spectronaut_tsv}"
            " --output {output.pg_matrix}"


# ==============================================================================
# SCP RULES
# ==============================================================================

rule scp_qc_filter:
    input:
        pg_matrix = INPUT_PG_MATRIX,
        manifest  = INPUT_MANIFEST,
    output:
        filtered_pg       = "scp/qc/filtered.pg_matrix.tsv",
        filtered_manifest = "scp/qc/filtered.manifest.tsv",
        qc_report         = "scp/qc/scp_qc_report.json",
    params:
        exclude_regex = _c("scp", "exclude_run_regex", default=r"^(library|lib|5cell|10cell|100cell)"),
        min_detected  = _c("scp", "min_protein_ids",   default=500),
        bottom_n      = _c("scp", "qc_bottom_n",       default=5),
        multiplier    = _c("scp", "qc_multiplier",     default=1.7),
    shell:
        "python {PY}/scp_qc_filter_cells_pg_matrix.py"
        " --pg-matrix                {input.pg_matrix}"
        " --manifest                 {input.manifest}"
        " --outdir                   scp/qc"
        " --run-col                  sample_id"
        " --exclude-run-regex        \"{params.exclude_regex}\""
        " --min-detected             {params.min_detected}"
        " --bottom-n                 {params.bottom_n}"
        " --multiplier               {params.multiplier}"
        " --out-filtered-pg-matrix   {output.filtered_pg}"
        " --out-filtered-manifest    {output.filtered_manifest}"

rule scp_pivot:
    input:
        pg_matrix  = "scp/qc/filtered.pg_matrix.tsv",
        annotation = "scp/qc/filtered.manifest.tsv",
    output:
        pivot  = "scp/pivot_pack.tsv",
        shifts = "scp/pivot_shifts.tsv",
        report = "scp/pivot_report.json",
    shell:
        "python {PY}/scp_make_pivot_from_pg_matrix.py"
        " --pg-matrix   {input.pg_matrix}"
        " --manifest    {input.annotation}"
        " --run-col     sample_id"
        " --out-pivot   {output.pivot}"
        " --out-shifts  {output.shifts}"
        " --out-report  {output.report}"

rule scp_joint_embedding:
    input:
        pg_matrix  = "scp/pivot_pack.tsv",
        annotation = "scp/qc/filtered.manifest.tsv",
    output:
        scp_annotation = "scp/clustering/scp_annotation.tsv",
        embeddings     = "scp/clustering/scp_cluster_assignments.tsv",
        report         = "scp/clustering/scp_clustering_report.json",
    params:
        n_pcs_int   = _c("scp", "joint_embedding", "n_pcs",             default=20),
        n_neighbors = _c("scp", "joint_embedding", "n_neighbors",       default=15),
        resolution  = _c("scp", "joint_embedding", "leiden_resolution", default=0.8),
        seed        = _c("scp", "joint_embedding", "seed",              default=0),
    shell:
        "python {PY}/scp_joint_embedding_leiden.py"
        " --pg-matrix         {input.pg_matrix}"
        " --manifest          {input.annotation}"
        " --run-col           sample_id"
        " --outdir            scp/clustering"
        " --n_pcs_int         {params.n_pcs_int}"
        " --n_neighbors       {params.n_neighbors}"
        " --leiden_resolution {params.resolution}"
        " --seed              {params.seed}"

rule scp_detection_markers:
    input:
        pg_matrix  = "scp/pivot_pack.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
    output:
        markers = "scp/markers/detection_markers.tsv",
    shell:
        "python {PY}/scp_markers_detection.py"
        " --matrix     {input.pg_matrix}"
        " --annotation {input.annotation}"
        " --out        {output.markers}"

rule scp_intensity_markers:
    input:
        pg_matrix  = "scp/pivot_pack.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
    output:
        markers = "scp/markers/intensity_markers_detected_only.tsv",
    shell:
        "python {PY}/scp_markers_intensity_detected_only.py"
        " --matrix     {input.pg_matrix}"
        " --annotation {input.annotation}"
        " --out        {output.markers}"

rule scp_aucell:
    input:
        pivot      = "scp/pivot_pack.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
    output:
        scores = "scp/aucell/tables/aucell_scores.tsv",
    params:
        collections = _c("scp", "aucell", "genesets",    default="msigdb_hallmark"),
        top_n       = _c("scp", "aucell", "n_top_genes", default=50),
    shell:
        "Rscript {RSC}/run_auc_from_pivot.R"
        " --pivot       {input.pivot}"
        " --annotation  {input.annotation}"
        " --outdir      scp/aucell"
        " --collections    {params.collections}"
        " --top-n-heatmap  {params.top_n}"

rule scplainer_fit:
    input:
        pg_matrix  = "scp/pivot_pack.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
    output:
        rds                = "scp/scplainer/sce_scplainer_fit.rds",
        variance_explained = "scp/scplainer/scplainer_variance_explained.tsv",
    shell:
        "Rscript {RSC}/run_scplainer_from_pg_matrix.R"
        " --pg-matrix  {input.pg_matrix}"
        " --annotation {input.annotation}"
        " --outdir     scp/scplainer"

rule scplainer_contrasts:
    input:
        rds = "scp/scplainer/sce_scplainer_fit.rds",
    output:
        da = "scp/scplainer/scplainer_cluster_DA.tsv",
    shell:
        "Rscript {RSC}/scplainer_contrasts_from_rds.R"
        " --rds {input.rds}"
        " --out {output.da}"

rule scplainer_batch_corrected:
    input:
        rds = "scp/scplainer/sce_scplainer_fit.rds",
    output:
        corrected = "scp/scplainer/batch_corrected_expression.tsv",
    shell:
        "Rscript {RSC}/scplainer_batch_corrected.R"
        " --rds {input.rds}"
        " --out {output.corrected}"

rule scplainer_volcano:
    input:
        da = "scp/scplainer/scplainer_cluster_DA.tsv",
    output:
        done = sentinel("scp/scplainer/plots"),
    params:
        fmt = _fmt(),
    shell:
        "python {PY}/plot_scplainer_volcano.py"
        " --input  {input.da}"
        " --outdir scp/scplainer/plots"
        " --format {params.fmt}"
        " && echo done > {output.done}"

rule scplainer_adapter:
    input:
        da = "scp/scplainer/scplainer_cluster_DA.tsv",
    output:
        markers     = "scp/markers/scplainer_intensity_markers.tsv",
        sig_markers = "scp/markers/scplainer_intensity_markers_significant.tsv",
        top_markers = "scp/markers/scplainer_intensity_markers_topN.tsv",
    shell:
        "Rscript {RSC}/adapter_scplainer_to_markers.R"
        " --scplainer-da {input.da}"
        " --outdir       scp/markers"

rule scp_consensus_markers:
    input:
        detection  = "scp/markers/detection_markers.tsv",
        intensity  = "scp/markers/intensity_markers_detected_only.tsv",
        scplainer  = "scp/markers/scplainer_intensity_markers.tsv",
    output:
        consensus = "scp/markers/consensus_markers.tsv",
    shell:
        "python {PY}/make_consensus_markers.py"
        " --detection  {input.detection}"
        " --intensity  {input.intensity}"
        " --scplainer  {input.scplainer}"
        " --out        {output.consensus}"

rule scp_cluster_summary:
    input:
        annotation  = "scp/clustering/scp_annotation.tsv",
        assignments = "scp/clustering/scp_cluster_assignments.tsv",
        det_markers = "scp/markers/detection_markers.tsv",
        int_markers = "scp/markers/intensity_markers_detected_only.tsv",
        scp_markers = "scp/markers/scplainer_intensity_markers.tsv",
        aucell      = "scp/aucell/tables/aucell_scores.tsv",
        consensus   = "scp/markers/consensus_markers.tsv",
    output:
        summary = "scp/llm/cluster_summary.tsv",
        prompt  = "scp/llm/cluster_llm_prompt.md",
    params:
        custom_proteins = _c("scp", "custom_proteins", default=""),
    run:
        cp_arg = f' --custom-proteins "{params.custom_proteins}"' if params.custom_proteins else ""
        shell(
            f"python {PY}/make_cluster_summary.py"
            f" --annotation        {{input.annotation}}"
            f" --assignments       {{input.assignments}}"
            f" --detection-markers {{input.det_markers}}"
            f" --intensity-markers {{input.int_markers}}"
            f" --scplainer-markers {{input.scp_markers}}"
            f" --aucell-scores     {{input.aucell}}"
            f" --consensus-markers {{input.consensus}}"
            f" --out               {{output.summary}}"
            f" --llm-export        {{output.prompt}}"
            f" --top-n-markers     15"
            f" --config            config/caspa.json"
            f"{cp_arg}"
        )

rule scp_llm_annotation:
    input:
        prompt     = "scp/llm/cluster_llm_prompt.md",
        pivot      = "scp/pivot_pack.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
    output:
        annotations = "scp/llm/cluster_cell_type_annotations.tsv",
    params:
        species     = _c("project", "species_label", default="human"),
        context     = _c("project", "description", default=""),
        model       = _c("scp", "llm", "model", default="gpt-4o"),
        base_url    = _c("scp", "llm", "base_url", default=""),
        api_key     = _c("scp", "llm", "api_key", default=""),
        condition_b      = _c("scp", "llm", "condition_b",     default=False),
        provider         = _c("scp", "llm", "provider",         default="openai"),
        thinking_budget  = _c("scp", "llm", "thinking_budget",  default=0),
    run:
        if params.api_key:
            os.environ["ELM_API_KEY"] = params.api_key
        cond_b_flag    = " --condition-b" if params.condition_b else ""
        thinking_flag  = f" --thinking-budget {params.thinking_budget}" if params.thinking_budget else ""
        shell(
            f"python {PY}/annotate_clusters_llm.py"
            f" --prompt-md    {{input.prompt}}"
            f" --pivot        {{input.pivot}}"
            f" --annotation   {{input.annotation}}"
            f" --out-tsv      {{output.annotations}}"
            f" --species      {{params.species}}"
            f" --experiment-context \"{{params.context}}\""
            f" --model        {{params.model}}"
            f" --base-url     \"{{params.base_url}}\""
            f" --provider     {{params.provider}}"
            f" --out-dir      scp/llm"
            f"{cond_b_flag}"
            f"{thinking_flag}"
        )

rule scp_llm_umap:
    input:
        assignments = "scp/clustering/scp_cluster_assignments.tsv",
        pivot       = "scp/pivot_pack.tsv",
        shifts      = "scp/pivot_shifts.tsv",
        corrected   = "scp/scplainer/batch_corrected_expression.tsv",
        annotation  = "scp/clustering/scp_annotation.tsv",
        aucell      = "scp/aucell/tables/aucell_scores.tsv",
        det_markers = "scp/markers/detection_markers.tsv",
        cell_types  = "scp/llm/cluster_cell_type_annotations.tsv",
    output:
        umap = "scp/llm/plots/umap_cell_types.pdf",
    params:
        proteins        = "scp/llm/scp_recommended_markers_flat.txt",
        custom_proteins = _c("scp", "custom_proteins", default=""),
    run:
        proteins_arg = ""
        if os.path.exists(params.proteins):
            proteins_arg = f" --proteins {params.proteins}"
        cp_arg = f' --custom-proteins "{params.custom_proteins}"' if params.custom_proteins else ""
        shell(
            f"python {PY}/plot_umap_overlays.py"
            f" --assignments            {{input.assignments}}"
            f" --pivot                  {{input.pivot}}"
            f" --shifts                 {{input.shifts}}"
            f" --corrected-pivot        {{input.corrected}}"
            f" --annotation             {{input.annotation}}"
            f" --aucell                 {{input.aucell}}"
            f" --detection-markers      {{input.det_markers}}"
            f" --cell-type-annotations  {{input.cell_types}}"
            f"{proteins_arg}"
            f"{cp_arg}"
            f" --outdir                 scp/llm"
        )

rule scp_harmony_diagnostics:
    input:
        assignments = "scp/clustering/scp_cluster_assignments.tsv",
        report      = "scp/clustering/scp_clustering_report.json",
    output:
        plot = "scp/viz/plots/harmony_before_after.pdf",
    params:
        fmt = _fmt(),
    shell:
        "python {PY}/plot_harmony_diagnostics.py"
        " --assignments {input.assignments}"
        " --report      {input.report}"
        " --outdir      scp/viz"
        " --format      {params.fmt}"

rule scp_normalisation_qc:
    input:
        pivot      = "scp/pivot_pack.tsv",
        shifts     = "scp/pivot_shifts.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
    output:
        done = sentinel("scp/viz/plots/normalisation"),
    params:
        fmt = _fmt(),
    shell:
        "python {PY}/plot_scp_qc_normalisation.py"
        " --pivot      {input.pivot}"
        " --shifts     {input.shifts}"
        " --annotation {input.annotation}"
        " --outdir     scp/viz"
        " --format     {params.fmt}"
        " && echo done > {output.done}"

rule scp_marker_volcanos:
    input:
        det_markers = "scp/markers/detection_markers.tsv",
        int_markers = "scp/markers/intensity_markers_detected_only.tsv",
        scp_markers = "scp/markers/scplainer_intensity_markers.tsv",
    output:
        done = sentinel("scp/viz/plots/marker_volcanos"),
    params:
        fmt = _fmt(),
    shell:
        "python {PY}/plot_scp_marker_volcanos.py"
        " --detection-markers {input.det_markers}"
        " --intensity-markers {input.int_markers}"
        " --scplainer-markers {input.scp_markers}"
        " --outdir            scp/viz"
        " --format            {params.fmt}"
        " && echo done > {output.done}"

rule scp_centroid_similarity:
    input:
        pivot      = "scp/pivot_pack.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
    output:
        corr = "scp/viz/tables/cluster_centroid_correlation.tsv",
    params:
        fmt = _fmt(),
    shell:
        "python {PY}/plot_cluster_centroid_similarity.py"
        " --pivot      {input.pivot}"
        " --annotation {input.annotation}"
        " --outdir     scp/viz"
        " --format     {params.fmt}"

rule scp_detection_matrix:
    input:
        pivot      = "scp/pivot_pack.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
        det_markers = "scp/markers/detection_markers.tsv",
    output:
        done = sentinel("scp/viz/plots/detection_matrix"),
    params:
        fmt = _fmt(),
    shell:
        "python {PY}/plot_detection_matrix.py"
        " --pivot             {input.pivot}"
        " --annotation        {input.annotation}"
        " --detection-markers {input.det_markers}"
        " --outdir            scp/viz"
        " --format            {params.fmt}"
        " && echo done > {output.done}"

rule scp_marker_dotplot:
    input:
        det_markers = "scp/markers/detection_markers.tsv",
        int_markers = "scp/markers/intensity_markers_detected_only.tsv",
    output:
        dotplot = "scp/markers/marker_dotplot_matrix.pdf",
    shell:
        "Rscript {RSC}/plot_scp_marker_dotplot_matrix.R"
        " --detection {input.det_markers}"
        " --intensity {input.int_markers}"
        " --out       {output.dotplot}"
        " --top-n     30"

rule scp_enrichment:
    input:
        pivot      = "scp/pivot_pack.tsv",
        annotation = "scp/clustering/scp_annotation.tsv",
        det_markers = "scp/markers/detection_markers.tsv",
        int_markers = "scp/markers/intensity_markers_detected_only.tsv",
        aucell      = "scp/aucell/tables/aucell_scores.tsv",
    output:
        done = sentinel("scp/enrichment"),
    shell:
        "Rscript {RSC}/run_scp_enrichment_from_markers.R"
        " --pivot      {input.pivot}"
        " --annotation {input.annotation}"
        " --detection  {input.det_markers}"
        " --intensity  {input.int_markers}"
        " --aucell     {input.aucell}"
        " --outdir     scp/enrichment"
        " --species    auto"
        " && echo done > {output.done}"
