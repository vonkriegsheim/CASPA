# CASPA — Context-Aware Single Cell Proteomic Analysis

Standalone Snakemake pipeline for single-cell proteomics analysis. Takes a
pre-exported protein-group matrix (DIA-NN, FragPipe) **or** a raw Spectronaut
long-format TSV. No bulk pipeline dependency.

---

## Quick start

### 1. Install

```bash
conda env create -f environment.yml
conda activate caspa
```

### 2. Scaffold a new experiment

```bash
# From a DIA-NN / FragPipe pg_matrix
python caspa/init.py \
    --workdir /path/to/MyExperiment \
    --pg-matrix /path/to/report.pg_matrix.tsv \
    --species human \
    --name "My Experiment"

# From a Spectronaut long-format export
python caspa/init.py \
    --workdir /path/to/MyExperiment \
    --spectronaut-tsv /path/to/spectronaut_export.tsv \
    --species human \
    --name "My Experiment"
```

This creates:
- `config/caspa.json` — pre-filled with your input path and species
- `config/ms_inputs.tsv` — sample sheet (auto-populated from pg_matrix column headers)

### 3. Edit config

```
config/caspa.json         — set scp.llm.api_key and review parameters
config/ms_inputs.tsv      — confirm sample_id, sample_file, batch columns
```

### 4. Run

```bash
python caspa/run.py --workdir /path/to/MyExperiment --cores 30
```

Dry-run to preview rules:

```bash
python caspa/run.py --workdir /path/to/MyExperiment --dry-run
```

Run to a specific target:

```bash
python caspa/run.py --workdir /path/to/MyExperiment --target scp_llm_annotation
```

---

## Input formats

### `config/ms_inputs.tsv`

| column | description |
|--------|-------------|
| `sample_id` | Must exactly match pg_matrix column name |
| `sample_file` | Raw file path (used for reference; can match sample_id) |
| `batch` | Batch label for Harmony batch correction (integer or string) |

### `config/caspa.json` — key fields

```json
{
  "project": {
    "name": "MyExperiment",
    "species_label": "human",
    "description": "One-sentence context for LLM cell type annotation"
  },
  "input": {
    "pg_matrix": "/path/to/report.pg_matrix.tsv",
    "spectronaut_tsv": null
  },
  "scp": {
    "custom_proteins": "INS,GCG,SST",
    "llm": {
      "api_key": "sk-...",
      "model": "gpt-4o"
    }
  }
}
```

---

## Output directory map

```
scp/
├── qc/
│   ├── filtered.pg_matrix.tsv
│   ├── filtered.manifest.tsv
│   └── scp_qc_report.json
├── pivot_pack.tsv                         ← canonical cell × protein matrix
├── pivot_shifts.tsv
├── clustering/
│   ├── scp_annotation.tsv
│   ├── scp_cluster_assignments.tsv        ← UMAP coordinates + cluster labels
│   └── scp_clustering_report.json
├── markers/
│   ├── detection_markers.tsv
│   ├── intensity_markers_detected_only.tsv
│   ├── scplainer_intensity_markers.tsv
│   ├── scplainer_intensity_markers_significant.tsv
│   ├── scplainer_intensity_markers_topN.tsv
│   ├── consensus_markers.tsv
│   └── marker_dotplot_matrix.pdf
├── scplainer/
│   ├── sce_scplainer_fit.rds
│   ├── scplainer_variance_explained.tsv
│   ├── scplainer_cluster_DA.tsv
│   ├── batch_corrected_expression.tsv
│   └── plots/
├── aucell/
│   └── tables/aucell_scores.tsv
├── llm/
│   ├── cluster_summary.tsv
│   ├── cluster_llm_prompt.md
│   ├── cluster_cell_type_annotations.tsv  ← LLM cell type calls
│   ├── scp_recommended_markers_flat.txt
│   └── plots/umap_cell_types.pdf
├── viz/
│   ├── plots/
│   │   ├── harmony_before_after.pdf
│   │   ├── normalisation/
│   │   ├── marker_volcanos/
│   │   └── detection_matrix/
│   └── tables/cluster_centroid_correlation.tsv
└── enrichment/
```

---

## Configuration reference

All defaults live in `caspa_defaults.json` (bundled with the repo).
Override any key in `config/caspa.json`.

### SCP parameters (`scp.*`)

| key | default | description |
|-----|---------|-------------|
| `min_protein_ids` | 500 | Min proteins per cell to pass QC |
| `exclude_run_regex` | `^(library\|lib\|...)` | Regex to exclude library/dilution runs |
| `joint_embedding.leiden_resolution` | 0.8 | Leiden clustering resolution |
| `joint_embedding.n_neighbors` | 15 | UMAP / kNN neighbours |
| `joint_embedding.harmony_batch_key` | `Batch` | Column in ms_inputs.tsv for batch correction |
| `custom_proteins` | `""` | Comma-separated gene names for custom UMAP overlays |
| `llm.model` | `gpt-4o` | OpenAI model for cell type annotation |

### Plot parameters (`plots.*`)

| key | default | description |
|-----|---------|-------------|
| `formats` | `["pdf","png"]` | Output formats |
| `adj_pval` | 0.05 | Adjusted p-value cutoff for volcano labels |

---

## PanglaoDB cross-validation

Cell type annotations are cross-validated against
`pipeline/assets/panglaodb_markers.tsv.gz` (8,286 markers, 178 cell types,
human + mouse) before the LLM call, providing an independent evidence layer.

---

## Citation

If you use CASPA, please cite the upstream tools:
- **scplainer**: Vanderaa & Gatto, 2023
- **AUCell**: Aibar et al., 2017
- **Harmony**: Korsunsky et al., 2019
- **Leiden**: Traag et al., 2019
