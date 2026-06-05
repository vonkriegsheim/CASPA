# CASPA — Context-Aware Single Cell Proteomic Analysis

Standalone Snakemake pipeline for single-cell proteomics analysis. Takes a
pre-exported protein-group matrix (DIA-NN, FragPipe) **or** a raw Spectronaut
long-format TSV. No bulk pipeline dependency.

---

## Quick start

### 1. Install

CASPA needs a Python stack **and** an R/Bioconductor stack.

**Linux / macOS** — one conda environment has everything:

```bash
bash install_unix.sh          # conda env create -f environment.yml + verify
conda activate caspa
```

**Windows (easiest)** — download
**[`caspa-setup.exe`](https://github.com/vonkriegsheim/CASPA/releases/download/v0.1.0/caspa-setup.exe)**
([all releases](https://github.com/vonkriegsheim/CASPA/releases)) and double-click.
It's a no-admin installer that sets up a private, version-pinned Python + R inside
its own folder (installs to `%LOCALAPPDATA%\CASPA`; ~20–40 min as it downloads
R + packages). When it finishes you get a Start-Menu **CASPA** folder with three
shortcuts (each is also a `.cmd` in the install folder):

| Shortcut | What to run it for |
|---|---|
| **CASPA Setup (GUI)** | Build an experiment by **form** — pick the `pg_matrix`, describe the experiment (condition/FACS for the LLM), auto-detect or edit batches, choose the LLM model + key, click **Generate**. It writes the config into a workdir and prints the exact run command. |
| **CASPA Console** | A command prompt with CASPA's Python + R already on PATH. **Run the analysis here:** `python caspa\run.py --workdir <workdir> --cores 8` |
| **CASPA Doctor** | Verifies the install (Python + R packages, and which R is used). Run this first if anything misbehaves. |

**End to end:** open **CASPA Setup (GUI)** → fill the form → **Generate** (note the
workdir it creates, e.g. `C:\Users\you\MyExperiment`) → open **CASPA Console** →
paste the run command it printed. Results land under `<workdir>\scp\` — cell-type
calls in `scp\llm\cluster_cell_type_annotations.tsv`, figures in `scp\llm\plots\`
and `scp\viz\plots\`. The annotation step needs an **LLM API key** (enter it in the
GUI): OpenAI/Anthropic/DeepSeek are paid; Gemini's free tier (`gemini-2.5-flash`,
key from aistudio.google.com) works.

**Windows (manual)** — bioconda has no win-64 Bioconductor builds, so install
Python with pip (miniforge) and R from CRAN. After installing
[miniforge](https://github.com/conda-forge/miniforge/releases) and
[R 4.6.0](https://cran.r-project.org/bin/windows/base/), from a
PowerShell prompt:

```powershell
powershell -ExecutionPolicy Bypass -File .\install_windows.ps1
```

That runs `pip install -r requirements-windows.txt`, installs the R packages via
`pipeline/scripts/R/install_r_packages.R`, and verifies the result. At run time
ensure `python`, `snakemake`, and `Rscript` are on PATH — the installer prints
the exact directories to add.

**Verify any install** (Python pkgs, R pkgs, and PATH resolution):

```bash
python caspa/doctor.py
```

**Docker** — no local install at all (good for Windows via Docker Desktop/WSL2,
and for reproducibility). The image bundles the full Python + R stack.

Pull the prebuilt image (nothing to build):

```bash
# run the pipeline
docker run --rm -v /host/MyExperiment:/work \
  ghcr.io/vonkriegsheim/caspa --workdir /work --cores 8

# …or launch the setup GUI at http://localhost:8501
docker run --rm -p 8501:8501 ghcr.io/vonkriegsheim/caspa gui

# …or check the install
docker run --rm ghcr.io/vonkriegsheim/caspa doctor
```

…or build it yourself:

```bash
docker build -t caspa .
docker run --rm -v /host/MyExperiment:/work caspa --workdir /work --cores 8
```

### Supported stack (pinned)

Versions are pinned for reproducibility — every install path lands on the same
tested stack rather than "whatever is newest". Bump deliberately and re-test.

| Component | Pinned to | Where |
|---|---|---|
| **R** (Windows native) | **4.6.0** (Bioconductor **3.23**) | the Windows installer, `install_windows.ps1`, `install_r_packages.R` — the combo tested **end-to-end** |
| **R** (conda / Docker) | **4.5** (Bioconductor **3.21**) | `environment.yml` (`r-base=4.5`) — its own known-good linux solve |
| **Python** | **3.11** | `environment.yml`; the Windows pip path runs on miniforge's 3.11–3.12 (both tested) |
| **Python packages** | exact versions | `requirements-windows.txt` (Windows) / `environment.yml` (conda/Docker) |
| **R / Bioconductor packages** | the Bioc release set | fixed once R is pinned |

If a package has no prebuilt Windows binary for your R/Bioc, `install_r_packages.R`
falls back to a source install (fine for pure-R/data packages such as `org.Hs.eg.db`,
`GO.db`). Both the installer and `caspa doctor` verify each R package by **loading**
it, not just by name — so a package that is installed but can't load (a missing
transitive dep like `GSEABase`/`GO.db`/`getopt`) is caught, not reported as green.

### GUI (optional)

Prefer a form to hand-editing JSON? Launch the setup GUI:

```bash
python caspa/gui.py          # opens in your browser (re-launches via Streamlit)
```

It builds `config/caspa.json` + `config/ms_inputs.tsv` for a workdir: point it at
your pg_matrix, describe the experiment (condition / FACS context for the LLM),
auto-detect or hand-edit batches in a grid, set parameters and the API key, then
click **Generate**. Needs `streamlit` (bundled in the conda env / Docker image;
`pip install streamlit` otherwise).

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

### 3. Assign batches (Evosep only)

If your data was acquired on an **Evosep One**, batch numbers can be detected
automatically from the file names. Skip this step for other instruments and fill
the `batch` column manually.

```bash
# Preview detected batches without writing
python caspa/assign_batches.py --workdir /path/to/MyExperiment --dry-run

# Write batch numbers into config/ms_inputs.tsv
python caspa/assign_batches.py --workdir /path/to/MyExperiment
```

See [Evosep batch detection](#evosep-batch-detection) for how batches are inferred
and how to tune the thresholds.

### 4. Edit config

```
config/caspa.json         — set scp.llm.api_key and review parameters
config/ms_inputs.tsv      — confirm sample_id, sample_file, batch columns
```

### 5. Run

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

> **Tip — restarting a partial run:**
> If the pipeline was interrupted, add `--rerun-incomplete` to avoid having to
> re-run already-completed steps.

```bash
python caspa/run.py --workdir /path/to/MyExperiment --cores 30 --rerun-incomplete
```

---

## Input formats

### `config/ms_inputs.tsv`

Auto-generated by `caspa init`. Three required columns:

| column | description |
|--------|-------------|
| `sample_id` | Bare filename stem (no path prefix, no extension). Must match the stem that `scp_qc_filter` derives from the pg_matrix column header. `caspa init` writes this correctly for all supported input types. |
| `sample_file` | Full pg_matrix column name (e.g. `input\filename.raw`). Written for reference; not used by the pipeline directly. |
| `batch` | Integer batch label for Harmony batch correction. All samples from the same acquisition session get the same number. Set to `1` for all samples if no batch correction is needed. |

> **DIA-NN / FragPipe note:** pg_matrix columns typically contain a path prefix
> and file extension (e.g. `input\run42.raw`). `caspa init` automatically strips
> these to produce the bare stem (`run42`) required in `sample_id`.

> **Batch note:** if all samples were acquired in a single uninterrupted session,
> set all batches to `1` and Harmony will be skipped automatically.

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

**LLM models, providers & token limits.** `scp.llm.model` accepts any OpenAI,
Anthropic (`claude-*`), or OpenAI-compatible model. For Gemini or DeepSeek, set
`scp.llm.provider` to `openai` and point `scp.llm.base_url` at their endpoint
(the GUI presets do this for you). `scp.llm.max_tokens` (default 16000) caps the
annotation output; you rarely need to touch it — CASPA **auto-reduces it to the
model's own limit** if that's lower (e.g. DeepSeek's 8192) and **warns** if a
response is truncated or a cluster is left unannotated. If you see such a warning
on a very large dataset, raise `max_tokens`.

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

## Evosep batch detection

`caspa/assign_batches.py` reads `config/ms_inputs.tsv` and auto-assigns batch
numbers from the Evosep file naming convention, then overwrites the `batch` column.

### Naming convention parsed

Evosep filenames follow the pattern:

```
{prefix}_{S<N>}-{well}_{something}_{run_number}
```

| field | example | meaning |
|-------|---------|---------|
| `S<N>` | `S4`, `s1` | Physical plate loaded on the Evosep |
| `{well}` | `A1`–`H12` | Position in the 96-well plate |
| `{run_number}` | `11154` | Absolute injection counter (global, monotonically increasing) |

Examples:
```
SC_01_s1-a1_1_7515
20250305_JI_AF_hexSCP_plate1_singlecell_pop1_3_S4-A9_1_11155
```

### How batches are detected

A **batch** is a single contiguous acquisition session — the set of injections
that ran without the plate being removed and reloaded. Two detectors run in order:

#### Primary — run-number gap (default threshold: 20)

For consecutive runs from the same plate (sorted by run number), a new batch
starts when the gap between run numbers exceeds `gap_threshold`.

- Within a session: cells are injected sequentially, gap = 1–6 (a few QC
  injections may sit between cell rows).
- Between sessions (different days or plate reloads): gap is typically > 50.

#### Secondary — well-position drop (default threshold: 80)

When consecutive same-plate runs have a small run-number gap (< `gap_threshold`)
but the 96-well position drops by more than `well_drop_threshold`, a new batch is
also triggered. This catches back-to-back full-plate reloads where the plate
completes (well H12) and is immediately reloaded (well A1 again), with only a 1–2
run-number gap between the two sessions.

| scenario | well drop | detected? |
|----------|-----------|-----------|
| H12 → A1 (full-plate reload) | 95 | ✓ new batch |
| A10 → A1 (test injections before main run) | 9 | ✗ same batch |
| A8 → A1 (startplate blank then cells) | 7 | ✗ same batch |

### Usage

```bash
# Preview — does not modify the manifest
python caspa/assign_batches.py --workdir /path/to/MyExperiment --dry-run

# Apply — overwrites config/ms_inputs.tsv batch column
python caspa/assign_batches.py --workdir /path/to/MyExperiment

# Override thresholds at the command line
python caspa/assign_batches.py --workdir /path/to/MyExperiment \
    --gap-threshold 50 \
    --well-drop-threshold 70
```

### Tuning thresholds

Both thresholds are stored in `caspa/evosep_batch_patterns.json` and can be
edited permanently there, or overridden per-run with `--gap-threshold` /
`--well-drop-threshold`.

```json
{
  "gap_threshold": 20,
  "well_drop_threshold": 80,
  "patterns": [...]
}
```

**When to increase `gap_threshold`:** if your lab runs many QC injections between
cell rows (more than 20 between any two cells of the same plate in the same
session), raise this value to avoid false batch splits.

**When to lower `well_drop_threshold`:** if plates are routinely loaded partially
(e.g. only to row E) and then reloaded immediately, the drop from E12 to A1 is 61.
Lower the threshold below 61 to catch these.

### What about non-Evosep data?

Leave the `batch` column in `ms_inputs.tsv` as-is (all `1` from `caspa init`) or
fill it manually. `assign_batches.py` only processes samples that match the
`evosep_standard` regex; unmatched samples are assigned batch `0` with a warning.

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
