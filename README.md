# CASPA — Context-Aware Single Cell Proteomic Analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24058602.svg)](https://doi.org/10.5281/zenodo.24058602)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Standalone Snakemake pipeline for single-cell proteomics analysis. Takes a
pre-exported protein-group matrix (DIA-NN, FragPipe) **or** a raw Spectronaut
long-format TSV. No bulk pipeline dependency.

**New to CASPA?** [Follow the tutorial](docs/TUTORIAL.md) — a complete first run
on a bundled example dataset, so you can confirm your install works (and see
what a normal run looks like) before pointing CASPA at your own data.

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
**[`caspa-setup.exe`](https://github.com/vonkriegsheim/CASPA/releases/latest/download/caspa-setup.exe)**
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

**Spectronaut users:** export your report using the bundled
[`directLFQfragmentOutput.rs`](directLFQfragmentOutput.rs) report scheme
(**Report → Import Scheme** in Spectronaut) to guarantee the export has the exact
columns CASPA's converter (`pipeline/scripts/python/spectronaut_to_diann_pg_matrix.py`)
requires (`PG.ProteinGroups`, `PEP.StrippedSequence`, `FG.Quantity`, etc. — the
full list is enforced at conversion time, with a clear error naming any missing
column). Building a custom export scheme by hand is error-prone; importing this
one avoids that entirely.

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
| `joint_embedding.detection_embedding` | `pca` | How the 0/1 detection block is embedded: `pca` (Gaussian PCA on centred calls, as published), `glmpca` (Bernoulli GLM-PCA with a per-protein intercept) or `glmpca_cellint` (GLM-PCA with per-protein and per-cell intercepts). See below. |
| `custom_proteins` | `""` | Comma-separated gene names for custom UMAP overlays |
| `llm.model` | `gpt-5.2` | Model for cell type annotation (any provider, see the provider table) |
| `llm.condition_b` | `true` | Three-round annotation architecture (Round 0 context constraints, then two annotation passes). This is the published default; `false` gives the plain two-pass prompt. |
| `llm.skip_without_key` | `false` | A missing API key is an **error**. Set `true` to deliberately skip annotation and write an empty table (e.g. for a QC-only run). |
| `llm.allow_incomplete` | `false` | An annotation table with missing or empty cluster calls is an **error**. Set `true` to write the partial table anyway. |
| `llm.marker_modalities` | `["detection","intensity","scplainer","aucell"]` | Which marker-mining modalities to show the LLM annotator (see below) |

### Plot parameters (`plots.*`)

| key | default | description |
|-----|---------|-------------|
| `formats` | `["pdf","png"]` | Output formats |
| `adj_pval` | 0.05 | Adjusted p-value cutoff for volcano labels |

---

### Annotation defaults and failure behaviour

The three-round architecture described in the paper is the default
(`llm.condition_b: true`): Round 0 asks the model to derive dataset-specific
constraints from the experiment description alone, Round 1 annotates the
clusters with the marker evidence, and Round 2 supplies the detection rates of
the supplementary markers the model requested in Round 1 and asks for a revised
table. Every rendered prompt, the Round 0
constraints and the raw model output are written to `scp/llm/` together with a
JSON sidecar recording model, prompt hash and settings.

Workdirs scaffolded before this version carry an explicit
`"condition_b": false` in `config/caspa.json` (the old template default), which
overrides the bundled default: set it to `true` or delete the key to get the
three-round architecture in an existing workdir.

The annotation step **fails loudly** rather than passing silently:

- no API key: the rule exits with an error naming the variable/config key to set
  (`llm.skip_without_key: true` opts out and writes an empty table);
- an empty model response: error;
- a table in which any cluster announced in the prompt has no row, or an empty
  cell-type call: error naming the clusters (`llm.allow_incomplete: true` opts out).
  The raw pass 1/pass 2 outputs are kept in `scp/llm/` for inspection either way.

### `joint_embedding.detection_embedding`: Gaussian PCA or Bernoulli GLM-PCA

The published default embeds the binary detection calls with a Gaussian PCA on
centred 0/1 values, which treats presence/absence as a continuous variable. The
error family the data call for is Bernoulli; `detection_embedding: "glmpca"`
fits a logistic (Bernoulli GLM-)PCA with a per-protein intercept, and
`"glmpca_cellint"` adds a per-cell intercept (a Rasch-type depth term, the
logistic analogue of a size factor). The fit is deterministic
(majorisation-minimisation; de Leeuw 2006, Landgraf & Lee 2020) given the
iteration cap and tolerance (`--glmpca_max_iter`, `--glmpca_tol`); the report
`scp/clustering/scp_clustering_report.json` records the deviance explained and
whether the tolerance was reached. Deviance explained is measured against the
per-protein-intercept null in both modes, so for `glmpca_cellint` it includes
the contribution of the cell intercept. In the
embedding ablation reported with the paper, the Bernoulli block did not improve
coarse lineage recovery against FACS ground truth, so Gaussian PCA remains the
default; the option is there so the choice can be tested on your own data.

### `llm.marker_modalities`: a diagnostic tool, not a "best modality" setting

CASPA's marker-mining step produces three independent evidence layers per
cluster — detection rate, intensity fold-change, and scplainer's
batch-corrected contrasts — and by default shows the LLM annotator all
three. `scp.llm.marker_modalities` lets you restrict which of
`"detection"`, `"intensity"`, `"scplainer"` are included in the prompt.

This exists as a **sensitivity check**, not a recommended alternative
default. A leave-one-out ablation across four datasets found that
restricting to a single modality can flip the call on an individual
hard/ambiguous cluster, but the *direction* of the effect is
dataset-specific — detection-only data resolved one dataset's hardest
cluster, while for a different dataset the same restriction broke a call
that the full evidence set got right. There is no modality subset that is
universally better; see `paperReview/response/figures/regenerate_figS4_leaveoneout.py`
for the full analysis. If an annotation looks wrong on a specific
cluster and you suspect it's being swayed by conflicting evidence across
modalities, re-running `scp_cluster_summary` with a restricted
`marker_modalities` list is a way to check that hypothesis on that one
cluster — not something to leave configured for a full production run.

We recommend running this **modality sensitivity check** as part of every
annotation rather than as an optional extra: see
[docs/TUTORIAL.md, section 6](docs/TUTORIAL.md#6-modality-sensitivity-check-which-calls-can-you-trust)
for the step-by-step procedure and how to read the result.

---

## PanglaoDB cross-validation

Cell type annotations are cross-validated against
`pipeline/assets/panglaodb_markers.tsv.gz` (8,286 markers, 178 cell types,
human + mouse) before the LLM call, providing an independent evidence layer.

---

## Tests

`tests/` holds unit and end-to-end tests (plain `unittest`, also pytest-compatible):
the Bernoulli GLM-PCA fit, the embedding script in all three detection-embedding
modes (skipped unless `CASPA_TEST_PIVOT` and `CASPA_TEST_MANIFEST` point at a
`pivot_pack.tsv` and `filtered.manifest.tsv` from a pipeline run), and the LLM
annotation step driven through a mock OpenAI-compatible server covering every
failure path (missing key, empty response, truncated or incomplete tables).

```bash
python -m unittest discover -s tests -v
```

---

## Citation and licence

If you use CASPA, please cite the paper and the archived software version:

- Salomó Coll C, Makar AN, Brenes AJ, Inns J, Trost M, Rajan N, Wilkinson S,
  von Kriegsheim A. A Context-Aware Single-Cell Proteomics Analysis pipeline.
  *Nature Communications* (2026), in press.
- CASPA software v0.2.0, the version used in the paper, archived in Zenodo:
  [doi:10.5281/zenodo.24058604](https://doi.org/10.5281/zenodo.24058604). The concept DOI
  [doi:10.5281/zenodo.24058602](https://doi.org/10.5281/zenodo.24058602) always resolves to the
  latest archived release; `CITATION.cff` carries the same metadata.

CASPA is released under the [MIT licence](LICENSE). It builds on the following
tools, each used under its own open-source licence (notices are retained in
their respective packages):
- **scplainer**: Vanderaa & Gatto, Genome Biol 26, 237 (2025)
- **AUCell** (via decoupler): Aibar et al., Nat Methods 14, 1083 (2017)
- **Harmony** (harmonypy): Korsunsky et al., Nat Methods 16, 1289 (2019)
- **Leiden** (leidenalg): Traag et al., Sci Rep 9, 5233 (2019)
- **scanpy**: Wolf et al., Genome Biol 19, 15 (2018)
- **Snakemake**: Mölder et al., F1000Res 10, 33 (2021)
- **PanglaoDB**: Franzén et al., Database 2019, baz046 (2019)
