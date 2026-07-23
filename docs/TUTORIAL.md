# CASPA tutorial: your first run

This walks through a complete CASPA run on a real, small example dataset bundled
with the repository, so you can confirm your install works end to end — and know
what a normal, successful run looks like — before pointing CASPA at your own data.

Expect the whole thing (install verification through inspecting the output) to
take about 15–20 minutes, most of which is the pipeline itself running.

## 0. Verify your install

```bash
python caspa/doctor.py
```

This checks Python + R packages and PATH resolution and prints exactly what's
missing if anything is wrong. Fix any reported problems before continuing — a
broken install produces confusing failures much later in the pipeline instead of
a clear error up front.

## 1. The example dataset

We'll use `testset/SkinCancer-human-diann/`, a real label-free DIA-NN dataset
already bundled in the repository (skin tumour tissue, DIA-NN protein-group
matrix). It's the smallest of the bundled example datasets, so it's the fastest
to run and cheapest to annotate — a good first check that everything works.

```
testset/SkinCancer-human-diann/
├── data/
│   ├── pg_matrix.tsv       ← DIA-NN protein-group matrix
│   └── ms_inputs.tsv       ← sample sheet (sample_id, sample_file, batch — already filled in)
```

**A note on reproducibility:** this tutorial uses CASPA's *default* QC threshold
(`min_protein_ids = 400`). The dataset published in the accompanying paper used a
stricter, dataset-specific threshold (`min_protein_ids = 500`) tuned for that
analysis, which resolves a cleaner 7-cluster solution. Running this tutorial with
defaults will most likely give you **8 clusters**, not 7 — that's expected, not a
bug. The point of this tutorial is to confirm the pipeline runs correctly on your
machine, not to reproduce the paper's exact published numbers. (If you do want to
match the published analysis, see [Reproducing the paper's SkinCancer
figures](#reproducing-the-papers-skincancer-figures) at the end.)

## 2. Scaffold a new experiment

```bash
python caspa/init.py \
    --workdir /path/to/CASPA_tutorial_run \
    --pg-matrix testset/SkinCancer-human-diann/data/pg_matrix.tsv \
    --species human \
    --name "CASPA Tutorial"
```

This creates `config/caspa.json` and `config/ms_inputs.tsv` under your workdir,
pre-filled from the pg_matrix's column headers.

Since this example dataset already ships a hand-verified `ms_inputs.tsv` with the
correct `batch` values, copy it over the auto-generated one rather than
re-deriving batches from scratch:

```bash
cp testset/SkinCancer-human-diann/data/ms_inputs.tsv /path/to/CASPA_tutorial_run/config/ms_inputs.tsv
```

(For your own data, you'd instead run `caspa/assign_batches.py` — Evosep only —
or hand-fill the `batch` column; see the main [README](../README.md#3-assign-batches-evosep-only).)

## 3. Edit the config

Open `config/caspa.json` in the workdir. Two things matter most for a first run:

**The experiment description.** This is the single highest-leverage field in the
whole config — it's the context the LLM annotation step uses to reason about
what cell types are plausible, what markers to expect, and what counts as
contamination versus genuine biology. A one-line placeholder produces
noticeably weaker, less specific annotations than a real description. For this
dataset, a description close to what was actually used for the published
analysis is:

```json
{
  "project": {
    "name": "CASPA Tutorial",
    "species_label": "human",
    "description": "Human skin tumour (CYLD cutaneous syndrome / cylindroma) resection, FACS-sorted into CD45-/CD200+ (tumour keratinocytes) and CD45+/CD74+ (immune/myeloid) populations before single-cell capture. Expected cell types: tumour keratinocytes of hair-follicle origin, macrophages, dendritic cells, stromal cells."
  },
  ...
}
```

If you're annotating your own data, write the equivalent for your experiment:
species, tissue, any sort/gating strategy, sample handling, and what cell types
you'd expect to see. Treat it like the "materials and methods" a colleague would
need to make sense of your data — the more specific, the better the annotation.

**The LLM API key.** Set `scp.llm.api_key` to a real API key from a commercial
provider (OpenAI, Anthropic, Gemini, DeepSeek, or any OpenAI-compatible
endpoint — see the [README's model/provider table](../README.md#configcaspajson--key-fields)).
**Use a real API key here, not a proxy, a shared/rate-limited key, or anything
else that isn't a standard, directly-billed provider account** — the failure
modes in the [appendix below](#appendix-llm-annotation-failure-modes-and-how-to-read-them)
are much harder to diagnose behind an intermediary. If you don't have a paid key
handy, Gemini's free tier (`gemini-2.5-flash`, from aistudio.google.com) is
enough to complete this tutorial.

```json
"scp": {
  "llm": {
    "api_key": "sk-...",
    "model": "gpt-5.2"
  }
}
```

## 4. Run

```bash
python caspa/run.py --workdir /path/to/CASPA_tutorial_run --cores 8
```

You'll see Snakemake step through QC → normalisation → batch correction →
clustering → marker mining → LLM annotation → plotting. On this dataset,
expect the whole run to take a few minutes of compute plus roughly 1–3 minutes
for the LLM annotation step itself (it's a small prompt — one of the fastest of
any dataset we've tested).

To preview the steps without running anything:

```bash
python caspa/run.py --workdir /path/to/CASPA_tutorial_run --dry-run
```

If the run is interrupted partway through, resume with:

```bash
python caspa/run.py --workdir /path/to/CASPA_tutorial_run --cores 8 --rerun-incomplete
```

## 5. Check the output

If everything worked, you should have:

```
<workdir>/scp/
├── clustering/scp_cluster_assignments.tsv   ← every cell's cluster + UMAP coords
├── markers/consensus_markers.tsv            ← top markers per cluster
├── llm/cluster_cell_type_annotations.tsv    ← the cell-type calls — the main deliverable
└── llm/plots/umap_cell_types.pdf            ← UMAP coloured by annotated cell type
```

Open `scp/llm/cluster_cell_type_annotations.tsv`. You should see one row per
cluster (most likely 7 or 8 on this dataset with default settings) with columns
including `cell_type`, `key_evidence`, `confidence`, and `contradictions`. A
successful run looks roughly like this (illustrative — your exact wording will
differ run to run, since these are genuine LLM calls, not fixed lookup values):

| Cluster | cell_type | confidence |
|---|---|---|
| C0 | Tumour keratinocyte (differentiated) | High |
| C1 | Macrophage (tissue-resident) | High |
| C2 | Basal/outer-root-sheath keratinocyte | Medium-high |
| ... | ... | ... |

If instead the file is empty, has fewer rows than clusters, or every row says
something generic like "Contaminant" or "Unknown" with Low confidence across
the board, something went wrong upstream — see the appendix below before
re-running.

**You're now ready to point CASPA at your own dataset.** Repeat steps 2–5 with
your own `pg_matrix.tsv` (or `--spectronaut-tsv` for a Spectronaut export) and a
description tailored to your experiment.

---

## Appendix: LLM-annotation failure modes and how to read them

The LLM annotation step is the one part of the pipeline that calls out to a
third-party API, and it's where nearly all real-world run failures happen —
everything upstream (QC, clustering, marker mining) is deterministic local
computation. The following are genuine failure modes encountered and diagnosed
while validating this pipeline against dozens of models and provider
configurations; if your own run produces incomplete or empty annotations, check
this list before assuming the pipeline itself is broken.

**Truncated or empty output on large/complex datasets.** If `scp/llm/
cluster_cell_type_annotations.tsv` has noticeably fewer rows than clusters, or
the file is empty after a run that otherwise completed without an error, this
is very likely an output-token budget problem, not a crash. Symptom: the run
finishes with exit code 0 (Snakemake reports success) but the annotation table
is silently incomplete. This happened repeatedly with datasets that have many
clusters (20+) and/or with a model's reasoning effort turned up — a model can
spend most of its allotted output tokens on internal reasoning before it ever
writes the answer table, leaving too little budget for the visible output.
**Fix:** raise `scp.llm.max_tokens` in `config/caspa.json` (try doubling it) and
re-run. CASPA already auto-clamps `max_tokens` down to a model's own hard cap
when needed and warns on truncation (check the pipeline log for a `WARNING:
... response was truncated` line) — but it cannot invent a higher budget than
you request, so if you've enabled a high thinking/reasoning-effort setting for
your model, budget generously.

**A run that appears to hang.** Some model/provider combinations can take
several minutes per cluster-annotation call on a large prompt — this is normal
and not a hang. Before assuming something is stuck: check whether the process
is still consuming CPU/network, and whether any output files under
`scp/llm/prompt_pass*_*.txt` have appeared or grown recently (a sign the call is
progressing, not stalled). Genuine multi-hour stalls with zero progress
(confirmed during pipeline validation, on certain small local open-weight
models specifically) do happen and are a real provider/model-side issue rather
than something to fix in your config; if a call has shown literally zero
progress for an extended period, it's reasonable to cancel and retry rather
than wait indefinitely. `CASPA_LLM_TIMEOUT` (seconds, env var) controls how long
CASPA itself will wait before giving up and reporting a timeout error, in case
you'd rather fail fast than guess.

**A transient provider-side error** (e.g. an HTTP 5xx "service unavailable," or
occasionally a spurious "model not found" for a real, correctly-spelled model
ID). These are infrastructure blips on the provider's end, not a configuration
problem — simply re-running the same command usually resolves them without any
change.

**Streaming-related failures on very long Claude calls.** If you're using a
Claude model with an extended thinking budget (`scp.llm.thinking_budget` set to
a large value) and see an error mentioning a 10-minute duration limit, this is
a client-library requirement (the Anthropic SDK requires streaming for calls it
estimates may run long) — CASPA's Claude call path already streams by default,
so this should not occur in current versions; if you do hit it, it indicates an
outdated CASPA install and is worth reporting.

**Everything scored "Low confidence" or generically labelled.** This is
usually not a bug at all — it means the pipeline is working as intended and
telling you the evidence genuinely doesn't support a confident call for that
cluster (small cluster size, few detected proteins, conflicting markers). Check
`scp/llm/cluster_summary.tsv` for that cluster's actual marker evidence before
assuming something is wrong; a real biological ambiguity should look low-
confidence, and forcing false confidence would be the actual bug.

**When in doubt:** `scp/llm/prompt_pass1_user.txt` and `prompt_pass2_user.txt`
contain the exact rendered prompt sent to the model for that run, and
`scp/llm/cluster_annotation_raw_pass1.txt` / `..._pass2.txt` contain the raw,
unparsed model response. Reading these directly is the fastest way to tell
whether a problem is upstream (bad/missing marker evidence reaching the prompt)
or downstream (a well-evidenced prompt that the model failed to answer
completely) — the fix is different in each case.

---

## Reproducing the paper's SkinCancer figures

To match the exact clustering used in the published analysis rather than the
tutorial defaults above, override the QC threshold in `config/caspa.json`:

```json
"scp": {
  "min_protein_ids": 500
}
```

This resolves the same 7-cluster solution reported in the paper's Figure 5,
rather than the 8-cluster default-settings result you get in the steps above.
