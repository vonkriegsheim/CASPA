#!/usr/bin/env python3
"""caspa GUI — a guided front-end for building a CASPA experiment.

A thin Streamlit wrapper over the existing CLI scripts: it generates
config/caspa.json and config/ms_inputs.tsv (with batches) for a workdir, so a
non-bioinformatician can set up a run by filling in a form instead of editing
JSON. It reuses init.py (sample-sheet generation) and assign_batches.py (Evosep
batch detection) — no analysis logic is duplicated here.

Launch (either works):
    python caspa/gui.py            # re-launches itself under Streamlit
    streamlit run caspa/gui.py
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def _under_streamlit() -> bool:
    try:
        from streamlit.runtime.scriptrunner import get_script_run_ctx
        return get_script_run_ctx() is not None
    except Exception:
        return False


# When run as `python caspa/gui.py`, relaunch through Streamlit so users don't
# need to remember the `streamlit run` incantation.
if not _under_streamlit():
    import subprocess
    try:
        import streamlit  # noqa: F401
    except ImportError:
        sys.exit("Streamlit is not installed. Install it with:\n"
                 "    pip install streamlit\n"
                 "(or use the conda env / Docker image, which include it).")
    raise SystemExit(subprocess.call(
        [sys.executable, "-m", "streamlit", "run", os.path.abspath(__file__), *sys.argv[1:]]))


# ---------------------------------------------------------------------------
# From here down we are running inside Streamlit.
# ---------------------------------------------------------------------------
import json

import pandas as pd
import streamlit as st

sys.path.insert(0, HERE)                       # import sibling CLI modules
import init as caspa_init                       # noqa: E402
import assign_batches as caspa_batches          # noqa: E402

DEFAULTS_PATH = os.path.join(os.path.dirname(HERE), "caspa_defaults.json")


def load_defaults() -> dict:
    try:
        with open(DEFAULTS_PATH, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return {}


def build_sample_sheet_df(pg_path: str) -> pd.DataFrame:
    """(sample_id, sample_file, batch=1) from a pg_matrix header — reuses init.py."""
    pairs = caspa_init.infer_sample_ids_from_pg_matrix(pg_path)   # [(col, stem), ...]
    return pd.DataFrame(
        {"sample_id": [stem for _c, stem in pairs],
         "sample_file": [col for col, _s in pairs],
         "batch": [1] * len(pairs)}
    )


def detect_evosep_batches(sample_ids):
    """Evosep batch numbers for a list of sample_ids — reuses assign_batches.py.

    Returns (mapping sample_id->batch, n_batches, n_unmatched).
    """
    patterns, gap, well_drop = caspa_batches.load_config(caspa_batches.DEFAULT_PATTERNS)
    tagged, row_objs = [], []
    for sid in sample_ids:
        ro = [sid]                              # unique object so id() is stable
        _pat, plate, well, run = caspa_batches.parse_sample(sid, patterns)
        if plate is not None:
            tagged.append((ro, _pat, plate, well, run))
        row_objs.append((sid, ro))
    if not tagged:
        return {sid: 0 for sid in sample_ids}, 0, len(sample_ids)
    rid_to_batch = caspa_batches.compute_batches(tagged, gap, well_drop)
    out, unmatched = {}, 0
    for sid, ro in row_objs:
        b = rid_to_batch.get(id(ro))
        if b is None:
            out[sid] = 0
            unmatched += 1
        else:
            out[sid] = int(b)
    n_batches = len({v for v in out.values() if v > 0})
    return out, n_batches, unmatched


def write_experiment(workdir, project, input_block, scp_block, sheet_df):
    """Write config/caspa.json + config/ms_inputs.tsv into the workdir."""
    cfg_dir = os.path.join(workdir, "config")
    os.makedirs(cfg_dir, exist_ok=True)
    config = {"project": project, "input": input_block, "scp": scp_block}
    with open(os.path.join(cfg_dir, "caspa.json"), "w", encoding="utf-8") as fh:
        json.dump(config, fh, indent=2)
    sheet_df.to_csv(os.path.join(cfg_dir, "ms_inputs.tsv"), sep="\t", index=False)
    return os.path.join(cfg_dir, "caspa.json"), os.path.join(cfg_dir, "ms_inputs.tsv")


# ===========================================================================
# UI
# ===========================================================================
st.set_page_config(page_title="CASPA experiment setup", page_icon="🧫", layout="wide")
D = load_defaults()
emb = D.get("scp", {}).get("joint_embedding", {})

st.title("🧫 CASPA — experiment setup")
st.caption("Build `config/caspa.json` + `config/ms_inputs.tsv` for a run. "
           "No coding required — fill in the form and click *Generate*.")

# ---- 1. Input data --------------------------------------------------------
st.header("1 · Input data")
input_kind = st.radio("Input type", ["DIA-NN / FragPipe pg_matrix", "Spectronaut long-format TSV"],
                      horizontal=True)
data_path = st.text_input("Path to the input file",
                          placeholder=r"D:\data\report.pg_matrix.tsv")
is_spectronaut = input_kind.startswith("Spectronaut")

if st.button("Load samples", type="primary", disabled=not data_path):
    if not os.path.isfile(data_path):
        st.error(f"File not found: {data_path}")
    elif is_spectronaut:
        # Spectronaut is converted to a pg_matrix by the pipeline at run time;
        # the sample sheet is filled in after conversion. Keep it minimal here.
        st.session_state.sheet = pd.DataFrame({"sample_id": [], "sample_file": [], "batch": []})
        st.session_state.is_spectronaut = True
        st.info("Spectronaut input recorded. The sample sheet is auto-populated "
                "after conversion at run time; batches can be set then.")
    else:
        try:
            st.session_state.sheet = build_sample_sheet_df(data_path)
            st.session_state.is_spectronaut = False
            st.success(f"Loaded {len(st.session_state.sheet)} samples from the pg_matrix header.")
        except Exception as e:  # noqa: BLE001
            st.error(f"Could not read pg_matrix header: {e}")

# ---- 2. Experiment context ------------------------------------------------
st.header("2 · Experiment")
c1, c2 = st.columns(2)
with c1:
    exp_name = st.text_input("Experiment name", value="My experiment")
    species = st.selectbox("Species", ["human", "mouse"])
with c2:
    custom_proteins = st.text_input("Custom proteins for UMAP overlays (comma-separated)",
                                    value=D.get("scp", {}).get("custom_proteins", ""))
description = st.text_area(
    "Experiment context (condition, FACS gating, tissue, developmental stage…)",
    placeholder="Single-cell proteomics of human ovarian fallopian-tube tissue: "
                "STIC lesions, tumour and stroma. CD45-/EpCAM+ sorted.",
    height=110,
    help="Free text fed to the LLM for context-aware annotation. The richer and more "
         "specific (cell types expected, sort gates, ambient signals), the better.")

# ---- 3. Sample sheet & batches -------------------------------------------
st.header("3 · Sample sheet & batches")
if st.session_state.get("sheet") is None:
    st.info("Load an input file above to populate the sample sheet.")
elif st.session_state.get("is_spectronaut"):
    st.info("Spectronaut sample sheet is generated at run time.")
else:
    b1, b2 = st.columns([1, 3])
    with b1:
        if st.button("⚙️ Auto-detect Evosep batches"):
            ids = list(st.session_state.sheet["sample_id"])
            bmap, n_b, n_un = detect_evosep_batches(ids)
            df = st.session_state.sheet.copy()
            df["batch"] = [bmap.get(s, 0) for s in df["sample_id"]]
            st.session_state.sheet = df
            st.session_state.batch_msg = (
                f"Detected {n_b} batch(es)."
                + (f" {n_un} sample(s) didn't match the Evosep pattern → batch 0 "
                   "(set these manually)." if n_un else ""))
            st.rerun()
    with b2:
        if st.session_state.get("batch_msg"):
            (st.warning if "batch 0" in st.session_state.batch_msg else st.success)(
                st.session_state.batch_msg)
        else:
            st.caption("Evosep timsTOF data → click auto-detect. Otherwise edit the "
                       "`batch` column directly (use `1` for everything if no batches).")

    edited = st.data_editor(
        st.session_state.sheet, num_rows="fixed", width="stretch", height=320,
        column_config={
            "sample_id": st.column_config.TextColumn("sample_id", disabled=True),
            "sample_file": st.column_config.TextColumn("sample_file", disabled=True),
            "batch": st.column_config.NumberColumn("batch", min_value=0, step=1),
        })
    st.session_state.sheet = edited
    nb = sorted(set(int(b) for b in edited["batch"]))
    st.caption(f"{len(edited)} samples · {len([b for b in nb if b > 0])} distinct batch(es): {nb}")

# ---- 4. Advanced settings -------------------------------------------------
st.header("4 · Settings")
with st.expander("Advanced clustering / QC parameters (defaults are sensible)"):
    a1, a2, a3 = st.columns(3)
    with a1:
        min_prot = st.number_input("Min proteins per cell (QC floor)",
                                   value=int(D.get("scp", {}).get("min_protein_ids", 400)), step=50)
        seed = st.number_input("Random seed", value=int(emb.get("seed", 0)), step=1)
    with a2:
        n_pcs = st.number_input("PCA components (intensity)", value=int(emb.get("n_pcs", 20)), step=1)
        n_neighbors = st.number_input("kNN / UMAP neighbours", value=int(emb.get("n_neighbors", 15)), step=1)
    with a3:
        leiden_res = st.number_input("Leiden resolution", value=float(emb.get("leiden_resolution", 0.8)),
                                     step=0.1, format="%.2f")
    st.caption("Adaptive Harmony (Shannon-entropy batch mixing) runs automatically; "
               "its targets are pipeline defaults.")

# ---- 5. LLM annotation ----------------------------------------------------
st.header("5 · LLM annotation")
llm_d = D.get("scp", {}).get("llm", {})
l1, l2 = st.columns(2)
with l1:
    provider = st.selectbox("Provider", ["openai", "claude"],
                            index=0 if llm_d.get("provider", "openai") == "openai" else 1)
    model = st.text_input("Model", value="gpt-5.2",
                          help="e.g. gpt-5.2 (recommended), gpt-4o, or a claude-* model.")
with l2:
    api_key = st.text_input("API key", type="password",
                            help="Stored only in this workdir's config/caspa.json (which is gitignored). "
                                 "Leave blank if using Claude Code to submit.")
    base_url = st.text_input("Custom base URL (optional)", value=llm_d.get("base_url", ""))
condition_b = st.checkbox("Use 3-round (Round-0 context) annotation",
                          value=bool(llm_d.get("condition_b", False)))

# ---- 6. Generate ----------------------------------------------------------
st.header("6 · Generate workdir")
workdir = st.text_input("Workdir (where config/ + outputs will live)",
                        placeholder=r"D:\CASPA_runs\MyExperiment")

if st.button("✅ Generate config + sample sheet", type="primary",
             disabled=not (workdir and data_path)):
    try:
        project = {"name": exp_name, "species_label": species, "description": description.strip()}
        input_block = {
            "pg_matrix": (None if is_spectronaut else os.path.abspath(data_path)),
            "spectronaut_tsv": (os.path.abspath(data_path) if is_spectronaut else None),
            "sample_sheet": "config/ms_inputs.tsv",
        }
        scp_block = {
            "min_protein_ids": int(min_prot),
            "custom_proteins": custom_proteins,
            "joint_embedding": {"n_pcs": int(n_pcs), "n_neighbors": int(n_neighbors),
                                "leiden_resolution": float(leiden_res), "seed": int(seed)},
            "llm": {"provider": provider, "model": model, "api_key": api_key,
                    "base_url": base_url, "condition_b": bool(condition_b)},
        }
        sheet = st.session_state.get("sheet")
        if sheet is None:
            sheet = pd.DataFrame({"sample_id": [], "sample_file": [], "batch": []})
        cfg_p, sheet_p = write_experiment(os.path.abspath(workdir), project, input_block,
                                          scp_block, sheet)
        st.success("Workdir is ready.")
        st.code(f"{cfg_p}\n{sheet_p}", language="text")
        st.markdown("**Next step — run the pipeline:**")
        st.code(f"python caspa/run.py --workdir \"{os.path.abspath(workdir)}\" --cores 8", language="bash")
        if not is_spectronaut and sheet is not None and len(sheet) and \
                len({int(b) for b in sheet["batch"] if int(b) > 0}) <= 1:
            st.info("All samples are in one batch — Harmony batch correction will be skipped.")
    except Exception as e:  # noqa: BLE001
        st.error(f"Could not write the workdir: {e}")
