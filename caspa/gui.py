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

# LLM services. Gemini and DeepSeek expose OpenAI-compatible APIs, so they ride
# the existing provider="openai" code path with just a base_url + model.
LLM_SERVICES = {
    "OpenAI":                     dict(provider="openai", base_url="",
                                       model="gpt-5.2"),
    "Anthropic (Claude)":         dict(provider="claude", base_url="",
                                       model="claude-sonnet-4-6"),
    "Google Gemini":              dict(provider="openai",
                                       base_url="https://generativelanguage.googleapis.com/v1beta/openai/",
                                       model="gemini-2.5-flash"),  # flash is on the free tier; 2.5-pro is paid-only
    "DeepSeek":                   dict(provider="openai",
                                       base_url="https://api.deepseek.com",
                                       model="deepseek-chat"),
    "Custom (OpenAI-compatible)": dict(provider="openai", base_url="", model=""),
}


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


def native_pick(kind="file"):
    """Open an OS file/folder dialog and return the chosen path ('' if cancelled
    or unavailable). Works for the local-desktop use case; no-ops headless."""
    try:
        import tkinter as tk
        from tkinter import filedialog
        root = tk.Tk()
        root.withdraw()
        root.wm_attributes("-topmost", 1)
        if kind == "dir":
            path = filedialog.askdirectory(parent=root)
        else:
            path = filedialog.askopenfilename(
                parent=root,
                filetypes=[("Protein matrix / TSV", "*.tsv *.txt *.csv"), ("All files", "*.*")])
        root.destroy()
        return os.path.normpath(path) if path else ""
    except Exception:
        st.session_state["_browse_unavailable"] = True
        return ""


def _browse(state_key, kind="file"):
    """on_click callback: open a dialog and write the result into session_state."""
    p = native_pick(kind)
    if p:
        st.session_state[state_key] = p
        if state_key == "data_path":            # default the workdir to the data folder
            st.session_state["workdir"] = os.path.dirname(p)


def num_field(label, default, key, cast=float, help=None):
    """A plain text box for a number (no +/- steppers). Falls back to default."""
    raw = st.text_input(label, value=str(default), key=key, help=help)
    try:
        return cast(raw)
    except (TypeError, ValueError):
        st.caption(f"⚠️ '{raw}' isn't a number — using {default}.")
        return cast(default)


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

# session defaults
for k, v in {"data_path": "", "workdir": "",
             "llm_service": "OpenAI",
             "llm_model": LLM_SERVICES["OpenAI"]["model"],
             "llm_base_url": LLM_SERVICES["OpenAI"]["base_url"]}.items():
    st.session_state.setdefault(k, v)

st.title("🧫 CASPA — experiment setup")
st.caption("Build `config/caspa.json` + `config/ms_inputs.tsv` for a run. "
           "No coding required — fill in the form and click *Generate*.")

# ---- 1. Input data --------------------------------------------------------
st.header("1 · Input data")
input_kind = st.radio("Input type", ["DIA-NN / FragPipe pg_matrix", "Spectronaut long-format TSV"],
                      horizontal=True)
is_spectronaut = input_kind.startswith("Spectronaut")
pc1, pc2 = st.columns([5, 1])
with pc1:
    st.text_input("Path to the input file", key="data_path",
                  placeholder=r"D:\data\report.pg_matrix.tsv")
with pc2:
    st.write("")
    st.button("📂 Browse…", on_click=_browse, args=("data_path", "file"),
              use_container_width=True)
if st.session_state.get("_browse_unavailable"):
    st.caption("File dialog unavailable in this environment — type/paste the path instead.")

data_path = st.session_state.get("data_path", "")
if st.button("Load samples", type="primary", disabled=not data_path):
    if not os.path.isfile(data_path):
        st.error(f"File not found: {data_path}")
    elif is_spectronaut:
        st.session_state.sheet = pd.DataFrame({"sample_id": [], "sample_file": [], "batch": []})
        st.session_state.is_spectronaut = True
        st.info("Spectronaut input recorded. The sample sheet is auto-populated "
                "after conversion at run time; batches can be set then.")
    else:
        try:
            st.session_state.sheet = build_sample_sheet_df(data_path)
            st.session_state.is_spectronaut = False
            if not st.session_state.get("workdir"):
                st.session_state["workdir"] = os.path.dirname(os.path.abspath(data_path))
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
                   "(set these by hand)." if n_un else ""))
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
        min_prot = num_field("Min proteins per cell (QC floor)",
                             int(D.get("scp", {}).get("min_protein_ids", 400)), "f_minprot", int)
        seed = num_field("Random seed", int(emb.get("seed", 0)), "f_seed", int)
    with a2:
        n_pcs = num_field("PCA components (intensity)", int(emb.get("n_pcs", 20)), "f_npcs", int)
        n_neighbors = num_field("kNN / UMAP neighbours", int(emb.get("n_neighbors", 15)), "f_nn", int)
    with a3:
        leiden_res = num_field("Leiden resolution", float(emb.get("leiden_resolution", 0.8)),
                               "f_leiden", float)
    st.markdown("**Adaptive Harmony batch correction**")
    h1, h2, h3 = st.columns(3)
    with h1:
        h_entropy = num_field("Shannon batch-mixing target (0–1)",
                              float(emb.get("harmony_entropy_target", 0.6)), "f_entropy", float,
                              help="Harmony's diversity penalty is raised until mean batch-mixing "
                                   "entropy across Leiden clusters reaches this target.")
    with h2:
        h_theta = num_field("Harmony θ (start)", float(emb.get("harmony_theta", 2.0)), "f_theta", float)
    with h3:
        h_theta_max = num_field("Harmony θ (max)", float(emb.get("harmony_theta_max", 5.0)),
                                "f_thetamax", float)

# ---- 5. LLM annotation ----------------------------------------------------
st.header("5 · LLM annotation")


def _on_service_change():
    p = LLM_SERVICES[st.session_state["llm_service"]]
    st.session_state["llm_model"] = p["model"]
    st.session_state["llm_base_url"] = p["base_url"]


l1, l2 = st.columns(2)
with l1:
    st.selectbox("Service", list(LLM_SERVICES), key="llm_service", on_change=_on_service_change)
    st.text_input("Model", key="llm_model",
                  help="e.g. gpt-5.2 (recommended), claude-sonnet-4-6, gemini-2.5-pro, deepseek-chat.")
with l2:
    api_key = st.text_input("API key", type="password",
                            help="Stored only in this workdir's config/caspa.json (gitignored). "
                                 "Use the key for the selected service. Leave blank if submitting "
                                 "via Claude Code.")
    st.text_input("Base URL (OpenAI-compatible)", key="llm_base_url",
                  help="Auto-filled for Gemini / DeepSeek; leave blank for OpenAI / Anthropic.")
condition_b = st.checkbox("Use 3-round (Round-0 context) annotation",
                          value=bool(D.get("scp", {}).get("llm", {}).get("condition_b", False)))
with st.expander("Advanced LLM options"):
    max_tokens = num_field(
        "Max output tokens",
        int(D.get("scp", {}).get("llm", {}).get("max_tokens", 16000)), "f_maxtok", int,
        help="Cap on the annotation output. Auto-reduced to the model's own limit if that's "
             "lower (e.g. DeepSeek 8192), and a run warns if output is truncated — raise this "
             "only if you see such a warning.")
provider = LLM_SERVICES[st.session_state["llm_service"]]["provider"]
model = st.session_state["llm_model"]
base_url = st.session_state["llm_base_url"]

# ---- 6. Generate ----------------------------------------------------------
st.header("6 · Generate workdir")
wc1, wc2 = st.columns([5, 1])
with wc1:
    st.text_input("Workdir (where config/ + outputs will live)", key="workdir",
                  placeholder="defaults to the folder of your input file")
with wc2:
    st.write("")
    st.button("📂 Browse…", on_click=_browse, args=("workdir", "dir"),
              use_container_width=True, key="browse_wd")
workdir = st.session_state.get("workdir", "")

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
                                "leiden_resolution": float(leiden_res), "seed": int(seed),
                                "harmony_entropy_target": float(h_entropy),
                                "harmony_theta": float(h_theta),
                                "harmony_theta_max": float(h_theta_max)},
            "llm": {"provider": provider, "model": model, "api_key": api_key,
                    "base_url": base_url, "condition_b": bool(condition_b),
                    "max_tokens": int(max_tokens)},
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
        if not is_spectronaut and len(sheet) and \
                len({int(b) for b in sheet["batch"] if int(b) > 0}) <= 1:
            st.info("All samples are in one batch — Harmony batch correction will be skipped.")
    except Exception as e:  # noqa: BLE001
        st.error(f"Could not write the workdir: {e}")
