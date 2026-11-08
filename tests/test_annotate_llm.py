#!/usr/bin/env python3
"""
Failure-handling tests for pipeline/scripts/python/annotate_clusters_llm.py.

The real script is run as a subprocess against a mock OpenAI-compatible HTTP
server (POST /v1/chat/completions and POST /v1/responses) that returns scripted
responses in sequence. Each test asserts on the exit code, the presence or
absence of the output TSV, and key stderr text.

Run with either of:
    python -m unittest tests/test_annotate_llm.py -v
    python -m pytest tests/test_annotate_llm.py -v

Skips cleanly when the `openai` package (needed by the script) is missing.
"""
from __future__ import annotations

import contextlib
import csv
import io
import json
import os
import re
import subprocess
import sys
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT_DIR = REPO / "pipeline" / "scripts" / "python"
SCRIPT = SCRIPT_DIR / "annotate_clusters_llm.py"
SNAKEFILE = REPO / "Snakefile"

try:
    import openai  # noqa: F401
    HAVE_OPENAI = True
except ImportError:  # pragma: no cover
    HAVE_OPENAI = False

CLUSTERS = ["C0", "C1", "C2"]
CELL_TYPES = {
    "C0": "Acinar cell",
    "C1": "Ductal cell",
    "C2": "Macrophage",
}


# ---------------------------------------------------------------------------
# Mock OpenAI-compatible server
# ---------------------------------------------------------------------------

class _Handler(BaseHTTPRequestHandler):
    """Serves scripted responses in order; records every request."""
    server_version = "MockLLM/1.0"

    def log_message(self, *args):  # silence default stderr logging
        pass

    def _send_json(self, status: int, payload: dict) -> None:
        data = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(data)))
        self.end_headers()
        self.wfile.write(data)

    def do_POST(self):  # noqa: N802
        n = int(self.headers.get("Content-Length") or 0)
        raw = self.rfile.read(n) if n else b""
        try:
            body = json.loads(raw or b"{}")
        except ValueError:
            body = {"_raw": raw.decode("utf-8", "replace")}
        mock = self.server.mock  # type: ignore[attr-defined]
        with mock.lock:
            idx = len(mock.requests)
            mock.requests.append({"path": self.path, "body": body})
            step = mock.script[idx] if idx < len(mock.script) else None

        if step is None:
            # Unexpected extra call. 400 is not retried by the openai client,
            # so the script fails fast instead of waiting out retry back-off.
            self._send_json(400, {"error": {"message": "mock: no scripted response "
                                                       f"for call #{idx + 1}",
                                            "type": "invalid_request_error"}})
            return
        if "status" in step and step["status"] != 200:
            self._send_json(step["status"], {"error": {"message": step.get("text", "mock error"),
                                                       "type": "mock_error"}})
            return

        text = step.get("text", "")
        finish = step.get("finish", "stop")
        model = body.get("model", "test-model")
        if self.path.rstrip("/").endswith("/chat/completions"):
            payload = {
                "id": f"chatcmpl-mock{idx}", "object": "chat.completion", "created": 0,
                "model": model,
                "choices": [{"index": 0, "finish_reason": finish,
                             "message": {"role": "assistant", "content": text}}],
                "usage": {"prompt_tokens": 1, "completion_tokens": 1, "total_tokens": 2},
            }
        elif self.path.rstrip("/").endswith("/responses"):
            # Responses API: status "incomplete" is the truncation signal.
            status = "incomplete" if finish in ("length", "incomplete") else "completed"
            payload = {
                "id": f"resp_mock{idx}", "object": "response", "created_at": 0,
                "status": status, "model": model,
                "output": [{"type": "message", "id": f"msg_mock{idx}", "status": "completed",
                            "role": "assistant",
                            "content": [{"type": "output_text", "text": text,
                                         "annotations": []}]}],
                "usage": {"input_tokens": 1, "output_tokens": 1, "total_tokens": 2},
            }
        else:
            self._send_json(404, {"error": {"message": f"mock: unknown path {self.path}"}})
            return
        self._send_json(200, payload)


class MockLLM:
    def __init__(self):
        self.lock = threading.Lock()
        self.script: list[dict] = []
        self.requests: list[dict] = []
        self.httpd = ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
        self.httpd.mock = self  # type: ignore[attr-defined]
        self.port = self.httpd.server_address[1]
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)

    def start(self):
        self.thread.start()

    def stop(self):
        self.httpd.shutdown()
        self.httpd.server_close()

    def reset(self, script: list[dict]) -> None:
        with self.lock:
            self.script = list(script)
            self.requests = []

    @property
    def base_url(self) -> str:
        return f"http://127.0.0.1:{self.port}/v1"

    def system_prompts(self) -> list[str]:
        out = []
        for r in self.requests:
            b = r["body"]
            if "messages" in b:
                out.append(next((m["content"] for m in b["messages"] if m["role"] == "system"), ""))
            else:
                out.append(b.get("instructions", ""))
        return out

    def user_prompts(self) -> list[str]:
        out = []
        for r in self.requests:
            b = r["body"]
            if "messages" in b:
                out.append(next((m["content"] for m in b["messages"] if m["role"] == "user"), ""))
            else:
                out.append(b.get("input", ""))
        return out


# ---------------------------------------------------------------------------
# Fixtures: prompt markdown, pivot, annotation (make_cluster_summary.py format)
# ---------------------------------------------------------------------------

GENES = ["PRSS1", "CPA1", "KRT19", "SOX9", "CD68", "PTPRC", "INS", "GCG", "VIM", "ACTB"]
# per-cluster "expressed" genes used to fill the pivot
EXPR = {
    "C0": {"PRSS1", "CPA1", "ACTB", "VIM"},
    "C1": {"KRT19", "SOX9", "ACTB", "VIM"},
    "C2": {"CD68", "PTPRC", "ACTB", "VIM", "INS"},
}
N_CELLS_PER_CLUSTER = 4


def runs_for(cluster: str) -> list[str]:
    return [f"20250101_mock_plate1_{cluster}_cell{i}_S1-A{i}_1_{1000 + i}"
            for i in range(1, N_CELLS_PER_CLUSTER + 1)]


def write_inputs(tmp: Path) -> dict[str, Path]:
    """Write cluster_llm_prompt.md, pivot_pack.tsv, scp_annotation.tsv."""
    # --- prompt (mirrors the sections make_cluster_summary.py writes) ---
    L = []
    L += ["# Cluster marker summary for cell type annotation", ""]
    L += ["## Experiment context", "",
          "- **Experiment**: mock_pancreas",
          "- **Species**: human",
          "- **Description**: Human pancreas, single-cell LFQ proteomics (test fixture)", ""]
    L += ["## Run metadata (auto-populated at runtime)", "",
          "- **Prompt version**: v2",
          "- **Model**: {model_name}",
          "- **Temperature**: 0",
          "- **Run timestamp**: {iso_timestamp}",
          "- **Prompt hash (SHA256)**: {sha256_of_prompt_text}", "",
          "*These fields are populated by the calling script and stored alongside the",
          "raw LLM output for audit trail purposes.*", ""]
    L += ["You are a single-cell proteomics expert. Below is a per-cluster summary from a "
          "label-free single-cell proteomics (SCP) experiment. For each cluster you are "
          "given statistical evidence across three complementary modalities.", "",
          "**Your task for each cluster:**", "",
          "1. Propose the most likely **cell type label**.", "",
          "The researcher can decide post-hoc whether to filter these clusters.", "",
          "---", ""]
    for cl in CLUSTERS:
        L += [f"## Cluster {cl}", "",
              f"- **n_cells**: {N_CELLS_PER_CLUSTER}",
              "- **median_n_detected**: 1200",
              "- **consensus markers**:", "",
              "  | Gene | score | n_modalities |",
              "  |------|-------|--------------|"]
        for g in sorted(EXPR[cl]):
            L.append(f"  | {g} | 0.90 | 3/3 |")
        L += ["", "- **detection markers**:", "",
              "  | Gene | det_rate_in | det_rate_out | qvalue |",
              "  |------|-------------|--------------|--------|"]
        for g in sorted(EXPR[cl]):
            L.append(f"  | {g} | 1.00 | 0.10 | 1.0e-03 |")
        L.append("")
    prompt = tmp / "cluster_llm_prompt.md"
    prompt.write_text("\n".join(L), encoding="utf-8", newline="\n")

    # --- annotation: Run / Condition / BioRep / Batch ---
    ann_rows = []
    for cl in CLUSTERS:
        for i, run in enumerate(runs_for(cl), start=1):
            ann_rows.append((run, cl, str(i), "1"))
    annotation = tmp / "scp_annotation.tsv"
    with open(annotation, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Run", "Condition", "BioRep", "Batch"])
        w.writerows(ann_rows)

    # --- pivot: proteins x cells, sample columns are raw-file paths (stem = Run) ---
    sample_cols = [f"input\\{run}.d" for cl in CLUSTERS for run in runs_for(cl)]
    pivot = tmp / "pivot_pack.tsv"
    with open(pivot, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["Protein.Group", "Protein.Names", "Genes", "First.Protein.Description"]
                   + sample_cols)
        for k, g in enumerate(GENES):
            vals = []
            for cl in CLUSTERS:
                for _ in runs_for(cl):
                    vals.append(f"{12.0 + k * 0.1:.3f}" if g in EXPR[cl] else "")
            w.writerow([f"P{10000 + k}", f"{g}_HUMAN", g, f"{g} protein"] + vals)
    return {"prompt": prompt, "pivot": pivot, "annotation": annotation}


# ---------------------------------------------------------------------------
# Scripted model responses
# ---------------------------------------------------------------------------

TABLE_HEADER = ("| Cluster | Suggested cell type | Key supporting markers | "
                "Contradictions / caveats | Confidence | Resolving markers (if Low/ambiguous) |\n"
                "|---|---|---|---|---|---|\n")


def annotation_table(cell_types: dict[str, str], id_fmt: str = "{}") -> str:
    rows = []
    for cl, ct in cell_types.items():
        rows.append(f"| {id_fmt.format(cl)} | {ct} | ACTB, VIM | none | High | - |")
    return TABLE_HEADER + "\n".join(rows) + "\n"


def pass1_response(cell_types=None, id_fmt="{}", requested="INS, GCG") -> str:
    cell_types = cell_types or CELL_TYPES
    return (annotation_table(cell_types, id_fmt)
            + "\n### Narrative\n\n"
            + "\n".join(f"{cl}: {ct} (ACTB det_rate_in=1.00 vs 0.10, q=0.0010)"
                        for cl, ct in cell_types.items())
            + f"\n\nREQUESTED_MARKERS: {requested}\n")


def pass2_response(cell_types=None, id_fmt="{}") -> str:
    cell_types = cell_types or CELL_TYPES
    return (annotation_table(cell_types, id_fmt)
            + "\n### Rationale\n\nNo change from Pass 1 for all clusters.\n\n"
            "RECOMMENDED_PANEL:\n"
            "| Cell type | Recommended markers | Rationale |\n"
            "| --- | --- | --- |\n"
            "| Acinar cell | PRSS1, CPA1 | Confirmed detectable |\n"
            "| Ductal cell | KRT19, SOX9 | Confirmed detectable |\n"
            "| Macrophage | CD68, PTPRC | Confirmed detectable |\n")


ROUND0_RESPONSE = ("### 1. Developmental / tissue vocabulary constraints\nAdult pancreas.\n"
                   "### 2. Expected ambient signals\nPRSS1, keratins.\n"
                   "### 3. Non-self protein acquisition\nMacrophages phagocytose.\n"
                   "### 4. Expected artefacts\nDissociation stress.\n")


def step(text: str, finish: str = "stop") -> dict:
    return {"text": text, "finish": finish}


# ---------------------------------------------------------------------------
# Test case base
# ---------------------------------------------------------------------------

@unittest.skipUnless(HAVE_OPENAI, "openai package not installed")
class AnnotateLLMEndToEnd(unittest.TestCase):
    mock: MockLLM

    @classmethod
    def setUpClass(cls):
        cls.mock = MockLLM()
        cls.mock.start()

    @classmethod
    def tearDownClass(cls):
        cls.mock.stop()

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory(prefix="annot_llm_")
        self.tmp = Path(self._tmp.name)
        self.inputs = write_inputs(self.tmp)
        self.out_dir = self.tmp / "llm"
        self.out_tsv = self.out_dir / "cluster_cell_type_annotations.tsv"

    def tearDown(self):
        self._tmp.cleanup()

    # -- helpers ------------------------------------------------------------
    def run_script(self, script: list[dict] | None, extra: list[str] = (),
                   with_key: bool = True, model: str = "test-model"):
        if script is not None:
            self.mock.reset(script)
        env = os.environ.copy()
        for k in ("ELM_API_KEY", "OPENAI_API_KEY", "ELM_BASE_URL", "ANTHROPIC_API_KEY"):
            env.pop(k, None)   # ELM_BASE_URL would override --base-url
        if with_key:
            env["ELM_API_KEY"] = "test"
        env["CASPA_LLM_TIMEOUT"] = "60"
        env["PYTHONIOENCODING"] = "utf-8"
        cmd = [sys.executable, str(SCRIPT),
               "--prompt-md", str(self.inputs["prompt"]),
               "--pivot", str(self.inputs["pivot"]),
               "--annotation", str(self.inputs["annotation"]),
               "--out-tsv", str(self.out_tsv),
               "--out-dir", str(self.out_dir),
               "--model", model, "--provider", "openai",
               "--base-url", self.mock.base_url,
               "--max-tokens", "2000",
               # point at a non-existent PanglaoDB file so coverage scoring is skipped
               "--panglaodb", str(self.tmp / "no_panglaodb.tsv.gz"),
               *extra]
        return subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8",
                              errors="replace", timeout=180, env=env, cwd=str(self.tmp))

    def read_tsv(self) -> list[dict]:
        with open(self.out_tsv, encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    def assert_round0_then_passes(self, n_expected: int, condition_b: bool = True):
        sys_prompts = self.mock.system_prompts()
        self.assertEqual(len(sys_prompts), n_expected,
                         f"expected {n_expected} LLM calls, got {len(sys_prompts)}")
        i = 0
        if condition_b:
            self.assertIn("preparing to annotate", sys_prompts[0], "first call is Round 0")
            i = 1
        self.assertIn("REQUESTED_MARKERS", sys_prompts[i], "pass 1 system prompt")
        if len(sys_prompts) > i + 1:
            self.assertIn("RECOMMENDED_PANEL", sys_prompts[i + 1], "pass 2 system prompt")

    # -- 1. happy path --------------------------------------------------------
    def test_01_happy_path_condition_b_default(self):
        r = self.run_script([step(ROUND0_RESPONSE), step(pass1_response()), step(pass2_response())])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(self.out_tsv.exists())
        rows = self.read_tsv()
        self.assertEqual([x["Cluster"] for x in rows], CLUSTERS)
        self.assertEqual({x["Cluster"]: x["cell_type"] for x in rows}, CELL_TYPES)
        self.assert_round0_then_passes(3, condition_b=True)
        # Round 0 only sees the experiment context, never the cluster data
        self.assertNotIn("## Cluster C0", self.mock.user_prompts()[0])
        self.assertIn("Experiment context", self.mock.user_prompts()[0])
        # pass 2 received the supplemental pivot lookup for INS/GCG
        p2_user = self.mock.user_prompts()[2]
        self.assertIn("Supplemental: raw per-cluster detection rates", p2_user)
        self.assertRegex(p2_user, r"\| INS \| 0\.00 \(0/4\) \| 0\.00 \(0/4\) \| 1\.00 \(4/4\) \|")
        self.assertNotIn("WARNING", r.stderr)
        self.assertTrue((self.out_dir / "round0_constraints.md").exists())
        self.assertTrue((self.out_dir / "scp_recommended_markers_flat.txt").exists())
        meta = json.loads((self.out_dir / "annotation_run_metadata.json").read_text())
        self.assertTrue(meta["condition_b"])
        self.assertEqual(meta["n_annotations_parsed"], 3)

    def test_01b_happy_path_no_condition_b(self):
        r = self.run_script([step(pass1_response()), step(pass2_response())],
                            extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual({x["Cluster"]: x["cell_type"] for x in self.read_tsv()}, CELL_TYPES)
        self.assert_round0_then_passes(2, condition_b=False)
        self.assertFalse((self.out_dir / "round0_constraints.md").exists())

    def test_01c_happy_path_responses_api_gpt5(self):
        """A gpt-5* model name must be routed to POST /v1/responses."""
        r = self.run_script([step(pass1_response()), step(pass2_response())],
                            extra=["--no-condition-b"], model="gpt-5-mock")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual({x["Cluster"]: x["cell_type"] for x in self.read_tsv()}, CELL_TYPES)
        self.assertEqual([q["path"] for q in self.mock.requests],
                         ["/v1/responses", "/v1/responses"])

    # -- 2. pass 2 table missing a cluster ----------------------------------
    def test_02_pass2_missing_cluster_exits_2(self):
        two = {"C0": CELL_TYPES["C0"], "C1": CELL_TYPES["C1"]}
        r = self.run_script([step(pass1_response(two)), step(pass2_response(two))],
                            extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertFalse(self.out_tsv.exists(), "no annotation TSV must be written")
        self.assertIn("ERROR: incomplete LLM annotation after pass 2", r.stderr)
        self.assertRegex(r.stderr, r"1 of 3 clusters have no row: C2")
        # raw model output is kept on disk for inspection
        self.assertTrue((self.out_dir / "cluster_annotation_raw_pass2.txt").exists())

    def test_02b_cluster_missing_from_pass2_only_is_filled_from_pass1(self):
        """A truncated pass 2 must not discard a complete pass 1 call: the cluster
        keeps its pass 1 row, the run succeeds with a WARNING and the substitution
        is recorded in the metadata sidecar."""
        p2 = pass2_response({"C0": CELL_TYPES["C0"], "C1": CELL_TYPES["C1"]})
        r = self.run_script([step(pass1_response()), step(p2)], extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("pass 2 returned no row for 1 cluster(s) (C2); their pass 1 "
                      "annotation was kept", r.stderr)
        rows = self.read_tsv()
        self.assertEqual([x["Cluster"] for x in rows], ["C0", "C1", "C2"])
        self.assertEqual(rows[2]["cell_type"], CELL_TYPES["C2"])
        meta = json.loads((self.out_dir / "annotation_run_metadata.json").read_text(encoding="utf-8"))
        self.assertEqual(meta["clusters_filled_from_pass1"], ["C2"])

    def test_02c_unknown_is_a_legitimate_call(self):
        """'Unknown' is a valid low-confidence answer, not an empty cell."""
        ct = dict(CELL_TYPES, C1="Unknown")
        r = self.run_script([step(pass1_response(ct)), step(pass2_response(ct))],
                            extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.read_tsv()[1]["cell_type"], "Unknown")

    # -- 3. placeholder cell type ------------------------------------------
    def test_03_pass2_placeholder_cell_type_exits_2(self):
        ct = dict(CELL_TYPES, C1="NA")
        r = self.run_script([step(pass1_response()), step(pass2_response(ct))],
                            extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertFalse(self.out_tsv.exists())
        self.assertRegex(r.stderr, r"1 cluster\(s\) have an empty cell-type call: C1")

    # -- 4. empty pass 2 response -------------------------------------------
    def test_04_pass2_empty_response_is_clean_error(self):
        r = self.run_script([step(pass1_response()), step("")], extra=["--no-condition-b"])
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse(self.out_tsv.exists())
        self.assertIn("returned an empty response", r.stderr)
        # must be a clean, actionable error (no Python traceback), exit 2 like
        # the other "annotation could not be produced" failures
        self.assertNotIn("Traceback", r.stderr)
        self.assertIn("ERROR:", r.stderr)
        self.assertEqual(r.returncode, 2, r.stderr)

    # -- 5. pass 2 prose, no table -> fallback to pass 1 ---------------------
    def test_05_pass2_prose_falls_back_to_pass1(self):
        prose = "I confirm my previous annotations without changes.\n\nRECOMMENDED_PANEL:\n" \
                "| Cell type | Recommended markers | Rationale |\n| --- | --- | --- |\n" \
                "| Acinar cell | PRSS1 | ok |\n"
        r = self.run_script([step(pass1_response()), step(prose)], extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("WARNING: pass 2 returned no parsable table; falling back to the pass 1",
                      r.stderr)
        self.assertEqual({x["Cluster"]: x["cell_type"] for x in self.read_tsv()}, CELL_TYPES)
        meta = json.loads((self.out_dir / "annotation_run_metadata.json").read_text())
        self.assertFalse(meta["pass2_used"])

    # -- 6. pass 1 incomplete, pass 2 complete ------------------------------
    def test_06_pass1_incomplete_pass2_complete_succeeds(self):
        """Pass 2 re-annotates every cluster from the full data, so a cluster
        missing after pass 1 can still be filled in. Pass 1 should only warn;
        the final table is what is enforced."""
        p1 = pass1_response({"C0": CELL_TYPES["C0"], "C2": CELL_TYPES["C2"]})
        r = self.run_script([step(p1), step(pass2_response())], extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual({x["Cluster"]: x["cell_type"] for x in self.read_tsv()}, CELL_TYPES)
        self.assertIn("WARNING: incomplete LLM annotation after pass 1", r.stderr)
        self.assertIn("C1", r.stderr)
        self.assertEqual(len(self.mock.requests), 2, "pass 2 must still be called")

    def test_06b_pass1_and_pass2_incomplete_exits_2(self):
        p1 = pass1_response({"C0": CELL_TYPES["C0"], "C2": CELL_TYPES["C2"]})
        p2 = pass2_response({"C0": CELL_TYPES["C0"], "C2": CELL_TYPES["C2"]})
        r = self.run_script([step(p1), step(p2)], extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertFalse(self.out_tsv.exists())
        self.assertIn("ERROR: incomplete LLM annotation after pass 2", r.stderr)
        self.assertRegex(r.stderr, r"no row: C1")

    # -- 7. --allow-incomplete ----------------------------------------------
    def test_07_allow_incomplete_writes_partial_table(self):
        two = {"C0": CELL_TYPES["C0"], "C1": CELL_TYPES["C1"]}
        r = self.run_script([step(pass1_response(two)), step(pass2_response(two))],
                            extra=["--no-condition-b", "--allow-incomplete"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("WARNING: incomplete LLM annotation after pass 2", r.stderr)
        self.assertRegex(r.stderr, r"no row: C2")
        rows = self.read_tsv()
        self.assertEqual([x["Cluster"] for x in rows], ["C0", "C1"])

    # -- 8. missing API key -------------------------------------------------
    def test_08_no_api_key_exits_1(self):
        r = self.run_script([], with_key=False)
        self.assertEqual(r.returncode, 1, r.stderr)
        self.assertFalse(self.out_tsv.exists())
        self.assertIn("ERROR: no API key found", r.stderr)
        self.assertEqual(len(self.mock.requests), 0)

    def test_08b_no_api_key_allow_skip_writes_header_only(self):
        r = self.run_script([], with_key=False, extra=["--allow-skip-without-key"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("WARNING: no API key found", r.stderr)
        self.assertTrue(self.out_tsv.exists())
        text = self.out_tsv.read_text(encoding="utf-8")
        self.assertEqual(text.splitlines(), ["Cluster\tcell_type\tkey_evidence\tcontradictions"
                                             "\tconfidence\tresolving_markers"])
        self.assertEqual(self.read_tsv(), [])
        self.assertEqual(len(self.mock.requests), 0)

    # -- 9. truncation -------------------------------------------------------
    def test_09_truncated_but_complete_warns_and_succeeds(self):
        r = self.run_script([step(pass1_response()), step(pass2_response(), finish="length")],
                            extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("WARNING: test-model response was truncated (finish=length)", r.stderr)
        self.assertEqual(len(self.read_tsv()), 3)

    def test_09b_truncated_and_incomplete_exits_2(self):
        one = {"C0": CELL_TYPES["C0"]}
        r = self.run_script([step(pass1_response(one)), step(pass2_response(one), finish="length")],
                            extra=["--no-condition-b"])
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertFalse(self.out_tsv.exists())
        self.assertIn("response was truncated (finish=length)", r.stderr)
        self.assertRegex(r.stderr, r"2 of 3 clusters have no row: C1, C2")

    def test_09c_truncated_responses_api_status_incomplete(self):
        r = self.run_script([step(pass1_response()), step(pass2_response(), finish="incomplete")],
                            extra=["--no-condition-b"], model="gpt-5-mock")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("truncated (finish=incomplete)", r.stderr)

    # -- 10. cluster id formatting -------------------------------------------
    def test_10_cluster_id_variants_do_not_fail_completeness(self):
        for i, fmt in enumerate(("**{}**", "Cluster {}", "{} ", "`{}`", "**Cluster {}**")):
            if i:
                self.tearDown()
                self.setUp()   # fresh tempdir per variant
            with self.subTest(id_fmt=fmt):
                r = self.run_script([step(pass1_response(id_fmt=fmt)),
                                     step(pass2_response(id_fmt=fmt))],
                                    extra=["--no-condition-b"])
                self.assertEqual(r.returncode, 0, f"id format {fmt!r}:\n{r.stderr}")
                self.assertEqual([x["Cluster"] for x in self.read_tsv()], CLUSTERS)


# ---------------------------------------------------------------------------
# Unit tests on the parser (import the script as a module)
# ---------------------------------------------------------------------------

class ParseAnnotationTableUnit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(SCRIPT_DIR))
        import annotate_clusters_llm as mod  # noqa: E402
        cls.mod = mod

    def test_cluster_id_normalisation(self):
        cases = {
            "C1": "C1", "**C1**": "C1", "Cluster C1": "C1", "C1 ": "C1", "`C1`": "C1",
            "**Cluster C1**": "C1", "cluster C1": "C1", "*C1*": "C1",
            # conservative: do not touch IDs that merely contain the word
            "Cluster_7": "Cluster_7", "Clusters C1": "Clusters C1",
        }
        for raw, want in cases.items():
            with self.subTest(raw=raw):
                rows = self.mod.parse_annotation_table(
                    annotation_table({raw: "Acinar cell"}))
                self.assertEqual(len(rows), 1)
                self.assertEqual(rows[0]["Cluster"], want)
                self.assertEqual(rows[0]["cell_type"], "Acinar cell")

    def test_bold_cell_type_is_unwrapped(self):
        rows = self.mod.parse_annotation_table(annotation_table({"C0": "**Acinar cell**"}))
        self.assertEqual(rows[0]["cell_type"], "Acinar cell")

    def test_expected_ids_from_prompt(self):
        with tempfile.TemporaryDirectory() as d:
            p = write_inputs(Path(d))["prompt"].read_text(encoding="utf-8")
        self.assertEqual(self.mod.expected_cluster_ids(p), CLUSTERS)

    def test_completeness_check(self):
        rows = self.mod.parse_annotation_table(annotation_table({"C0": "A", "C1": "NA"}))
        err = io.StringIO()
        with contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as cm:
            self.mod.check_annotation_completeness(rows, CLUSTERS, stage="pass 2")
        self.assertEqual(cm.exception.code, 2)
        self.assertIn("ERROR: incomplete LLM annotation after pass 2", err.getvalue())
        self.assertIn("no row: C2", err.getvalue())
        self.assertIn("empty cell-type call: C1", err.getvalue())
        # --allow-incomplete -> warning, problems returned
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            problems = self.mod.check_annotation_completeness(rows, CLUSTERS, allow_incomplete=True)
        self.assertEqual(len(problems), 2)
        self.assertTrue(err.getvalue().startswith("WARNING: "))
        # enforce=False (pass 1) -> warning even without --allow-incomplete
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            problems = self.mod.check_annotation_completeness(rows, CLUSTERS, stage="pass 1",
                                                              enforce=False)
        self.assertEqual(len(problems), 2)
        self.assertIn("WARNING: incomplete LLM annotation after pass 1", err.getvalue())
        self.assertEqual(self.mod.check_annotation_completeness(
            self.mod.parse_annotation_table(annotation_table(CELL_TYPES)), CLUSTERS), [])

    def test_empty_response_error_is_runtime_error(self):
        self.assertTrue(issubclass(self.mod.LLMResponseError, RuntimeError))


# ---------------------------------------------------------------------------
# 11. Snakefile wiring (static)
# ---------------------------------------------------------------------------

class SnakefileWiring(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        src = SNAKEFILE.read_text(encoding="utf-8")
        m = re.search(r"^rule scp_llm_annotation:\n(.*?)(?=^rule |\Z)", src, re.S | re.M)
        assert m, "rule scp_llm_annotation not found in Snakefile"
        cls.rule = m.group(1)

    def test_params_read_the_config_keys(self):
        self.assertRegex(self.rule, r'condition_b\s*=\s*_c\("scp",\s*"llm",\s*"condition_b"')
        self.assertRegex(self.rule, r'skip_no_key\s*=\s*_c\("scp",\s*"llm",\s*"skip_without_key"')
        self.assertRegex(self.rule, r'allow_incomplete\s*=\s*_c\("scp",\s*"llm",\s*"allow_incomplete"')

    def test_flags_follow_config(self):
        self.assertRegex(self.rule, r'cond_b_flag\s*=\s*" --condition-b" if params\.condition_b '
                                    r'else " --no-condition-b"')
        self.assertRegex(self.rule, r'skip_flag\s*=\s*" --allow-skip-without-key" '
                                    r'if params\.skip_no_key else ""')
        self.assertRegex(self.rule, r'incomplete_flag\s*=\s*" --allow-incomplete" '
                                    r'if params\.allow_incomplete else ""')
        shell_call = self.rule[self.rule.index("shell("):]
        for var in ("{cond_b_flag}", "{skip_flag}", "{incomplete_flag}"):
            self.assertIn(var, shell_call, f"{var} is not interpolated into the shell command")

    def test_script_accepts_the_flags(self):
        src = SCRIPT.read_text(encoding="utf-8")
        for flag in ("--no-condition-b", "--condition-b", "--allow-skip-without-key",
                     "--allow-incomplete"):
            self.assertIn(f'"{flag}"', src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
