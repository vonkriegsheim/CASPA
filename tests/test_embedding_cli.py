"""
End-to-end tests for scp_joint_embedding_leiden.py --det_embedding {pca,glmpca,glmpca_cellint}.

These run the full script as a subprocess on a real (small) SCP test set and are
slow (a few minutes in total). They are skipped when the input files are absent.
Override the input locations with CASPA_TEST_PIVOT / CASPA_TEST_MANIFEST.

Run with either
    python -m unittest tests.test_embedding_cli -v
    python -m pytest tests/test_embedding_cli.py -v
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "pipeline" / "scripts" / "python" / "scp_joint_embedding_leiden.py"
PRE_CHANGE_COMMIT = "84292d5"

_DEFAULT_DIR = Path(r"C:\Users\Alex\AppData\Local\Temp\claude\a4dfa14bdcb4ac46"
                    r"\01001967-73e9-4235-b4b9-6d4c84b21701\scratchpad\sk_glm\scp")
PIVOT = Path(os.environ.get("CASPA_TEST_PIVOT", _DEFAULT_DIR / "pivot_pack.tsv"))
MANIFEST = Path(os.environ.get("CASPA_TEST_MANIFEST", _DEFAULT_DIR / "qc" / "filtered.manifest.tsv"))

MODES = ("pca", "glmpca", "glmpca_cellint")
REQUIRED_OUTPUTS = ("scp_cluster_assignments.tsv", "scp_annotation.tsv", "scp_clustering_report.json")

inputs_present = PIVOT.is_file() and MANIFEST.is_file()


def _run_script(script: Path, outdir: Path, extra_args=()):
    cmd = [sys.executable, str(script),
           "--pg-matrix", str(PIVOT), "--manifest", str(MANIFEST),
           "--run-col", "sample_id", "--outdir", str(outdir),
           "--harmony_theta_max", "0", *extra_args]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc, time.time() - t0


@unittest.skipUnless(inputs_present, f"test inputs missing: {PIVOT} / {MANIFEST} (slow end-to-end test)")
class TestDetectionEmbeddingModes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="caspa_glm_cli_"))
        cls.runs = {}
        for mode in MODES:
            outdir = cls.tmp / mode
            proc, secs = _run_script(SCRIPT, outdir, ["--det_embedding", mode])
            cls.runs[mode] = {"proc": proc, "outdir": outdir, "seconds": secs}
        cls.manifest = pd.read_csv(MANIFEST, sep="\t", dtype=str)

    @classmethod
    def tearDownClass(cls):
        if os.environ.get("CASPA_TEST_KEEP"):
            print(f"\nkeeping outputs in {cls.tmp}")
        else:
            shutil.rmtree(cls.tmp, ignore_errors=True)

    def _assert_ok(self, mode):
        r = self.runs[mode]
        proc = r["proc"]
        self.assertEqual(proc.returncode, 0,
                         f"{mode}: exit {proc.returncode}\nSTDOUT:\n{proc.stdout[-3000:]}\n"
                         f"STDERR:\n{proc.stderr[-3000:]}")
        for name in REQUIRED_OUTPUTS:
            self.assertTrue((r["outdir"] / name).is_file(), f"{mode}: missing {name}")
        return r

    def _report(self, mode):
        return json.loads((self.runs[mode]["outdir"] / "scp_clustering_report.json").read_text(encoding="utf-8"))

    def test_all_modes_exit_zero_and_write_outputs(self):
        for mode in MODES:
            with self.subTest(mode=mode):
                self._assert_ok(mode)

    def test_report_records_mode(self):
        for mode in MODES:
            with self.subTest(mode=mode):
                self._assert_ok(mode)
                rep = self._report(mode)
                self.assertEqual(rep["params"]["detection_embedding"], mode)
                self.assertEqual(rep["status"], "OK")

    def test_glmpca_fit_block(self):
        self._assert_ok("pca")
        self.assertIsNone(self._report("pca")["glmpca_fit"])
        for mode in ("glmpca", "glmpca_cellint"):
            with self.subTest(mode=mode):
                self._assert_ok(mode)
                fit = self._report(mode)["glmpca_fit"]
                self.assertIsInstance(fit, dict)
                for key in ("n_iter", "deviance", "null_deviance", "deviance_explained",
                            "cell_intercept", "rank"):
                    self.assertIn(key, fit)
                self.assertGreater(fit["deviance_explained"], 0.0)
                self.assertLess(fit["deviance_explained"], 1.0)
                self.assertEqual(fit["cell_intercept"], mode == "glmpca_cellint")
                self.assertEqual(fit["rank"], self._report(mode)["params"]["n_pcs_det"])

    def test_every_manifest_cell_gets_a_cluster(self):
        wanted = set(self.manifest["sample_id"].astype(str).str.strip())
        wanted.discard("")
        for mode in MODES:
            with self.subTest(mode=mode):
                r = self._assert_ok(mode)
                assign = pd.read_csv(r["outdir"] / "scp_cluster_assignments.tsv", sep="\t", dtype=str)
                ann = pd.read_csv(r["outdir"] / "scp_annotation.tsv", sep="\t", dtype=str)
                self.assertEqual(set(assign["Run"]), wanted)
                self.assertEqual(len(assign), len(wanted))
                self.assertFalse(assign["cluster_id"].isna().any())
                self.assertTrue((assign["cluster_id"].str.len() > 0).all())
                self.assertEqual(set(ann["Run"]), wanted)
                self.assertFalse(ann["Condition"].isna().any())
                self.assertGreaterEqual(assign["cluster_id"].nunique(), 1)

    def test_glmpca_modes_differ_from_pca(self):
        # The detection block must actually change with the mode (otherwise the
        # flag would be a no-op). Compare UMAP coordinates, which are a
        # deterministic function of the embedding.
        frames = {}
        for mode in MODES:
            r = self._assert_ok(mode)
            frames[mode] = pd.read_csv(r["outdir"] / "scp_cluster_assignments.tsv", sep="\t")
        for mode in ("glmpca", "glmpca_cellint"):
            with self.subTest(mode=mode):
                same = (frames[mode][["umap_1", "umap_2"]].to_numpy()
                        == frames["pca"][["umap_1", "umap_2"]].to_numpy()).all()
                self.assertFalse(same, f"{mode} produced identical UMAP to pca")

    def test_pca_mode_matches_pre_change_script(self):
        """Regression: default (pca) mode must reproduce the pre-GLM-PCA script exactly."""
        if shutil.which("git") is None:
            self.skipTest("git not available")
        old_script = self.tmp / "old_script.py"
        proc = subprocess.run(
            ["git", "-C", str(REPO), "show",
             f"{PRE_CHANGE_COMMIT}:pipeline/scripts/python/scp_joint_embedding_leiden.py"],
            capture_output=True, text=True)
        if proc.returncode != 0:
            self.skipTest(f"git show {PRE_CHANGE_COMMIT} failed: {proc.stderr.strip()}")
        old_script.write_bytes(proc.stdout.encode("utf-8"))

        old_out = self.tmp / "pre_change"
        proc_old, _ = _run_script(old_script, old_out)
        self.assertEqual(proc_old.returncode, 0,
                         f"pre-change script failed:\n{proc_old.stderr[-3000:]}")

        new = self._assert_ok("pca")
        a_new = pd.read_csv(new["outdir"] / "scp_cluster_assignments.tsv", sep="\t", dtype=str)
        a_old = pd.read_csv(old_out / "scp_cluster_assignments.tsv", sep="\t", dtype=str)
        self.assertEqual(a_new["Run"].tolist(), a_old["Run"].tolist())
        self.assertEqual(a_new["cluster_id"].tolist(), a_old["cluster_id"].tolist())
        self.assertEqual(a_new["cluster_id_raw"].tolist(), a_old["cluster_id_raw"].tolist())
        ann_new = pd.read_csv(new["outdir"] / "scp_annotation.tsv", sep="\t", dtype=str)
        ann_old = pd.read_csv(old_out / "scp_annotation.tsv", sep="\t", dtype=str)
        pd.testing.assert_frame_equal(ann_new, ann_old)


if __name__ == "__main__":
    unittest.main()
