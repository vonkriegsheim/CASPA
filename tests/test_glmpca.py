"""
Unit tests for bernoulli_glmpca() in pipeline/scripts/python/scp_joint_embedding_leiden.py.

Run with either
    python -m unittest discover -s tests -v
    python -m pytest tests -v
"""

import contextlib
import io
import sys
import unittest
import warnings
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "pipeline" / "scripts" / "python"
sys.path.insert(0, str(SCRIPT_DIR))

from scp_joint_embedding_leiden import bernoulli_glmpca  # noqa: E402

INFO_KEYS = {"n_iter", "deviance", "null_deviance", "deviance_explained",
             "cell_intercept", "rank", "converged"}


def _sigmoid(e):
    return 1.0 / (1.0 + np.exp(-e))


def simulate(n_cells=200, n_prot=300, seed=0, depth_sd=0.0):
    """
    Binary matrix drawn from logit P = a_i + b_j + (U V^T)_ij with rank-2 latent
    structure: dimension 1 is a two-group indicator, dimension 2 a continuous
    nuisance axis. depth_sd > 0 adds a per-cell intercept (detection depth).
    Returns (X01, group_labels).
    """
    rng = np.random.default_rng(seed)
    groups = np.repeat([0, 1], [n_cells // 2, n_cells - n_cells // 2])
    U = np.column_stack([np.where(groups == 0, -1.5, 1.5),
                         rng.standard_normal(n_cells)])
    V = rng.standard_normal((n_prot, 2)) * 0.8
    b = rng.normal(-0.5, 1.0, size=n_prot)
    a = rng.normal(0.0, depth_sd, size=n_cells) if depth_sd > 0 else np.zeros(n_cells)
    eta = a[:, None] + b[None, :] + U @ V.T
    X = (rng.random((n_cells, n_prot)) < _sigmoid(eta)).astype(float)
    return X, groups


def nearest_centroid_accuracy(scores, groups):
    cents = np.stack([scores[groups == g].mean(axis=0) for g in (0, 1)])
    d = ((scores[:, None, :] - cents[None, :, :]) ** 2).sum(axis=2)
    pred = d.argmin(axis=1)
    return (pred == groups).mean()


class TestRecoversStructure(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.X, cls.groups = simulate()
        cls.scores, cls.info = bernoulli_glmpca(cls.X, 2, verbose=False)

    def test_shapes_and_finite(self):
        self.assertEqual(self.scores.shape, (200, 2))
        self.assertTrue(np.isfinite(self.scores).all())

    def test_groups_separated(self):
        acc = nearest_centroid_accuracy(self.scores, self.groups)
        self.assertGreater(acc, 0.95, f"nearest-centroid accuracy {acc:.3f}")
        # Group means must be far apart relative to the within-group spread along
        # the axis that joins them (one-dimensional separability check).
        cents = np.stack([self.scores[self.groups == g].mean(axis=0) for g in (0, 1)])
        axis = cents[1] - cents[0]
        axis /= np.linalg.norm(axis)
        proj = self.scores @ axis
        within = np.sqrt(np.mean([proj[self.groups == g].var() for g in (0, 1)]))
        gap = proj[self.groups == 1].mean() - proj[self.groups == 0].mean()
        self.assertGreater(gap / within, 3.0)

    def test_deviance_explained_positive(self):
        self.assertGreater(self.info["deviance_explained"], 0.0)
        self.assertLess(self.info["deviance_explained"], 1.0)
        self.assertLess(self.info["deviance"], self.info["null_deviance"])

    def test_not_worse_than_gaussian_pca_on_separation(self):
        # Guard against a silently broken fit: the logistic embedding should
        # separate the groups about as well as a plain SVD of the centred calls.
        Xc = self.X - self.X.mean(axis=0, keepdims=True)
        U, S, _ = np.linalg.svd(Xc, full_matrices=False)
        acc_pca = nearest_centroid_accuracy(U[:, :2] * S[:2], self.groups)
        acc_glm = nearest_centroid_accuracy(self.scores, self.groups)
        self.assertGreaterEqual(acc_glm, acc_pca - 0.02)


class TestMonotoneDeviance(unittest.TestCase):
    def test_deviance_non_increasing_in_iterations(self):
        X, _ = simulate(n_cells=120, n_prot=150, seed=1)
        for cell_intercept in (False, True):
            devs = []
            for m in (1, 2, 5, 20, 60):
                _, info = bernoulli_glmpca(X, 3, cell_intercept=cell_intercept,
                                           max_iter=m, tol=0.0, verbose=False)
                self.assertEqual(info["n_iter"], m)
                devs.append(info["deviance"])
            for d_prev, d_next in zip(devs, devs[1:]):
                self.assertLessEqual(
                    d_next, d_prev * (1 + 1e-12),
                    f"deviance increased: {devs} (cell_intercept={cell_intercept})")
            self.assertLess(devs[-1], devs[0])

    def test_first_iteration_improves_on_null(self):
        X, _ = simulate(n_cells=120, n_prot=150, seed=2)
        _, info = bernoulli_glmpca(X, 2, max_iter=1, tol=0.0, verbose=False)
        self.assertLess(info["deviance"], info["null_deviance"])

    def test_tol_stops_early_and_flags_convergence(self):
        X, _ = simulate(n_cells=80, n_prot=100, seed=3)
        _, info = bernoulli_glmpca(X, 2, max_iter=500, tol=1e-2, verbose=False)
        self.assertLess(info["n_iter"], 500)
        self.assertTrue(info["converged"])
        _, info = bernoulli_glmpca(X, 2, max_iter=3, tol=0.0, verbose=False)
        self.assertEqual(info["n_iter"], 3)
        self.assertFalse(info["converged"])

    def test_unconverged_fit_prints_note(self):
        X, _ = simulate(n_cells=40, n_prot=50, seed=6)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            bernoulli_glmpca(X, 2, max_iter=2, tol=0.0, verbose=True)
        self.assertIn("max_iter", buf.getvalue())


class TestDeterminism(unittest.TestCase):
    def test_two_calls_identical(self):
        X, _ = simulate(n_cells=100, n_prot=120, seed=4)
        s1, i1 = bernoulli_glmpca(X, 4, max_iter=40, verbose=False)
        s2, i2 = bernoulli_glmpca(X, 4, max_iter=40, verbose=False)
        self.assertTrue(np.array_equal(s1, s2))
        self.assertEqual(i1, i2)

    def test_independent_of_global_rng_state(self):
        X, _ = simulate(n_cells=100, n_prot=120, seed=5)
        np.random.seed(1)
        s1, _ = bernoulli_glmpca(X, 3, max_iter=30, verbose=False)
        np.random.seed(999)
        np.random.random(10)
        s2, _ = bernoulli_glmpca(X, 3, max_iter=30, verbose=False)
        self.assertTrue(np.array_equal(s1, s2))


class TestEdgeCases(unittest.TestCase):
    def _fit_ok(self, X, k, expected_k=None, **kw):
        X = np.asarray(X)
        n, p = X.shape
        if expected_k is None:
            expected_k = max(1, min(k, n - 1, p - 1))
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            scores, info = bernoulli_glmpca(X, k, max_iter=60, verbose=False, **kw)
        self.assertEqual(scores.shape, (n, expected_k))
        self.assertTrue(np.isfinite(scores).all(), "non-finite scores")
        for key in ("deviance", "null_deviance"):
            self.assertTrue(np.isfinite(info[key]), f"non-finite {key}")
        self.assertEqual(info["rank"], expected_k)
        return scores, info

    def test_protein_detected_everywhere_and_nowhere(self):
        rng = np.random.default_rng(10)
        X = (rng.random((60, 50)) < 0.4).astype(float)
        X[:, 0] = 1.0
        X[:, 1] = 0.0
        for ci in (False, True):
            self._fit_ok(X, 5, cell_intercept=ci)

    def test_cell_with_all_zeros_or_all_ones(self):
        rng = np.random.default_rng(11)
        X = (rng.random((60, 50)) < 0.4).astype(float)
        X[0, :] = 0.0
        X[1, :] = 1.0
        for ci in (False, True):
            self._fit_ok(X, 5, cell_intercept=ci)

    def test_n_components_exceeds_dimensions(self):
        rng = np.random.default_rng(12)
        X = (rng.random((10, 8)) < 0.5).astype(float)
        self._fit_ok(X, 50, expected_k=7)
        X = (rng.random((8, 10)) < 0.5).astype(float)
        self._fit_ok(X, 50, expected_k=7)
        self._fit_ok(X, 0, expected_k=1)

    def test_tiny_3x3(self):
        X = np.array([[1, 0, 1], [0, 1, 1], [1, 1, 0]])
        self._fit_ok(X, 2, expected_k=2)
        self._fit_ok(X, 5, expected_k=2, cell_intercept=True)

    def test_bool_and_int_dtypes_match_float(self):
        rng = np.random.default_rng(13)
        Xb = rng.random((40, 30)) < 0.5
        sf, _ = self._fit_ok(Xb.astype(float), 3)
        sb, _ = self._fit_ok(Xb, 3)
        si, _ = self._fit_ok(Xb.astype(int), 3)
        self.assertTrue(np.array_equal(sf, sb))
        self.assertTrue(np.array_equal(sf, si))

    def test_constant_matrices(self):
        self._fit_ok(np.zeros((10, 10)), 2)
        self._fit_ok(np.ones((10, 10)), 2)


class TestCellIntercept(unittest.TestCase):
    def test_cell_intercept_runs_and_changes_fit(self):
        X, groups = simulate(n_cells=150, n_prot=200, seed=20, depth_sd=1.5)
        s0, i0 = bernoulli_glmpca(X, 2, cell_intercept=False, max_iter=60, verbose=False)
        s1, i1 = bernoulli_glmpca(X, 2, cell_intercept=True, max_iter=60, verbose=False)
        self.assertFalse(i0["cell_intercept"])
        self.assertTrue(i1["cell_intercept"])
        self.assertEqual(s0.shape, s1.shape)
        self.assertFalse(np.allclose(s0, s1))
        self.assertNotAlmostEqual(i0["deviance"], i1["deviance"])
        # Absorbing the depth term into an intercept should improve the fit ...
        self.assertLess(i1["deviance"], i0["deviance"])
        # ... and the group structure must still be recovered.
        self.assertGreater(nearest_centroid_accuracy(s1, groups), 0.95)


class TestInfoDict(unittest.TestCase):
    def test_keys_and_rank(self):
        X, _ = simulate(n_cells=50, n_prot=70, seed=30)
        for k, expected in ((3, 3), (100, 49)):
            scores, info = bernoulli_glmpca(X, k, max_iter=20, verbose=False)
            self.assertTrue(INFO_KEYS <= set(info), f"missing keys: {INFO_KEYS - set(info)}")
            self.assertEqual(info["rank"], expected)
            self.assertEqual(info["rank"], scores.shape[1])
            self.assertIsInstance(info["n_iter"], int)
            self.assertGreaterEqual(info["n_iter"], 1)
            self.assertLessEqual(info["n_iter"], 20)
            self.assertIsInstance(info["cell_intercept"], bool)
            self.assertIsInstance(info["converged"], bool)
            self.assertAlmostEqual(info["deviance_explained"],
                                   1.0 - info["deviance"] / info["null_deviance"])

    def test_verbose_prints(self):
        X, _ = simulate(n_cells=30, n_prot=40, seed=31)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            bernoulli_glmpca(X, 2, max_iter=5, verbose=True)
        self.assertIn("GLM-PCA", buf.getvalue())


if __name__ == "__main__":
    unittest.main()
