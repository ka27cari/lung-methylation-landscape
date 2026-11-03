"""Frozen synthetic regression fixture; contains no clinical observations."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pfl.bayesian import nig_fit


class BayesianRegressionTests(unittest.TestCase):
    def test_weighted_model_preserves_pre_refactor_coefficients(self):
        x = np.column_stack([
            np.ones(7), np.arange(7) / 7.0,
            np.array([1, 0, 1, 0, 0, 1, 1], dtype=float),
        ])
        y = np.array([0.3, -0.2, 0.7, 1.0, 0.0, 0.5, 1.5])
        w = np.array([0.8, 0.2, 0.5, 0.6, 0.1, 0.9, 0.6])
        fit = nig_fit(x, y, w, n_draw=9, prior_variance=10.0,
                      a0=2.0, b0=1.0, rng=np.random.default_rng(840))
        np.testing.assert_allclose(fit["mn"], [0.32027062, 0.78769383, 0.00671151], atol=1e-8)
        np.testing.assert_allclose(fit["qlo"], [-0.84963191, -0.74668657, -1.14891692], atol=1e-8)
        np.testing.assert_allclose(fit["qhi"], [1.49017316, 2.32207422, 1.16233995], atol=1e-8)
        self.assertAlmostEqual(fit["effective_n"], 5.542510121457488, places=10)
        np.testing.assert_allclose(fit["draws"][0], [0.41780146, 0.90493669, -0.53903142], atol=1e-8)


if __name__ == "__main__":
    unittest.main()
