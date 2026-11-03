"""Synthetic-data tests; no participant data included."""
from pathlib import Path
import sys
import unittest
import numpy as np

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(BASE/'src'))
from pfl.external_projection import align_controls, terminal_mass, frozen_projection, state_posterior


class TestProjection(unittest.TestCase):
    def test_alignment(self):
        source = np.array([[1., 10.], [2., 20.], [3., 30.]])
        ext_control = np.array([[11., 110.], [12., 120.], [13., 130.]])
        ext_all = np.vstack([ext_control, [14., 140.]])
        got = align_controls(ext_all, source, ext_control)
        np.testing.assert_allclose(got[:3], source)
        np.testing.assert_allclose(got[3], [4.,40.])

    def test_zero_iqr_is_rejected(self):
        with self.assertRaises(ValueError):
            align_controls(np.ones((3, 2)),np.ones((3, 2)),np.ones((3, 2)))

    def test_same_frozen_state_posterior_for_raw_and_aligned(self):
        source = np.array([[0., 0.], [1., 1.], [2., 2.]])
        external = source * 10 + 100
        params = (np.array([[0., 0.], [2., 2.]]), np.ones(2), 0.2, [1])
        result = frozen_projection(external, source, [True] * 3, [True] * 3, *params)
        np.testing.assert_allclose(result["scores"], source)
        np.testing.assert_allclose(result["terminal_mass"], terminal_mass(source, *params))
        np.testing.assert_allclose(result["posterior"], state_posterior(source, *params[:3]))

    def test_mass_is_probability(self):
        x = np.array([[0.],[1.],[2.]])
        cent = np.array([[0.],[2.]])
        mass = terminal_mass(x,cent,np.ones(2),0.2,[1])
        self.assertGreater(mass[-1], mass[0])
        self.assertTrue(((mass >= 0) & (mass <= 1)).all())


if __name__ == '__main__':
    unittest.main()
