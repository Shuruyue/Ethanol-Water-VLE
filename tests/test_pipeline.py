"""Tests for the VLE pipeline."""

from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from analysis import compute_isothermal_excess_curves
from parameter_store import load_system_parameters
from pipeline import run_analysis


class TestPipeline(unittest.TestCase):
    """Smoke-test major outputs."""

    def test_run_without_plots(self):
        result = run_analysis(points=20, save_plots=False)
        self.assertEqual(len(result.nrtl.x), 20)
        self.assertEqual(len(result.nrtl.y), 20)
        self.assertEqual(len(result.nrtl.T_celsius), 20)
        self.assertEqual(len(result.output_files), 0)

    def test_excess_energy_sign(self):
        result = run_analysis(points=30, save_plots=False)
        self.assertGreater(result.GE_j_mol.max(), 0.0)

    def test_isothermal_cpe_curve(self):
        system = load_system_parameters("ethanol_water_1atm")
        x = np.linspace(0.05, 0.95, 25)
        ge, he, cpe = compute_isothermal_excess_curves(x, 78.2, system.nrtl)
        self.assertEqual(len(ge), 25)
        self.assertEqual(len(he), 25)
        self.assertEqual(len(cpe), 25)
        self.assertTrue(np.isfinite(cpe).all())


if __name__ == "__main__":
    unittest.main()
