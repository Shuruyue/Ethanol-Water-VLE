"""Tests for excess property calculations and azeotrope detection."""

import unittest

import numpy as np

from analysis import (
    AzeotropePoint,
    compute_excess_curves,
    compute_gamma_curves,
    compute_isothermal_excess_curves,
    excess_gibbs_energy,
    excess_enthalpy,
    excess_heat_capacity,
    find_azeotrope,
)
from models import gamma_nrtl
from parameter_store import load_system_parameters
from solver import compute_txy


class TestExcessGibbsEnergy(unittest.TestCase):
    """Validate excess Gibbs energy calculations."""

    def setUp(self):
        system = load_system_parameters("ethanol_water_1atm")
        self.nrtl = system.nrtl
        self.T_K = 351.38  # ~78.2 deg C

    def test_boundary_near_zero(self):
        """GE should be near zero at composition boundaries."""
        ge_low = excess_gibbs_energy(0.001, self.T_K, self.nrtl)
        ge_high = excess_gibbs_energy(0.999, self.T_K, self.nrtl)
        self.assertLess(abs(ge_low), 50.0)
        self.assertLess(abs(ge_high), 50.0)

    def test_positive_at_midpoint(self):
        """GE must be positive at x1=0.5 for ethanol-water (positive deviation)."""
        ge = excess_gibbs_energy(0.5, self.T_K, self.nrtl)
        self.assertGreater(ge, 0.0)

    def test_formula_consistency(self):
        """GE = RT * sum(xi * ln(gamma_i))."""
        from math import log
        x1 = 0.4
        g1, g2 = gamma_nrtl(x1, self.T_K, self.nrtl)
        ge_expected = 8.314462618 * self.T_K * (x1 * log(g1) + (1 - x1) * log(g2))
        ge = excess_gibbs_energy(x1, self.T_K, self.nrtl)
        self.assertAlmostEqual(ge, ge_expected, places=6)


class TestExcessEnthalpy(unittest.TestCase):
    """Validate excess enthalpy calculations."""

    def setUp(self):
        system = load_system_parameters("ethanol_water_1atm")
        self.nrtl = system.nrtl
        self.T_K = 351.38

    def test_finite_at_midpoint(self):
        """HE must be a finite number."""
        he = excess_enthalpy(0.5, self.T_K, self.nrtl)
        self.assertTrue(np.isfinite(he))

    def test_stable_across_delta(self):
        """HE should be stable across different numerical derivative step sizes."""
        he_tight = excess_enthalpy(0.5, self.T_K, self.nrtl, delta_t=1e-4)
        he_default = excess_enthalpy(0.5, self.T_K, self.nrtl, delta_t=1e-3)
        he_loose = excess_enthalpy(0.5, self.T_K, self.nrtl, delta_t=1e-2)
        # All should agree within 1%
        self.assertAlmostEqual(he_tight, he_default, delta=abs(he_default) * 0.01 + 1.0)
        self.assertAlmostEqual(he_default, he_loose, delta=abs(he_default) * 0.01 + 1.0)


class TestExcessHeatCapacity(unittest.TestCase):
    """Validate excess heat capacity calculations."""

    def setUp(self):
        system = load_system_parameters("ethanol_water_1atm")
        self.nrtl = system.nrtl
        self.T_K = 351.38

    def test_finite_across_composition(self):
        """CPE must be finite for all compositions."""
        for x1 in np.linspace(0.05, 0.95, 20):
            cpe = excess_heat_capacity(float(x1), self.T_K, self.nrtl)
            self.assertTrue(np.isfinite(cpe), f"CPE not finite at x1={x1}")


class TestIsothermalCurves(unittest.TestCase):
    """Validate isothermal excess curve arrays."""

    def test_array_lengths(self):
        system = load_system_parameters("ethanol_water_1atm")
        x = np.linspace(0.05, 0.95, 30)
        ge, he, cpe = compute_isothermal_excess_curves(x, 78.2, system.nrtl)
        self.assertEqual(len(ge), 30)
        self.assertEqual(len(he), 30)
        self.assertEqual(len(cpe), 30)


class TestGammaCurves(unittest.TestCase):
    """Validate gamma curve arrays."""

    def test_array_lengths(self):
        system = load_system_parameters("ethanol_water_1atm")
        x = np.linspace(0.05, 0.95, 25)
        g1, g2 = compute_gamma_curves(x, 78.2, system.nrtl)
        self.assertEqual(len(g1), 25)
        self.assertEqual(len(g2), 25)
        self.assertTrue(np.all(g1 > 0))
        self.assertTrue(np.all(g2 > 0))


class TestAzeotropeDetection(unittest.TestCase):
    """Validate azeotrope detection against literature."""

    def setUp(self):
        self.system = load_system_parameters("ethanol_water_1atm")

    def test_azeotrope_exists_at_1atm(self):
        """Ethanol-water should have an azeotrope at 1 atm."""
        x_grid = np.linspace(0.01, 0.99, 200)
        x, y, T = compute_txy(
            pressure_kpa=101.325,
            component_1=self.system.component_1.antoine,
            component_2=self.system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=self.system.nrtl,
            x_grid=x_grid,
        )
        azeo = find_azeotrope(x, y, T, tolerance=0.01)
        self.assertIsNotNone(azeo)

    def test_azeotrope_composition(self):
        """Azeotrope should be near x1 = 0.894 (literature value)."""
        x_grid = np.linspace(0.01, 0.99, 200)
        x, y, T = compute_txy(
            pressure_kpa=101.325,
            component_1=self.system.component_1.antoine,
            component_2=self.system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=self.system.nrtl,
            x_grid=x_grid,
        )
        azeo = find_azeotrope(x, y, T, tolerance=0.01)
        self.assertAlmostEqual(azeo.x1, 0.894, delta=0.03)

    def test_azeotrope_temperature(self):
        """Azeotrope temperature should be near 78.15 degC (literature)."""
        x_grid = np.linspace(0.01, 0.99, 200)
        x, y, T = compute_txy(
            pressure_kpa=101.325,
            component_1=self.system.component_1.antoine,
            component_2=self.system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=self.system.nrtl,
            x_grid=x_grid,
        )
        azeo = find_azeotrope(x, y, T, tolerance=0.01)
        self.assertAlmostEqual(azeo.temperature_c, 78.15, delta=1.0)

    def test_no_azeotrope_below_tolerance(self):
        """With very tight tolerance, we might not find azeotrope on coarse grid."""
        x = np.array([0.1, 0.5, 0.9])
        y = np.array([0.3, 0.7, 0.95])
        T = np.array([90.0, 82.0, 79.0])
        azeo = find_azeotrope(x, y, T, tolerance=0.001)
        self.assertIsNone(azeo)

    def test_interpolated_azeotrope_crossing(self):
        """Crossing between grid points should return an interpolated azeotrope."""
        x = np.array([0.40, 0.60])
        y = np.array([0.45, 0.55])  # y-x: +0.05 -> -0.05, crossing at x=0.5
        T = np.array([81.0, 79.0])
        azeo = find_azeotrope(x, y, T, tolerance=1e-6)
        self.assertIsNotNone(azeo)
        self.assertAlmostEqual(azeo.x1, 0.5, places=8)
        self.assertAlmostEqual(azeo.y1, 0.5, places=8)
        self.assertAlmostEqual(azeo.temperature_c, 80.0, places=8)
        self.assertAlmostEqual(azeo.abs_error, 0.0, places=12)


if __name__ == "__main__":
    unittest.main()
