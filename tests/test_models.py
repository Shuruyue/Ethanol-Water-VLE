"""Precision tests for thermodynamic models (NRTL, Van Laar, Antoine)."""

import unittest
from math import exp, log

import numpy as np

from models import (
    AntoineParams,
    NRTLParams,
    VanLaarParams,
    dpsat_dT_kpa_per_celsius,
    gamma_ideal,
    gamma_nrtl,
    gamma_nrtl_binary,
    gamma_van_laar,
    psat_kpa,
    tau_value,
)
from parameter_store import load_system_parameters


class TestAntoineEquation(unittest.TestCase):
    """Validate Antoine saturation pressure calculations."""

    def setUp(self):
        system = load_system_parameters("ethanol_water_1atm")
        self.ethanol = system.component_1.antoine
        self.water = system.component_2.antoine

    def test_ethanol_boiling_point_psat(self):
        """Ethanol Psat at 78.37 degC should be ~101.3 kPa."""
        psat = psat_kpa(78.37, self.ethanol)
        self.assertAlmostEqual(psat, 101.325, delta=3.0)

    def test_water_boiling_point_psat(self):
        """Water Psat at 100.0 degC should be ~101.3 kPa."""
        psat = psat_kpa(100.0, self.water)
        self.assertAlmostEqual(psat, 101.325, delta=2.0)

    def test_psat_monotonically_increasing(self):
        """Psat must increase with temperature for both components."""
        for comp in [self.ethanol, self.water]:
            temperatures = [50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
            pressures = [psat_kpa(t, comp) for t in temperatures]
            for i in range(1, len(pressures)):
                self.assertGreater(pressures[i], pressures[i - 1])

    def test_psat_mmhg_to_kpa_conversion(self):
        """Verify mmHg to kPa conversion factor (1 mmHg = 0.133322 kPa)."""
        params = AntoineParams(name="test", A=7.0, B=1000.0, C=200.0, unit="mmHg")
        result_mmhg = psat_kpa(50.0, params)
        # Should be 10^(7 - 1000/250) * 0.133322 = 10^3 * 0.133322
        expected = (10.0 ** (7.0 - 1000.0 / 250.0)) * 0.133322
        self.assertAlmostEqual(result_mmhg, expected, places=6)

    def test_psat_temperature_derivative_matches_finite_difference(self):
        """Antoine dPsat/dT should match central finite difference."""
        T = 78.0
        delta = 1e-4
        dpsat_analytic = dpsat_dT_kpa_per_celsius(T, self.ethanol)
        dpsat_fd = (psat_kpa(T + delta, self.ethanol) - psat_kpa(T - delta, self.ethanol)) / (2.0 * delta)
        self.assertAlmostEqual(dpsat_analytic, dpsat_fd, delta=abs(dpsat_fd) * 1e-6 + 1e-8)


class TestNRTLModel(unittest.TestCase):
    """Validate NRTL activity coefficient calculations."""

    def setUp(self):
        system = load_system_parameters("ethanol_water_1atm")
        self.nrtl = system.nrtl

    def test_hand_calculation_at_x05(self):
        """NRTL gamma at x1=0.5, T=351.15 K must match hand calculation."""
        T_K = 351.15
        x1 = 0.5
        g1, g2 = gamma_nrtl(x1, T_K, self.nrtl)

        # Hand-calculate
        x2 = 1.0 - x1
        tau12 = self.nrtl.a12 + self.nrtl.b12 / T_K
        tau21 = self.nrtl.a21 + self.nrtl.b21 / T_K
        G12 = exp(-self.nrtl.alpha12 * tau12)
        G21 = exp(-self.nrtl.alpha21 * tau21)
        D1 = x1 + x2 * G21
        D2 = x2 + x1 * G12
        ln_g1_hand = x2**2 * (tau21 * (G21 / D1) ** 2 + tau12 * G12 / D2**2)
        ln_g2_hand = x1**2 * (tau12 * (G12 / D2) ** 2 + tau21 * G21 / D1**2)

        self.assertAlmostEqual(g1, exp(ln_g1_hand), places=10)
        self.assertAlmostEqual(g2, exp(ln_g2_hand), places=10)

    def test_infinite_dilution_gamma1(self):
        """At x1->0, gamma2 should approach 1.0."""
        T_K = 351.15
        _, g2 = gamma_nrtl(0.001, T_K, self.nrtl)
        self.assertAlmostEqual(g2, 1.0, delta=0.01)

    def test_infinite_dilution_gamma2(self):
        """At x1->1, gamma1 should approach 1.0."""
        T_K = 351.15
        g1, _ = gamma_nrtl(0.999, T_K, self.nrtl)
        self.assertAlmostEqual(g1, 1.0, delta=0.01)

    def test_gamma_positive(self):
        """Activity coefficients must always be positive."""
        T_K = 351.15
        for x1 in np.linspace(0.01, 0.99, 50):
            g1, g2 = gamma_nrtl(float(x1), T_K, self.nrtl)
            self.assertGreater(g1, 0.0)
            self.assertGreater(g2, 0.0)

    def test_tau_value(self):
        """tau(T) = a + b/T."""
        self.assertAlmostEqual(tau_value(1.0, 100.0, 350.0), 1.0 + 100.0 / 350.0, places=10)

    def test_vectorized_nrtl_matches_scalar(self):
        """Vectorized NRTL implementation should match scalar outputs."""
        T_K = 351.15
        x = np.linspace(0.05, 0.95, 12)
        g1_vec, g2_vec = gamma_nrtl_binary(x, T_K, self.nrtl)

        for i, x1 in enumerate(x):
            g1, g2 = gamma_nrtl(float(x1), T_K, self.nrtl)
            self.assertAlmostEqual(float(g1_vec[i]), g1, places=12)
            self.assertAlmostEqual(float(g2_vec[i]), g2, places=12)


class TestVanLaarModel(unittest.TestCase):
    """Validate Van Laar activity coefficient calculations."""

    def setUp(self):
        self.params = VanLaarParams(A12=1.65, A21=0.95)

    def test_infinite_dilution_limits(self):
        """At x1->0, gamma1 -> exp(A12); at x1->1, gamma2 -> exp(A21)."""
        T_K = 351.15
        g1_inf, _ = gamma_van_laar(0.001, T_K, self.params)
        _, g2_inf = gamma_van_laar(0.999, T_K, self.params)
        self.assertAlmostEqual(g1_inf, exp(self.params.A12), delta=exp(self.params.A12) * 0.01)
        self.assertAlmostEqual(g2_inf, exp(self.params.A21), delta=exp(self.params.A21) * 0.01)

    def test_hand_calculation_at_x05(self):
        """Van Laar at x1=0.5 must match hand calculation."""
        x1, x2 = 0.5, 0.5
        denom = (self.params.A12 * x1 + self.params.A21 * x2) ** 2
        ln_g1 = self.params.A12 * (self.params.A21 * x2) ** 2 / denom
        ln_g2 = self.params.A21 * (self.params.A12 * x1) ** 2 / denom

        g1, g2 = gamma_van_laar(0.5, 351.15, self.params)
        self.assertAlmostEqual(g1, exp(ln_g1), places=10)
        self.assertAlmostEqual(g2, exp(ln_g2), places=10)


class TestIdealModel(unittest.TestCase):
    """Validate ideal solution model."""

    def test_always_unity(self):
        """Ideal gammas are always (1.0, 1.0)."""
        for x1 in [0.0, 0.25, 0.5, 0.75, 1.0]:
            g1, g2 = gamma_ideal(x1, 300.0, None)
            self.assertEqual(g1, 1.0)
            self.assertEqual(g2, 1.0)


if __name__ == "__main__":
    unittest.main()
