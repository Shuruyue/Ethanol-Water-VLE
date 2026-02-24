"""Tests for bubble-point solver and T-x-y curve generation."""

import unittest

import numpy as np

from models import gamma_ideal, gamma_nrtl
from parameter_store import load_system_parameters
from solver import compute_txy, solve_bubble_temperature


class TestBubblePoint(unittest.TestCase):
    """Validate bubble-point temperature solver."""

    def setUp(self):
        self.system = load_system_parameters("ethanol_water_1atm")
        self.eth = self.system.component_1.antoine
        self.wat = self.system.component_2.antoine

    def test_pure_ethanol_boiling_point(self):
        """x1=1 (pure ethanol) should boil at ~78.37 degC."""
        T = solve_bubble_temperature(
            x1=1.0, pressure_kpa=101.325,
            component_1=self.eth, component_2=self.wat,
            gamma_fn=gamma_ideal, gamma_params=None,
        )
        self.assertAlmostEqual(T, 78.37, delta=0.5)

    def test_pure_water_boiling_point(self):
        """x1=0 (pure water) should boil at ~100.0 degC."""
        T = solve_bubble_temperature(
            x1=0.0, pressure_kpa=101.325,
            component_1=self.eth, component_2=self.wat,
            gamma_fn=gamma_ideal, gamma_params=None,
        )
        self.assertAlmostEqual(T, 100.0, delta=0.5)

    def test_nrtl_mixture_bubble_point(self):
        """T at x1=0.5, NRTL should be between pure boiling points (or at azeotrope)."""
        T = solve_bubble_temperature(
            x1=0.5, pressure_kpa=101.325,
            component_1=self.eth, component_2=self.wat,
            gamma_fn=gamma_nrtl, gamma_params=self.system.nrtl,
        )
        self.assertGreater(T, 70.0)
        self.assertLess(T, 105.0)


class TestTxyGeneration(unittest.TestCase):
    """Validate T-x-y curve generation."""

    def setUp(self):
        self.system = load_system_parameters("ethanol_water_1atm")

    def test_txy_array_lengths(self):
        """Output arrays must have same length as input grid."""
        x_grid = np.linspace(0.01, 0.99, 40)
        x, y, T = compute_txy(
            pressure_kpa=101.325,
            component_1=self.system.component_1.antoine,
            component_2=self.system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=self.system.nrtl,
            x_grid=x_grid,
        )
        self.assertEqual(len(x), 40)
        self.assertEqual(len(y), 40)
        self.assertEqual(len(T), 40)

    def test_y_bounded_zero_one(self):
        """y1 values must be between 0 and 1."""
        x_grid = np.linspace(0.01, 0.99, 40)
        _, y, _ = compute_txy(
            pressure_kpa=101.325,
            component_1=self.system.component_1.antoine,
            component_2=self.system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=self.system.nrtl,
            x_grid=x_grid,
        )
        self.assertTrue(np.all(y >= 0.0))
        self.assertTrue(np.all(y <= 1.0))

    def test_y_greater_than_x_at_low_x(self):
        """For positive deviation, y1 > x1 at low x1 (ethanol more volatile)."""
        x_grid = np.linspace(0.05, 0.3, 10)
        x, y, _ = compute_txy(
            pressure_kpa=101.325,
            component_1=self.system.component_1.antoine,
            component_2=self.system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=self.system.nrtl,
            x_grid=x_grid,
        )
        # At low ethanol concentrations, y1 should exceed x1
        self.assertTrue(np.all(y > x), "y1 should be > x1 at low ethanol fractions")

    def test_temperature_monotonically_decreasing_ideal(self):
        """For ideal solution, T should decrease as ethanol fraction increases."""
        x_grid = np.linspace(0.01, 0.99, 50)
        _, _, T = compute_txy(
            pressure_kpa=101.325,
            component_1=self.system.component_1.antoine,
            component_2=self.system.component_2.antoine,
            gamma_fn=gamma_ideal,
            gamma_params=None,
            x_grid=x_grid,
        )
        # Temperature should generally decrease as we add more ethanol
        self.assertGreater(T[0], T[-1])


if __name__ == "__main__":
    unittest.main()
