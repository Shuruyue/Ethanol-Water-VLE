"""Tests for plotting helpers with algorithmic content."""

import unittest

import numpy as np

from models import gamma_nrtl
from parameter_store import load_system_parameters
from plotting import _mccabe_thiele_steps
from solver import compute_txy


class TestMcCabeThieleStepping(unittest.TestCase):
    """Validate McCabe-Thiele stepping direction and feasibility."""

    def test_steps_move_left_toward_bottoms(self):
        system = load_system_parameters("ethanol_water_1atm")
        x, y, _ = compute_txy(
            pressure_kpa=system.pressure_kpa,
            component_1=system.component_1.antoine,
            component_2=system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=system.nrtl,
            x_grid=np.linspace(0.01, 0.99, 120),
        )

        stages_x, stages_y, n_stages, *_ = _mccabe_thiele_steps(
            x_nrtl=x,
            y_nrtl=y,
            xD=0.85,
            xB=0.05,
            xF=0.30,
            q=1.0,
            R_reflux=2.0,
            max_stages=50,
        )

        self.assertGreater(n_stages, 0)

        # Horizontal segments are every 4 points: (i+1 -> i+2). They should move left.
        for i in range(0, len(stages_x) - 3, 4):
            x_start = stages_x[i + 1]
            x_end = stages_x[i + 2]
            self.assertLessEqual(x_end, x_start + 1e-9)
            self.assertAlmostEqual(stages_y[i + 1], stages_y[i + 2], places=12)

        # Overall should progress from distillate toward bottoms.
        self.assertLess(stages_x[-1], 0.85)
        self.assertLess(stages_x[-1], 0.35)


if __name__ == "__main__":
    unittest.main()
