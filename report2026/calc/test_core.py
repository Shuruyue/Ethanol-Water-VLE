"""Consistency checks for vle_core (run with: pytest -q calc/test_core.py)."""
import numpy as np
import pytest

from vle_core import NRTL_LIT, R, Virial, azeotrope, bubble_T, tsat

M = NRTL_LIT


@pytest.mark.parametrize("x", [0.1, 0.4, 0.8])
@pytest.mark.parametrize("T", [298.15, 351.35])
def test_analytic_hE_matches_finite_difference(x, T):
    h = 1e-3
    num = -R * T**2 * (M.gE_RT(x, T + h) - M.gE_RT(x, T - h)) / (2 * h)
    assert abs(num - M.hE(x, T)) < 1e-4 * max(1.0, abs(num))


@pytest.mark.parametrize("x", [0.05, 0.5, 0.95])
def test_gE_equals_sum_x_lngamma(x):
    l1, l2 = M.lngamma(x, 350.0)
    assert abs(x * l1 + (1 - x) * l2 - M.gE_RT(x, 350.0)) < 1e-12


def test_normal_boiling_points():
    assert tsat(0, 101.325) == pytest.approx(351.45, abs=0.05)  # ethanol, 78.30 degC
    assert tsat(1, 101.325) == pytest.approx(373.124, abs=0.01)  # water (IAPWS-95)


def test_pure_component_limits():
    T0, y0 = bubble_T(0.0, 101.325, M)
    T1, y1 = bubble_T(1.0, 101.325, M)
    assert y0 == pytest.approx(0.0, abs=1e-12) and y1 == pytest.approx(1.0, abs=1e-12)
    assert T0 == pytest.approx(tsat(1, 101.325), abs=1e-6)
    assert T1 == pytest.approx(tsat(0, 101.325), abs=1e-6)


def test_virial_correction_is_small_at_1atm():
    T, y = bubble_T(0.5, 101.325, M)
    T2, y2 = bubble_T(0.5, 101.325, M, Virial())
    assert abs(T2 - T) < 0.1 and abs(y2 - y) < 0.005


def test_azeotrope_matches_experiment():
    xa, Ta = azeotrope(M, 101.325)
    assert xa == pytest.approx(0.8946, abs=0.01)  # Kamihama et al. (2012), interpolated
    assert Ta - 273.15 == pytest.approx(78.11, abs=0.1)
    assert np.isfinite(Ta)
