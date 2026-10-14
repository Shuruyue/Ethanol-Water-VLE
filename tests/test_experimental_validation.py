"""Validate the NRTL parameter set against measured isobaric VLE data.

Data: Kamihama et al., J. Chem. Eng. Data 57 (2012) 339-344, 101.3 kPa,
21 points, with TRC expanded uncertainties.
"""

import numpy as np
import pytest

from experimental_data import load_experimental_txy
from models import gamma_ideal, gamma_nrtl
from parameter_store import load_system_parameters
from solver import compute_txy


@pytest.fixture(scope="module")
def data():
    return load_experimental_txy()


@pytest.fixture(scope="module")
def system():
    return load_system_parameters()


def _predict(system, data, gamma_fn, params):
    _, y, T = compute_txy(
        data.pressure_kpa, system.component_1.antoine, system.component_2.antoine, gamma_fn, params, data.x1
    )
    return y, T


def test_data_file_is_complete(data):
    assert data.x1.size == 21
    assert data.pressure_kpa == pytest.approx(101.3)
    assert np.all(np.diff(data.x1) > 0)


def test_nrtl_reproduces_data_within_uncertainty(system, data):
    y, T = _predict(system, data, gamma_nrtl, system.nrtl)
    assert np.mean(np.abs(T - data.T_celsius)) < 0.15
    assert np.mean(np.abs(y - data.y1)) < 0.005
    assert np.all(np.abs(T - data.T_celsius) <= data.U_T)
    assert np.all(np.abs(y - data.y1) <= data.U_y1)


def test_ideal_solution_is_rejected_by_data(system, data):
    y, T = _predict(system, data, gamma_ideal, None)
    assert np.mean(np.abs(T - data.T_celsius)) > 3.0
    # the largest error is on the dilute-ethanol side, not mid-range
    assert data.x1[np.argmax(np.abs(T - data.T_celsius))] < 0.3
