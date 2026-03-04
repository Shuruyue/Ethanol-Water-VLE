"""Reusable profile builders for VLE post-processing workflows."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from analysis import (
    compute_gamma_curves,
    compute_isothermal_excess_curves,
    find_azeotrope,
)
from models import gamma_nrtl, psat_kpa
from parameter_store import SystemParameters
from solver import compute_txy


@dataclass(frozen=True)
class IsothermalProfile:
    """Isothermal composition profile used by UI and exports."""

    x: np.ndarray
    GE_j_mol: np.ndarray
    HE_j_mol: np.ndarray
    CPE_j_molK: np.ndarray
    gamma1: np.ndarray
    gamma2: np.ndarray
    alpha12: np.ndarray
    temperature_c: float


@dataclass(frozen=True)
class PressureSensitivityProfile:
    """Azeotrope movement across an isobaric pressure sweep."""

    pressure_kpa: np.ndarray
    azeotrope_x1: np.ndarray
    azeotrope_y1: np.ndarray
    azeotrope_temperature_c: np.ndarray


def build_isothermal_profile(
    system: SystemParameters,
    temperature_c: float,
    points: int,
) -> IsothermalProfile:
    """Build GE/HE/CPE/gamma/alpha profile at one temperature."""
    x = np.linspace(0.01, 0.99, points)

    GE, HE, CPE = compute_isothermal_excess_curves(x, temperature_c, system.nrtl)
    gamma1, gamma2 = compute_gamma_curves(x, temperature_c, system.nrtl)

    psat1 = psat_kpa(temperature_c, system.component_1.antoine)
    psat2 = psat_kpa(temperature_c, system.component_2.antoine)
    alpha12 = (gamma1 * psat1) / np.maximum(gamma2 * psat2, 1e-30)

    return IsothermalProfile(
        x=x,
        GE_j_mol=GE,
        HE_j_mol=HE,
        CPE_j_molK=CPE,
        gamma1=gamma1,
        gamma2=gamma2,
        alpha12=alpha12,
        temperature_c=float(temperature_c),
    )


def build_pressure_sensitivity_profile(
    system: SystemParameters,
    pressure_values_kpa: np.ndarray,
    points: int,
    azeotrope_tolerance: float = 0.01,
) -> PressureSensitivityProfile:
    """Build pressure-sweep azeotrope profile."""
    x_grid = np.linspace(0.01, 0.99, points)

    azeotrope_x: list[float] = []
    azeotrope_y: list[float] = []
    azeotrope_t: list[float] = []

    for pressure in pressure_values_kpa:
        x_n, y_n, T_n = compute_txy(
            pressure_kpa=float(pressure),
            component_1=system.component_1.antoine,
            component_2=system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=system.nrtl,
            x_grid=x_grid,
        )

        azeo = find_azeotrope(x_n, y_n, T_n, tolerance=azeotrope_tolerance)
        if azeo is None:
            azeotrope_x.append(np.nan)
            azeotrope_y.append(np.nan)
            azeotrope_t.append(np.nan)
        else:
            azeotrope_x.append(azeo.x1)
            azeotrope_y.append(azeo.y1)
            azeotrope_t.append(azeo.temperature_c)

    return PressureSensitivityProfile(
        pressure_kpa=np.asarray(pressure_values_kpa, dtype=float),
        azeotrope_x1=np.asarray(azeotrope_x, dtype=float),
        azeotrope_y1=np.asarray(azeotrope_y, dtype=float),
        azeotrope_temperature_c=np.asarray(azeotrope_t, dtype=float),
    )
