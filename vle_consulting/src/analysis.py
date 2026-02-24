"""Post-processing for excess properties and azeotrope detection."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from models import NRTLParams, R_GAS, gamma_nrtl, gamma_nrtl_with_derivatives


@dataclass(frozen=True)
class AzeotropePoint:
    """Azeotrope descriptor."""

    x1: float
    y1: float
    temperature_c: float
    abs_error: float


def excess_gibbs_energy(x1: float, T_kelvin: float, nrtl_params: NRTLParams) -> float:
    """Molar excess Gibbs energy, J/mol."""
    gamma1, gamma2 = gamma_nrtl(x1, T_kelvin, nrtl_params)
    return R_GAS * T_kelvin * (x1 * np.log(gamma1) + (1.0 - x1) * np.log(gamma2))


def excess_enthalpy(x1: float, T_kelvin: float, nrtl_params: NRTLParams, delta_t: float = 1e-3) -> float:
    r"""Molar excess enthalpy, J/mol.

    Uses analytical derivatives: HE = -RT^2 * sum(xi * d(ln gamma_i)/dT).
    The ``delta_t`` parameter is kept for backward compatibility but is ignored
    when analytical derivatives are available.
    """
    _ = delta_t  # kept for API compatibility
    g1, g2, dlng1_dT, dlng2_dT = gamma_nrtl_with_derivatives(x1, T_kelvin, nrtl_params)
    return -R_GAS * (T_kelvin**2) * (x1 * dlng1_dT + (1.0 - x1) * dlng2_dT)


def excess_heat_capacity(
    x1: float,
    T_kelvin: float,
    nrtl_params: NRTLParams,
    delta_t: float = 0.01,
) -> float:
    """Molar excess heat capacity, J/(mol*K), as d(HE)/dT.

    Uses numerical derivative of the analytical HE for second-order accuracy.
    """
    h_plus = excess_enthalpy(x1, T_kelvin + delta_t, nrtl_params)
    h_minus = excess_enthalpy(x1, T_kelvin - delta_t, nrtl_params)
    return (h_plus - h_minus) / (2.0 * delta_t)


def compute_excess_curves(x: np.ndarray, T_celsius: np.ndarray, nrtl_params: NRTLParams) -> tuple[np.ndarray, np.ndarray]:
    """Return GE and HE arrays along a bubble line."""
    T_kelvin = T_celsius + 273.15
    GE = np.array([excess_gibbs_energy(x1, T, nrtl_params) for x1, T in zip(x, T_kelvin, strict=False)])
    HE = np.array([excess_enthalpy(x1, T, nrtl_params) for x1, T in zip(x, T_kelvin, strict=False)])
    return GE, HE


def compute_isothermal_excess_curves(
    x: np.ndarray,
    T_celsius: float,
    nrtl_params: NRTLParams,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return GE, HE, and CPE at one temperature."""
    T_kelvin = T_celsius + 273.15
    GE = np.array([excess_gibbs_energy(x1, T_kelvin, nrtl_params) for x1 in x])
    HE = np.array([excess_enthalpy(x1, T_kelvin, nrtl_params) for x1 in x])
    CPE = np.array([excess_heat_capacity(x1, T_kelvin, nrtl_params) for x1 in x])
    return GE, HE, CPE


def compute_gamma_curves(x: np.ndarray, T_celsius: float, nrtl_params: NRTLParams) -> tuple[np.ndarray, np.ndarray]:
    """Return gamma1 and gamma2 arrays at one temperature."""
    T_kelvin = T_celsius + 273.15
    gamma1: list[float] = []
    gamma2: list[float] = []
    for x1 in x:
        g1, g2 = gamma_nrtl(float(x1), T_kelvin, nrtl_params)
        gamma1.append(g1)
        gamma2.append(g2)
    return np.asarray(gamma1), np.asarray(gamma2)


def find_azeotrope(x: np.ndarray, y: np.ndarray, T_celsius: np.ndarray, tolerance: float = 0.01) -> AzeotropePoint | None:
    """Return azeotrope point if |x-y| is within tolerance."""
    idx = int(np.argmin(np.abs(y - x)))
    error = float(abs(y[idx] - x[idx]))
    if error > tolerance:
        return None
    return AzeotropePoint(
        x1=float(x[idx]),
        y1=float(y[idx]),
        temperature_c=float(T_celsius[idx]),
        abs_error=error,
    )
