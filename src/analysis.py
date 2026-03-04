"""Post-processing for excess properties and azeotrope detection."""

from __future__ import annotations

from dataclasses import dataclass
from math import log as _log

import numpy as np

from models import NRTLParams, R_GAS, gamma_nrtl, gamma_nrtl_binary, gamma_nrtl_with_derivatives


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
    return R_GAS * T_kelvin * (x1 * _log(gamma1) + (1.0 - x1) * _log(gamma2))


def excess_enthalpy(x1: float, T_kelvin: float, nrtl_params: NRTLParams, delta_t: float = 1e-3) -> float:
    r"""Molar excess enthalpy, J/mol.

    Uses analytical derivatives: HE = -RT^2 * sum(xi * d(ln gamma_i)/dT).
    The ``delta_t`` parameter is kept for backward compatibility but is ignored
    when analytical derivatives are available.
    """
    _ = delta_t  # kept for API compatibility
    g1, g2, dlng1_dT, dlng2_dT = gamma_nrtl_with_derivatives(x1, T_kelvin, nrtl_params)
    return -R_GAS * (T_kelvin**2) * (x1 * dlng1_dT + (1.0 - x1) * dlng2_dT)


def _excess_enthalpy_array(x: np.ndarray, T_kelvin: np.ndarray | float, nrtl_params: NRTLParams) -> np.ndarray:
    """Vectorized HE helper."""
    _, _, dlng1_dT, dlng2_dT = gamma_nrtl_binary(x, T_kelvin, nrtl_params, return_dln_dT=True)
    return -R_GAS * (np.asarray(T_kelvin, dtype=float) ** 2) * (x * dlng1_dT + (1.0 - x) * dlng2_dT)


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
    x_arr = np.asarray(x, dtype=float)
    T_kelvin = np.asarray(T_celsius, dtype=float) + 273.15
    gamma1, gamma2 = gamma_nrtl_binary(x_arr, T_kelvin, nrtl_params, return_dln_dT=False)
    GE = R_GAS * T_kelvin * (x_arr * np.log(gamma1) + (1.0 - x_arr) * np.log(gamma2))
    HE = _excess_enthalpy_array(x_arr, T_kelvin, nrtl_params)
    return GE, HE


def compute_isothermal_excess_curves(
    x: np.ndarray,
    T_celsius: float,
    nrtl_params: NRTLParams,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return GE, HE, and CPE at one temperature."""
    x_arr = np.asarray(x, dtype=float)
    T_kelvin = T_celsius + 273.15
    gamma1, gamma2 = gamma_nrtl_binary(x_arr, T_kelvin, nrtl_params, return_dln_dT=False)
    GE = R_GAS * T_kelvin * (x_arr * np.log(gamma1) + (1.0 - x_arr) * np.log(gamma2))
    HE = _excess_enthalpy_array(x_arr, T_kelvin, nrtl_params)
    delta_t = 0.01
    HE_plus = _excess_enthalpy_array(x_arr, T_kelvin + delta_t, nrtl_params)
    HE_minus = _excess_enthalpy_array(x_arr, T_kelvin - delta_t, nrtl_params)
    CPE = (HE_plus - HE_minus) / (2.0 * delta_t)
    return GE, HE, CPE


def compute_gamma_curves(x: np.ndarray, T_celsius: float, nrtl_params: NRTLParams) -> tuple[np.ndarray, np.ndarray]:
    """Return gamma1 and gamma2 arrays at one temperature."""
    T_kelvin = T_celsius + 273.15
    x_arr = np.asarray(x, dtype=float)
    gamma1, gamma2 = gamma_nrtl_binary(x_arr, T_kelvin, nrtl_params, return_dln_dT=False)
    return np.asarray(gamma1), np.asarray(gamma2)


def find_azeotrope(x: np.ndarray, y: np.ndarray, T_celsius: np.ndarray, tolerance: float = 0.01) -> AzeotropePoint | None:
    """Return azeotrope point if |x-y| is within tolerance."""
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    T_arr = np.asarray(T_celsius, dtype=float)
    diff = y_arr - x_arr

    crossing_idx = np.where(diff[:-1] * diff[1:] <= 0.0)[0]
    if crossing_idx.size > 0:
        idx = int(crossing_idx[np.argmin(np.abs(diff[crossing_idx]))])
        x0, x1 = x_arr[idx], x_arr[idx + 1]
        y0, y1 = y_arr[idx], y_arr[idx + 1]
        T0, T1 = T_arr[idx], T_arr[idx + 1]
        d0, d1 = diff[idx], diff[idx + 1]

        if abs(d1 - d0) > 1e-15:
            frac = -d0 / (d1 - d0)
            frac = float(np.clip(frac, 0.0, 1.0))
        else:
            frac = 0.0

        x_star = x0 + frac * (x1 - x0)
        y_star = y0 + frac * (y1 - y0)
        T_star = T0 + frac * (T1 - T0)
        error = abs(y_star - x_star)
        if error <= tolerance:
            return AzeotropePoint(
                x1=float(x_star),
                y1=float(y_star),
                temperature_c=float(T_star),
                abs_error=float(error),
            )

    idx = int(np.argmin(np.abs(diff)))
    error = float(abs(diff[idx]))
    if error > tolerance:
        return None
    return AzeotropePoint(
        x1=float(x_arr[idx]),
        y1=float(y_arr[idx]),
        temperature_c=float(T_arr[idx]),
        abs_error=error,
    )
