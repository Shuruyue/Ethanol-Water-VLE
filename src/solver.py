"""Bubble-point solving and T-x-y generation."""

from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np
from scipy.optimize import brentq

from models import AntoineParams, psat_kpa

GammaFunction = Callable[[float, float, object], tuple[float, float]]


def bubble_residual(
    T_celsius: float,
    x1: float,
    pressure_kpa: float,
    component_1: AntoineParams,
    component_2: AntoineParams,
    gamma_fn: GammaFunction,
    gamma_params: object,
) -> float:
    """Residual for the bubble-point equation."""
    T_kelvin = T_celsius + 273.15
    gamma1, gamma2 = gamma_fn(x1, T_kelvin, gamma_params)
    p_calc = (
        x1 * gamma1 * psat_kpa(T_celsius, component_1)
        + (1.0 - x1) * gamma2 * psat_kpa(T_celsius, component_2)
    )
    return p_calc - pressure_kpa


def solve_bubble_temperature(
    x1: float,
    pressure_kpa: float,
    component_1: AntoineParams,
    component_2: AntoineParams,
    gamma_fn: GammaFunction,
    gamma_params: object,
    t_low: float = 40.0,
    t_high: float = 120.0,
) -> float:
    """Solve bubble temperature in deg C at fixed pressure and composition."""
    return float(
        brentq(
            bubble_residual,
            t_low,
            t_high,
            args=(x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params),
        )
    )


def compute_txy(
    pressure_kpa: float,
    component_1: AntoineParams,
    component_2: AntoineParams,
    gamma_fn: GammaFunction,
    gamma_params: object,
    x_grid: Iterable[float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute x, y, and T arrays for a binary system."""
    x_array = np.asarray(list(x_grid), dtype=float)
    t_list: list[float] = []
    y_list: list[float] = []

    for x1 in x_array:
        T_celsius = solve_bubble_temperature(
            x1=x1,
            pressure_kpa=pressure_kpa,
            component_1=component_1,
            component_2=component_2,
            gamma_fn=gamma_fn,
            gamma_params=gamma_params,
        )
        T_kelvin = T_celsius + 273.15
        gamma1, gamma2 = gamma_fn(x1, T_kelvin, gamma_params)

        y1 = x1 * gamma1 * psat_kpa(T_celsius, component_1) / pressure_kpa
        y2 = (1.0 - x1) * gamma2 * psat_kpa(T_celsius, component_2) / pressure_kpa
        y1 /= max(y1 + y2, 1e-30)

        t_list.append(T_celsius)
        y_list.append(float(y1))

    return x_array, np.asarray(y_list), np.asarray(t_list)
