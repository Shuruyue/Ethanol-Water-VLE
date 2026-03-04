"""Bubble-point solving and T-x-y generation."""

from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np
from scipy.optimize import brentq

from models import (
    AntoineParams,
    NRTLParams,
    VanLaarParams,
    _FLOAT_FLOOR,
    dpsat_dT_kpa_per_celsius,
    gamma_with_temperature_derivatives,
    psat_kpa,
)

GammaParams = NRTLParams | VanLaarParams | None
GammaFunction = Callable[[float, float, GammaParams], tuple[float, float]]


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


def bubble_residual_and_derivative(
    T_celsius: float,
    x1: float,
    pressure_kpa: float,
    component_1: AntoineParams,
    component_2: AntoineParams,
    gamma_fn: GammaFunction,
    gamma_params: object,
) -> tuple[float, float]:
    """Return bubble residual and temperature derivative."""
    T_kelvin = T_celsius + 273.15
    gamma1, gamma2, dln_g1_dT, dln_g2_dT = gamma_with_temperature_derivatives(
        x1=x1,
        T_kelvin=T_kelvin,
        gamma_fn=gamma_fn,
        gamma_params=gamma_params,
    )
    x2 = 1.0 - x1

    psat1 = psat_kpa(T_celsius, component_1)
    psat2 = psat_kpa(T_celsius, component_2)
    dpsat1_dT = dpsat_dT_kpa_per_celsius(T_celsius, component_1)
    dpsat2_dT = dpsat_dT_kpa_per_celsius(T_celsius, component_2)

    dgamma1_dT = gamma1 * dln_g1_dT
    dgamma2_dT = gamma2 * dln_g2_dT

    p_calc = x1 * gamma1 * psat1 + x2 * gamma2 * psat2
    dp_dT = x1 * (dgamma1_dT * psat1 + gamma1 * dpsat1_dT) + x2 * (dgamma2_dT * psat2 + gamma2 * dpsat2_dT)
    return p_calc - pressure_kpa, float(dp_dT)


def _is_bracketed(f_low: float, f_high: float) -> bool:
    return f_low == 0.0 or f_high == 0.0 or (f_low < 0.0 < f_high) or (f_high < 0.0 < f_low)


def _find_bracket(
    x1: float,
    pressure_kpa: float,
    component_1: AntoineParams,
    component_2: AntoineParams,
    gamma_fn: GammaFunction,
    gamma_params: object,
    t_low: float,
    t_high: float,
    T_guess: float | None,
) -> tuple[float, float, float, float]:
    f_low = bubble_residual(t_low, x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params)
    f_high = bubble_residual(t_high, x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params)
    if _is_bracketed(f_low, f_high):
        return t_low, t_high, f_low, f_high

    center = T_guess if T_guess is not None else 0.5 * (t_low + t_high)
    span = max(5.0, 0.5 * (t_high - t_low))
    min_temp = -20.0
    max_temp = 180.0

    for _ in range(12):
        lo = max(min_temp, center - span)
        hi = min(max_temp, center + span)
        f_lo = bubble_residual(lo, x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params)
        f_hi = bubble_residual(hi, x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params)
        if _is_bracketed(f_lo, f_hi):
            return lo, hi, f_lo, f_hi
        span *= 1.5

    raise ValueError(
        f"Unable to bracket bubble-point root at x1={x1:.4f}, P={pressure_kpa:.2f} kPa "
        f"from initial range [{t_low}, {t_high}] deg C."
    )


def solve_bubble_temperature(
    x1: float,
    pressure_kpa: float,
    component_1: AntoineParams,
    component_2: AntoineParams,
    gamma_fn: GammaFunction,
    gamma_params: object,
    t_low: float = 40.0,
    t_high: float = 120.0,
    T_guess: float | None = None,
    max_newton_steps: int = 0,
) -> float:
    """Solve bubble temperature in deg C at fixed pressure and composition."""
    try:
        lo = t_low
        hi = t_high

        if max_newton_steps > 0:
            lo, hi, f_lo, _ = _find_bracket(
                x1=x1,
                pressure_kpa=pressure_kpa,
                component_1=component_1,
                component_2=component_2,
                gamma_fn=gamma_fn,
                gamma_params=gamma_params,
                t_low=t_low,
                t_high=t_high,
                T_guess=T_guess,
            )

            T_current = float(np.clip(T_guess, lo, hi)) if T_guess is not None else 0.5 * (lo + hi)
            for _ in range(max_newton_steps):
                f_curr, df_curr = bubble_residual_and_derivative(
                    T_current, x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params
                )
                if abs(f_curr) < 1e-10:
                    return float(T_current)
                if abs(df_curr) < 1e-12:
                    break

                T_next = T_current - f_curr / df_curr
                if not (lo < T_next < hi):
                    break

                f_next = bubble_residual(T_next, x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params)
                if _is_bracketed(f_lo, f_next):
                    hi = T_next
                else:
                    lo, f_lo = T_next, f_next
                T_current = float(T_next)

        return float(
            brentq(
                bubble_residual,
                lo,
                hi,
                args=(x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params),
            )
        )
    except ValueError:
        lo, hi, _, _ = _find_bracket(
            x1=x1,
            pressure_kpa=pressure_kpa,
            component_1=component_1,
            component_2=component_2,
            gamma_fn=gamma_fn,
            gamma_params=gamma_params,
            t_low=t_low,
            t_high=t_high,
            T_guess=T_guess,
        )
        try:
            return float(
                brentq(
                    bubble_residual,
                    lo,
                    hi,
                    args=(x1, pressure_kpa, component_1, component_2, gamma_fn, gamma_params),
                )
            )
        except ValueError as exc:
            raise ValueError(
                f"Bubble-point solver failed at x1={x1:.4f}, P={pressure_kpa:.2f} kPa, "
                f"search range=[{t_low}, {t_high}] deg C. "
                f"Check Antoine correlation validity. Original: {exc}"
            ) from exc


def compute_txy(
    pressure_kpa: float,
    component_1: AntoineParams,
    component_2: AntoineParams,
    gamma_fn: GammaFunction,
    gamma_params: object,
    x_grid: Iterable[float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute x, y, and T arrays for a binary system."""
    x_array = np.asarray(x_grid, dtype=float)
    t_list: list[float] = []
    y_list: list[float] = []
    previous_temperature: float | None = None

    for x1 in x_array:
        T_celsius = solve_bubble_temperature(
            x1=x1,
            pressure_kpa=pressure_kpa,
            component_1=component_1,
            component_2=component_2,
            gamma_fn=gamma_fn,
            gamma_params=gamma_params,
            t_low=40.0,
            t_high=120.0,
            T_guess=previous_temperature,
        )

        T_kelvin = T_celsius + 273.15
        gamma1, gamma2 = gamma_fn(x1, T_kelvin, gamma_params)

        y1 = x1 * gamma1 * psat_kpa(T_celsius, component_1) / pressure_kpa
        y2 = (1.0 - x1) * gamma2 * psat_kpa(T_celsius, component_2) / pressure_kpa
        y1 /= max(y1 + y2, _FLOAT_FLOOR)

        t_list.append(T_celsius)
        y_list.append(float(y1))
        previous_temperature = T_celsius

    return x_array, np.asarray(y_list), np.asarray(t_list)
