"""Core thermodynamic models used by the VLE workflow."""

from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Callable
from math import exp as _exp, log as _log

import numpy as np

R_GAS = 8.314462618  # J/(mol*K)
_FLOAT_FLOOR = 1e-30  # guard for division-by-zero


@dataclass(frozen=True)
class AntoineParams:
    """Antoine equation parameters."""

    name: str
    A: float
    B: float
    C: float
    unit: str = "mmHg"


@dataclass(frozen=True)
class VanLaarParams:
    """Van Laar model parameters."""

    A12: float
    A21: float


@dataclass(frozen=True)
class NRTLParams:
    """NRTL parameter set with explicit constant and temperature terms."""

    alpha12: float
    alpha21: float
    a12: float
    b12: float
    a21: float
    b21: float


def psat_kpa(T_celsius: float, params: AntoineParams) -> float:
    r"""Return saturation pressure in kPa from Antoine correlation.

    .. math::
        \log_{10} P^{\mathrm{sat}} = A - \frac{B}{T + C}

    Parameters
    ----------
    T_celsius : float
        Temperature in degrees Celsius.
    params : AntoineParams
        Antoine equation parameters (A, B, C, unit).

    Returns
    -------
    float
        Saturation pressure in kPa.
    """
    log10_psat = params.A - params.B / (T_celsius + params.C)
    psat = 10.0 ** log10_psat
    if params.unit.lower() == "mmhg":
        return psat * 0.133322
    return psat


def dpsat_dT_kpa_per_celsius(T_celsius: float, params: AntoineParams) -> float:
    r"""Temperature derivative of Antoine saturation pressure, kPa/degC.

    For :math:`\log_{10}(P) = A - B/(T+C)`, derivative is:
        dP/dT = P * ln(10) * B / (T + C)^2
    """
    psat = psat_kpa(T_celsius, params)
    return psat * np.log(10.0) * params.B / (T_celsius + params.C) ** 2


def gamma_ideal(x1: float, T_kelvin: float, _params: object) -> tuple[float, float]:
    """Ideal-solution activity coefficients."""
    _ = x1, T_kelvin
    return 1.0, 1.0


def gamma_van_laar(x1: float, T_kelvin: float, params: VanLaarParams) -> tuple[float, float]:
    r"""Binary Van Laar activity coefficient model.

    .. math::
        \ln \gamma_1 = \frac{A_{12} (A_{21} x_2)^2}{(A_{12} x_1 + A_{21} x_2)^2}

    .. math::
        \ln \gamma_2 = \frac{A_{21} (A_{12} x_1)^2}{(A_{12} x_1 + A_{21} x_2)^2}
    """
    _ = T_kelvin
    x2 = 1.0 - x1
    denom = max((params.A12 * x1 + params.A21 * x2) ** 2, _FLOAT_FLOOR)
    ln_gamma1 = params.A12 * (params.A21 * x2) ** 2 / denom
    ln_gamma2 = params.A21 * (params.A12 * x1) ** 2 / denom
    return _exp(ln_gamma1), _exp(ln_gamma2)


def tau_value(a_term: float, b_term: float, T_kelvin: float) -> float:
    """NRTL tau(T) = a + b/T."""
    return a_term + b_term / T_kelvin


def gamma_nrtl_binary(
    x1: np.ndarray | float,
    T_kelvin: np.ndarray | float,
    params: NRTLParams,
    return_dln_dT: bool = False,
) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized binary NRTL model.

    Supports scalar or array ``x1`` and ``T_kelvin`` with NumPy broadcasting.
    """
    x1_arr = np.asarray(x1, dtype=float)
    T_arr = np.asarray(T_kelvin, dtype=float)
    x2_arr = 1.0 - x1_arr

    tau12 = params.a12 + params.b12 / T_arr
    tau21 = params.a21 + params.b21 / T_arr

    G12 = np.exp(-params.alpha12 * tau12)
    G21 = np.exp(-params.alpha21 * tau21)

    D1 = np.maximum(x1_arr + x2_arr * G21, _FLOAT_FLOOR)
    D2 = np.maximum(x2_arr + x1_arr * G12, _FLOAT_FLOOR)

    ln_gamma1 = (x2_arr**2) * (tau21 * (G21 / D1) ** 2 + (tau12 * G12) / (D2**2))
    ln_gamma2 = (x1_arr**2) * (tau12 * (G12 / D2) ** 2 + (tau21 * G21) / (D1**2))

    gamma1 = np.exp(ln_gamma1)
    gamma2 = np.exp(ln_gamma2)

    if not return_dln_dT:
        return gamma1, gamma2

    T2 = T_arr * T_arr
    dtau12_dT = -params.b12 / T2
    dtau21_dT = -params.b21 / T2

    dG12_dT = -params.alpha12 * dtau12_dT * G12
    dG21_dT = -params.alpha21 * dtau21_dT * G21

    dD1_dT = x2_arr * dG21_dT
    dD2_dT = x1_arr * dG12_dT

    A = tau21 * G21**2 / D1**2
    B = tau12 * G12 / D2**2
    C = tau12 * G12**2 / D2**2
    D_val = tau21 * G21 / D1**2

    dA_dT = (dtau21_dT * G21**2 + 2.0 * tau21 * G21 * dG21_dT) / D1**2 \
            - 2.0 * tau21 * G21**2 * dD1_dT / D1**3
    dB_dT = (dtau12_dT * G12 + tau12 * dG12_dT) / D2**2 \
            - 2.0 * tau12 * G12 * dD2_dT / D2**3
    dC_dT = (dtau12_dT * G12**2 + 2.0 * tau12 * G12 * dG12_dT) / D2**2 \
            - 2.0 * tau12 * G12**2 * dD2_dT / D2**3
    dD_dT = (dtau21_dT * G21 + tau21 * dG21_dT) / D1**2 \
            - 2.0 * tau21 * G21 * dD1_dT / D1**3

    dln_g1_dT = x2_arr**2 * (dA_dT + dB_dT)
    dln_g2_dT = x1_arr**2 * (dC_dT + dD_dT)
    return gamma1, gamma2, dln_g1_dT, dln_g2_dT


def gamma_nrtl(x1: float, T_kelvin: float, params: NRTLParams) -> tuple[float, float]:
    r"""Binary NRTL activity coefficient model.

    .. math::
        \ln \gamma_1 = x_2^2 \left[ \tau_{21} \left(\frac{G_{21}}{D_1}\right)^2
        + \frac{\tau_{12} G_{12}}{D_2^2} \right]

    where :math:`\tau_{ij} = a_{ij} + b_{ij}/T`,
    :math:`G_{ij} = \exp(-\alpha_{ij} \tau_{ij})`,
    :math:`D_1 = x_1 + x_2 G_{21}`, :math:`D_2 = x_2 + x_1 G_{12}`.
    """
    x2 = 1.0 - x1

    tau12 = tau_value(params.a12, params.b12, T_kelvin)
    tau21 = tau_value(params.a21, params.b21, T_kelvin)

    G12 = _exp(-params.alpha12 * tau12)
    G21 = _exp(-params.alpha21 * tau21)

    D1 = max(x1 + x2 * G21, _FLOAT_FLOOR)
    D2 = max(x2 + x1 * G12, _FLOAT_FLOOR)

    ln_gamma1 = (x2**2) * (tau21 * (G21 / D1) ** 2 + (tau12 * G12) / (D2**2))
    ln_gamma2 = (x1**2) * (tau12 * (G12 / D2) ** 2 + (tau21 * G21) / (D1**2))
    return _exp(ln_gamma1), _exp(ln_gamma2)


def gamma_nrtl_with_derivatives(
    x1: float, T_kelvin: float, params: NRTLParams,
) -> tuple[float, float, float, float]:
    r"""Binary NRTL model with analytical temperature derivatives.

    Returns (gamma1, gamma2, d_ln_gamma1_dT, d_ln_gamma2_dT).

    Analytical derivatives use:
        d(tau)/dT = -b / T^2
        d(G)/dT  = -alpha * d(tau)/dT * G
    """
    gamma1, gamma2, dln_g1_dT, dln_g2_dT = gamma_nrtl_binary(x1, T_kelvin, params, return_dln_dT=True)
    return float(gamma1), float(gamma2), float(dln_g1_dT), float(dln_g2_dT)


GammaFunction = Callable[[float, float, object], tuple[float, float]]


def gamma_with_temperature_derivatives(
    x1: float,
    T_kelvin: float,
    gamma_fn: GammaFunction,
    gamma_params: object,
    delta_t: float = 1e-3,
) -> tuple[float, float, float, float]:
    """Return gamma and d(ln gamma)/dT for any supported activity model."""
    if gamma_fn is gamma_nrtl:
        return gamma_nrtl_with_derivatives(x1, T_kelvin, gamma_params)

    gamma1, gamma2 = gamma_fn(x1, T_kelvin, gamma_params)

    # Ideal and Van Laar are T-independent in this project parameterization.
    if gamma_fn is gamma_ideal or gamma_fn is gamma_van_laar:
        return gamma1, gamma2, 0.0, 0.0

    g1_plus, g2_plus = gamma_fn(x1, T_kelvin + delta_t, gamma_params)
    g1_minus, g2_minus = gamma_fn(x1, T_kelvin - delta_t, gamma_params)
    dln_g1_dT = (_log(max(g1_plus, _FLOAT_FLOOR)) - _log(max(g1_minus, _FLOAT_FLOOR))) / (2.0 * delta_t)
    dln_g2_dT = (_log(max(g2_plus, _FLOAT_FLOOR)) - _log(max(g2_minus, _FLOAT_FLOOR))) / (2.0 * delta_t)
    return gamma1, gamma2, float(dln_g1_dT), float(dln_g2_dT)
