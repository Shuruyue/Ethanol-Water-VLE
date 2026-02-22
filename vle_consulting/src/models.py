"""Core thermodynamic models used by the VLE workflow."""

from dataclasses import dataclass

import numpy as np

R_GAS = 8.314462618  # J/(mol*K)


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
    """Return saturation pressure in kPa from Antoine correlation."""
    log10_psat = params.A - params.B / (T_celsius + params.C)
    psat = 10.0 ** log10_psat
    if params.unit.lower() == "mmhg":
        return psat * 0.133322
    return psat


def gamma_ideal(x1: float, T_kelvin: float, _params: object) -> tuple[float, float]:
    """Ideal-solution activity coefficients."""
    _ = x1, T_kelvin
    return 1.0, 1.0


def gamma_van_laar(x1: float, T_kelvin: float, params: VanLaarParams) -> tuple[float, float]:
    """Binary Van Laar model."""
    _ = T_kelvin
    x2 = 1.0 - x1
    denom = max((params.A12 * x1 + params.A21 * x2) ** 2, 1e-30)
    ln_gamma1 = params.A12 * (params.A21 * x2) ** 2 / denom
    ln_gamma2 = params.A21 * (params.A12 * x1) ** 2 / denom
    return float(np.exp(ln_gamma1)), float(np.exp(ln_gamma2))


def tau_value(a_term: float, b_term: float, T_kelvin: float) -> float:
    """NRTL tau(T) = a + b/T."""
    return a_term + b_term / T_kelvin


def gamma_nrtl(x1: float, T_kelvin: float, params: NRTLParams) -> tuple[float, float]:
    """Binary NRTL model."""
    x2 = 1.0 - x1

    tau12 = tau_value(params.a12, params.b12, T_kelvin)
    tau21 = tau_value(params.a21, params.b21, T_kelvin)

    G12 = np.exp(-params.alpha12 * tau12)
    G21 = np.exp(-params.alpha21 * tau21)

    D1 = max(x1 + x2 * G21, 1e-30)
    D2 = max(x2 + x1 * G12, 1e-30)

    ln_gamma1 = (x2**2) * (tau21 * (G21 / D1) ** 2 + (tau12 * G12) / (D2**2))
    ln_gamma2 = (x1**2) * (tau12 * (G12 / D2) ** 2 + (tau21 * G21) / (D1**2))

    return float(np.exp(ln_gamma1)), float(np.exp(ln_gamma2))
