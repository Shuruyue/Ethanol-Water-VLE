"""Core thermodynamics for the ethanol(1)-water(2) consulting report.

Gamma-phi VLE with NIST Antoine correlations, NRTL / Van Laar liquid models,
Tsonopoulos second-virial vapour correction, excess functions and simple
McCabe-Thiele (constant-molar-overflow) distillation tools.

Units: T in K, P in kPa, energies in J/mol unless stated otherwise.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

R = 8.314462618  # J mol^-1 K^-1
MW = (46.068, 18.015)  # g/mol, ethanol, water

# ---------------------------------------------------------------------------
# Pure-component vapour pressure
#   ethanol : NIST WebBook Antoine set, Ambrose & Sprake (1970),
#             log10(P/bar) = A - B/(T/K + C), fitted range 292.77-366.63 K
#             (checked against the Schroeder et al. 2014 reference EOS via CoolProp:
#              |dev| < 0.8 % for 293-373 K, ~1.2 % at 393 K; see run_analysis.py)
#   water   : IAPWS-95 saturation pressure (Wagner & Pruss 2002)
# ---------------------------------------------------------------------------
ANTOINE_ETOH = (5.24677, 1598.673, -46.424)  # bar, K


def psat(i: int, T: float) -> float:
    """Vapour pressure of component i (0 = ethanol, 1 = water) in kPa."""
    if i == 0:
        A, B, C = ANTOINE_ETOH
        return 100.0 * 10.0 ** (A - B / (T + C))
    from chemicals.iapws import iapws95_Psat
    return iapws95_Psat(T) / 1000.0


def tsat(i: int, P: float) -> float:
    return brentq(lambda T: psat(i, T) - P, 250.0, 500.0)


# Liquid molar volumes (cm3/mol) for the Poynting factor, near 350 K
VL = (61.0, 18.8)
# Enthalpy of vaporisation at the normal boiling point (NIST WebBook), J/mol
DHVAP_NBP = (38560.0, 40650.0)

# ---------------------------------------------------------------------------
# Activity-coefficient models
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class NRTL:
    """Binary NRTL: tau_ij = a_ij + b_ij/T, G_ij = exp(-alpha tau_ij)."""

    a12: float
    b12: float
    a21: float
    b21: float
    alpha: float = 0.3
    name: str = "NRTL"

    def taus(self, T):
        return self.a12 + self.b12 / T, self.a21 + self.b21 / T

    def lngamma(self, x1, T):
        x1 = np.asarray(x1, dtype=float)
        x2 = 1.0 - x1
        t12, t21 = self.taus(T)
        G12, G21 = np.exp(-self.alpha * t12), np.exp(-self.alpha * t21)
        D1, D2 = x1 + x2 * G21, x2 + x1 * G12
        l1 = x2**2 * (t21 * (G21 / D1) ** 2 + t12 * G12 / D2**2)
        l2 = x1**2 * (t12 * (G12 / D2) ** 2 + t21 * G21 / D1**2)
        return l1, l2

    def gE_RT(self, x1, T):
        x1 = np.asarray(x1, dtype=float)
        x2 = 1.0 - x1
        t12, t21 = self.taus(T)
        G12, G21 = np.exp(-self.alpha * t12), np.exp(-self.alpha * t21)
        return x1 * x2 * (t21 * G21 / (x1 + x2 * G21) + t12 * G12 / (x2 + x1 * G12))

    def hE(self, x1, T):
        """Analytical Gibbs-Helmholtz: H^E = -R T^2 d(g^E/RT)/dT."""
        x1 = np.asarray(x1, dtype=float)
        x2 = 1.0 - x1
        t12, t21 = self.taus(T)
        dt12, dt21 = -self.b12 / T**2, -self.b21 / T**2
        G12, G21 = np.exp(-self.alpha * t12), np.exp(-self.alpha * t21)
        dG12, dG21 = -self.alpha * dt12 * G12, -self.alpha * dt21 * G21

        def term(t, dt, G, dG, num_x, den_other_x):
            # f = t G / (num_x + den_other_x G)
            D = num_x + den_other_x * G
            return (dt * G + t * dG) / D - t * G * den_other_x * dG / D**2

        dg = x1 * x2 * (term(t21, dt21, G21, dG21, x1, x2) + term(t12, dt12, G12, dG12, x2, x1))
        return -R * T**2 * dg

    def cpE(self, x1, T, h=0.01):
        return (self.hE(x1, T + h) - self.hE(x1, T - h)) / (2 * h)


@dataclass(frozen=True)
class VanLaar:
    """ln g1 = A12 (A21 x2)^2/(A12 x1 + A21 x2)^2 (T-independent)."""

    A12: float
    A21: float
    name: str = "Van Laar"

    def lngamma(self, x1, T):
        x1 = np.asarray(x1, dtype=float)
        x2 = 1.0 - x1
        D = (self.A12 * x1 + self.A21 * x2) ** 2
        return self.A12 * (self.A21 * x2) ** 2 / D, self.A21 * (self.A12 * x1) ** 2 / D

    def gE_RT(self, x1, T):
        x1 = np.asarray(x1, dtype=float)
        x2 = 1.0 - x1
        return self.A12 * self.A21 * x1 * x2 / (self.A12 * x1 + self.A21 * x2)


@dataclass(frozen=True)
class Ideal:
    name: str = "Ideal (Raoult)"

    def lngamma(self, x1, T):
        z = np.zeros_like(np.asarray(x1, dtype=float))
        return z, z

    def gE_RT(self, x1, T):
        return np.zeros_like(np.asarray(x1, dtype=float))


# Literature NRTL set (Aspen Plus databank, regressed from Dortmund Data Bank;
# tabulated in Beneke et al., 2013, Appendix B).  alpha = 0.3.
NRTL_LIT = NRTL(a12=-0.8009, b12=246.18, a21=3.4578, b21=-586.0809, alpha=0.3,
                name="NRTL (Aspen/DDB)")

# ---------------------------------------------------------------------------
# Vapour phase: Tsonopoulos second virial coefficients (cm3/mol)
# ---------------------------------------------------------------------------


@dataclass
class Virial:
    enabled: bool = True
    _cache: dict = field(default_factory=dict)

    def B(self, T):
        key = round(T, 3)
        if key in self._cache:
            return self._cache[key]
        from chemicals.virial import BVirial_Tsonopoulos_extended as BT
        Tc = (514.71, 647.096)
        Pc = (6.268e6, 22.064e6)
        om = (0.646, 0.3443)
        mu = (1.69, 1.85)
        B11 = BT(T, Tc[0], Pc[0], om[0], species_type="alcohol", dipole=mu[0]) * 1e6
        B22 = BT(T, Tc[1], Pc[1], om[1], species_type="water", dipole=mu[1]) * 1e6
        # cross term: Lorentz-Berthelot type combining, no polar term (Tsonopoulos 1974)
        Tc12 = np.sqrt(Tc[0] * Tc[1])
        Vc = (167.1, 55.95)
        Zc = (0.241, 0.229)
        Vc12 = ((Vc[0] ** (1 / 3) + Vc[1] ** (1 / 3)) / 2) ** 3
        Zc12 = 0.5 * (Zc[0] + Zc[1])
        Pc12 = Zc12 * 8.314462618 * Tc12 / (Vc12 * 1e-6)
        om12 = 0.5 * (om[0] + om[1])
        B12 = BT(T, Tc12, Pc12, om12) * 1e6
        out = (B11, B22, B12)
        self._cache[key] = out
        return out

    def Phi(self, y1, T, P):
        """Phi_i = phi_hat_i/phi_i^sat * exp(-V_i^L (P - P_i^sat)/RT);  y_i Phi_i P = x_i gamma_i P_i^sat."""
        if not self.enabled:
            return 1.0, 1.0
        B11, B22, B12 = self.B(T)
        d12 = 2 * B12 - B11 - B22
        y2 = 1.0 - y1
        RT = 8.314462618 * T  # J/mol ; B in cm3/mol, P in kPa -> B*P [cm3 kPa/mol = mJ/mol]
        p1s, p2s = psat(0, T), psat(1, T)
        ln_phi1 = (B11 * P + P * y2**2 * d12) * 1e-3 / RT
        ln_phi2 = (B22 * P + P * y1**2 * d12) * 1e-3 / RT
        ln_phi1s = B11 * p1s * 1e-3 / RT
        ln_phi2s = B22 * p2s * 1e-3 / RT
        poy1 = -VL[0] * (P - p1s) * 1e-3 / RT
        poy2 = -VL[1] * (P - p2s) * 1e-3 / RT
        return np.exp(ln_phi1 - ln_phi1s + poy1), np.exp(ln_phi2 - ln_phi2s + poy2)


NO_VIRIAL = Virial(enabled=False)

# ---------------------------------------------------------------------------
# Bubble / dew point solvers
# ---------------------------------------------------------------------------


def bubble_T(x1: float, P: float, model, vap: Virial = NO_VIRIAL, T0: float | None = None):
    """Return (T, y1) at bubble point.  Brent on the summation residual,
    successive substitution on the vapour composition for the Phi loop."""
    x2 = 1.0 - x1
    y1 = x1
    for _ in range(50):
        def f(T):
            l1, l2 = model.lngamma(x1, T)
            P1, P2 = vap.Phi(y1, T, P)
            return (x1 * np.exp(l1) * psat(0, T) / P1 + x2 * np.exp(l2) * psat(1, T) / P2) / P - 1.0
        T = brentq(f, 250.0, 420.0, xtol=1e-10)
        l1, l2 = model.lngamma(x1, T)
        P1, P2 = vap.Phi(y1, T, P)
        k1 = np.exp(l1) * psat(0, T) / (P1 * P)
        k2 = np.exp(l2) * psat(1, T) / (P2 * P)
        y_new = x1 * k1 / (x1 * k1 + x2 * k2)
        if abs(y_new - y1) < 1e-11 or not vap.enabled:
            y1 = y_new
            break
        y1 = y_new
    return float(T), float(y1)


def dew_T(y1: float, P: float, model, vap: Virial = NO_VIRIAL):
    """Return (T, x1) at dew point by nested solution in x."""
    def resid_x(x1):
        T, y = bubble_T(x1, P, model, vap)
        return y - y1
    # y(x) is monotonic below the azeotrope for this system; search bracket
    xs = np.linspace(1e-6, 1 - 1e-6, 400)
    ys = np.array([bubble_T(x, P, model, vap)[1] for x in xs]) - y1
    idx = np.where(np.sign(ys[:-1]) != np.sign(ys[1:]))[0]
    if idx.size == 0:
        raise ValueError("dew point not bracketed")
    i = idx[0]
    x1 = brentq(resid_x, xs[i], xs[i + 1], xtol=1e-12)
    T, _ = bubble_T(x1, P, model, vap)
    return T, x1


def txy(model, P, n=401, vap: Virial = NO_VIRIAL, xgrid=None):
    xs = np.linspace(0.0, 1.0, n) if xgrid is None else np.asarray(xgrid)
    res = np.array([bubble_T(x, P, model, vap) for x in xs])
    return xs, res[:, 1], res[:, 0]


def azeotrope(model, P, vap: Virial = NO_VIRIAL, lo=0.5, hi=0.999999):
    """Locate x=y; return (x_az, T_az) or None."""
    g = lambda x: bubble_T(x, P, model, vap)[1] - x
    xs = np.linspace(lo, hi, 300)
    vals = np.array([g(x) for x in xs])
    idx = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    if idx.size == 0:
        return None
    i = idx[-1]
    xa = brentq(g, xs[i], xs[i + 1], xtol=1e-12)
    return xa, bubble_T(xa, P, model, vap)[0]


def mole_to_mass(x1):
    return x1 * MW[0] / (x1 * MW[0] + (1 - x1) * MW[1])


# ---------------------------------------------------------------------------
# McCabe-Thiele (CMO) tools
# ---------------------------------------------------------------------------


class EqCurve:
    """Tabulated y(x) and inverse x(y) from a VLE model at fixed P."""

    def __init__(self, model, P, n=4001, vap: Virial = NO_VIRIAL):
        x = np.concatenate([np.linspace(0, 0.02, 400), np.linspace(0.02, 1.0, n)[1:]])
        r = np.array([bubble_T(xx, P, model, vap) for xx in x])
        self.x, self.y, self.T = x, r[:, 1], r[:, 0]
        d = self.y - self.x
        cross = np.where((d[1:-1] > 0) & (d[2:] <= 0))[0]
        self.x_az = float(self.x[cross[0] + 1]) if cross.size else 1.0
        m = self.x <= self.x_az
        self._xm, self._ym = self.x[m], self.y[m]

    def y_of_x(self, x):
        return np.interp(x, self.x, self.y)

    def T_of_x(self, x):
        return np.interp(x, self.x, self.T)

    def x_of_y(self, y):
        """Inverse on the branch below the azeotrope (y > x)."""
        return np.interp(y, self._ym, self._xm)


def rmin_cmo(eq: EqCurve, xF, xD, q=1.0):
    """Minimum reflux for saturated-liquid feed (q=1) accounting for tangent
    pinches: the rectifying line must stay below y_eq on [xF, xD)."""
    xs = np.linspace(xF, xD - 1e-6, 20000)
    s = (xD - eq.y_of_x(xs)) / (xD - xs)
    k = int(np.argmax(s))
    smin = s[k]
    return smin / (1 - smin), xs[k]


def step_stages(eq: EqCurve, xF, xD, xB, Rr, q=1.0, max_stages=400):
    """Step off stages from the top (total condenser). Returns list of
    (x_n, y_n) per equilibrium stage (reboiler counted as a stage) and the
    optimal feed stage."""
    D_F = (xF - xB) / (xD - xB)
    L_D = Rr
    # per mole of feed
    D = D_F
    L = Rr * D
    V = L + D
    Lbar = L + q * 1.0
    Vbar = V - (1 - q) * 1.0
    B = 1.0 - D
    stages = []
    y = xD
    feed_stage = None
    rect = True
    for n in range(1, max_stages + 1):
        x = float(eq.x_of_y(y))
        stages.append((x, y))
        if x <= xB:
            break
        if rect and x <= xF:
            rect = False
            feed_stage = n
        if rect:
            y = (L / V) * x + D * xD / V
        else:
            y = (Lbar / Vbar) * x - B * xB / Vbar
        if y <= x + 1e-9:  # pinch: cannot proceed
            return stages, feed_stage, False
    return stages, feed_stage, True


def nmin_total_reflux(eq: EqCurve, xD, xB, max_stages=400):
    y = xD
    for n in range(1, max_stages + 1):
        x = float(eq.x_of_y(y))
        if x <= xB:
            return n
        y = x
    return np.inf


def rate_column(eq_true: EqCurve, N, NF, Rr, D_F, xF, q=1.0):
    """CMO rating: given N stages (incl. reboiler, total condenser not counted),
    feed stage NF, reflux ratio and distillate rate, find xD (and xB) that close
    the balance under the *true* equilibrium curve."""
    def bottom_mismatch(xD):
        B = 1.0 - D_F
        xB = (xF - D_F * xD) / B
        L, V = Rr * D_F, (Rr + 1) * D_F
        Lbar, Vbar = L + q, V - (1 - q)
        y = xD
        x = None
        for n in range(1, N + 1):
            x = float(eq_true.x_of_y(y))
            if n < NF:
                y = (L / V) * x + D_F * xD / V
            else:
                y = (Lbar / Vbar) * x - B * xB / Vbar
        return x - xB
    # scan for root
    xDs = np.linspace(xF + 1e-4, min(0.9999, xF / D_F - 1e-6), 600)
    vals = np.array([bottom_mismatch(v) for v in xDs])
    idx = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    if idx.size == 0:
        return None
    i = idx[-1]
    xD = brentq(bottom_mismatch, xDs[i], xDs[i + 1], xtol=1e-12)
    xB = (xF - D_F * xD) / (1 - D_F)
    return xD, xB
