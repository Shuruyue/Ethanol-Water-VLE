#!/usr/bin/env python3
"""Phase 1: Scientific Correctness Verification Script.

Validates all VLE formulas against known reference values:
- NRTL activity coefficients vs upstream thermo library
- Van Laar activity coefficients vs hand-calculated values
- Pure-component boiling points from Antoine equation
- Excess properties (GE, HE, CPE) sign and magnitude
- Azeotrope prediction vs literature
"""

from __future__ import annotations

import sys
from pathlib import Path
from math import exp, log

import numpy as np

# Ensure UTF-8 output on Windows (cp950/cp1252 cannot encode ≈, →, etc.)
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from models import (
    AntoineParams,
    NRTLParams,
    VanLaarParams,
    gamma_ideal,
    gamma_nrtl,
    gamma_van_laar,
    psat_kpa,
    tau_value,
)
from solver import compute_txy, solve_bubble_temperature
from analysis import (
    excess_gibbs_energy,
    excess_enthalpy,
    excess_heat_capacity,
    find_azeotrope,
)
from parameter_store import load_system_parameters

# ── Reference Data ──────────────────────────────────────────────────────────

# Ethanol Antoine (mmHg): A=8.20417, B=1642.89, C=230.3
# Water Antoine (mmHg):   A=8.07131, B=1730.63, C=233.426

# Literature boiling points at 1 atm (101.325 kPa):
#   Ethanol: 78.37 °C
#   Water:   100.0 °C

# Literature azeotrope at 1 atm:
#   x_ethanol ≈ 0.8943, T ≈ 78.15 °C

# DDBST Example P05.01b: ethanol-water at 70°C (343.15 K), x=[0.252, 0.748]
# with tauB=[[0, -61.025], [673.236, 0]], alphaC=[[0, 0.2974],[0.2974, 0]]
# Result: gamma = [1.936, 1.154]

PASS_COUNT = 0
FAIL_COUNT = 0


def check(name: str, computed, expected, rtol: float = 0.01):
    """Check if computed value matches expected within relative tolerance."""
    global PASS_COUNT, FAIL_COUNT
    if isinstance(expected, bool):
        ok = bool(computed) == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            FAIL_COUNT += 1
        else:
            PASS_COUNT += 1
        print(f"  [{status}] {name}: computed={computed}, expected={expected}")
    elif isinstance(expected, (int, float)):
        if abs(expected) < 1e-15:
            error = abs(computed - expected)
            ok = error < 1e-10
        else:
            error = abs(computed - expected) / abs(expected)
            ok = error < rtol
        status = "PASS" if ok else "FAIL"
        if not ok:
            FAIL_COUNT += 1
        else:
            PASS_COUNT += 1
        print(f"  [{status}] {name}: computed={computed:.6g}, expected={expected:.6g}, err={error:.2e}")
    else:
        ok = computed == expected
        status = "PASS" if ok else "FAIL"
        if not ok:
            FAIL_COUNT += 1
        else:
            PASS_COUNT += 1
        print(f"  [{status}] {name}: {computed}")


def section(title: str):
    print(f"\n{'='*70}")
    print(f"  {title}")
    print(f"{'='*70}")


# ── 1.1 NRTL Formula Verification ──────────────────────────────────────────

def verify_nrtl_formula():
    section("1.1 NRTL Activity Coefficients")

    # Use project's NRTL parameters
    system = load_system_parameters("ethanol_water_1atm")
    nrtl = system.nrtl

    # Test point: x1=0.5, T=351.15 K (78°C)
    T_K = 351.15
    x1 = 0.5
    g1, g2 = gamma_nrtl(x1, T_K, nrtl)

    # Hand-calculate for comparison
    x2 = 1.0 - x1
    tau12 = nrtl.a12 + nrtl.b12 / T_K
    tau21 = nrtl.a21 + nrtl.b21 / T_K
    G12 = exp(-nrtl.alpha12 * tau12)
    G21 = exp(-nrtl.alpha21 * tau21)
    D1 = x1 + x2 * G21
    D2 = x2 + x1 * G12
    ln_g1_hand = x2**2 * (tau21 * (G21 / D1)**2 + tau12 * G12 / D2**2)
    ln_g2_hand = x1**2 * (tau12 * (G12 / D2)**2 + tau21 * G21 / D1**2)
    g1_hand = exp(ln_g1_hand)
    g2_hand = exp(ln_g2_hand)

    print(f"  Parameters: alpha12={nrtl.alpha12}, a12={nrtl.a12}, b12={nrtl.b12}")
    print(f"  Parameters: alpha21={nrtl.alpha21}, a21={nrtl.a21}, b21={nrtl.b21}")
    print(f"  Test point: x1={x1}, T={T_K} K")
    print(f"  tau12={tau12:.6f}, tau21={tau21:.6f}")
    print(f"  G12={G12:.6f}, G21={G21:.6f}")

    check("gamma1 vs hand calc", g1, g1_hand, rtol=1e-10)
    check("gamma2 vs hand calc", g2, g2_hand, rtol=1e-10)

    # Boundary tests
    print("\n  Boundary tests:")
    g1_0, g2_0 = gamma_nrtl(0.001, T_K, nrtl)
    g1_1, g2_1 = gamma_nrtl(0.999, T_K, nrtl)
    print(f"  x1→0: gamma1={g1_0:.4f} (should be gamma1_inf), gamma2={g2_0:.4f} (should→1)")
    print(f"  x1→1: gamma1={g1_1:.4f} (should→1), gamma2={g2_1:.4f} (should be gamma2_inf)")
    check("gamma2 → 1 as x1→0", g2_0, 1.0, rtol=0.01)
    check("gamma1 → 1 as x1→1", g1_1, 1.0, rtol=0.01)

    # Cross-validate with upstream thermo (using same tau/alpha)
    print("\n  Cross-validation with upstream NRTL_gammas_binaries:")
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
        from thermo.nrtl import NRTL_gammas_binaries
        upstream_gammas = NRTL_gammas_binaries(
            [x1, x2], tau12, tau21, nrtl.alpha12, nrtl.alpha21
        )
        check("gamma1 vs upstream thermo", g1, upstream_gammas[0], rtol=1e-8)
        check("gamma2 vs upstream thermo", g2, upstream_gammas[1], rtol=1e-8)
    except ImportError:
        print("  [SKIP] Could not import thermo.nrtl for cross-validation")


# ── 1.2 Van Laar Verification ──────────────────────────────────────────────

def verify_van_laar():
    section("1.2 Van Laar Activity Coefficients")

    params = VanLaarParams(A12=1.65, A21=0.95)
    T_K = 351.15  # not used by Van Laar

    # x1=0.5 hand calculation
    x1, x2 = 0.5, 0.5
    denom = (params.A12 * x1 + params.A21 * x2) ** 2
    ln_g1 = params.A12 * (params.A21 * x2) ** 2 / denom
    ln_g2 = params.A21 * (params.A12 * x1) ** 2 / denom
    g1_hand = exp(ln_g1)
    g2_hand = exp(ln_g2)

    g1, g2 = gamma_van_laar(x1, T_K, params)
    check("gamma1 (x1=0.5) vs hand", g1, g1_hand, rtol=1e-10)
    check("gamma2 (x1=0.5) vs hand", g2, g2_hand, rtol=1e-10)

    # Infinite dilution: x1→0, gamma1 → exp(A12)
    g1_inf, _ = gamma_van_laar(0.001, T_K, params)
    check("gamma1_inf vs exp(A12)", g1_inf, exp(params.A12), rtol=0.01)

    _, g2_inf = gamma_van_laar(0.999, T_K, params)
    check("gamma2_inf vs exp(A21)", g2_inf, exp(params.A21), rtol=0.01)


# ── 1.3 Antoine & Bubble Point Verification ────────────────────────────────

def verify_bubble_point():
    section("1.3 Antoine & Bubble Point")

    system = load_system_parameters("ethanol_water_1atm")
    eth = system.component_1.antoine
    wat = system.component_2.antoine

    # Pure ethanol boiling point at 1 atm
    psat_eth_78 = psat_kpa(78.37, eth)
    print(f"  Ethanol Psat at 78.37°C = {psat_eth_78:.3f} kPa (should ≈ 101.325)")
    check("Ethanol Psat@78.37°C ≈ 101.3 kPa", psat_eth_78, 101.325, rtol=0.02)

    # Pure water boiling point at 1 atm
    psat_wat_100 = psat_kpa(100.0, wat)
    print(f"  Water Psat at 100.0°C = {psat_wat_100:.3f} kPa (should ≈ 101.325)")
    check("Water Psat@100°C ≈ 101.3 kPa", psat_wat_100, 101.325, rtol=0.02)

    # Solve pure-component bubble points
    T_eth = solve_bubble_temperature(
        x1=1.0, pressure_kpa=101.325,
        component_1=eth, component_2=wat,
        gamma_fn=gamma_ideal, gamma_params=None,
    )
    check("Pure ethanol Tb", T_eth, 78.37, rtol=0.005)

    T_wat = solve_bubble_temperature(
        x1=0.0, pressure_kpa=101.325,
        component_1=eth, component_2=wat,
        gamma_fn=gamma_ideal, gamma_params=None,
    )
    check("Pure water Tb", T_wat, 100.0, rtol=0.005)

    # y normalization check: verify y1+y2 ≈ 1 without forced renormalization
    print("\n  y-normalization check at x1=0.5 (NRTL):")
    T_bub = solve_bubble_temperature(
        x1=0.5, pressure_kpa=101.325,
        component_1=eth, component_2=wat,
        gamma_fn=gamma_nrtl, gamma_params=system.nrtl,
    )
    T_K = T_bub + 273.15
    g1, g2 = gamma_nrtl(0.5, T_K, system.nrtl)
    y1_raw = 0.5 * g1 * psat_kpa(T_bub, eth) / 101.325
    y2_raw = 0.5 * g2 * psat_kpa(T_bub, wat) / 101.325
    sum_y = y1_raw + y2_raw
    print(f"  y1_raw={y1_raw:.6f}, y2_raw={y2_raw:.6f}, sum={sum_y:.6f}")
    check("y1+y2 ≈ 1 (raw, no renormalization)", sum_y, 1.0, rtol=0.001)


# ── 1.4 Excess Properties Verification ─────────────────────────────────────

def verify_excess_properties():
    section("1.4 Excess Properties (GE, HE, CPE)")

    system = load_system_parameters("ethanol_water_1atm")
    nrtl = system.nrtl

    # GE at boundaries must be 0
    T_K = 351.38  # ~78.2 deg C
    ge_0 = excess_gibbs_energy(0.001, T_K, nrtl)  # near 0, avoiding 0*ln(1) indeterminate
    ge_1 = excess_gibbs_energy(0.999, T_K, nrtl)
    print(f"  GE(x1=0.001) = {ge_0:.4f} J/mol (should be near 0)")
    print(f"  GE(x1=0.999) = {ge_1:.4f} J/mol (should be near 0)")
    check("GE(x1->0) near 0", abs(ge_0) < 50, True)
    check("GE(x1->1) near 0", abs(ge_1) < 50, True)

    # GE at midpoint: should be positive for ethanol-water
    ge_05 = excess_gibbs_energy(0.5, T_K, nrtl)
    print(f"  GE(x1=0.5, T=78.2°C) = {ge_05:.2f} J/mol")
    check("GE(x1=0.5) > 0", ge_05 > 0, True)

    # GE hand-check: GE = RT * (x1*ln(g1) + x2*ln(g2))
    g1, g2 = gamma_nrtl(0.5, T_K, nrtl)
    ge_hand = 8.314462618 * T_K * (0.5 * log(g1) + 0.5 * log(g2))
    check("GE vs RT(x1*lnγ1+x2*lnγ2)", ge_05, ge_hand, rtol=1e-10)

    # HE sign: ethanol-water mixing is generally exothermic at dilute,
    # but can be endothermic at moderate concentrations
    he_05 = excess_enthalpy(0.5, T_K, nrtl)
    print(f"  HE(x1=0.5, T=78.2°C) = {he_05:.2f} J/mol")
    # Literature: HE for ethanol-water at equimolar ~25°C is ~-700 to -900 J/mol (exothermic)
    # At 78°C the behavior may differ, but should be finite
    check("HE is finite", np.isfinite(he_05), True)

    # CPE check: should be finite
    cpe_05 = excess_heat_capacity(0.5, T_K, nrtl)
    print(f"  CPE(x1=0.5, T=78.2°C) = {cpe_05:.4f} J/(mol·K)")
    check("CPE is finite", np.isfinite(cpe_05), True)

    # Numerical derivative precision test: compare HE with tighter delta
    he_tight = excess_enthalpy(0.5, T_K, nrtl, delta_t=1e-4)
    he_loose = excess_enthalpy(0.5, T_K, nrtl, delta_t=1e-2)
    print(f"  HE(delta=1e-4)={he_tight:.4f}, HE(delta=1e-3)={he_05:.4f}, HE(delta=1e-2)={he_loose:.4f}")
    check("HE stable across delta_t", he_tight, he_05, rtol=0.01)


# ── 1.5 Azeotrope Verification ─────────────────────────────────────────────

def verify_azeotrope():
    section("1.5 Azeotrope Detection")

    system = load_system_parameters("ethanol_water_1atm")
    x_grid = np.linspace(0.01, 0.99, 200)

    x, y, T = compute_txy(
        pressure_kpa=101.325,
        component_1=system.component_1.antoine,
        component_2=system.component_2.antoine,
        gamma_fn=gamma_nrtl,
        gamma_params=system.nrtl,
        x_grid=x_grid,
    )

    azeo = find_azeotrope(x, y, T, tolerance=0.01)
    if azeo:
        print(f"  Detected azeotrope: x1={azeo.x1:.4f}, T={azeo.temperature_c:.2f}°C, |x-y|={azeo.abs_error:.5f}")
        # Literature: x_eth ≈ 0.8943, T ≈ 78.15°C
        check("Azeotrope x1 ≈ 0.894", azeo.x1, 0.894, rtol=0.02)
        check("Azeotrope T ≈ 78.15°C", azeo.temperature_c, 78.15, rtol=0.01)
    else:
        global FAIL_COUNT
        FAIL_COUNT += 1
        print("  [FAIL ❌] No azeotrope detected (should exist!)")


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print("╔══════════════════════════════════════════════════════════════════════╗")
    print("║  Phase 1: Scientific Correctness Verification                      ║")
    print("║  Ethanol-Water VLE Project                                         ║")
    print("╚══════════════════════════════════════════════════════════════════════╝")

    verify_nrtl_formula()
    verify_van_laar()
    verify_bubble_point()
    verify_excess_properties()
    verify_azeotrope()

    print(f"\n{'='*70}")
    print(f"  SUMMARY: {PASS_COUNT} passed, {FAIL_COUNT} failed")
    print(f"{'='*70}")

    if FAIL_COUNT > 0:
        print("\n  ⚠️  Some checks failed! Review the details above.")
        sys.exit(1)
    else:
        print("\n  ✅ All scientific checks passed!")


if __name__ == "__main__":
    main()
