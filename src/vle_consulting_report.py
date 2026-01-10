#!/usr/bin/env python3
"""
Ethanol-Water VLE Consulting Report - Thermodynamic Analysis

This script performs vapor-liquid equilibrium calculations for the
Ethanol(1) - Water(2) binary system at 1 atm (101.3 kPa).

Models implemented:
    - Raoult's Law (Ideal solution)
    - Van Laar activity coefficient model
    - NRTL (Non-Random Two-Liquid) activity coefficient model

Outputs (PNG) - saved to current directory:
    1. T-x-y diagram: TxY_Ethanol_Water_1atm.png
    2. y-x diagram: Yx_Ethanol_Water_1atm.png
    3. G^E vs x (NRTL): GE_vs_x_Ethanol_Water_1atm.png
    4. ΔHmix vs x (NRTL): Hmix_vs_x_Ethanol_Water_1atm.png

Usage:
    cd report/
    python vle_consulting_report.py

Reference: Sandler, Chemical, Biochemical, and Engineering Thermodynamics
           Chapters 6-10

Author: VLE Consulting Team
Date: 2026
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass
from scipy.optimize import brentq

# Universal gas constant
R = 8.314462618  # J/mol/K

# Output directory for figures
FIGURES_DIR = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)


# ============================================================
# Section 1: Antoine Vapor Pressure Correlation
# ============================================================

@dataclass(frozen=True)
class AntoineParams:
    """Antoine equation parameters for vapor pressure calculation.
    
    Equation: log10(Psat) = A - B / (T + C)
    where T is in °C and Psat is in the specified unit.
    """
    A: float
    B: float
    C: float
    psat_unit: str = "mmHg"


def psat_antoine_kpa(T_C: float, p: AntoineParams) -> float:
    """Calculate saturation pressure using Antoine equation.
    
    Args:
        T_C: Temperature in degrees Celsius
        p: AntoineParams dataclass with coefficients
        
    Returns:
        Saturation pressure in kPa
    """
    log10_psat = p.A - p.B / (T_C + p.C)
    Psat = 10.0 ** log10_psat
    # Convert from mmHg to kPa if necessary
    return Psat * 0.133322 if p.psat_unit.lower() == "mmhg" else Psat


# ============================================================
# Section 2: Activity Coefficient Models
# ============================================================

def gamma_ideal(x1: float, T_K: float, params=None):
    """Calculate activity coefficients for ideal solution (Raoult's Law).
    
    Args:
        x1: Mole fraction of component 1 in liquid phase
        T_K: Temperature in Kelvin (unused for ideal case)
        params: Unused parameter placeholder
        
    Returns:
        Tuple of (gamma1, gamma2) = (1.0, 1.0)
    """
    return 1.0, 1.0


@dataclass(frozen=True)
class VanLaarParams:
    """Van Laar model parameters for binary mixture."""
    A12: float
    A21: float


def gamma_vanlaar(x1: float, T_K: float, p: VanLaarParams):
    """Calculate Van Laar activity coefficients for a binary mixture.
    
    Args:
        x1: Mole fraction of component 1
        T_K: Temperature in Kelvin (unused, Van Laar is T-independent)
        p: VanLaarParams dataclass
        
    Returns:
        Tuple of (gamma1, gamma2)
    """
    x2 = 1.0 - x1
    denom = (p.A12 * x1 + p.A21 * x2) ** 2
    denom = max(denom, 1e-30)  # Prevent division by zero
    
    ln_g1 = p.A12 * (p.A21 * x2) ** 2 / denom
    ln_g2 = p.A21 * (p.A12 * x1) ** 2 / denom
    
    return np.exp(ln_g1), np.exp(ln_g2)


@dataclass(frozen=True)
class NRTLParams:
    """NRTL model parameters for binary mixture.
    
    The temperature-dependent tau parameters are calculated as:
        tau_ij(T) = a_ij + b_ij / T
    where T is in Kelvin.
    """
    alpha12: float  # Non-randomness factor for 1-2 interaction
    alpha21: float  # Non-randomness factor for 2-1 interaction
    a12: float      # Constant term for tau12
    b12: float      # Temperature-dependent term for tau12 (K)
    a21: float      # Constant term for tau21
    b21: float      # Temperature-dependent term for tau21 (K)


def tau(a: float, b: float, T: float) -> float:
    """Calculate NRTL tau parameter.
    
    Args:
        a: Constant term
        b: Temperature-dependent term (K)
        T: Temperature in Kelvin
        
    Returns:
        tau value
    """
    return a + b / T


def gamma_nrtl(x1: float, T_K: float, p: NRTLParams):
    """Calculate NRTL activity coefficients for a binary mixture.
    
    Based on the Non-Random Two-Liquid model (Renon & Prausnitz, 1968).
    
    Args:
        x1: Mole fraction of component 1
        T_K: Temperature in Kelvin
        p: NRTLParams dataclass
        
    Returns:
        Tuple of (gamma1, gamma2)
    """
    x2 = 1.0 - x1
    
    # Calculate temperature-dependent tau parameters
    tau12 = tau(p.a12, p.b12, T_K)
    tau21 = tau(p.a21, p.b21, T_K)
    
    # Calculate G parameters
    G12 = np.exp(-p.alpha12 * tau12)
    G21 = np.exp(-p.alpha21 * tau21)
    
    # Denominators (with protection against division by zero)
    D1 = max(x1 + x2 * G21, 1e-30)
    D2 = max(x2 + x1 * G12, 1e-30)
    
    # Calculate ln(gamma) values
    ln_g1 = (x2**2) * (tau21 * (G21 / D1)**2 + (tau12 * G12) / (D2**2))
    ln_g2 = (x1**2) * (tau12 * (G12 / D2)**2 + (tau21 * G21) / (D1**2))
    
    return np.exp(ln_g1), np.exp(ln_g2)


# ============================================================
# Section 3: Bubble-Point VLE Calculations
# ============================================================

def bubble_residual(T_C, x1, P_kPa, psat1, psat2, gamma_model, gamma_params):
    """Calculate the bubble point residual for VLE.
    
    The bubble point condition is:
        sum(x_i * gamma_i * Psat_i(T)) = P
    
    Args:
        T_C: Temperature in degrees Celsius
        x1: Liquid mole fraction of component 1
        P_kPa: System pressure in kPa
        psat1: Antoine parameters for component 1
        psat2: Antoine parameters for component 2
        gamma_model: Activity coefficient function
        gamma_params: Parameters for the activity coefficient model
        
    Returns:
        Residual value (should be zero at bubble point)
    """
    T_K = T_C + 273.15
    g1, g2 = gamma_model(x1, T_K, gamma_params)
    
    P_calc = (x1 * g1 * psat_antoine_kpa(T_C, psat1) +
              (1 - x1) * g2 * psat_antoine_kpa(T_C, psat2))
    
    return P_calc - P_kPa


def solve_bubble_T(x1, P_kPa, psat1, psat2, gamma_model, gamma_params):
    """Solve for bubble point temperature at given composition and pressure.
    
    Args:
        x1: Liquid mole fraction of component 1
        P_kPa: System pressure in kPa
        psat1: Antoine parameters for component 1
        psat2: Antoine parameters for component 2
        gamma_model: Activity coefficient function
        gamma_params: Parameters for the activity coefficient model
        
    Returns:
        Bubble point temperature in degrees Celsius
    """
    return brentq(
        bubble_residual, 40.0, 120.0,
        args=(x1, P_kPa, psat1, psat2, gamma_model, gamma_params)
    )


def compute_Txy(P_kPa, psat1, psat2, gamma_model, gamma_params, x_grid):
    """Compute T-x-y data for VLE at fixed pressure.
    
    Args:
        P_kPa: System pressure in kPa
        psat1: Antoine parameters for component 1
        psat2: Antoine parameters for component 2
        gamma_model: Activity coefficient function
        gamma_params: Parameters for the activity coefficient model
        x_grid: Array of liquid compositions to evaluate
        
    Returns:
        Tuple of (x_grid, y_array, T_array)
    """
    T_list, y_list = [], []
    
    for x1 in x_grid:
        # Solve for bubble temperature
        T_C = solve_bubble_T(x1, P_kPa, psat1, psat2, gamma_model, gamma_params)
        T_K = T_C + 273.15
        
        # Calculate activity coefficients at bubble point
        g1, g2 = gamma_model(x1, T_K, gamma_params)
        
        # Calculate vapor phase compositions
        y1 = x1 * g1 * psat_antoine_kpa(T_C, psat1) / P_kPa
        y2 = (1 - x1) * g2 * psat_antoine_kpa(T_C, psat2) / P_kPa
        
        # Normalize to ensure y1 + y2 = 1
        y1 /= (y1 + y2)
        
        T_list.append(T_C)
        y_list.append(y1)
    
    return x_grid, np.array(y_list), np.array(T_list)


# ============================================================
# Section 4: Excess Property Calculations (NRTL)
# ============================================================

def excess_gibbs_GE(x1: float, T_K: float, p: NRTLParams) -> float:
    """Calculate molar excess Gibbs energy using NRTL model.
    
    G^E = RT * sum(x_i * ln(gamma_i))
    
    Args:
        x1: Mole fraction of component 1
        T_K: Temperature in Kelvin
        p: NRTLParams dataclass
        
    Returns:
        Excess Gibbs energy in J/mol
    """
    g1, g2 = gamma_nrtl(x1, T_K, p)
    return R * T_K * (x1 * np.log(g1) + (1 - x1) * np.log(g2))


def excess_enthalpy_HE(x1: float, T_K: float, p: NRTLParams, dT: float = 1e-3) -> float:
    """Calculate molar excess enthalpy using numerical differentiation.
    
    H^E = -R * T^2 * sum(x_i * d(ln(gamma_i))/dT)
    
    This approximates the heat of mixing (ΔH_mix) when the ideal mixing
    enthalpy is negligible.
    
    Args:
        x1: Mole fraction of component 1
        T_K: Temperature in Kelvin
        p: NRTLParams dataclass
        dT: Temperature step for numerical differentiation
        
    Returns:
        Excess enthalpy in J/mol
    """
    # Central difference for derivative
    g1p, g2p = gamma_nrtl(x1, T_K + dT, p)
    g1m, g2m = gamma_nrtl(x1, T_K - dT, p)
    
    dln1 = (np.log(g1p) - np.log(g1m)) / (2 * dT)
    dln2 = (np.log(g2p) - np.log(g2m)) / (2 * dT)
    
    return -R * T_K**2 * (x1 * dln1 + (1 - x1) * dln2)


# ============================================================
# Section 5: Main Execution
# ============================================================

def main():
    """Main function to generate all VLE plots and data."""
    
    print("=" * 60)
    print("Ethanol-Water VLE Consulting Report")
    print("Thermodynamic Analysis at 1 atm (101.3 kPa)")
    print("=" * 60)
    
    # System pressure
    P_kPa = 101.3  # 1 atm
    
    # Antoine constants (valid for the boiling range)
    # Source: NIST Chemistry WebBook / Perry's Chemical Engineers' Handbook
    ethanol = AntoineParams(A=8.20417, B=1642.89, C=230.3, psat_unit="mmHg")
    water = AntoineParams(A=8.07131, B=1730.63, C=233.426, psat_unit="mmHg")
    
    print("\nAntoine Parameters:")
    print(f"  Ethanol: A={ethanol.A}, B={ethanol.B}, C={ethanol.C}")
    print(f"  Water:   A={water.A}, B={water.B}, C={water.C}")
    
    # Van Laar parameters (illustrative values for Ethanol-Water)
    vanlaar_p = VanLaarParams(A12=1.65, A21=0.95)
    print(f"\nVan Laar Parameters: A12={vanlaar_p.A12}, A21={vanlaar_p.A21}")
    
    # NRTL parameters - ChemSep Database (Validated)
    # Source: thermo/Interaction Parameters/ChemSep/nrtl.json
    # System: Ethanol (64-17-5) - Water (7732-18-5)
    # 
    # ChemSep NRTL format: tau_ij = b_ij / T  (T in Kelvin, no constant term)
    # This format ensures temperature-dependent activity coefficients,
    # enabling non-zero excess enthalpy (H^E) calculations.
    # 
    # Validated by: tests/test_nrtl.py::test_NRTL_chemsep
    # At T=343.15K, xs=[0.252, 0.748]: gammas ≈ [1.985, 1.146]
    nrtl_p = NRTLParams(
        alpha12=0.2937,   # ChemSep alphaij (symmetric for this system)
        alpha21=0.2937,
        a12=0.0,          # ChemSep format: tau = b/T (no constant term)
        b12=-29.167,      # Ethanol→Water interaction parameter (K)
        a21=0.0,          # ChemSep format: tau = b/T (no constant term)
        b21=624.868       # Water→Ethanol interaction parameter (K)
    )
    print(f"\nNRTL Parameters (DECHEMA/Gmehling):")
    print(f"  alpha12={nrtl_p.alpha12}, alpha21={nrtl_p.alpha21}")
    print(f"  a12={nrtl_p.a12}, b12={nrtl_p.b12} K")
    print(f"  a21={nrtl_p.a21}, b21={nrtl_p.b21} K")
    
    # Composition grid
    x = np.linspace(0.01, 0.99, 60)
    
    # ========================
    # VLE Calculations
    # ========================
    print("\nCalculating VLE data...")
    
    # Ideal model (Raoult's Law)
    xi, yi, Ti = compute_Txy(P_kPa, ethanol, water, gamma_ideal, None, x)
    print("  Ideal (Raoult's Law): Complete")
    
    # Van Laar model
    xv, yv, Tv = compute_Txy(P_kPa, ethanol, water, gamma_vanlaar, vanlaar_p, x)
    print("  Van Laar model: Complete")
    
    # NRTL model
    xn, yn, Tn = compute_Txy(P_kPa, ethanol, water, gamma_nrtl, nrtl_p, x)
    print("  NRTL model: Complete")
    
    # ========================
    # Plot 1: T-x-y Diagram
    # ========================
    print("\nGenerating plots...")
    
    plt.figure(figsize=(10, 7))
    
    # Liquid lines (solid)
    plt.plot(xi, Ti, 'b-', linewidth=1.5, label="Ideal (Raoult): Bubble")
    plt.plot(xv, Tv, 'g-', linewidth=1.5, label="Van Laar: Bubble")
    plt.plot(xn, Tn, 'r-', linewidth=2, label="NRTL: Bubble")
    
    # Vapor lines (dashed)
    plt.plot(yi, Ti, 'b--', linewidth=1.5, label="Ideal (Raoult): Dew")
    plt.plot(yv, Tv, 'g--', linewidth=1.5, label="Van Laar: Dew")
    plt.plot(yn, Tn, 'r--', linewidth=2, label="NRTL: Dew")
    
    plt.xlabel("Mole Fraction of Ethanol ($x_1$ or $y_1$)", fontsize=12)
    plt.ylabel("Temperature (°C)", fontsize=12)
    plt.title("Ethanol-Water T-x-y Diagram at 1 atm (101.3 kPa)", fontsize=14)
    plt.legend(loc='best', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "TxY_Ethanol_Water_1atm.png"), dpi=300, bbox_inches='tight')
    print("  Saved: TxY_Ethanol_Water_1atm.png")
    plt.close()
    
    # ========================
    # Plot 2: y-x Diagram
    # ========================
    plt.figure(figsize=(8, 8))
    
    # Reference line y = x
    plt.plot([0, 1], [0, 1], 'k:', linewidth=1, label="y = x (Reference)")
    
    # Model predictions
    plt.plot(xi, yi, 'b-', linewidth=1.5, label="Ideal (Raoult)")
    plt.plot(xv, yv, 'g-', linewidth=1.5, label="Van Laar")
    plt.plot(xn, yn, 'r-', linewidth=2, label="NRTL")
    
    plt.xlabel("Liquid Mole Fraction of Ethanol, $x_1$", fontsize=12)
    plt.ylabel("Vapor Mole Fraction of Ethanol, $y_1$", fontsize=12)
    plt.title("Ethanol-Water y-x Diagram at 1 atm (101.3 kPa)", fontsize=14)
    plt.legend(loc='lower right', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.gca().set_aspect('equal')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "Yx_Ethanol_Water_1atm.png"), dpi=300, bbox_inches='tight')
    print("  Saved: Yx_Ethanol_Water_1atm.png")
    plt.close()
    
    # ========================
    # Excess Properties (NRTL)
    # ========================
    print("\nCalculating excess properties (NRTL model)...")
    
    TnK = Tn + 273.15  # Convert to Kelvin
    GE = np.array([excess_gibbs_GE(x1, T_K, nrtl_p) for x1, T_K in zip(xn, TnK)])
    HE = np.array([excess_enthalpy_HE(x1, T_K, nrtl_p) for x1, T_K in zip(xn, TnK)])
    
    print(f"  Max G^E = {np.max(GE):.2f} J/mol at x1 = {xn[np.argmax(GE)]:.3f}")
    print(f"  Max H^E = {np.max(HE):.2f} J/mol at x1 = {xn[np.argmax(HE)]:.3f}")
    
    # ========================
    # Plot 3: G^E vs x
    # ========================
    plt.figure(figsize=(10, 6))
    plt.plot(xn, GE / 1000.0, 'b-', linewidth=2)
    plt.xlabel("Liquid Mole Fraction of Ethanol, $x_1$", fontsize=12)
    plt.ylabel("Excess Gibbs Energy, $G^E$ (kJ/mol)", fontsize=12)
    plt.title("NRTL Excess Gibbs Energy vs Composition\n(Evaluated along bubble line at 1 atm)", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.axhline(0, color='k', linewidth=0.5)
    plt.xlim(0, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "GE_vs_x_Ethanol_Water_1atm.png"), dpi=300, bbox_inches='tight')
    print("  Saved: GE_vs_x_Ethanol_Water_1atm.png")
    plt.close()
    
    # ========================
    # Plot 4: ΔH_mix vs x
    # ========================
    plt.figure(figsize=(10, 6))
    plt.plot(xn, HE / 1000.0, 'm-', linewidth=2)
    plt.xlabel("Liquid Mole Fraction of Ethanol, $x_1$", fontsize=12)
    plt.ylabel("Excess Enthalpy, $H^E \\approx \\Delta H_{mix}$ (kJ/mol)", fontsize=12)
    plt.title("NRTL Excess Enthalpy vs Composition\n(Evaluated along bubble line at 1 atm)", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.axhline(0, color='k', linewidth=0.5)
    plt.xlim(0, 1)
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, "Hmix_vs_x_Ethanol_Water_1atm.png"), dpi=300, bbox_inches='tight')
    print("  Saved: Hmix_vs_x_Ethanol_Water_1atm.png")
    plt.close()
    
    # ========================
    # Summary
    # ========================
    print("\n" + "=" * 60)
    print("CALCULATION SUMMARY")
    print("=" * 60)
    print(f"System: Ethanol (1) - Water (2)")
    print(f"Pressure: {P_kPa} kPa (1 atm)")
    print(f"\nBoiling Points (pure components):")
    print(f"  Ethanol: {Ti[0]:.1f}°C (at x1 → 1)")
    print(f"  Water:   {Ti[-1]:.1f}°C (at x1 → 0)")
    
    # Find azeotrope (if y1 = x1 approximately)
    diff = np.abs(yn - xn)
    azeotrope_idx = np.argmin(diff)
    if diff[azeotrope_idx] < 0.01:  # Close to azeotrope
        print(f"\nAzeotrope detected (NRTL):")
        print(f"  x1 = y1 ≈ {xn[azeotrope_idx]:.3f}")
        print(f"  T ≈ {Tn[azeotrope_idx]:.1f}°C")
    
    print("\n" + "=" * 60)
    print("All calculations and plots completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
