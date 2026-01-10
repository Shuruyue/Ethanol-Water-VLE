"""
Diffusion Calculations Module
擴散計算模組

Contains functions for calculating diffusion coefficients and fluxes
using Fuller-Schettler-Giddings correlation, Fick's Law, and Stefan diffusion.
"""

import math
from typing import Tuple, Optional

from .properties import (
    ChemicalProperties, 
    R_GAS, 
    AIR_MW, 
    AIR_DIFFUSION_VOLUME,
    antoine_vapor_pressure,
    ETHANOL,
)


# ============================================================
# Diffusion Coefficient Calculations
# ============================================================

def diffusion_coefficient_fsg(
    T_celsius: float,
    P_atm: float,
    props_A: ChemicalProperties,
    mw_B: float = AIR_MW,
    diffusion_vol_B: float = AIR_DIFFUSION_VOLUME
) -> float:
    """
    Calculate binary gas diffusion coefficient using Fuller-Schettler-Giddings.
    使用 Fuller-Schettler-Giddings 關係式計算二元氣相擴散係數
    
    The FSG correlation:
        D_AB = 0.001 * T^1.75 * sqrt(1/M_A + 1/M_B) / 
               [P * (Σv_A^(1/3) + Σv_B^(1/3))^2]
    
    Args:
        T_celsius: Temperature (°C)
        P_atm: Pressure (atm)
        props_A: Properties of diffusing species (e.g., ethanol)
        mw_B: Molecular weight of carrier gas (default: air)
        diffusion_vol_B: Diffusion volume of carrier gas (default: air)
        
    Returns:
        Binary diffusion coefficient D_AB (m²/s)
        
    Example:
        >>> D_ab = diffusion_coefficient_fsg(25.0, 1.0, ETHANOL)
        >>> print(f"D_ethanol-air = {D_ab:.2e} m²/s")
        D_ethanol-air = 1.18e-05 m²/s
    
    Reference:
        Fuller, E.N., Schettler, P.D., Giddings, J.C. (1966).
        Industrial & Engineering Chemistry, 58(5), 18-27.
    """
    T_K = T_celsius + 273.15
    
    M_A = props_A.mw
    M_B = mw_B
    v_A = props_A.diffusion_volume
    v_B = diffusion_vol_B
    
    # FSG equation
    numerator = 0.001 * (T_K ** 1.75) * math.sqrt(1/M_A + 1/M_B)
    denominator = P_atm * (v_A**(1/3) + v_B**(1/3))**2
    
    D_AB = numerator / denominator  # cm²/s
    
    # Convert cm²/s to m²/s
    D_AB_m2s = D_AB * 1e-4
    
    return D_AB_m2s


def diffusion_coefficient_wilke_lee(
    T_celsius: float,
    P_atm: float,
    props_A: ChemicalProperties,
    mw_B: float = AIR_MW,
    sigma_AB: Optional[float] = None,
    omega_D: Optional[float] = 1.0
) -> float:
    """
    Calculate diffusion coefficient using Wilke-Lee correlation.
    使用 Wilke-Lee 關係式計算擴散係數（替代方法）
    
    Simplified form:
        D_AB = 1.084e-4 * T^1.5 * sqrt(1/M_A + 1/M_B) / (P * σ_AB² * Ω_D)
    
    Args:
        T_celsius: Temperature (°C)
        P_atm: Pressure (atm)
        props_A: Properties of species A
        mw_B: Molecular weight of species B
        sigma_AB: Collision diameter (Å), estimated if None
        omega_D: Collision integral, default 1.0
        
    Returns:
        Diffusion coefficient (m²/s)
    """
    T_K = T_celsius + 273.15
    M_A = props_A.mw
    M_B = mw_B
    
    # Estimate sigma_AB if not provided (rough estimate)
    if sigma_AB is None:
        sigma_AB = 4.0  # Approximate for organic-air
    
    numerator = 1.084e-4 * (T_K ** 1.5) * math.sqrt(1/M_A + 1/M_B)
    denominator = P_atm * (sigma_AB ** 2) * omega_D
    
    D_AB = numerator / denominator  # cm²/s
    D_AB_m2s = D_AB * 1e-4
    
    return D_AB_m2s


# ============================================================
# Diffusion Flux Calculations
# ============================================================

def fick_diffusion_flux(
    D_AB: float,
    dC_dz: float
) -> float:
    """
    Calculate diffusion flux using Fick's First Law.
    使用 Fick 第一定律計算擴散通量
    
    J_A = -D_AB * (dC_A/dz)
    
    Args:
        D_AB: Binary diffusion coefficient (m²/s)
        dC_dz: Concentration gradient (mol/m³/m = mol/m⁴)
        
    Returns:
        Molar flux J_A (mol/(m²·s))
    """
    J_A = -D_AB * dC_dz
    return J_A


def stefan_diffusion_flux(
    T_celsius: float,
    P_total_kPa: float,
    D_AB: float,
    z_film: float,
    P_A_surface_kPa: float,
    P_A_bulk_kPa: float = 0.0
) -> float:
    """
    Calculate molar flux for Stefan diffusion (evaporation through stagnant film).
    計算 Stefan 擴散的莫耳通量（經過靜止薄膜的蒸發）
    
    For unimolecular diffusion (A diffusing through stagnant B):
        N_A = (D_AB * P) / (R * T * z) * ln[(P - P_A∞) / (P - P_A0)]
    
    Args:
        T_celsius: Temperature (°C)
        P_total_kPa: Total pressure (kPa)
        D_AB: Diffusion coefficient (m²/s)
        z_film: Film thickness / diffusion path length (m)
        P_A_surface_kPa: Partial pressure of A at surface (saturation pressure)
        P_A_bulk_kPa: Partial pressure of A in bulk (usually ~0 for evaporation)
        
    Returns:
        Molar flux N_A (mol/(m²·s))
        
    Example:
        >>> # Ethanol evaporation at 25°C
        >>> D_ab = 1.18e-5  # m²/s
        >>> P_sat = 7.87  # kPa for ethanol at 25°C
        >>> N_A = stefan_diffusion_flux(25.0, 101.325, D_ab, 0.01, P_sat)
        >>> print(f"N_A = {N_A:.4e} mol/(m²·s)")
    """
    T_K = T_celsius + 273.15
    P = P_total_kPa * 1000  # Convert to Pa
    P_A0 = P_A_surface_kPa * 1000  # Pa
    P_A_inf = P_A_bulk_kPa * 1000  # Pa
    
    # Stefan equation
    # N_A = (D_AB * P) / (R * T * z) * ln[(P - P_A∞) / (P - P_A0)]
    
    # Check for valid pressure difference
    if P_A0 >= P:
        raise ValueError("Surface pressure cannot exceed total pressure")
    
    log_term = math.log((P - P_A_inf) / (P - P_A0))
    N_A = (D_AB * P) / (R_GAS * T_K * z_film) * log_term
    
    return N_A


def mass_flux_from_molar(N_A: float, mw: float) -> float:
    """
    Convert molar flux to mass flux.
    將莫耳通量轉換為質量通量
    
    Args:
        N_A: Molar flux (mol/(m²·s))
        mw: Molecular weight (g/mol)
        
    Returns:
        Mass flux (g/(m²·s))
    """
    return N_A * mw


# ============================================================
# Convenience Functions
# ============================================================

def calculate_ethanol_diffusion(T_celsius: float = 25.0, P_atm: float = 1.0) -> dict:
    """
    Calculate all diffusion parameters for ethanol in air.
    計算乙醇在空氣中的所有擴散參數
    
    Args:
        T_celsius: Temperature (°C)
        P_atm: Pressure (atm)
        
    Returns:
        Dictionary containing diffusion coefficient, vapor pressure, and flux estimates
    """
    from .properties import antoine_vapor_pressure, saturation_concentration
    
    # Diffusion coefficient
    D_AB = diffusion_coefficient_fsg(T_celsius, P_atm, ETHANOL)
    
    # Vapor pressure
    P_sat = antoine_vapor_pressure(T_celsius, ETHANOL)
    
    # Saturation concentration
    C_sat = saturation_concentration(T_celsius, ETHANOL)
    
    # Estimate Stefan flux with 1 cm film thickness
    z_film = 0.01  # m
    P_total = P_atm * 101.325  # kPa
    N_A = stefan_diffusion_flux(T_celsius, P_total, D_AB, z_film, P_sat, 0.0)
    
    return {
        "temperature_C": T_celsius,
        "pressure_atm": P_atm,
        "D_AB_m2s": D_AB,
        "P_sat_kPa": P_sat,
        "C_sat_mol_m3": C_sat,
        "N_A_mol_m2s": N_A,
        "mass_flux_g_m2s": mass_flux_from_molar(N_A, ETHANOL.mw),
    }


def print_diffusion_summary(T_celsius: float = 25.0, P_atm: float = 1.0):
    """
    Print summary of diffusion calculations for ethanol.
    印出乙醇擴散計算摘要
    """
    results = calculate_ethanol_diffusion(T_celsius, P_atm)
    
    print(f"\n{'='*55}")
    print(f"Ethanol-Air Diffusion Summary at {T_celsius}°C, {P_atm} atm")
    print(f"{'='*55}")
    print(f"Diffusion Coefficient D_AB: {results['D_AB_m2s']:.3e} m²/s")
    print(f"Vapor Pressure P_sat:       {results['P_sat_kPa']:.3f} kPa")
    print(f"Saturation Concentration:   {results['C_sat_mol_m3']:.3f} mol/m³")
    print(f"Stefan Flux (z=1cm):        {results['N_A_mol_m2s']:.4e} mol/(m²·s)")
    print(f"Mass Flux:                  {results['mass_flux_g_m2s']:.4e} g/(m²·s)")
    print(f"{'='*55}\n")


if __name__ == "__main__":
    # Demo calculations
    print_diffusion_summary(25.0, 1.0)
    print_diffusion_summary(30.0, 1.0)
