"""
Mass Transfer Coefficient Module
質傳係數模組

Calculate mass transfer coefficients for natural and forced convection,
evaporation rates, and total evaporation time.
"""

import math
from typing import Tuple, Optional
from dataclasses import dataclass

from .properties import (
    ChemicalProperties,
    R_GAS,
    ETHANOL,
    antoine_vapor_pressure,
    saturation_concentration,
)
from .diffusion import diffusion_coefficient_fsg


# ============================================================
# Dimensionless Numbers
# ============================================================

def schmidt_number(nu: float, D_AB: float) -> float:
    """
    Calculate Schmidt number.
    計算 Schmidt 數
    
    Sc = ν / D_AB = momentum diffusivity / mass diffusivity
    
    Args:
        nu: Kinematic viscosity of gas (m²/s)
        D_AB: Binary diffusion coefficient (m²/s)
        
    Returns:
        Schmidt number (dimensionless)
    """
    return nu / D_AB


def reynolds_number(u: float, L: float, nu: float) -> float:
    """
    Calculate Reynolds number.
    計算 Reynolds 數
    
    Re = u * L / ν
    
    Args:
        u: Velocity (m/s)
        L: Characteristic length (m)
        nu: Kinematic viscosity (m²/s)
        
    Returns:
        Reynolds number (dimensionless)
    """
    return u * L / nu


def sherwood_number_from_km(k_m: float, L: float, D_AB: float) -> float:
    """
    Calculate Sherwood number from mass transfer coefficient.
    從質傳係數計算 Sherwood 數
    
    Sh = k_m * L / D_AB
    
    Args:
        k_m: Mass transfer coefficient (m/s)
        L: Characteristic length (m)
        D_AB: Diffusion coefficient (m²/s)
        
    Returns:
        Sherwood number (dimensionless)
    """
    return k_m * L / D_AB


# ============================================================
# Mass Transfer Coefficients
# ============================================================

NU_AIR_25C = 1.56e-5  # Kinematic viscosity of air at 25°C (m²/s)


def mass_transfer_natural(
    T_celsius: float,
    L_char: float,
    props: ChemicalProperties = ETHANOL,
    nu_air: float = NU_AIR_25C,
    g: float = 9.81
) -> Tuple[float, float]:
    """
    Calculate mass transfer coefficient for natural convection.
    計算自然對流質傳係數
    
    For horizontal surface facing up (evaporation):
        Sh = 0.54 * (Gr * Sc)^0.25   for 10^4 < Gr*Sc < 10^7
        Sh = 0.15 * (Gr * Sc)^0.33   for 10^7 < Gr*Sc < 10^11
    
    Args:
        T_celsius: Temperature (°C)
        L_char: Characteristic length (m), typically plate width or diameter
        props: Chemical properties of evaporating species
        nu_air: Kinematic viscosity of air (m²/s)
        g: Gravitational acceleration (m/s²)
        
    Returns:
        Tuple of (k_m in m/s, Sherwood number)
        
    Example:
        >>> k_m, Sh = mass_transfer_natural(25.0, 0.08, ETHANOL)
        >>> print(f"k_m = {k_m:.4e} m/s, Sh = {Sh:.2f}")
    """
    # Calculate diffusion coefficient
    D_AB = diffusion_coefficient_fsg(T_celsius, 1.0, props)
    
    # Schmidt number
    Sc = schmidt_number(nu_air, D_AB)
    
    # Estimate Grashof number for mass transfer
    # Gr = g * L³ * Δρ/ρ / ν²
    # For dilute vapor, approximate Δρ/ρ from concentration difference
    P_sat = antoine_vapor_pressure(T_celsius, props)
    P_total = 101.325  # kPa
    y_A = P_sat / P_total  # Mole fraction at surface
    
    # Density ratio approximation (vapor lighter than air)
    delta_rho_ratio = y_A * (1 - props.mw / 28.97)
    delta_rho_ratio = abs(delta_rho_ratio)
    
    Gr = g * (L_char ** 3) * delta_rho_ratio / (nu_air ** 2)
    
    # Rayleigh-like number for mass transfer
    Gr_Sc = Gr * Sc
    
    # Choose correlation based on Gr*Sc
    if Gr_Sc < 1e4:
        # Very low natural convection - use limiting Sh
        Sh = 2.0
    elif Gr_Sc < 1e7:
        # Laminar natural convection
        Sh = 0.54 * (Gr_Sc ** 0.25)
    else:
        # Turbulent natural convection
        Sh = 0.15 * (Gr_Sc ** 0.33)
    
    # Mass transfer coefficient
    k_m = Sh * D_AB / L_char
    
    return k_m, Sh


def mass_transfer_forced(
    T_celsius: float,
    u_air: float,
    L_char: float,
    props: ChemicalProperties = ETHANOL,
    nu_air: float = NU_AIR_25C
) -> Tuple[float, float]:
    """
    Calculate mass transfer coefficient for forced convection.
    計算強制對流質傳係數
    
    For flow over flat plate:
        Laminar (Re < 5×10^5):  Sh = 0.664 * Re^0.5 * Sc^0.33
        Turbulent:              Sh = 0.037 * Re^0.8 * Sc^0.33
    
    Args:
        T_celsius: Temperature (°C)
        u_air: Air velocity (m/s)
        L_char: Characteristic length (m)
        props: Chemical properties of evaporating species
        nu_air: Kinematic viscosity of air (m²/s)
        
    Returns:
        Tuple of (k_m in m/s, Sherwood number)
        
    Example:
        >>> k_m, Sh = mass_transfer_forced(25.0, 1.0, 0.08, ETHANOL)
        >>> print(f"k_m = {k_m:.4e} m/s, Sh = {Sh:.2f}")
    """
    # Calculate diffusion coefficient
    D_AB = diffusion_coefficient_fsg(T_celsius, 1.0, props)
    
    # Dimensionless numbers
    Sc = schmidt_number(nu_air, D_AB)
    Re = reynolds_number(u_air, L_char, nu_air)
    
    # Sherwood number correlation
    Re_crit = 5e5
    
    if Re < Re_crit:
        # Laminar flow
        Sh = 0.664 * (Re ** 0.5) * (Sc ** (1/3))
    else:
        # Turbulent flow (combined laminar-turbulent)
        Sh = 0.037 * (Re ** 0.8) * (Sc ** (1/3))
    
    # Mass transfer coefficient
    k_m = Sh * D_AB / L_char
    
    return k_m, Sh


def mass_transfer_from_experiment(
    mass_loss_rate: float,
    surface_area: float,
    T_celsius: float,
    props: ChemicalProperties = ETHANOL
) -> float:
    """
    Calculate mass transfer coefficient from experimental data.
    從實驗數據計算質傳係數
    
    k_m = ṁ / (A * M * (C_s - C_∞))
    
    Assuming C_∞ ≈ 0 (well-ventilated conditions)
    
    Args:
        mass_loss_rate: Mass loss rate (g/s)
        surface_area: Evaporation surface area (m²)
        T_celsius: Temperature (°C)
        props: Chemical properties
        
    Returns:
        Mass transfer coefficient (m/s)
    """
    C_sat = saturation_concentration(T_celsius, props)
    C_bulk = 0.0  # Assume zero concentration in bulk
    
    # Convert mass rate to molar rate
    molar_rate = mass_loss_rate / props.mw  # mol/s
    
    # Molar flux
    N_A = molar_rate / surface_area  # mol/(m²·s)
    
    # Mass transfer coefficient
    k_m = N_A / (C_sat - C_bulk)  # m/s
    
    return k_m


# ============================================================
# Evaporation Rate and Time
# ============================================================

def evaporation_rate(
    k_m: float,
    surface_area: float,
    T_celsius: float,
    props: ChemicalProperties = ETHANOL,
    C_bulk: float = 0.0
) -> float:
    """
    Calculate evaporation rate.
    計算蒸發速率
    
    ṁ = k_m * A * M * (C_s - C_∞)
    
    Args:
        k_m: Mass transfer coefficient (m/s)
        surface_area: Evaporation area (m²)
        T_celsius: Temperature (°C)
        props: Chemical properties
        C_bulk: Bulk concentration (mol/m³), default 0
        
    Returns:
        Mass evaporation rate (g/s)
    """
    C_sat = saturation_concentration(T_celsius, props)
    
    # Molar flux
    N_A = k_m * (C_sat - C_bulk)  # mol/(m²·s)
    
    # Mass rate
    m_dot = N_A * surface_area * props.mw  # g/s
    
    return m_dot


def evaporation_time(
    initial_mass: float,
    k_m: float,
    surface_area: float,
    T_celsius: float,
    props: ChemicalProperties = ETHANOL
) -> float:
    """
    Estimate total evaporation time (assuming constant rate).
    估算完全蒸發時間（假設恆定速率）
    
    t = m_0 / ṁ
    
    Note: This is a simplified estimate. Actual evaporation slows
    as liquid level drops (increasing diffusion path).
    
    Args:
        initial_mass: Initial liquid mass (g)
        k_m: Mass transfer coefficient (m/s)
        surface_area: Evaporation area (m²)
        T_celsius: Temperature (°C)
        props: Chemical properties
        
    Returns:
        Estimated evaporation time (seconds)
    """
    m_dot = evaporation_rate(k_m, surface_area, T_celsius, props)
    
    if m_dot <= 0:
        return float('inf')
    
    t_total = initial_mass / m_dot
    
    return t_total


# ============================================================
# Comparison Functions
# ============================================================

@dataclass
class ConvectionComparison:
    """Results of natural vs forced convection comparison."""
    T_celsius: float
    L_char: float
    u_air: float
    
    # Natural convection
    km_natural: float
    Sh_natural: float
    rate_natural: float
    time_natural: float
    
    # Forced convection  
    km_forced: float
    Sh_forced: float
    rate_forced: float
    time_forced: float
    
    # Ratio
    enhancement_factor: float


def compare_convection_modes(
    T_celsius: float = 25.0,
    L_char: float = 0.08,
    u_air: float = 1.0,
    surface_area: float = 0.005,
    initial_mass: float = 10.0,
    props: ChemicalProperties = ETHANOL
) -> ConvectionComparison:
    """
    Compare natural and forced convection mass transfer.
    比較自然對流與強制對流的質傳
    
    Args:
        T_celsius: Temperature (°C)
        L_char: Characteristic length (m)
        u_air: Air velocity for forced convection (m/s)
        surface_area: Evaporation surface area (m²)
        initial_mass: Initial liquid mass (g)
        props: Chemical properties
        
    Returns:
        ConvectionComparison dataclass with all results
    """
    # Natural convection
    km_nat, Sh_nat = mass_transfer_natural(T_celsius, L_char, props)
    rate_nat = evaporation_rate(km_nat, surface_area, T_celsius, props)
    time_nat = evaporation_time(initial_mass, km_nat, surface_area, T_celsius, props)
    
    # Forced convection
    km_for, Sh_for = mass_transfer_forced(T_celsius, u_air, L_char, props)
    rate_for = evaporation_rate(km_for, surface_area, T_celsius, props)
    time_for = evaporation_time(initial_mass, km_for, surface_area, T_celsius, props)
    
    # Enhancement factor
    enhancement = km_for / km_nat if km_nat > 0 else float('inf')
    
    return ConvectionComparison(
        T_celsius=T_celsius,
        L_char=L_char,
        u_air=u_air,
        km_natural=km_nat,
        Sh_natural=Sh_nat,
        rate_natural=rate_nat,
        time_natural=time_nat,
        km_forced=km_for,
        Sh_forced=Sh_for,
        rate_forced=rate_for,
        time_forced=time_for,
        enhancement_factor=enhancement,
    )


def print_comparison_summary(result: ConvectionComparison):
    """
    Print comparison summary.
    印出比較摘要
    """
    print(f"\n{'='*60}")
    print(f"Natural vs Forced Convection Comparison")
    print(f"T = {result.T_celsius}°C, L = {result.L_char*100:.1f} cm, u = {result.u_air} m/s")
    print(f"{'='*60}")
    print(f"\n{'Parameter':<25} {'Natural':<15} {'Forced':<15}")
    print(f"{'-'*55}")
    print(f"{'k_m (m/s)':<25} {result.km_natural:.4e}    {result.km_forced:.4e}")
    print(f"{'Sherwood number':<25} {result.Sh_natural:.2f}          {result.Sh_forced:.2f}")
    print(f"{'Evap. rate (g/min)':<25} {result.rate_natural*60:.4f}       {result.rate_forced*60:.4f}")
    print(f"{'Evap. time (min)':<25} {result.time_natural/60:.1f}          {result.time_forced/60:.1f}")
    print(f"\nEnhancement factor: {result.enhancement_factor:.2f}x")
    print(f"{'='*60}\n")


if __name__ == "__main__":
    # Demo: Compare natural vs forced convection
    result = compare_convection_modes(
        T_celsius=25.0,
        L_char=0.08,  # 8 cm dish
        u_air=1.5,    # 1.5 m/s fan
        surface_area=0.005,  # ~8 cm diameter
        initial_mass=15.0,   # 15 g
    )
    print_comparison_summary(result)
