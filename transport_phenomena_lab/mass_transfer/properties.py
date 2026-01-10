"""
Physical Properties Module
物理性質模組

Contains chemical properties, Antoine equation for vapor pressure,
and other thermodynamic data needed for mass transfer calculations.
"""

from dataclasses import dataclass
from typing import Optional
import math


# ============================================================
# Constants
# ============================================================

R_GAS = 8.314  # J/(mol·K) - Universal gas constant
R_GAS_ATM = 0.08206  # L·atm/(mol·K)


# ============================================================
# Data Classes
# ============================================================

@dataclass
class ChemicalProperties:
    """
    Chemical properties data class.
    化學物質性質資料類別
    
    Attributes:
        name: Chemical name
        formula: Chemical formula
        mw: Molecular weight (g/mol)
        density: Liquid density at 25°C (kg/m³)
        diffusion_volume: Fuller diffusion volume (dimensionless)
        antoine_A: Antoine equation A parameter (log10, mmHg)
        antoine_B: Antoine equation B parameter
        antoine_C: Antoine equation C parameter
    """
    name: str
    formula: str
    mw: float  # g/mol
    density: float  # kg/m³
    diffusion_volume: float  # Fuller volume
    antoine_A: float
    antoine_B: float
    antoine_C: float
    

@dataclass
class ExperimentConditions:
    """
    Experiment conditions data class.
    實驗條件資料類別
    
    Attributes:
        temperature: Ambient temperature (°C)
        pressure: Atmospheric pressure (kPa)
        humidity: Relative humidity (%)
        surface_area: Evaporation surface area (m²)
        initial_mass: Initial liquid mass (g)
    """
    temperature: float  # °C
    pressure: float = 101.325  # kPa
    humidity: float = 50.0  # %
    surface_area: float = 0.005  # m² (default ~8cm diameter dish)
    initial_mass: float = 10.0  # g


# ============================================================
# Chemical Property Database
# ============================================================

# Ethanol (乙醇) properties
ETHANOL = ChemicalProperties(
    name="Ethanol",
    formula="C2H5OH",
    mw=46.07,  # g/mol
    density=789.0,  # kg/m³ at 25°C
    diffusion_volume=51.77,  # Fuller volume
    # Antoine coefficients (T in °C, P in mmHg)
    # log10(P) = A - B/(T + C)
    antoine_A=8.20417,
    antoine_B=1642.89,
    antoine_C=230.300,
)

# Water (水) properties
WATER = ChemicalProperties(
    name="Water",
    formula="H2O",
    mw=18.015,  # g/mol
    density=997.0,  # kg/m³ at 25°C
    diffusion_volume=12.7,  # Fuller volume
    # Antoine coefficients
    antoine_A=8.07131,
    antoine_B=1730.63,
    antoine_C=233.426,
)

# Air (空氣) properties for diffusion calculations
AIR_MW = 28.97  # g/mol
AIR_DIFFUSION_VOLUME = 20.1  # Fuller volume


# ============================================================
# Property Functions
# ============================================================

def get_ethanol_properties() -> ChemicalProperties:
    """Return ethanol properties. 取得乙醇性質"""
    return ETHANOL


def get_water_properties() -> ChemicalProperties:
    """Return water properties. 取得水性質"""
    return WATER


def antoine_vapor_pressure(T_celsius: float, props: ChemicalProperties) -> float:
    """
    Calculate vapor pressure using Antoine equation.
    使用 Antoine 方程計算蒸氣壓
    
    Args:
        T_celsius: Temperature in degrees Celsius (°C)
        props: ChemicalProperties with Antoine coefficients
        
    Returns:
        Vapor pressure in kPa
        
    Example:
        >>> props = get_ethanol_properties()
        >>> P_sat = antoine_vapor_pressure(25.0, props)
        >>> print(f"Ethanol P_sat at 25°C: {P_sat:.2f} kPa")
    """
    # Antoine equation: log10(P_mmHg) = A - B/(T + C)
    log_P_mmHg = props.antoine_A - props.antoine_B / (T_celsius + props.antoine_C)
    P_mmHg = 10 ** log_P_mmHg
    
    # Convert mmHg to kPa (1 mmHg = 0.133322 kPa)
    P_kPa = P_mmHg * 0.133322
    
    return P_kPa


def saturation_concentration(T_celsius: float, props: ChemicalProperties, 
                            P_total_kPa: float = 101.325) -> float:
    """
    Calculate saturation concentration of vapor at liquid surface.
    計算液體表面的飽和濃度
    
    Args:
        T_celsius: Temperature (°C)
        props: Chemical properties
        P_total_kPa: Total pressure (kPa)
        
    Returns:
        Saturation concentration (mol/m³)
    """
    T_K = T_celsius + 273.15
    P_sat_kPa = antoine_vapor_pressure(T_celsius, props)
    
    # C = P_sat / (R * T)
    # Convert kPa to Pa for SI units
    P_sat_Pa = P_sat_kPa * 1000
    C_sat = P_sat_Pa / (R_GAS * T_K)  # mol/m³
    
    return C_sat


def mole_fraction_from_humidity(T_celsius: float, humidity_percent: float,
                                P_total_kPa: float = 101.325) -> float:
    """
    Calculate water vapor mole fraction from relative humidity.
    從相對濕度計算水蒸氣莫耳分率
    
    Args:
        T_celsius: Temperature (°C)
        humidity_percent: Relative humidity (%)
        P_total_kPa: Total atmospheric pressure (kPa)
        
    Returns:
        Water vapor mole fraction in air
    """
    P_sat = antoine_vapor_pressure(T_celsius, WATER)
    P_water = (humidity_percent / 100.0) * P_sat
    y_water = P_water / P_total_kPa
    
    return y_water


def mixture_density(ethanol_mass_fraction: float, T_celsius: float = 25.0) -> float:
    """
    Estimate density of ethanol-water mixture.
    估算乙醇-水混合液密度
    
    Args:
        ethanol_mass_fraction: Mass fraction of ethanol (0-1)
        T_celsius: Temperature (°C)
        
    Returns:
        Mixture density (kg/m³)
        
    Note:
        Uses simple linear mixing rule; actual mixture shows volume contraction.
    """
    w1 = ethanol_mass_fraction
    w2 = 1 - w1
    
    rho1 = ETHANOL.density
    rho2 = WATER.density
    
    # Simple mixing rule (approximate)
    rho_mix = 1 / (w1/rho1 + w2/rho2)
    
    return rho_mix


# ============================================================
# Convenience Functions
# ============================================================

def print_properties_summary(props: ChemicalProperties, T_celsius: float = 25.0):
    """
    Print a summary of chemical properties.
    印出化學性質摘要
    """
    P_sat = antoine_vapor_pressure(T_celsius, props)
    C_sat = saturation_concentration(T_celsius, props)
    
    print(f"\n{'='*50}")
    print(f"Chemical: {props.name} ({props.formula})")
    print(f"{'='*50}")
    print(f"Molecular Weight:    {props.mw:.2f} g/mol")
    print(f"Liquid Density:      {props.density:.1f} kg/m³")
    print(f"Diffusion Volume:    {props.diffusion_volume:.2f}")
    print(f"At T = {T_celsius}°C:")
    print(f"  Vapor Pressure:    {P_sat:.3f} kPa ({P_sat*7.5:.1f} mmHg)")
    print(f"  Sat. Conc.:        {C_sat:.3f} mol/m³")
    print(f"{'='*50}\n")


if __name__ == "__main__":
    # Demo: Print properties for ethanol and water
    print_properties_summary(ETHANOL, 25.0)
    print_properties_summary(WATER, 25.0)
    
    # Calculate for 75% ethanol solution
    print("75% Ethanol Solution:")
    print(f"  Mixture density: {mixture_density(0.75):.1f} kg/m³")
