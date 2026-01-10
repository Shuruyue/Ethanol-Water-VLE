"""
Experimental Data Analysis Module
實驗數據分析模組

Functions for loading experimental data, calculating mass loss rates,
and generating publication-quality plots.
"""

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple
from pathlib import Path

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False
    
try:
    import matplotlib.pyplot as plt
    import matplotlib.ticker as ticker
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ============================================================
# Data Classes
# ============================================================

@dataclass
class ExperimentData:
    """
    Container for evaporation experiment data.
    
    Attributes:
        name: Experiment identifier (e.g., "natural_convection_run1")
        condition: "natural" or "forced"
        time_min: List of time points (minutes)
        mass_g: List of mass measurements (grams)
        temperature_C: Ambient temperature (°C)
        humidity_percent: Relative humidity (%)
        surface_area_m2: Evaporation surface area (m²)
        air_velocity_ms: Air velocity for forced convection (m/s), None for natural
        notes: Additional notes
    """
    name: str
    condition: str  # "natural" or "forced"
    time_min: List[float]
    mass_g: List[float]
    temperature_C: float
    humidity_percent: float
    surface_area_m2: float
    air_velocity_ms: Optional[float] = None
    notes: str = ""
    
    @property
    def initial_mass(self) -> float:
        """Return initial mass (g)."""
        return self.mass_g[0] if self.mass_g else 0.0
    
    @property
    def final_mass(self) -> float:
        """Return final mass (g)."""
        return self.mass_g[-1] if self.mass_g else 0.0
    
    @property
    def total_mass_loss(self) -> float:
        """Return total mass loss (g)."""
        return self.initial_mass - self.final_mass
    
    @property
    def duration_min(self) -> float:
        """Return experiment duration (minutes)."""
        return self.time_min[-1] - self.time_min[0] if len(self.time_min) > 1 else 0.0


@dataclass
class AnalysisResults:
    """
    Results from experimental data analysis.
    
    Attributes:
        experiment_name: Source experiment name
        avg_mass_loss_rate_g_min: Average mass loss rate (g/min)
        avg_mass_loss_rate_g_s: Average mass loss rate (g/s)
        calculated_km_ms: Calculated mass transfer coefficient (m/s)
        linear_fit_slope: Slope from linear regression (g/min)
        linear_fit_r2: R² coefficient of determination
    """
    experiment_name: str
    avg_mass_loss_rate_g_min: float
    avg_mass_loss_rate_g_s: float
    calculated_km_ms: float
    linear_fit_slope: Optional[float] = None
    linear_fit_r2: Optional[float] = None


# ============================================================
# Data Loading Functions
# ============================================================

def load_experiment_data(filepath: str) -> ExperimentData:
    """
    Load experiment data from CSV file.
    
    Expected CSV format:
        # name: experiment_name
        # condition: natural
        # temperature_C: 25.0
        # humidity_percent: 50.0
        # surface_area_m2: 0.005
        # air_velocity_ms: 1.5  (optional, for forced convection)
        time_min,mass_g
        0,15.00
        5,14.85
        ...
    
    Args:
        filepath: Path to CSV file
        
    Returns:
        ExperimentData object
    """
    path = Path(filepath)
    
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {filepath}")
    
    # Parse file
    metadata = {}
    time_data = []
    mass_data = []
    
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
                
            # Parse metadata comments
            if line.startswith('#'):
                parts = line[1:].split(':', 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip()
                    metadata[key] = value
                continue
            
            # Skip header
            if line.startswith('time'):
                continue
            
            # Parse data
            parts = line.split(',')
            if len(parts) >= 2:
                try:
                    time_data.append(float(parts[0]))
                    mass_data.append(float(parts[1]))
                except ValueError:
                    continue
    
    # Create ExperimentData object
    return ExperimentData(
        name=metadata.get('name', path.stem),
        condition=metadata.get('condition', 'unknown'),
        time_min=time_data,
        mass_g=mass_data,
        temperature_C=float(metadata.get('temperature_C', 25.0)),
        humidity_percent=float(metadata.get('humidity_percent', 50.0)),
        surface_area_m2=float(metadata.get('surface_area_m2', 0.005)),
        air_velocity_ms=float(metadata.get('air_velocity_ms')) if 'air_velocity_ms' in metadata else None,
        notes=metadata.get('notes', ''),
    )


def create_sample_data(condition: str = "natural") -> ExperimentData:
    """
    Create sample experiment data for testing.
    
    Args:
        condition: "natural" or "forced"
        
    Returns:
        ExperimentData with simulated values
    """
    # Simulated data - replace with actual measurements
    if condition == "natural":
        time = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
        mass = [15.00, 14.92, 14.84, 14.76, 14.68, 14.60, 14.52, 
                14.44, 14.36, 14.28, 14.20, 14.12, 14.04]
        velocity = None
    else:
        time = [0, 5, 10, 15, 20, 25, 30]
        mass = [15.00, 14.80, 14.60, 14.40, 14.20, 14.00, 13.80]
        velocity = 1.5
    
    return ExperimentData(
        name=f"sample_{condition}",
        condition=condition,
        time_min=time,
        mass_g=mass,
        temperature_C=25.0,
        humidity_percent=50.0,
        surface_area_m2=0.005,
        air_velocity_ms=velocity,
        notes="Sample data for testing",
    )


# ============================================================
# Analysis Functions
# ============================================================

def calculate_mass_loss_rate(data: ExperimentData, method: str = "average") -> float:
    """
    Calculate mass loss rate from experimental data.
    
    Args:
        data: ExperimentData object
        method: "average" for simple average, "linear_fit" for regression
        
    Returns:
        Mass loss rate in g/min
    """
    if len(data.time_min) < 2:
        return 0.0
    
    if method == "average":
        # Simple average rate
        delta_m = data.initial_mass - data.final_mass
        delta_t = data.duration_min
        return delta_m / delta_t if delta_t > 0 else 0.0
    
    elif method == "linear_fit" and HAS_NUMPY:
        # Linear regression
        t = np.array(data.time_min)
        m = np.array(data.mass_g)
        
        # Fit: m = m0 - rate * t
        coeffs = np.polyfit(t, m, 1)
        rate = -coeffs[0]  # Negative slope = positive loss rate
        return rate
    
    else:
        return calculate_mass_loss_rate(data, method="average")


def linear_regression(x: List[float], y: List[float]) -> Tuple[float, float, float]:
    """
    Perform simple linear regression without numpy.
    
    Args:
        x: Independent variable values
        y: Dependent variable values
        
    Returns:
        Tuple of (slope, intercept, R²)
    """
    n = len(x)
    if n < 2:
        return 0.0, 0.0, 0.0
    
    # Calculate means
    x_mean = sum(x) / n
    y_mean = sum(y) / n
    
    # Calculate slope and intercept
    numerator = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y))
    denominator = sum((xi - x_mean) ** 2 for xi in x)
    
    if denominator == 0:
        return 0.0, y_mean, 0.0
    
    slope = numerator / denominator
    intercept = y_mean - slope * x_mean
    
    # Calculate R²
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    ss_res = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
    
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0
    
    return slope, intercept, r_squared


def analyze_experiment(data: ExperimentData) -> AnalysisResults:
    """
    Perform complete analysis of experiment data.
    
    Args:
        data: ExperimentData object
        
    Returns:
        AnalysisResults with calculated parameters
    """
    from .mass_transfer_coeff import mass_transfer_from_experiment
    from .properties import ETHANOL
    
    # Calculate mass loss rate
    avg_rate_g_min = calculate_mass_loss_rate(data, method="average")
    avg_rate_g_s = avg_rate_g_min / 60.0
    
    # Linear regression
    slope, intercept, r2 = linear_regression(data.time_min, data.mass_g)
    
    # Calculate mass transfer coefficient
    km = mass_transfer_from_experiment(
        mass_loss_rate=avg_rate_g_s,
        surface_area=data.surface_area_m2,
        T_celsius=data.temperature_C,
        props=ETHANOL,
    )
    
    return AnalysisResults(
        experiment_name=data.name,
        avg_mass_loss_rate_g_min=avg_rate_g_min,
        avg_mass_loss_rate_g_s=avg_rate_g_s,
        calculated_km_ms=km,
        linear_fit_slope=-slope,  # Convert to positive loss rate
        linear_fit_r2=r2,
    )


# ============================================================
# Plotting Functions
# ============================================================

def plot_mass_vs_time(
    data: ExperimentData,
    output_path: Optional[str] = None,
    show_fit: bool = True,
    figsize: Tuple[float, float] = (8, 6)
) -> Optional[object]:
    """
    Plot mass vs time for a single experiment.
    
    Args:
        data: ExperimentData object
        output_path: Path to save figure (None to display)
        show_fit: Whether to show linear fit
        figsize: Figure size (width, height) in inches
        
    Returns:
        Matplotlib figure object if available, None otherwise
    """
    if not HAS_MATPLOTLIB:
        print("Warning: matplotlib not available, skipping plot")
        return None
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Plot data points
    ax.plot(data.time_min, data.mass_g, 'o-', 
            markersize=8, linewidth=1.5, 
            label='Experimental data')
    
    # Add linear fit
    if show_fit and len(data.time_min) >= 2:
        slope, intercept, r2 = linear_regression(data.time_min, data.mass_g)
        t_fit = [data.time_min[0], data.time_min[-1]]
        m_fit = [slope * t + intercept for t in t_fit]
        ax.plot(t_fit, m_fit, 'r--', linewidth=1.5,
                label=f'Linear fit (R² = {r2:.4f})')
    
    # Labels and formatting
    ax.set_xlabel('Time (min)', fontsize=12)
    ax.set_ylabel('Mass (g)', fontsize=12)
    ax.set_title(f'Mass Loss vs Time - {data.name}', fontsize=14)
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)
    
    # Add annotations
    rate = calculate_mass_loss_rate(data)
    ax.annotate(f'Evaporation rate: {rate:.4f} g/min',
                xy=(0.05, 0.95), xycoords='axes fraction',
                fontsize=10, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {output_path}")
    
    return fig


def compare_conditions(
    natural_data: ExperimentData,
    forced_data: ExperimentData,
    output_path: Optional[str] = None,
    figsize: Tuple[float, float] = (10, 6)
) -> Optional[object]:
    """
    Compare natural and forced convection experiments.
    
    Args:
        natural_data: Natural convection experiment data
        forced_data: Forced convection experiment data
        output_path: Path to save figure
        figsize: Figure size
        
    Returns:
        Matplotlib figure object if available
    """
    if not HAS_MATPLOTLIB:
        print("Warning: matplotlib not available, skipping plot")
        return None
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    
    # Left plot: Mass vs Time comparison
    ax1.plot(natural_data.time_min, natural_data.mass_g, 'o-',
             markersize=6, label='Natural convection')
    ax1.plot(forced_data.time_min, forced_data.mass_g, 's-',
             markersize=6, label='Forced convection')
    
    ax1.set_xlabel('Time (min)', fontsize=12)
    ax1.set_ylabel('Mass (g)', fontsize=12)
    ax1.set_title('Mass Loss Comparison', fontsize=14)
    ax1.legend(loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Right plot: Evaporation rate bar chart
    rate_natural = calculate_mass_loss_rate(natural_data)
    rate_forced = calculate_mass_loss_rate(forced_data)
    
    bars = ax2.bar(['Natural\nConvection', 'Forced\nConvection'],
                   [rate_natural, rate_forced],
                   color=['steelblue', 'coral'],
                   edgecolor='black')
    
    ax2.set_ylabel('Evaporation Rate (g/min)', fontsize=12)
    ax2.set_title('Evaporation Rate Comparison', fontsize=14)
    
    # Add value labels on bars
    for bar, rate in zip(bars, [rate_natural, rate_forced]):
        height = bar.get_height()
        ax2.annotate(f'{rate:.4f}',
                     xy=(bar.get_x() + bar.get_width() / 2, height),
                     ha='center', va='bottom', fontsize=11)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {output_path}")
    
    return fig


def print_analysis_summary(results: AnalysisResults):
    """Print analysis results summary."""
    print(f"\n{'='*55}")
    print(f"Analysis Results: {results.experiment_name}")
    print(f"{'='*55}")
    print(f"Average evaporation rate: {results.avg_mass_loss_rate_g_min:.4f} g/min")
    print(f"                        : {results.avg_mass_loss_rate_g_s:.4e} g/s")
    print(f"Mass transfer coeff k_m: {results.calculated_km_ms:.4e} m/s")
    if results.linear_fit_r2 is not None:
        print(f"Linear fit R²:          {results.linear_fit_r2:.4f}")
    print(f"{'='*55}\n")


# ============================================================
# Main Demo
# ============================================================

if __name__ == "__main__":
    # Create sample data
    natural = create_sample_data("natural")
    forced = create_sample_data("forced")
    
    # Analyze
    results_nat = analyze_experiment(natural)
    results_for = analyze_experiment(forced)
    
    print_analysis_summary(results_nat)
    print_analysis_summary(results_for)
    
    # Plot if matplotlib available
    if HAS_MATPLOTLIB:
        compare_conditions(natural, forced)
        plt.show()
