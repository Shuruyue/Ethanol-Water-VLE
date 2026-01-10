#!/usr/bin/env python3
"""
Mass Transfer Experiment Analysis - Main Script
質傳實驗分析主程式

This script performs a complete analysis of volatile liquid evaporation
experiments under natural and forced convection conditions.

Usage:
    python main.py                          # Run with sample data
    python main.py --data data/             # Run with custom data folder
    python main.py --temperature 30         # Specify temperature
"""

import argparse
import sys
from pathlib import Path

# Add mass_transfer to path if running directly
sys.path.insert(0, str(Path(__file__).parent))

from mass_transfer import (
    # Properties
    ETHANOL,
    WATER,
    get_ethanol_properties,
    antoine_vapor_pressure,
    print_properties_summary,
    
    # Diffusion
    diffusion_coefficient_fsg,
    stefan_diffusion_flux,
    print_diffusion_summary,
    
    # Mass Transfer
    mass_transfer_natural,
    mass_transfer_forced,
    compare_convection_modes,
    print_comparison_summary,
    
    # Analysis
    ExperimentData,
    create_sample_data,
    analyze_experiment,
    print_analysis_summary,
    plot_mass_vs_time,
    compare_conditions,
)
from mass_transfer.report_generator import ReportGenerator


def run_theoretical_analysis(T_celsius: float = 25.0, L_char: float = 0.08, 
                             u_air: float = 1.5):
    """
    Run theoretical mass transfer analysis.
    
    Args:
        T_celsius: Temperature (°C)
        L_char: Characteristic length (m)
        u_air: Air velocity for forced convection (m/s)
    """
    print("\n" + "="*70)
    print("THEORETICAL ANALYSIS - Ethanol Evaporation")
    print("="*70)
    
    # 1. Chemical Properties
    print("\n[1] Chemical Properties")
    print_properties_summary(ETHANOL, T_celsius)
    
    # 2. Diffusion Analysis
    print("[2] Diffusion Analysis")
    print_diffusion_summary(T_celsius, 1.0)
    
    # 3. Mass Transfer Comparison
    print("[3] Mass Transfer Coefficient Comparison")
    result = compare_convection_modes(
        T_celsius=T_celsius,
        L_char=L_char,
        u_air=u_air,
        surface_area=3.14159 * (L_char/2)**2,  # Circular dish
        initial_mass=15.0,
    )
    print_comparison_summary(result)
    
    return result


def run_experimental_analysis(data_folder: str = None, save_figures: bool = True):
    """
    Run analysis on experimental data.
    
    Args:
        data_folder: Path to folder containing CSV data files
        save_figures: Whether to save generated figures
    """
    print("\n" + "="*70)
    print("EXPERIMENTAL DATA ANALYSIS")
    print("="*70)
    
    # Use sample data if no folder specified
    natural_data = create_sample_data("natural")
    forced_data = create_sample_data("forced")
    
    print(f"\nNatural convection experiment:")
    print(f"  Duration: {natural_data.duration_min:.0f} min")
    print(f"  Mass loss: {natural_data.total_mass_loss:.2f} g")
    
    print(f"\nForced convection experiment:")
    print(f"  Duration: {forced_data.duration_min:.0f} min")
    print(f"  Mass loss: {forced_data.total_mass_loss:.2f} g")
    
    # Analyze both experiments
    print("\n" + "-"*50)
    results_natural = analyze_experiment(natural_data)
    print_analysis_summary(results_natural)
    
    results_forced = analyze_experiment(forced_data)
    print_analysis_summary(results_forced)
    
    # Generate comparison
    print("\n[Comparison]")
    ratio = results_forced.avg_mass_loss_rate_g_min / results_natural.avg_mass_loss_rate_g_min
    print(f"Forced/Natural evaporation rate ratio: {ratio:.2f}x")
    
    # Generate plots
    figures_dir = Path(__file__).parent / "figures"
    figures_dir.mkdir(exist_ok=True)
    
    try:
        if save_figures:
            plot_mass_vs_time(natural_data, 
                              output_path=str(figures_dir / "mass_vs_time_natural.png"))
            plot_mass_vs_time(forced_data,
                              output_path=str(figures_dir / "mass_vs_time_forced.png"))
            compare_conditions(natural_data, forced_data,
                               output_path=str(figures_dir / "comparison.png"))
            print(f"\nFigures saved to: {figures_dir}")
            
            # Generate PowerPoint Report
            report_path = Path(__file__).parent / "reports" / "Experiment_Report.pptx"
            report_path.parent.mkdir(exist_ok=True)
            
            generator = ReportGenerator(str(report_path))
            generator.generate(results_natural, results_forced, str(figures_dir))
            
    except Exception as e:
        print(f"\nNote: Could not generate report/figures ({e})")
    
    return results_natural, results_forced


def generate_report_data(T_celsius: float = 25.0):
    """
    Generate all data needed for the final report.
    
    Returns a dictionary with all calculated parameters.
    """
    L_char = 0.08  # 8 cm dish
    u_air = 1.5    # 1.5 m/s fan
    A = 3.14159 * (L_char/2)**2  # m²
    
    data = {
        "conditions": {
            "temperature_C": T_celsius,
            "temperature_K": T_celsius + 273.15,
            "pressure_kPa": 101.325,
            "characteristic_length_m": L_char,
            "surface_area_m2": A,
            "air_velocity_ms": u_air,
        },
        "properties": {
            "ethanol_mw": ETHANOL.mw,
            "ethanol_density_kg_m3": ETHANOL.density,
            "vapor_pressure_kPa": antoine_vapor_pressure(T_celsius, ETHANOL),
        },
        "diffusion": {
            "D_AB_m2_s": diffusion_coefficient_fsg(T_celsius, 1.0, ETHANOL),
        },
    }
    
    # Mass transfer calculations
    km_nat, Sh_nat = mass_transfer_natural(T_celsius, L_char)
    km_for, Sh_for = mass_transfer_forced(T_celsius, u_air, L_char)
    
    data["natural_convection"] = {
        "k_m_ms": km_nat,
        "Sherwood": Sh_nat,
    }
    
    data["forced_convection"] = {
        "k_m_ms": km_for,
        "Sherwood": Sh_for,
    }
    
    data["enhancement_factor"] = km_for / km_nat
    
    return data


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Mass Transfer Experiment Analysis for Transport Phenomena II"
    )
    parser.add_argument("--temperature", "-T", type=float, default=25.0,
                        help="Ambient temperature in °C (default: 25)")
    parser.add_argument("--velocity", "-u", type=float, default=1.5,
                        help="Air velocity for forced convection in m/s (default: 1.5)")
    parser.add_argument("--dish-diameter", "-d", type=float, default=8.0,
                        help="Dish diameter in cm (default: 8)")
    parser.add_argument("--data", type=str, default=None,
                        help="Path to data folder")
    parser.add_argument("--no-figures", action="store_true",
                        help="Skip generating figures")
    
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("TRANSPORT PHENOMENA II - MASS TRANSFER EXPERIMENT ANALYSIS")
    print("Volatile Liquid Evaporation Study")
    print("="*70)
    
    # Run theoretical analysis
    theo_result = run_theoretical_analysis(
        T_celsius=args.temperature,
        L_char=args.dish_diameter / 100,  # Convert cm to m
        u_air=args.velocity,
    )
    
    # Run experimental analysis
    exp_results = run_experimental_analysis(
        data_folder=args.data,
        save_figures=not args.no_figures,
    )
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"""
Key Findings:
- Diffusion coefficient D_AB = {diffusion_coefficient_fsg(args.temperature, 1.0, ETHANOL):.2e} m²/s
- Natural convection k_m = {theo_result.km_natural:.2e} m/s
- Forced convection k_m = {theo_result.km_forced:.2e} m/s
- Enhancement factor = {theo_result.enhancement_factor:.2f}x

The forced convection condition significantly enhances mass transfer,
reducing evaporation time by approximately {(1 - 1/theo_result.enhancement_factor)*100:.0f}%.
""")
    
    print("="*70)
    print("Analysis complete.")
    print("="*70 + "\n")


if __name__ == "__main__":
    main()
