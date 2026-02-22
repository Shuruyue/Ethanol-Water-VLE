#!/usr/bin/env python3
"""Generate a practical output image pack for engineering use."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from analysis import compute_gamma_curves, compute_isothermal_excess_curves, find_azeotrope
from models import gamma_nrtl, psat_kpa
from parameter_store import load_system_parameters
from pipeline import RunResult, run_analysis
from plotting import (
    plot_activity_coefficients,
    plot_azeotrope_pressure_sensitivity,
    plot_excess_heat_capacity,
    plot_isothermal_excess_bundle,
    plot_relative_volatility,
)
from solver import compute_txy


def _temp_tag(T_celsius: float) -> str:
    token = f"{T_celsius:.1f}".replace(".", "p")
    return f"t{token}c"


def _pressure_sensitivity(system_id: str, pressures_kpa: np.ndarray, points: int) -> tuple[np.ndarray, np.ndarray]:
    system = load_system_parameters(system_id)
    x_grid = np.linspace(0.01, 0.99, points)

    azeotrope_x: list[float] = []
    azeotrope_t: list[float] = []

    for pressure in pressures_kpa:
        x_n, y_n, T_n = compute_txy(
            pressure_kpa=float(pressure),
            component_1=system.component_1.antoine,
            component_2=system.component_2.antoine,
            gamma_fn=gamma_nrtl,
            gamma_params=system.nrtl,
            x_grid=x_grid,
        )
        azeo = find_azeotrope(x_n, y_n, T_n, tolerance=0.01)
        if azeo is None:
            azeotrope_x.append(np.nan)
            azeotrope_t.append(np.nan)
        else:
            azeotrope_x.append(azeo.x1)
            azeotrope_t.append(azeo.temperature_c)

    return np.asarray(azeotrope_x), np.asarray(azeotrope_t)


def generate_output_pack(
    system_id: str,
    output_dir: Path,
    points: int,
    t_ref_celsius: float,
    include_baseline: bool = True,
    baseline_result: RunResult | None = None,
) -> list[Path]:
    """Generate practical engineering plots."""
    system = load_system_parameters(system_id)
    output_dir.mkdir(parents=True, exist_ok=True)

    files: list[Path] = []
    if include_baseline:
        baseline = baseline_result or run_analysis(
            system_id=system_id,
            points=points,
            save_plots=True,
            output_dir=output_dir,
        )
        files.extend(list(baseline.output_files))

    x_grid = np.linspace(0.01, 0.99, points)
    GE_ref, HE_ref, CPE_ref = compute_isothermal_excess_curves(x_grid, t_ref_celsius, system.nrtl)
    gamma1_ref, gamma2_ref = compute_gamma_curves(x_grid, t_ref_celsius, system.nrtl)

    psat1 = psat_kpa(t_ref_celsius, system.component_1.antoine)
    psat2 = psat_kpa(t_ref_celsius, system.component_2.antoine)
    alpha12 = (gamma1_ref * psat1) / np.maximum(gamma2_ref * psat2, 1e-30)

    tag = _temp_tag(t_ref_celsius)
    files.append(plot_isothermal_excess_bundle(output_dir, x_grid, GE_ref, HE_ref, CPE_ref, tag=tag))
    files.append(plot_excess_heat_capacity(output_dir, x_grid, CPE_ref, tag=tag))
    files.append(plot_activity_coefficients(output_dir, x_grid, gamma1_ref, gamma2_ref, tag=tag))
    files.append(plot_relative_volatility(output_dir, x_grid, alpha12, tag=tag))

    pressure_sweep = np.linspace(70.0, 130.0, 13)
    azeo_x, azeo_t = _pressure_sensitivity(system_id, pressure_sweep, points)
    files.append(plot_azeotrope_pressure_sensitivity(output_dir, pressure_sweep, azeo_x, azeo_t))

    return files


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser."""
    parser = argparse.ArgumentParser(description="Generate practical VLE output pack.")
    parser.add_argument("--system-id", default="ethanol_water_1atm")
    parser.add_argument("--output-dir", default="figures/practical_pack")
    parser.add_argument("--points", type=int, default=80)
    parser.add_argument("--t-ref", type=float, default=78.2, help="Reference isothermal temperature in deg C")
    parser.add_argument("--without-baseline", action="store_true", help="Skip baseline VLE plots and export only practical add-ons")
    return parser


def main() -> None:
    """CLI entrypoint."""
    args = build_parser().parse_args()
    out_dir = Path(args.output_dir)
    files = generate_output_pack(
        system_id=args.system_id,
        output_dir=out_dir,
        points=args.points,
        t_ref_celsius=args.t_ref,
        include_baseline=not args.without_baseline,
    )

    print("=" * 64)
    print("Practical output pack generated")
    print("=" * 64)
    for path in files:
        print(path)
    print("=" * 64)


if __name__ == "__main__":
    main()
