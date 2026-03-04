#!/usr/bin/env python3
"""Generate a practical output image pack for engineering use."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from parameter_store import load_system_parameters
from pipeline import RunResult, run_analysis
from plotting import (
    plot_activity_coefficients,
    plot_azeotrope_pressure_sensitivity,
    plot_excess_heat_capacity,
    plot_isothermal_excess_bundle,
    plot_relative_volatility,
)
from profiles import build_isothermal_profile, build_pressure_sensitivity_profile


def _temp_tag(T_celsius: float) -> str:
    token = f"{T_celsius:.1f}".replace(".", "p")
    return f"t{token}c"


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

    iso = build_isothermal_profile(system=system, temperature_c=t_ref_celsius, points=points)

    tag = _temp_tag(t_ref_celsius)
    files.append(
        plot_isothermal_excess_bundle(
            output_dir,
            iso.x,
            iso.GE_j_mol,
            iso.HE_j_mol,
            iso.CPE_j_molK,
            tag=tag,
        )
    )
    files.append(plot_excess_heat_capacity(output_dir, iso.x, iso.CPE_j_molK, tag=tag))
    files.append(plot_activity_coefficients(output_dir, iso.x, iso.gamma1, iso.gamma2, tag=tag))
    files.append(plot_relative_volatility(output_dir, iso.x, iso.alpha12, tag=tag))

    pressure_sweep = np.linspace(70.0, 130.0, 13)
    sweep = build_pressure_sensitivity_profile(system=system, pressure_values_kpa=pressure_sweep, points=points)
    files.append(
        plot_azeotrope_pressure_sensitivity(
            output_dir,
            sweep.pressure_kpa,
            sweep.azeotrope_x1,
            sweep.azeotrope_temperature_c,
        )
    )

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
