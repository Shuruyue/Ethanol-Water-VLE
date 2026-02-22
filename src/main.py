#!/usr/bin/env python3
"""CLI entrypoint for VLE analysis."""

from __future__ import annotations

import argparse
from pathlib import Path

from pipeline import run_analysis


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser."""
    parser = argparse.ArgumentParser(description="Ethanol-water VLE analysis at fixed pressure.")
    parser.add_argument("--system-id", default="ethanol_water_1atm", help="System ID in parameter_db.json")
    parser.add_argument("--points", type=int, default=60, help="Number of composition points")
    parser.add_argument("--output-dir", default=None, help="Directory for generated plots")
    parser.add_argument("--no-plots", action="store_true", help="Skip plot generation")
    return parser


def print_summary(result) -> None:
    """Print a compact run summary."""
    print("=" * 60)
    print("VLE analysis summary")
    print("=" * 60)
    print(f"system        : {result.system.system_name}")
    print(f"pressure (kPa): {result.system.pressure_kpa:.3f}")
    print(f"nrtl source   : {result.system.nrtl_source.name}")
    print(f"a12, b12      : {result.system.nrtl.a12:.6g}, {result.system.nrtl.b12:.6g}")
    print(f"a21, b21      : {result.system.nrtl.a21:.6g}, {result.system.nrtl.b21:.6g}")
    print(f"max GE (J/mol): {result.GE_j_mol.max():.3f}")
    print(f"max HE (J/mol): {result.HE_j_mol.max():.3f}")
    print(f"min HE (J/mol): {result.HE_j_mol.min():.3f}")

    ethanol_bp = float(result.nrtl.T_celsius[-1])
    water_bp = float(result.nrtl.T_celsius[0])
    print(f"boiling point water   : {water_bp:.3f} deg C")
    print(f"boiling point ethanol : {ethanol_bp:.3f} deg C")

    if result.azeotrope:
        print("azeotrope      : detected")
        print(f"x1             : {result.azeotrope.x1:.4f}")
        print(f"y1             : {result.azeotrope.y1:.4f}")
        print(f"T (deg C)      : {result.azeotrope.temperature_c:.3f}")
        print(f"|x-y|          : {result.azeotrope.abs_error:.5f}")
    else:
        print("azeotrope      : not detected in current grid")

    if result.output_files:
        print("output files:")
        for item in result.output_files:
            print(f"  - {item}")
    print("=" * 60)


def main() -> None:
    """Run CLI."""
    parser = build_parser()
    args = parser.parse_args()

    output_dir = Path(args.output_dir) if args.output_dir else None
    result = run_analysis(
        system_id=args.system_id,
        points=args.points,
        save_plots=not args.no_plots,
        output_dir=output_dir,
    )
    print_summary(result)


if __name__ == "__main__":
    main()
