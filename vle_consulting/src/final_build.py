#!/usr/bin/env python3
"""Build final deliverable bundle (plots + data tables + summary)."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from export_data import export_full_data_bundle
from output_pack import generate_output_pack
from pipeline import run_analysis


def build_parser() -> argparse.ArgumentParser:
    """Build CLI parser."""
    parser = argparse.ArgumentParser(description="Build final VLE deliverable bundle.")
    parser.add_argument("--system-id", default="ethanol_water_1atm")
    parser.add_argument("--output-root", default="final_outputs")
    parser.add_argument("--points", type=int, default=80)
    parser.add_argument("--t-ref", type=float, default=78.2)
    parser.add_argument("--p-min", type=float, default=70.0)
    parser.add_argument("--p-max", type=float, default=130.0)
    parser.add_argument("--p-steps", type=int, default=13)
    return parser


def main() -> None:
    """Run full final-build pipeline."""
    args = build_parser().parse_args()
    output_root = Path(args.output_root)
    plots_dir = output_root / "plots"
    data_dir = output_root / "data"

    baseline = run_analysis(
        system_id=args.system_id,
        points=args.points,
        save_plots=True,
        output_dir=plots_dir / "baseline",
    )

    practical_plot_files = generate_output_pack(
        system_id=args.system_id,
        output_dir=plots_dir / "practical_pack",
        points=args.points,
        t_ref_celsius=args.t_ref,
        include_baseline=False,
        baseline_result=baseline,
    )

    pressure_values = np.linspace(args.p_min, args.p_max, args.p_steps)
    data_files = export_full_data_bundle(
        result=baseline,
        t_ref_celsius=args.t_ref,
        points=args.points,
        pressure_values_kpa=pressure_values,
        output_dir=data_dir,
    )

    print("=" * 72)
    print("Final build completed")
    print("=" * 72)
    print("Plot files:")
    for path in list(baseline.output_files) + practical_plot_files:
        print(f"  {path}")
    print("Data files:")
    for path in data_files:
        print(f"  {path}")
    print("=" * 72)


if __name__ == "__main__":
    main()
