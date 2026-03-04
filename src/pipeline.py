"""End-to-end VLE analysis pipeline."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from analysis import AzeotropePoint, compute_excess_curves, find_azeotrope
from models import gamma_ideal, gamma_nrtl, gamma_van_laar
from parameter_store import SystemParameters, load_system_parameters
from plotting import (
    plot_excess_enthalpy,
    plot_excess_gibbs,
    plot_mccabe_thiele,
    plot_txy,
    plot_txy_with_experiment,
    plot_yx,
)
from solver import compute_txy

CurveBuilderFn = Callable[[float, float, object], tuple[float, float]]


@dataclass(frozen=True)
class CurveResult:
    """One model curve set."""

    x: np.ndarray
    y: np.ndarray
    T_celsius: np.ndarray


@dataclass(frozen=True)
class RunResult:
    """Complete run output."""

    system: SystemParameters
    ideal: CurveResult
    van_laar: CurveResult
    nrtl: CurveResult
    GE_j_mol: np.ndarray
    HE_j_mol: np.ndarray
    azeotrope: AzeotropePoint | None
    output_files: tuple[Path, ...]


def _build_curve(
    system: SystemParameters,
    x_grid: np.ndarray,
    gamma_fn: CurveBuilderFn,
    gamma_params: object,
) -> CurveResult:
    x, y, T = compute_txy(
        pressure_kpa=system.pressure_kpa,
        component_1=system.component_1.antoine,
        component_2=system.component_2.antoine,
        gamma_fn=gamma_fn,
        gamma_params=gamma_params,
        x_grid=x_grid,
    )
    return CurveResult(x=x, y=y, T_celsius=T)


def run_analysis(
    system_id: str = "ethanol_water_1atm",
    points: int = 60,
    save_plots: bool = True,
    output_dir: Path | str | None = None,
) -> RunResult:
    """Run full VLE workflow for one system ID."""
    system = load_system_parameters(system_id=system_id)
    x_grid = np.linspace(0.01, 0.99, points)
    ideal_curve = _build_curve(system, x_grid, gamma_ideal, None)
    van_laar_curve = _build_curve(system, x_grid, gamma_van_laar, system.van_laar)
    nrtl_curve = _build_curve(system, x_grid, gamma_nrtl, system.nrtl)

    GE, HE = compute_excess_curves(nrtl_curve.x, nrtl_curve.T_celsius, system.nrtl)
    azeotrope = find_azeotrope(nrtl_curve.x, nrtl_curve.y, nrtl_curve.T_celsius, tolerance=0.01)

    output_files: list[Path] = []
    if save_plots:
        plot_dir = Path(output_dir) if output_dir else Path(__file__).resolve().parents[1] / "figures"
        output_files.extend(
            [
                plot_txy(
                    plot_dir,
                    ideal_curve.x,
                    ideal_curve.y,
                    ideal_curve.T_celsius,
                    van_laar_curve.x,
                    van_laar_curve.y,
                    van_laar_curve.T_celsius,
                    nrtl_curve.x,
                    nrtl_curve.y,
                    nrtl_curve.T_celsius,
                ),
                plot_yx(
                    plot_dir,
                    ideal_curve.x,
                    ideal_curve.y,
                    van_laar_curve.x,
                    van_laar_curve.y,
                    nrtl_curve.x,
                    nrtl_curve.y,
                ),
                plot_excess_gibbs(plot_dir, nrtl_curve.x, GE),
                plot_excess_enthalpy(plot_dir, nrtl_curve.x, HE),
                plot_mccabe_thiele(plot_dir, nrtl_curve.x, nrtl_curve.y, xD=0.85, xB=0.05, xF=0.3, R_reflux=2.0),
                plot_txy_with_experiment(plot_dir, nrtl_curve.x, nrtl_curve.y, nrtl_curve.T_celsius),
            ]
        )

    return RunResult(
        system=system,
        ideal=ideal_curve,
        van_laar=van_laar_curve,
        nrtl=nrtl_curve,
        GE_j_mol=GE,
        HE_j_mol=HE,
        azeotrope=azeotrope,
        output_files=tuple(output_files),
    )
