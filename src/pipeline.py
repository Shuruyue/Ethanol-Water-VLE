"""End-to-end VLE analysis pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from analysis import AzeotropePoint, compute_excess_curves, find_azeotrope
from models import gamma_ideal, gamma_nrtl, gamma_van_laar
from parameter_store import SystemParameters, load_system_parameters
from plotting import plot_excess_enthalpy, plot_excess_gibbs, plot_txy, plot_yx
from solver import compute_txy


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


def run_analysis(
    system_id: str = "ethanol_water_1atm",
    points: int = 60,
    save_plots: bool = True,
    output_dir: Path | str | None = None,
) -> RunResult:
    """Run full VLE workflow for one system ID."""
    system = load_system_parameters(system_id=system_id)
    x_grid = np.linspace(0.01, 0.99, points)

    x_i, y_i, T_i = compute_txy(
        pressure_kpa=system.pressure_kpa,
        component_1=system.component_1.antoine,
        component_2=system.component_2.antoine,
        gamma_fn=gamma_ideal,
        gamma_params=None,
        x_grid=x_grid,
    )

    x_v, y_v, T_v = compute_txy(
        pressure_kpa=system.pressure_kpa,
        component_1=system.component_1.antoine,
        component_2=system.component_2.antoine,
        gamma_fn=gamma_van_laar,
        gamma_params=system.van_laar,
        x_grid=x_grid,
    )

    x_n, y_n, T_n = compute_txy(
        pressure_kpa=system.pressure_kpa,
        component_1=system.component_1.antoine,
        component_2=system.component_2.antoine,
        gamma_fn=gamma_nrtl,
        gamma_params=system.nrtl,
        x_grid=x_grid,
    )

    GE, HE = compute_excess_curves(x_n, T_n, system.nrtl)
    azeotrope = find_azeotrope(x_n, y_n, T_n, tolerance=0.01)

    output_files: list[Path] = []
    if save_plots:
        plot_dir = Path(output_dir) if output_dir else Path(__file__).resolve().parents[1] / "figures"
        output_files.extend(
            [
                plot_txy(plot_dir, x_i, y_i, T_i, x_v, y_v, T_v, x_n, y_n, T_n),
                plot_yx(plot_dir, x_i, y_i, x_v, y_v, x_n, y_n),
                plot_excess_gibbs(plot_dir, x_n, GE),
                plot_excess_enthalpy(plot_dir, x_n, HE),
            ]
        )

    return RunResult(
        system=system,
        ideal=CurveResult(x=x_i, y=y_i, T_celsius=T_i),
        van_laar=CurveResult(x=x_v, y=y_v, T_celsius=T_v),
        nrtl=CurveResult(x=x_n, y=y_n, T_celsius=T_n),
        GE_j_mol=GE,
        HE_j_mol=HE,
        azeotrope=azeotrope,
        output_files=tuple(output_files),
    )
