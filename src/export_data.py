#!/usr/bin/env python3
"""Export numerical datasets for final-report appendices."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from parameter_store import SystemParameters, load_system_parameters
from pipeline import RunResult
from profiles import build_isothermal_profile, build_pressure_sensitivity_profile


def _write_csv(path: Path, header: list[str], rows: list[list[float | str]]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def export_baseline_curve_data(result: RunResult, output_dir: Path) -> dict[str, Path]:
    """Export baseline model curves and excess properties."""
    output_dir.mkdir(parents=True, exist_ok=True)

    model_rows: list[list[float | str]] = []
    for i in range(len(result.ideal.x)):
        model_rows.append(["ideal", float(result.ideal.x[i]), float(result.ideal.y[i]), float(result.ideal.T_celsius[i])])
        model_rows.append(["van_laar", float(result.van_laar.x[i]), float(result.van_laar.y[i]), float(result.van_laar.T_celsius[i])])
        model_rows.append(["nrtl", float(result.nrtl.x[i]), float(result.nrtl.y[i]), float(result.nrtl.T_celsius[i])])

    model_path = _write_csv(
        output_dir / "vle_model_curves.csv",
        ["model", "x1", "y1", "temperature_c"],
        model_rows,
    )

    excess_rows: list[list[float | str]] = []
    for i in range(len(result.nrtl.x)):
        excess_rows.append(
            [
                float(result.nrtl.x[i]),
                float(result.nrtl.T_celsius[i]),
                float(result.GE_j_mol[i]),
                float(result.HE_j_mol[i]),
            ]
        )

    excess_path = _write_csv(
        output_dir / "excess_properties_bubble_line.csv",
        ["x1", "temperature_c", "GE_j_mol", "HE_j_mol"],
        excess_rows,
    )

    return {
        "vle_model_curves": model_path,
        "excess_properties_bubble_line": excess_path,
    }


def export_isothermal_data(
    system: SystemParameters,
    temperature_c: float,
    points: int,
    output_dir: Path,
) -> dict[str, Path]:
    """Export GE/HE/CPE/gamma at a fixed temperature."""
    output_dir.mkdir(parents=True, exist_ok=True)
    iso = build_isothermal_profile(system=system, temperature_c=temperature_c, points=points)

    iso_rows: list[list[float | str]] = []
    for i in range(points):
        iso_rows.append(
            [
                float(iso.x[i]),
                float(iso.GE_j_mol[i]),
                float(iso.HE_j_mol[i]),
                float(iso.CPE_j_molK[i]),
                float(iso.gamma1[i]),
                float(iso.gamma2[i]),
                float(iso.alpha12[i]),
            ]
        )

    iso_path = _write_csv(
        output_dir / f"isothermal_properties_{temperature_c:.1f}C.csv",
        ["x1", "GE_j_mol", "HE_j_mol", "CPE_j_molK", "gamma1", "gamma2", "alpha12"],
        iso_rows,
    )

    return {"isothermal_properties": iso_path}


def export_pressure_sensitivity_data(
    system: SystemParameters,
    pressure_values_kpa: np.ndarray,
    points: int,
    output_dir: Path,
) -> dict[str, Path]:
    """Export pressure-vs-azeotrope table."""
    output_dir.mkdir(parents=True, exist_ok=True)
    sweep = build_pressure_sensitivity_profile(system=system, pressure_values_kpa=pressure_values_kpa, points=points)

    rows: list[list[float | str]] = []
    for i, pressure in enumerate(sweep.pressure_kpa):
        if np.isnan(sweep.azeotrope_x1[i]):
            rows.append([float(pressure), "", "", ""])
        else:
            rows.append(
                [
                    float(pressure),
                    float(sweep.azeotrope_x1[i]),
                    float(sweep.azeotrope_y1[i]),
                    float(sweep.azeotrope_temperature_c[i]),
                ]
            )

    sensitivity_path = _write_csv(
        output_dir / "azeotrope_pressure_sensitivity.csv",
        ["pressure_kpa", "azeotrope_x1", "azeotrope_y1", "azeotrope_temperature_c"],
        rows,
    )

    return {"azeotrope_pressure_sensitivity": sensitivity_path}


def export_run_summary_json(
    result: RunResult,
    t_ref_celsius: float,
    output_dir: Path,
) -> Path:
    """Export key metrics and source metadata."""
    output_dir.mkdir(parents=True, exist_ok=True)
    iso = build_isothermal_profile(result.system, temperature_c=t_ref_celsius, points=len(result.nrtl.x))

    payload = {
        "system": {
            "id": result.system.system_id,
            "name": result.system.system_name,
            "pressure_kpa": result.system.pressure_kpa,
        },
        "parameter_sources": {
            "component_1": {
                "name": result.system.component_1.source.name,
                "citation": result.system.component_1.source.citation,
                "url": result.system.component_1.source.url,
            },
            "component_2": {
                "name": result.system.component_2.source.name,
                "citation": result.system.component_2.source.citation,
                "url": result.system.component_2.source.url,
            },
            "nrtl": {
                "name": result.system.nrtl_source.name,
                "citation": result.system.nrtl_source.citation,
                "url": result.system.nrtl_source.url,
            },
        },
        "nrtl_parameters": {
            "alpha12": result.system.nrtl.alpha12,
            "alpha21": result.system.nrtl.alpha21,
            "a12": result.system.nrtl.a12,
            "b12": result.system.nrtl.b12,
            "a21": result.system.nrtl.a21,
            "b21": result.system.nrtl.b21,
        },
        "metrics": {
            "max_ge_j_mol": float(np.max(result.GE_j_mol)),
            "max_he_j_mol": float(np.max(result.HE_j_mol)),
            "min_he_j_mol": float(np.min(result.HE_j_mol)),
            "reference_temperature_c": float(t_ref_celsius),
            "max_relative_volatility_ref": float(np.max(iso.alpha12)),
            "min_relative_volatility_ref": float(np.min(iso.alpha12)),
        },
        "azeotrope": None,
    }

    if result.azeotrope is not None:
        payload["azeotrope"] = {
            "x1": float(result.azeotrope.x1),
            "y1": float(result.azeotrope.y1),
            "temperature_c": float(result.azeotrope.temperature_c),
            "abs_error": float(result.azeotrope.abs_error),
            "x1_vs_literature_0p8943": float(result.azeotrope.x1 - 0.8943),
            "temperature_vs_literature_78p15": float(result.azeotrope.temperature_c - 78.15),
        }

    summary_path = output_dir / "run_summary.json"
    with summary_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)

    return summary_path


def export_full_data_bundle(
    result: RunResult,
    t_ref_celsius: float,
    points: int,
    pressure_values_kpa: np.ndarray,
    output_dir: Path,
) -> list[Path]:
    """Export all CSV/JSON data files for final report."""
    files: list[Path] = []
    files.extend(export_baseline_curve_data(result, output_dir / "baseline").values())
    files.extend(export_isothermal_data(result.system, t_ref_celsius, points, output_dir / "isothermal").values())
    files.extend(
        export_pressure_sensitivity_data(result.system, pressure_values_kpa, points, output_dir / "sensitivity").values()
    )
    files.append(export_run_summary_json(result, t_ref_celsius, output_dir))
    return files
