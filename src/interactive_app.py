#!/usr/bin/env python3
"""Streamlit dashboard for Wolfram-style interactive VLE exploration."""

from __future__ import annotations

import csv
from io import StringIO
from dataclasses import replace

import numpy as np

from analysis import find_azeotrope
from models import NRTLParams, gamma_ideal, gamma_nrtl, gamma_van_laar
from parameter_store import load_system_parameters
from profiles import build_isothermal_profile
from solver import compute_txy

try:
    import matplotlib.pyplot as plt
    import streamlit as st
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "Missing dependency. Install with: pip install streamlit matplotlib"
    ) from exc


@st.cache_data(show_spinner=False)
def compute_dashboard(
    system_id: str,
    pressure_kpa: float,
    temp_ref_c: float,
    points: int,
    alpha12: float,
    alpha21: float,
    a12: float,
    b12: float,
    a21: float,
    b21: float,
):
    """Compute dashboard datasets."""
    system = load_system_parameters(system_id)
    params = NRTLParams(
        alpha12=alpha12,
        alpha21=alpha21,
        a12=a12,
        b12=b12,
        a21=a21,
        b21=b21,
    )

    x_grid = np.linspace(0.01, 0.99, points)

    x_i, y_i, T_i = compute_txy(
        pressure_kpa=pressure_kpa,
        component_1=system.component_1.antoine,
        component_2=system.component_2.antoine,
        gamma_fn=gamma_ideal,
        gamma_params=None,
        x_grid=x_grid,
    )

    x_n, y_n, T_n = compute_txy(
        pressure_kpa=pressure_kpa,
        component_1=system.component_1.antoine,
        component_2=system.component_2.antoine,
        gamma_fn=gamma_nrtl,
        gamma_params=params,
        x_grid=x_grid,
    )
    x_v, y_v, T_v = compute_txy(
        pressure_kpa=pressure_kpa,
        component_1=system.component_1.antoine,
        component_2=system.component_2.antoine,
        gamma_fn=gamma_van_laar,
        gamma_params=system.van_laar,
        x_grid=x_grid,
    )

    tuned_system = replace(system, pressure_kpa=pressure_kpa, nrtl=params)
    iso = build_isothermal_profile(system=tuned_system, temperature_c=temp_ref_c, points=points)

    azeo = find_azeotrope(x_n, y_n, T_n, tolerance=0.01)

    return {
        "system": system,
        "params": params,
        "x_i": x_i,
        "y_i": y_i,
        "T_i": T_i,
        "x_n": x_n,
        "y_n": y_n,
        "T_n": T_n,
        "x_v": x_v,
        "y_v": y_v,
        "T_v": T_v,
        "GE_ref": iso.GE_j_mol,
        "HE_ref": iso.HE_j_mol,
        "CPE_ref": iso.CPE_j_molK,
        "gamma1_ref": iso.gamma1,
        "gamma2_ref": iso.gamma2,
        "alpha_rel": iso.alpha12,
        "azeo": azeo,
    }


def _plot_txy(data):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(data["x_i"], data["T_i"], "b-", linewidth=1.4, label="ideal bubble")
    ax.plot(data["y_i"], data["T_i"], "b--", linewidth=1.4, label="ideal dew")
    ax.plot(data["x_v"], data["T_v"], color="green", linewidth=1.5, label="van laar bubble")
    ax.plot(data["y_v"], data["T_v"], color="green", linestyle="--", linewidth=1.5, label="van laar dew")
    ax.plot(data["x_n"], data["T_n"], "r-", linewidth=2.0, label="nrtl bubble")
    ax.plot(data["y_n"], data["T_n"], "r--", linewidth=2.0, label="nrtl dew")
    ax.scatter([0.8943], [78.15], marker="D", color="black", s=28, label="lit azeotrope")
    ax.set_xlabel("ethanol mole fraction")
    ax.set_ylabel("temperature (deg C)")
    ax.set_title("T-x-y")
    ax.set_xlim(0.0, 1.0)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    return fig


def _plot_yx(data):
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot([0.0, 1.0], [0.0, 1.0], "k:", linewidth=1.0, label="y=x")
    ax.plot(data["x_i"], data["y_i"], "b-", linewidth=1.4, label="ideal")
    ax.plot(data["x_v"], data["y_v"], color="green", linewidth=1.5, label="van laar")
    ax.plot(data["x_n"], data["y_n"], "r-", linewidth=2.0, label="nrtl")
    ax.set_xlabel("x1")
    ax.set_ylabel("y1")
    ax.set_title("y-x")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_aspect("equal")
    fig.tight_layout()
    return fig


def _plot_excess(data, temp_ref_c: float):
    fig, ax = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    ax[0].plot(data["x_n"], data["HE_ref"] / 1000.0, "m-", linewidth=2.0)
    ax[0].axhline(0.0, color="k", linewidth=0.8)
    ax[0].set_ylabel("HE (kJ/mol)")
    ax[0].set_title(f"Excess Enthalpy at {temp_ref_c:.1f} deg C")
    ax[0].grid(alpha=0.3)

    ax[1].plot(data["x_n"], data["CPE_ref"], color="darkorange", linewidth=2.0)
    ax[1].axhline(0.0, color="k", linewidth=0.8)
    ax[1].set_xlabel("x1")
    ax[1].set_ylabel("CPE (J/mol/K)")
    ax[1].set_title("Excess Heat Capacity")
    ax[1].grid(alpha=0.3)

    fig.tight_layout()
    return fig


def _plot_gamma_alpha(data, temp_ref_c: float):
    fig, ax = plt.subplots(2, 1, figsize=(7, 8), sharex=True)
    ax[0].plot(data["x_n"], data["gamma1_ref"], "r-", linewidth=2.0, label="gamma1")
    ax[0].plot(data["x_n"], data["gamma2_ref"], "b-", linewidth=2.0, label="gamma2")
    ax[0].set_ylabel("activity coefficient")
    ax[0].set_title(f"Activity Coefficients at {temp_ref_c:.1f} deg C")
    ax[0].grid(alpha=0.3)
    ax[0].legend(loc="best")

    ax[1].plot(data["x_n"], data["alpha_rel"], color="teal", linewidth=2.0)
    ax[1].axhline(1.0, color="k", linestyle="--", linewidth=0.9)
    ax[1].set_xlabel("x1")
    ax[1].set_ylabel("alpha12")
    ax[1].set_title("Relative Volatility")
    ax[1].grid(alpha=0.3)

    fig.tight_layout()
    return fig


def _build_dashboard_csv(data, temp_ref_c: float) -> str:
    """Build downloadable CSV for the dashboard state."""
    handle = StringIO()
    writer = csv.writer(handle)
    writer.writerow(
        [
            "x1",
            "y1_ideal",
            "T_ideal_c",
            "y1_van_laar",
            "T_van_laar_c",
            "y1_nrtl",
            "T_nrtl_c",
            f"GE_{temp_ref_c:.1f}C_j_mol",
            f"HE_{temp_ref_c:.1f}C_j_mol",
            f"CPE_{temp_ref_c:.1f}C_j_molK",
            f"gamma1_{temp_ref_c:.1f}C",
            f"gamma2_{temp_ref_c:.1f}C",
            f"alpha12_{temp_ref_c:.1f}C",
        ]
    )
    for i in range(len(data["x_n"])):
        writer.writerow(
            [
                float(data["x_n"][i]),
                float(data["y_i"][i]),
                float(data["T_i"][i]),
                float(data["y_v"][i]),
                float(data["T_v"][i]),
                float(data["y_n"][i]),
                float(data["T_n"][i]),
                float(data["GE_ref"][i]),
                float(data["HE_ref"][i]),
                float(data["CPE_ref"][i]),
                float(data["gamma1_ref"][i]),
                float(data["gamma2_ref"][i]),
                float(data["alpha_rel"][i]),
            ]
        )
    return handle.getvalue()


def main() -> None:
    """Render dashboard."""
    st.set_page_config(page_title="Ethanol-Water Interactive VLE", layout="wide")
    st.title("Ethanol-Water Interactive VLE Dashboard")
    st.caption("Wolfram-style exploration for VLE, HE, and CPE using NRTL")

    base_system = load_system_parameters("ethanol_water_1atm")

    with st.sidebar:
        st.header("Controls")
        pressure_kpa = st.slider(
            "Pressure (kPa)",
            min_value=70.0,
            max_value=130.0,
            value=float(base_system.pressure_kpa),
            step=0.5,
        )
        temp_ref_c = st.slider(
            "Reference Temperature (deg C)",
            min_value=20.0,
            max_value=110.0,
            value=78.2,
            step=0.1,
        )
        points = st.slider("Grid points", min_value=30, max_value=120, value=70, step=5)

        st.subheader("NRTL Parameters")
        alpha12 = st.slider("alpha12", min_value=0.10, max_value=0.60, value=float(base_system.nrtl.alpha12), step=0.01)
        alpha21 = st.slider("alpha21", min_value=0.10, max_value=0.60, value=float(base_system.nrtl.alpha21), step=0.01)
        a12 = st.slider("a12", min_value=-3.0, max_value=3.0, value=float(base_system.nrtl.a12), step=0.001)
        b12 = st.slider("b12", min_value=-1200.0, max_value=1200.0, value=float(base_system.nrtl.b12), step=1.0)
        a21 = st.slider("a21", min_value=-3.0, max_value=6.0, value=float(base_system.nrtl.a21), step=0.001)
        b21 = st.slider("b21", min_value=-1200.0, max_value=1200.0, value=float(base_system.nrtl.b21), step=1.0)

    if abs(a12) < 1e-12 and abs(a21) < 1e-12:
        st.error("a12 和 a21 不可同時為 0。請調整常數項參數。")
        st.stop()

    with st.spinner("Computing equilibrium and excess-property curves..."):
        data = compute_dashboard(
            system_id="ethanol_water_1atm",
            pressure_kpa=pressure_kpa,
            temp_ref_c=temp_ref_c,
            points=points,
            alpha12=alpha12,
            alpha21=alpha21,
            a12=a12,
            b12=b12,
            a21=a21,
            b21=b21,
        )

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("max GE (J/mol)", f"{data['GE_ref'].max():.1f}")
    col2.metric("max HE (J/mol)", f"{data['HE_ref'].max():.1f}")
    col3.metric("min HE (J/mol)", f"{data['HE_ref'].min():.1f}")
    if data["azeo"] is None:
        col4.metric("azeotrope", "not found")
    else:
        col4.metric("azeotrope x1", f"{data['azeo'].x1:.4f}", delta=f"{data['azeo'].x1 - 0.8943:+.4f} vs lit")

    tab_vle, tab_excess, tab_data = st.tabs(["VLE Maps", "Excess Thermodynamics", "Data & Diagnostics"])

    with tab_vle:
        row1_col1, row1_col2 = st.columns(2)
        row1_col1.pyplot(_plot_txy(data), clear_figure=True)
        row1_col2.pyplot(_plot_yx(data), clear_figure=True)

    with tab_excess:
        row2_col1, row2_col2 = st.columns(2)
        row2_col1.pyplot(_plot_excess(data, temp_ref_c), clear_figure=True)
        row2_col2.pyplot(_plot_gamma_alpha(data, temp_ref_c), clear_figure=True)

    with tab_data:
        csv_payload = _build_dashboard_csv(data, temp_ref_c)
        st.download_button(
            label="Download current dashboard data (CSV)",
            data=csv_payload,
            file_name=f"dashboard_state_{temp_ref_c:.1f}C_{pressure_kpa:.1f}kPa.csv",
            mime="text/csv",
        )
        if data["azeo"] is not None:
            st.write(
                f"Detected azeotrope: x1={data['azeo'].x1:.5f}, y1={data['azeo'].y1:.5f}, "
                f"T={data['azeo'].temperature_c:.3f} deg C, |x-y|={data['azeo'].abs_error:.6f}"
            )
        st.write(
            f"Relative volatility range at {temp_ref_c:.1f} deg C: "
            f"{float(np.min(data['alpha_rel'])):.3f} to {float(np.max(data['alpha_rel'])):.3f}"
        )

    st.subheader("Parameter Sources")
    st.write(f"NRTL: {base_system.nrtl_source.name}")
    st.write(base_system.nrtl_source.citation)
    st.write(f"Antoine: {base_system.component_1.source.name} / {base_system.component_2.source.name}")


if __name__ == "__main__":
    main()
