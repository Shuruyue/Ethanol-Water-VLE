#!/usr/bin/env python3
"""Streamlit dashboard for Wolfram-style interactive VLE exploration."""

from __future__ import annotations

import numpy as np

from analysis import compute_gamma_curves, compute_isothermal_excess_curves, find_azeotrope
from models import NRTLParams, gamma_ideal, gamma_nrtl, psat_kpa
from parameter_store import load_system_parameters
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

    GE_ref, HE_ref, CPE_ref = compute_isothermal_excess_curves(x_grid, temp_ref_c, params)
    gamma1_ref, gamma2_ref = compute_gamma_curves(x_grid, temp_ref_c, params)

    psat1 = psat_kpa(temp_ref_c, system.component_1.antoine)
    psat2 = psat_kpa(temp_ref_c, system.component_2.antoine)
    alpha_rel = (gamma1_ref * psat1) / np.maximum(gamma2_ref * psat2, 1e-30)

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
        "GE_ref": GE_ref,
        "HE_ref": HE_ref,
        "CPE_ref": CPE_ref,
        "gamma1_ref": gamma1_ref,
        "gamma2_ref": gamma2_ref,
        "alpha_rel": alpha_rel,
        "azeo": azeo,
    }


def _plot_txy(data):
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(data["x_i"], data["T_i"], "b-", linewidth=1.4, label="ideal bubble")
    ax.plot(data["y_i"], data["T_i"], "b--", linewidth=1.4, label="ideal dew")
    ax.plot(data["x_n"], data["T_n"], "r-", linewidth=2.0, label="nrtl bubble")
    ax.plot(data["y_n"], data["T_n"], "r--", linewidth=2.0, label="nrtl dew")
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
        col4.metric("azeotrope x1", f"{data['azeo'].x1:.4f}")

    row1_col1, row1_col2 = st.columns(2)
    row1_col1.pyplot(_plot_txy(data), clear_figure=True)
    row1_col2.pyplot(_plot_yx(data), clear_figure=True)

    row2_col1, row2_col2 = st.columns(2)
    row2_col1.pyplot(_plot_excess(data, temp_ref_c), clear_figure=True)
    row2_col2.pyplot(_plot_gamma_alpha(data, temp_ref_c), clear_figure=True)

    st.subheader("Parameter Sources")
    st.write(f"NRTL: {base_system.nrtl_source.name}")
    st.write(base_system.nrtl_source.citation)
    st.write(f"Antoine: {base_system.component_1.source.name} / {base_system.component_2.source.name}")


if __name__ == "__main__":
    main()
