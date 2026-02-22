"""Plot generation helpers."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def _ensure_output_dir(output_dir: Path | str) -> Path:
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def plot_txy(
    output_dir: Path | str,
    x_ideal: np.ndarray,
    y_ideal: np.ndarray,
    t_ideal: np.ndarray,
    x_van_laar: np.ndarray,
    y_van_laar: np.ndarray,
    t_van_laar: np.ndarray,
    x_nrtl: np.ndarray,
    y_nrtl: np.ndarray,
    t_nrtl: np.ndarray,
    dpi: int = 300,
) -> Path:
    """Save T-x-y comparison plot."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / "txy_ethanol_water_1atm.png"

    plt.figure(figsize=(10, 7))
    plt.plot(x_ideal, t_ideal, "b-", linewidth=1.5, label="ideal bubble")
    plt.plot(y_ideal, t_ideal, "b--", linewidth=1.5, label="ideal dew")
    plt.plot(x_van_laar, t_van_laar, "g-", linewidth=1.5, label="van_laar bubble")
    plt.plot(y_van_laar, t_van_laar, "g--", linewidth=1.5, label="van_laar dew")
    plt.plot(x_nrtl, t_nrtl, "r-", linewidth=2.0, label="nrtl bubble")
    plt.plot(y_nrtl, t_nrtl, "r--", linewidth=2.0, label="nrtl dew")
    plt.xlabel("ethanol mole fraction (x1 or y1)")
    plt.ylabel("temperature (deg C)")
    plt.title("ethanol-water T-x-y at 1 atm")
    plt.xlim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out_path


def plot_yx(
    output_dir: Path | str,
    x_ideal: np.ndarray,
    y_ideal: np.ndarray,
    x_van_laar: np.ndarray,
    y_van_laar: np.ndarray,
    x_nrtl: np.ndarray,
    y_nrtl: np.ndarray,
    dpi: int = 300,
) -> Path:
    """Save y-x comparison plot."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / "yx_ethanol_water_1atm.png"

    plt.figure(figsize=(8, 8))
    plt.plot([0.0, 1.0], [0.0, 1.0], "k:", linewidth=1.0, label="y=x")
    plt.plot(x_ideal, y_ideal, "b-", linewidth=1.5, label="ideal")
    plt.plot(x_van_laar, y_van_laar, "g-", linewidth=1.5, label="van_laar")
    plt.plot(x_nrtl, y_nrtl, "r-", linewidth=2.0, label="nrtl")
    plt.xlabel("liquid mole fraction x1")
    plt.ylabel("vapor mole fraction y1")
    plt.title("ethanol-water y-x at 1 atm")
    plt.xlim(0.0, 1.0)
    plt.ylim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.legend(loc="lower right")
    plt.gca().set_aspect("equal")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out_path


def plot_excess_gibbs(output_dir: Path | str, x_nrtl: np.ndarray, GE_j_mol: np.ndarray, dpi: int = 300) -> Path:
    """Save excess Gibbs energy plot."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / "ge_ethanol_water_1atm.png"

    plt.figure(figsize=(10, 6))
    plt.plot(x_nrtl, GE_j_mol / 1000.0, "b-", linewidth=2.0)
    plt.xlabel("liquid mole fraction x1")
    plt.ylabel("GE (kJ/mol)")
    plt.title("nrtl excess Gibbs energy")
    plt.axhline(0.0, color="k", linewidth=0.5)
    plt.xlim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out_path


def plot_excess_enthalpy(output_dir: Path | str, x_nrtl: np.ndarray, HE_j_mol: np.ndarray, dpi: int = 300) -> Path:
    """Save excess enthalpy plot."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / "he_ethanol_water_1atm.png"

    plt.figure(figsize=(10, 6))
    plt.plot(x_nrtl, HE_j_mol / 1000.0, "m-", linewidth=2.0)
    plt.xlabel("liquid mole fraction x1")
    plt.ylabel("HE (kJ/mol)")
    plt.title("nrtl excess enthalpy")
    plt.axhline(0.0, color="k", linewidth=0.5)
    plt.xlim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out_path


def plot_excess_heat_capacity(output_dir: Path | str, x_nrtl: np.ndarray, CPE_j_molK: np.ndarray, tag: str, dpi: int = 300) -> Path:
    """Save excess heat capacity plot."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / f"cpe_ethanol_water_{tag}.png"

    plt.figure(figsize=(10, 6))
    plt.plot(x_nrtl, CPE_j_molK, color="darkorange", linewidth=2.0)
    plt.xlabel("liquid mole fraction x1")
    plt.ylabel("CPE (J/mol/K)")
    plt.title("nrtl excess heat capacity")
    plt.axhline(0.0, color="k", linewidth=0.5)
    plt.xlim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out_path


def plot_activity_coefficients(
    output_dir: Path | str,
    x_nrtl: np.ndarray,
    gamma1: np.ndarray,
    gamma2: np.ndarray,
    tag: str,
    dpi: int = 300,
) -> Path:
    """Save activity-coefficient curves."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / f"gamma_ethanol_water_{tag}.png"

    plt.figure(figsize=(10, 6))
    plt.plot(x_nrtl, gamma1, "r-", linewidth=2.0, label="gamma1 ethanol")
    plt.plot(x_nrtl, gamma2, "b-", linewidth=2.0, label="gamma2 water")
    plt.xlabel("liquid mole fraction x1")
    plt.ylabel("activity coefficient")
    plt.title("nrtl activity coefficients")
    plt.xlim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out_path


def plot_relative_volatility(
    output_dir: Path | str,
    x_nrtl: np.ndarray,
    alpha12: np.ndarray,
    tag: str,
    dpi: int = 300,
) -> Path:
    """Save relative-volatility curve."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / f"relative_volatility_ethanol_water_{tag}.png"

    plt.figure(figsize=(10, 6))
    plt.plot(x_nrtl, alpha12, color="teal", linewidth=2.0)
    plt.axhline(1.0, color="k", linewidth=0.8, linestyle="--")
    plt.xlabel("liquid mole fraction x1")
    plt.ylabel("relative volatility alpha12")
    plt.title("relative volatility at fixed temperature")
    plt.xlim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out_path


def plot_azeotrope_pressure_sensitivity(
    output_dir: Path | str,
    pressures_kpa: np.ndarray,
    azeotrope_x: np.ndarray,
    azeotrope_t: np.ndarray,
    dpi: int = 300,
) -> Path:
    """Save azeotrope sensitivity to pressure."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / "azeotrope_pressure_sensitivity.png"

    fig, ax = plt.subplots(2, 1, figsize=(9, 9), sharex=True)

    ax[0].plot(pressures_kpa, azeotrope_x, "o-", color="purple", linewidth=1.5)
    ax[0].set_ylabel("azeotrope x1")
    ax[0].set_title("azeotrope composition vs pressure")
    ax[0].grid(alpha=0.3)

    ax[1].plot(pressures_kpa, azeotrope_t, "o-", color="brown", linewidth=1.5)
    ax[1].set_xlabel("pressure (kPa)")
    ax[1].set_ylabel("azeotrope temperature (deg C)")
    ax[1].set_title("azeotrope temperature vs pressure")
    ax[1].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_isothermal_excess_bundle(
    output_dir: Path | str,
    x_nrtl: np.ndarray,
    GE_j_mol: np.ndarray,
    HE_j_mol: np.ndarray,
    CPE_j_molK: np.ndarray,
    tag: str,
    dpi: int = 300,
) -> Path:
    """Save bundled isothermal excess-property plot."""
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / f"isothermal_excess_bundle_{tag}.png"

    fig, ax = plt.subplots(3, 1, figsize=(9, 12), sharex=True)

    ax[0].plot(x_nrtl, GE_j_mol / 1000.0, "b-", linewidth=2.0)
    ax[0].set_ylabel("GE (kJ/mol)")
    ax[0].grid(alpha=0.3)

    ax[1].plot(x_nrtl, HE_j_mol / 1000.0, "m-", linewidth=2.0)
    ax[1].set_ylabel("HE (kJ/mol)")
    ax[1].grid(alpha=0.3)

    ax[2].plot(x_nrtl, CPE_j_molK, color="darkorange", linewidth=2.0)
    ax[2].set_xlabel("liquid mole fraction x1")
    ax[2].set_ylabel("CPE (J/mol/K)")
    ax[2].grid(alpha=0.3)

    fig.suptitle("isothermal excess-property bundle")
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path
