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


def plot_mccabe_thiele(
    output_dir: Path | str,
    x_nrtl: np.ndarray,
    y_nrtl: np.ndarray,
    xD: float = 0.80,
    xB: float = 0.05,
    xF: float = 0.40,
    q: float = 1.0,
    R_reflux: float = 2.0,
    max_stages: int = 50,
    dpi: int = 300,
) -> Path:
    """Save McCabe-Thiele step diagram for binary distillation.

    Parameters
    ----------
    x_nrtl, y_nrtl : arrays
        VLE equilibrium curve (x, y for ethanol).
    xD : float
        Distillate composition (mol fraction ethanol).
    xB : float
        Bottoms composition (mol fraction ethanol).
    xF : float
        Feed composition.
    q : float
        Feed quality (1 = saturated liquid, 0 = saturated vapor).
    R_reflux : float
        Reflux ratio (L/D).
    max_stages : int
        Maximum stepping stages for safety.
    """
    from scipy.interpolate import interp1d

    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / "mccabe_thiele_ethanol_water.png"

    # Interpolation for y = f(x)
    y_eq = interp1d(x_nrtl, y_nrtl, kind="linear", fill_value="extrapolate")

    # Rectifying operating line: y = (R/(R+1))*x + xD/(R+1)
    slope_rect = R_reflux / (R_reflux + 1.0)
    intercept_rect = xD / (R_reflux + 1.0)

    # q-line intersection with rectifying line
    if abs(q - 1.0) < 1e-6:
        x_intersect = xF
    else:
        # q-line: y = (q/(q-1))*x - xF/(q-1)
        slope_q = q / (q - 1.0)
        intercept_q = -xF / (q - 1.0)
        x_intersect = (intercept_rect - intercept_q) / (slope_q - slope_rect)
    y_intersect = slope_rect * x_intersect + intercept_rect

    # Stripping operating line: through (xB, xB) and (x_intersect, y_intersect)
    slope_strip = (y_intersect - xB) / max(x_intersect - xB, 1e-10)
    intercept_strip = xB - slope_strip * xB

    # Stepping from xD down
    stages_x = [xD]
    stages_y = [xD]
    x_curr = xD
    n_stages = 0

    for _ in range(max_stages):
        # Horizontal line to equilibrium curve
        y_curr = float(y_eq(x_curr))
        if y_curr <= xB:
            break
        stages_x.extend([x_curr, x_curr])
        stages_y.extend([x_curr, y_curr])

        # Vertical line down to operating line
        if x_curr >= x_intersect:
            x_next = (y_curr - intercept_rect) / max(slope_rect, 1e-10)
        else:
            x_next = (y_curr - intercept_strip) / max(slope_strip, 1e-10)

        stages_x.extend([x_curr, x_next])
        stages_y.extend([y_curr, y_curr])
        x_curr = x_next
        n_stages += 1

        if x_curr <= xB:
            break

    # Plot
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.plot([0, 1], [0, 1], "k:", linewidth=0.8, label="y = x")
    ax.plot(x_nrtl, y_nrtl, "r-", linewidth=2.0, label="NRTL equilibrium")

    # Operating lines
    x_rect = np.linspace(x_intersect, xD, 50)
    ax.plot(x_rect, slope_rect * x_rect + intercept_rect, "b--", linewidth=1.5, label="rectifying")
    x_strip = np.linspace(xB, x_intersect, 50)
    ax.plot(x_strip, slope_strip * x_strip + intercept_strip, "g--", linewidth=1.5, label="stripping")

    # Steps
    ax.plot(stages_x, stages_y, color="orange", linewidth=1.0, alpha=0.8)

    # Annotations
    ax.axvline(xD, color="blue", linewidth=0.5, linestyle=":", alpha=0.5)
    ax.axvline(xB, color="green", linewidth=0.5, linestyle=":", alpha=0.5)
    ax.axvline(xF, color="gray", linewidth=0.5, linestyle=":", alpha=0.5)
    ax.text(xD, 0.02, f"xD={xD:.2f}", fontsize=8, ha="center")
    ax.text(xB, 0.02, f"xB={xB:.2f}", fontsize=8, ha="center")
    ax.text(xF, 0.02, f"xF={xF:.2f}", fontsize=8, ha="center")

    ax.set_xlabel("liquid mole fraction x1 (ethanol)")
    ax.set_ylabel("vapor mole fraction y1 (ethanol)")
    ax.set_title(f"McCabe-Thiele diagram ({n_stages} stages, R={R_reflux})")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_aspect("equal")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return out_path


def plot_txy_with_experiment(
    output_dir: Path | str,
    x_nrtl: np.ndarray,
    y_nrtl: np.ndarray,
    t_nrtl: np.ndarray,
    exp_x: np.ndarray | None = None,
    exp_y: np.ndarray | None = None,
    exp_t: np.ndarray | None = None,
    data_source: str = "experimental",
    add_literature_azeotrope: bool = True,
    dpi: int = 300,
) -> Path:
    """Save T-x-y plot with experimental data overlay.

    Parameters
    ----------
    exp_x, exp_y, exp_t : arrays, optional
        Experimental liquid composition, vapor composition, and temperature.
    data_source : str
        Label for the experimental data source.
    """
    out_dir = _ensure_output_dir(output_dir)
    out_path = out_dir / "txy_model_vs_experiment.png"

    plt.figure(figsize=(10, 7))
    plt.plot(x_nrtl, t_nrtl, "r-", linewidth=2.0, label="NRTL bubble")
    plt.plot(y_nrtl, t_nrtl, "r--", linewidth=2.0, label="NRTL dew")

    if exp_x is not None and exp_t is not None:
        plt.scatter(exp_x, exp_t, marker="o", c="blue", s=40, zorder=5,
                    label=f"{data_source} (liquid)", edgecolors="darkblue", linewidths=0.5)
    if exp_y is not None and exp_t is not None:
        plt.scatter(exp_y, exp_t, marker="^", c="cyan", s=40, zorder=5,
                    label=f"{data_source} (vapor)", edgecolors="darkblue", linewidths=0.5)
    if add_literature_azeotrope:
        plt.scatter(
            [0.8943],
            [78.15],
            marker="D",
            c="black",
            s=32,
            zorder=6,
            label="literature azeotrope",
        )

    plt.xlabel("ethanol mole fraction (x1 or y1)")
    plt.ylabel("temperature (deg C)")
    plt.title("ethanol-water T-x-y: model vs experiment")
    plt.xlim(0.0, 1.0)
    plt.grid(alpha=0.3)
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close()
    return out_path
