"""Experimental VLE data used to validate the models.

The default data set is the isobaric ethanol(1)-water(2) data of
Kamihama et al., J. Chem. Eng. Data 57 (2012) 339-344 (doi:10.1021/je2008704),
transcribed from the NIST TRC ThermoML archive, including the TRC expanded
(k = 2) uncertainties of T and y1.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

DEFAULT_FILE = "kamihama2012_ethanol_water_101p3kPa.csv"
DEFAULT_LABEL = "Kamihama et al. (2012)"


@dataclass(frozen=True)
class ExperimentalTxy:
    """Isobaric T-x-y data (T in deg C, P in kPa)."""

    x1: np.ndarray
    y1: np.ndarray
    T_celsius: np.ndarray
    pressure_kpa: float
    U_T: np.ndarray
    U_y1: np.ndarray
    label: str


def default_data_path() -> Path:
    return Path(__file__).resolve().parents[1] / "data" / "experimental" / DEFAULT_FILE


def load_experimental_txy(path: Path | None = None, label: str = DEFAULT_LABEL) -> ExperimentalTxy:
    """Load an isobaric data file (``#`` comment lines are skipped)."""
    data_path = path or default_data_path()
    with data_path.open("r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(line for line in handle if not line.lstrip().startswith("#")))
    if not rows:
        raise ValueError(f"No data rows in {data_path}")
    col = lambda key: np.array([float(r[key]) for r in rows])  # noqa: E731
    pressures = col("P_kPa")
    if np.ptp(pressures) > 1e-9:
        raise ValueError("Experimental file must be isobaric.")
    return ExperimentalTxy(
        x1=col("x1"),
        y1=col("y1"),
        T_celsius=col("T_K") - 273.15,
        pressure_kpa=float(pressures[0]),
        U_T=col("U_T_K"),
        U_y1=col("U_y1"),
        label=label,
    )
