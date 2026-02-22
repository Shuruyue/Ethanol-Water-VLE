"""Load parameter data with source metadata checks."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from models import AntoineParams, NRTLParams, VanLaarParams


@dataclass(frozen=True)
class SourceRecord:
    """Citation metadata for a parameter block."""

    name: str
    citation: str
    url: str = ""


@dataclass(frozen=True)
class ComponentParameters:
    """Component-level parameters."""

    name: str
    cas: str
    antoine: AntoineParams
    source: SourceRecord


@dataclass(frozen=True)
class SystemParameters:
    """All parameter blocks required by the calculation pipeline."""

    system_id: str
    system_name: str
    pressure_kpa: float
    component_1: ComponentParameters
    component_2: ComponentParameters
    nrtl: NRTLParams
    van_laar: VanLaarParams
    nrtl_source: SourceRecord
    van_laar_source: SourceRecord


def default_db_path() -> Path:
    """Return the default JSON parameter database path."""
    return Path(__file__).resolve().parents[1] / "data" / "parameter_db.json"


def load_parameter_db(path: Path | None = None) -> dict[str, Any]:
    """Load raw JSON database."""
    db_path = path or default_db_path()
    with db_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _read_source(payload: dict[str, Any], section_name: str) -> SourceRecord:
    source = payload.get("source", {})
    name = str(source.get("name", "")).strip()
    citation = str(source.get("citation", "")).strip()
    if not name or not citation:
        raise ValueError(f"Missing source metadata in section: {section_name}")
    return SourceRecord(name=name, citation=citation, url=str(source.get("url", "")).strip())


def _read_component(component_key: str, payload: dict[str, Any]) -> ComponentParameters:
    source = _read_source(payload, f"component:{component_key}")
    antoine = payload["antoine"]
    params = AntoineParams(
        name=str(payload["name"]),
        A=float(antoine["A"]),
        B=float(antoine["B"]),
        C=float(antoine["C"]),
        unit=str(antoine.get("unit", "mmHg")),
    )
    return ComponentParameters(
        name=str(payload["name"]),
        cas=str(payload["cas"]),
        antoine=params,
        source=source,
    )


def load_system_parameters(system_id: str = "ethanol_water_1atm", path: Path | None = None) -> SystemParameters:
    """
    Load one system and enforce non-empty source metadata.

    Also enforces that NRTL constant terms are explicitly present and not both zero.
    """
    db = load_parameter_db(path)
    systems = db.get("systems", {})
    if system_id not in systems:
        raise KeyError(f"System not found in parameter database: {system_id}")

    system = systems[system_id]
    components = system.get("components", {})
    component_1 = _read_component("component_1", components["component_1"])
    component_2 = _read_component("component_2", components["component_2"])

    nrtl_raw = system["nrtl"]
    nrtl_source = _read_source(nrtl_raw, "nrtl")
    nrtl = NRTLParams(
        alpha12=float(nrtl_raw["alpha12"]),
        alpha21=float(nrtl_raw["alpha21"]),
        a12=float(nrtl_raw["a12"]),
        b12=float(nrtl_raw["b12"]),
        a21=float(nrtl_raw["a21"]),
        b21=float(nrtl_raw["b21"]),
    )
    if nrtl.a12 == 0.0 and nrtl.a21 == 0.0:
        raise ValueError("NRTL constants a12 and a21 cannot both be zero for this project.")

    van_laar_raw = system["van_laar"]
    van_laar_source = _read_source(van_laar_raw, "van_laar")
    van_laar = VanLaarParams(
        A12=float(van_laar_raw["A12"]),
        A21=float(van_laar_raw["A21"]),
    )

    return SystemParameters(
        system_id=system_id,
        system_name=str(system["name"]),
        pressure_kpa=float(system["pressure_kpa"]),
        component_1=component_1,
        component_2=component_2,
        nrtl=nrtl,
        van_laar=van_laar,
        nrtl_source=nrtl_source,
        van_laar_source=van_laar_source,
    )
