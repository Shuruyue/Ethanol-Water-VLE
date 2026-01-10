# Transport Phenomena Laboratory
## Volatile Liquid Evaporation Analysis

A Python module for analyzing mass transfer in volatile liquid evaporation experiments,
supporting both natural and forced convection conditions.

---

## Module Structure

```
transport_phenomena_lab/
├── mass_transfer/           # Core calculation module
│   ├── __init__.py         # Module exports
│   ├── properties.py       # Chemical properties & Antoine equation
│   ├── diffusion.py        # Diffusion coefficients & fluxes
│   ├── mass_transfer_coeff.py  # Mass transfer correlations
│   └── analysis.py         # Data analysis & plotting
├── data/                   # Experimental data (CSV)
├── figures/                # Generated plots
├── reports/                # Final reports
└── README.md               # This file
```

---

## Quick Start

### Installation

```bash
pip install numpy matplotlib python-pptx
```

### Basic Usage

```python
from mass_transfer import (
    diffusion_coefficient_fsg,
    mass_transfer_natural,
    mass_transfer_forced,
    compare_convection_modes,
    ETHANOL,
)

# Calculate diffusion coefficient (ethanol-air at 25°C)
D_ab = diffusion_coefficient_fsg(T_celsius=25.0, P_atm=1.0, props_A=ETHANOL)
print(f"D_ab = {D_ab:.2e} m²/s")

# Compare natural vs forced convection
result = compare_convection_modes(
    T_celsius=25.0,
    L_char=0.08,       # 8 cm dish diameter
    u_air=1.5,         # 1.5 m/s fan velocity
    surface_area=0.005,
    initial_mass=15.0,
)
print(f"Enhancement factor: {result.enhancement_factor:.2f}x")
```

### Analyze Experimental Data & Generate Report

```python
# Run entire analysis and generate PowerPoint report
# python main.py
```

### Analyze Experimental Data

```python
from mass_transfer import load_experiment_data, analyze_experiment

# Load and analyze data
data = load_experiment_data("data/natural_convection.csv")
results = analyze_experiment(data)
print(f"Mass transfer coefficient: {results.calculated_km_ms:.4e} m/s")
```

---

## Key Equations

| Parameter | Equation | Unit |
|-----------|----------|------|
| Diffusion Coefficient | Fuller-Schettler-Giddings | m²/s |
| Vapor Pressure | Antoine Equation | kPa |
| Stefan Flux | $N_A = \frac{D_{AB}P}{RTz}\ln\frac{P-P_{A\infty}}{P-P_{A0}}$ | mol/(m²·s) |
| Mass Transfer Coeff. | Sherwood correlations | m/s |

---

## Data File Format

CSV files with metadata header:

```csv
# name: experiment_natural_run1
# condition: natural
# temperature_C: 25.0
# humidity_percent: 50.0
# surface_area_m2: 0.005
time_min,mass_g
0,15.00
5,14.92
10,14.84
...
```

---

## References

- Bird, R.B., Stewart, W.E., Lightfoot, E.N. *Transport Phenomena*, 2nd Ed.
- Fuller, E.N., Schettler, P.D., Giddings, J.C. (1966). *Ind. Eng. Chem.*, 58(5), 18-27.
- Perry's Chemical Engineers' Handbook, 9th Ed.

---

## Author

VLE Consulting Team - Transport Phenomena II Self-Study Assignment
