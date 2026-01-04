# VLE Consulting Project
## Ethanol-Water Non-Ideal Mixing System Analysis

A comprehensive thermodynamic analysis of the Ethanol-Water binary system using
NRTL (Non-Random Two-Liquid) activity coefficient model for vapor-liquid equilibrium
calculations.

## Project Structure

```
vle_consulting/
├── src/                    # Source code
│   ├── __init__.py
│   ├── vle_consulting_report.py  # Main calculation script
│   └── vle_calculations.py       # Alternative calculation module
├── docs/                   # Documentation
│   └── VLE_Consulting_Report.md  # Complete consulting report
├── figures/                # Generated plots
│   ├── TxY_Ethanol_Water_1atm.png
│   ├── Yx_Ethanol_Water_1atm.png
│   ├── GE_vs_x_Ethanol_Water_1atm.png
│   └── Hmix_vs_x_Ethanol_Water_1atm.png
├── data/                   # Data files (if any)
├── references/             # Reference materials
│   ├── thermo.pdf          # Original assignment
│   └── PDF詳細內容.txt      # Assignment details
├── README.md               # This file
├── requirements.txt        # Python dependencies
└── config.yaml             # Configuration parameters
```

## Quick Start

### Installation

```bash
pip install -r requirements.txt
```

### Generate VLE Plots

```bash
cd src
python vle_consulting_report.py
```

Output files will be saved to the `figures/` directory.

## System Parameters

| Parameter | Value | Source |
|-----------|-------|--------|
| Pressure | 101.3 kPa (1 atm) | Specification |
| Antoine (Ethanol) | A=8.20417, B=1642.89, C=230.3 | NIST |
| Antoine (Water) | A=8.07131, B=1730.63, C=233.426 | NIST |
| NRTL α₁₂ = α₂₁ | 0.30 | Literature |
| NRTL τ₁₂ | a₁₂=0.88, b₁₂=0 K | Literature |
| NRTL τ₂₁ | a₂₁=1.45, b₂₁=0 K | Literature |

## Models Implemented

1. **Raoult's Law** - Ideal solution assumption
2. **Van Laar Model** - Simple activity coefficient model
3. **NRTL Model** - Recommended for polar mixtures

## Key Findings

- **Azeotrope**: Detected at x₁ ≈ 0.89 (95.6 wt% ethanol), T ≈ 78.2°C
- **Positive Deviation**: G^E > 0 throughout composition range
- **Separation Limit**: Simple distillation cannot exceed 95.6 wt% ethanol

## Report Sections

| Part | Weight | Description |
|------|--------|-------------|
| I | 10% | Executive Summary |
| II | 20% | Problem Definition & Theory |
| III | 30% | NRTL Analysis & VLE Comparison |
| IV | 25% | Excess Properties & Energy |
| V | 15% | Engineering Recommendations |

## References

- Sandler, S.I., *Chemical, Biochemical, and Engineering Thermodynamics*, Ch. 6-10
- Renon, H. & Prausnitz, J.M. (1968). Local compositions in thermodynamic excess functions for liquid mixtures. *AIChE Journal*, 14(1), 135-144.

## Author

VLE Consulting Team - January 2026
