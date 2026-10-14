# Ethanol-Water VLE Final Project

This folder contains the computation code for the chemical thermodynamics final project.
The scope is ethanol(1)-water(2) with non-ideal liquid behavior.

## Current Capabilities

- VLE prediction at fixed pressure: `T-x-y`, `y-x`
- Excess properties from NRTL: `GE`, `HE`, `CPE`
- Model comparison: ideal, Van Laar, NRTL
- Interactive dashboard (Wolfram-style exploration)
- Practical image-pack export for engineering decisions
- Final build pipeline: plots + CSV tables + summary JSON
- Shared profile engine for CLI/export/UI consistency (`src/profiles.py`)
- Interpolated azeotrope detection (sub-grid estimate)
- Isothermal export includes `alpha12` (relative volatility) directly
- Dashboard supports one-click CSV download of current state

## Directory

```text
vle_consulting/
├── data/
│   └── parameter_db.json
├── figures/
├── references/
│   └── thermo.pdf
├── src/
│   ├── main.py
│   ├── interactive_app.py
│   ├── output_pack.py
│   ├── export_data.py
│   ├── final_build.py
│   ├── pipeline.py
│   ├── parameter_store.py
│   ├── models.py
│   ├── solver.py
│   ├── analysis.py
│   ├── profiles.py
│   └── plotting.py
└── tests/
    ├── conftest.py
    ├── test_models.py
    ├── test_solver.py
    ├── test_analysis.py
    ├── test_plotting.py
    ├── test_pipeline.py
    ├── test_parameter_store.py
    ├── test_output_pack.py
    ├── test_export_data.py
    └── verify_scientific_correctness.py
```

## Install

```bash
pip install -r requirements.txt
```

## Run Baseline CLI

```bash
python src/main.py
python src/main.py --points 80
python src/main.py --output-dir figures
python src/main.py --no-plots
```

## Run Interactive Dashboard

```bash
streamlit run src/interactive_app.py
```

Interactive controls:

- Pressure (kPa)
- Reference temperature (deg C)
- Grid size
- Full NRTL parameter set: `alpha12`, `alpha21`, `a12`, `b12`, `a21`, `b21`

## Generate Practical Output Pack

```bash
python src/output_pack.py --output-dir figures/practical_pack --t-ref 78.2 --points 80
```

Default output set (recommended for real-world discussion):

1. `txy_ethanol_water_1atm.png`: boiling/dew temperature design map
2. `yx_ethanol_water_1atm.png`: separation trajectory map
3. `ge_ethanol_water_1atm.png`: non-ideality intensity
4. `he_ethanol_water_1atm.png`: thermal effect of mixing
5. `isothermal_excess_bundle_tXXc.png`: `GE`, `HE`, `CPE` at chosen temperature
6. `cpe_ethanol_water_tXXc.png`: excess heat-capacity signal
7. `gamma_ethanol_water_tXXc.png`: activity coefficients vs composition
8. `relative_volatility_ethanol_water_tXXc.png`: separation driving force map
9. `azeotrope_pressure_sensitivity.png`: how azeotrope location shifts with pressure

## Build Final Deliverable Bundle

```bash
python src/final_build.py --output-root final_outputs --points 80 --t-ref 78.2
```

This generates:

- `final_outputs/plots/`:
  baseline and practical-pack figures
- `final_outputs/data/baseline/`:
  model curves and bubble-line excess properties CSV
- `final_outputs/data/isothermal/`:
  fixed-temperature `GE/HE/CPE/gamma/alpha12` CSV
- `final_outputs/data/sensitivity/`:
  azeotrope-pressure sensitivity CSV
- `final_outputs/data/run_summary.json`:
  key metrics + relative-volatility range + azeotrope-vs-literature offsets

## Optimization Notes (Current Version)

### 1) Architecture

- Unified repeated post-processing into `src/profiles.py`.
- `output_pack.py`, `export_data.py`, and `interactive_app.py` now use the same profile builders.
- Reduced divergence risk between CLI plots, exported CSV, and dashboard values.

### 2) Core Algorithms

- Added vectorized NRTL evaluator (`gamma_nrtl_binary`) and reused it in excess-property curve generation.
- Added Antoine derivative function `dPsat/dT` for solver diagnostics and derivative-ready extensions.
- Bubble-point solver now supports continuation-style bracketing (`T_guess`) to improve robustness.
- Azeotrope finder now interpolates across `y-x` sign changes for sub-grid accuracy.

### 3) UI and Output

- Dashboard now overlays Van Laar against Ideal/NRTL in both `T-x-y` and `y-x`.
- Added literature azeotrope marker for quick visual sanity check.
- Added dashboard CSV download for reporting/review handoff.
- Extended summary JSON with relative-volatility envelope and azeotrope deviation metrics.

### 4) Verification Loop

- Unit/regression tests: `pytest -q tests` (45 passed, incl. validation against experimental data)
- Scientific verification script:
  - Windows recommendation: `set PYTHONUTF8=1` then `python tests/verify_scientific_correctness.py`

## Data Rules

- Parameters are loaded only from `data/parameter_db.json`.
- Every parameter block must include source metadata.
- NRTL uses explicit constant and temperature terms:
  - `tau_ij(T) = a_ij + b_ij / T`
- Loader rejects datasets where both `a12` and `a21` are zero.

## Real-World Gaps To Add Next

To improve real-plant fidelity, add these datasets next:

- pressure-dependent vapor-phase non-ideality (`phi` correction)
- liquid density and heat-capacity data for full energy balance
- validated pressure-swing azeotrope points from literature
- uncertainty bounds for NRTL parameters (confidence intervals)
- pilot or plant composition-temperature samples for model calibration

## References Used In Database

- Antoine coefficients: common tabulated sets in mmHg / deg C (ethanol valid -57 to 80 deg C,
  water -20 to 100 deg C). These are **not** NIST WebBook sets (NIST uses bar / K);
  the 1 atm water-rich end (up to 100 deg C) is an extrapolation for ethanol.
- NRTL parameters: Aspen Plus databank values regressed from the Dortmund Data Bank,
  as tabulated in Beneke, Peters, Glasser & Hildebrandt, *Understanding Distillation
  Using Column Profile Maps*, Wiley (2013), Appendix B. (Earlier versions of this README
  attributed them to DECHEMA Vol. I, which is not correct.)
- Van Laar parameters: the illustrative values from the course handout code (not literature); kept only for model comparison.
- Experimental validation data: Kamihama et al., *J. Chem. Eng. Data* 57 (2012) 339-344,
  doi:10.1021/je2008704 (`data/experimental/`, via NIST TRC ThermoML).
- Course textbook: Sandler, *Chemical, Biochemical, and Engineering Thermodynamics*, 5th ed., Wiley (2017).

## 2026 Report Revision (`report2026/`)

`report2026/` is a self-contained, fully reproducible rewrite of the final report
(XeLaTeX + Python). Every number in the PDF is generated from `results/results.json`.

```bash
cd report2026
pip install numpy scipy pandas matplotlib chemicals CoolProp
python calc/run_analysis.py && python calc/make_figures.py
python calc/make_numbers.py && python calc/make_tables.py
cd tex && xelatex report.tex && xelatex report.tex
```

## Algorithm References

- Renon, H.; Prausnitz, J. M. (1968). *Local compositions in thermodynamic excess functions for liquid mixtures*. AIChE Journal, 14(1), 135-144. (NRTL origin)
- Michelsen, M. L. (1982). *The Isothermal Flash Problem. Part I. Stability*. Fluid Phase Equilibria, 9, 1-19.
- Michelsen, M. L. (1982). *The Isothermal Flash Problem. Part II. Phase-Split Calculation*. Fluid Phase Equilibria, 9, 21-40.
- SciPy documentation for Brent root solver (`scipy.optimize.brentq`) used for robust 1D root finding.
