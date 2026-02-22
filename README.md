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

## Directory

```text
vle_consulting/
├── data/
│   └── parameter_db.json
├── figures/
├── plans/
│   └── five_phase_plan.md
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
│   └── plotting.py
└── tests/
    ├── test_parameter_store.py
    ├── test_pipeline.py
    └── test_output_pack.py
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
  fixed-temperature `GE/HE/CPE/gamma` CSV
- `final_outputs/data/sensitivity/`:
  azeotrope-pressure sensitivity CSV
- `final_outputs/data/run_summary.json`:
  key metrics + parameter-source metadata

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

- NIST Chemistry WebBook for Antoine coefficients
- DECHEMA data collection (Gmehling et al.) for NRTL parameters
- Sandler textbook context for course setting
