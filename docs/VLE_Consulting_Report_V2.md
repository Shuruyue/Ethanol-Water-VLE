# Vapor-Liquid Equilibrium Consulting Report V2
## Non-Ideal Binary System: Ethanol-Water at 1 atm
## Comprehensive Thermodynamic Analysis

**Prepared for:** Management Team, Fine Chemicals & Biotechnology Materials Co.  
**Prepared by:** External Process Engineering Consulting Team  
**Date:** January 2026  
**Document Version:** 2.0  
**Reference Framework:** Sandler, *Chemical, Biochemical, and Engineering Thermodynamics*, 5th Ed.

---

# Table of Contents

1. [Executive Summary](#part-i-executive-summary)
2. [Thermodynamic Fundamentals](#part-ii-thermodynamic-fundamentals)
3. [Mathematical Framework](#part-iii-mathematical-framework)
4. [NRTL Model Implementation](#part-iv-nrtl-model-implementation)
5. [Phase Equilibrium Analysis](#part-v-phase-equilibrium-analysis)
6. [Excess Property Analysis](#part-vi-excess-property-analysis)
7. [Engineering Recommendations](#part-vii-engineering-recommendations)
8. [Technical Appendix](#appendix)

---

# Part I: Executive Summary

## Problem Statement

The client's distillation process for Ethanol-Water separation has experienced operational failures due to **incorrect ideal solution assumptions**. This report provides:

- Complete thermodynamic analysis using the NRTL activity coefficient model
- Validated parameters from ChemSep database
- Quantitative phase equilibrium predictions
- Engineering recommendations for process redesign

## Key Findings

| Parameter | Ideal Model | NRTL Model | Difference |
|-----------|-------------|------------|------------|
| Azeotrope prediction | None | 78.1°C, x₁=0.874 | Critical design limit |
| Max activity coeff. (γ₁) | 1.0 | ~5.7 | 470% deviation |
| Bubble point error | - | Up to 8°C | Temperature control failure |
| Excess Gibbs energy | 0 | 925 J/mol | Non-ideal behavior |
| Excess enthalpy | 0 | 663 J/mol | Heat effects ignored |

## Primary Recommendation

**Redesign required:** Implement pressure-swing or extractive distillation to break the azeotrope at x₁ ≈ 0.874 (95.6 wt% ethanol).

---

# Part II: Thermodynamic Fundamentals

## 2.1 Phase Equilibrium Criterion

At thermodynamic equilibrium between liquid (L) and vapor (V) phases, the fugacity of each component must be equal:

$$\hat{f}_i^{(L)} = \hat{f}_i^{(V)} \quad \text{for all components } i$$

## 2.2 Fugacity Expressions

### Liquid Phase Fugacity

For a component in a liquid mixture:

$$\hat{f}_i^{(L)} = x_i \gamma_i f_i^{0,L}$$

where:
- $x_i$ = liquid mole fraction of component $i$
- $\gamma_i$ = activity coefficient (accounts for non-ideal molecular interactions)
- $f_i^{0,L}$ = fugacity of pure liquid $i$ at system T and P

### Vapor Phase Fugacity

For low-to-moderate pressures (ideal gas assumption):

$$\hat{f}_i^{(V)} = y_i P$$

where:
- $y_i$ = vapor mole fraction
- $P$ = total system pressure

## 2.3 Modified Raoult's Law

Combining the above expressions at equilibrium:

$$y_i P = x_i \gamma_i P_i^{sat}(T)$$

This is the **modified Raoult's Law**, which reduces to the ideal Raoult's Law when $\gamma_i = 1$.

## 2.4 Why Ideal Assumptions Fail

For the Ethanol-Water system, activity coefficients deviate significantly from unity due to:

1. **Hydrogen bond disruption**: Mixing breaks water's tetrahedral H-bond network
2. **Molecular size asymmetry**: Ethanol (46.07 g/mol) vs. Water (18.02 g/mol)
3. **Polarity differences**: Both polar, but different dipole moments

Result: $\gamma_1$ ranges from **1.0 to 5.7** and $\gamma_2$ from **1.0 to 2.7**.

---

# Part III: Mathematical Framework

## 3.1 Excess Gibbs Energy

The excess Gibbs energy quantifies deviation from ideal solution behavior:

$$G^E = G^{mixture} - G^{ideal\ mixture} = RT \sum_i x_i \ln \gamma_i$$

For positive $G^E$: **Positive deviation** from Raoult's Law (γ > 1)  
For negative $G^E$: **Negative deviation** from Raoult's Law (γ < 1)

## 3.2 Relationship Between G^E and Activity Coefficients

From the fundamental excess property relation:

$$\ln \gamma_i = \frac{1}{RT} \left( \frac{\partial (nG^E)}{\partial n_i} \right)_{T,P,n_{j \neq i}}$$

This shows that $G^E$ contains all information needed to compute activity coefficients.

## 3.3 Excess Enthalpy (Heat of Mixing)

The excess enthalpy (approximately equal to heat of mixing for incompressible liquids):

$$H^E = G^E - T \left( \frac{\partial G^E}{\partial T} \right)_P = -RT^2 \sum_i x_i \frac{\partial \ln \gamma_i}{\partial T}$$

This derivative relationship is why **temperature-dependent** activity coefficient parameters are essential for energy calculations.

## 3.4 Antoine Equation for Saturation Pressure

The pure component saturation pressure is computed using the Antoine equation:

$$\log_{10} P_i^{sat} = A_i - \frac{B_i}{T + C_i}$$

where $T$ is in °C and $P^{sat}$ is in mmHg (converted to kPa by multiplying by 0.133322).

### Parameters Used (NIST Chemistry WebBook)

| Component | A | B | C | Valid Range |
|-----------|---|---|---|-------------|
| **Ethanol** | 8.20417 | 1642.89 | 230.3 | 14-93°C |
| **Water** | 8.07131 | 1730.63 | 233.426 | 60-150°C |

---

# Part IV: NRTL Model Implementation

## 4.1 The NRTL (Non-Random Two-Liquid) Model

Developed by **Renon and Prausnitz (1968)**, the NRTL model is based on the local composition concept:

> "The local concentration of molecules around a central molecule differs from the bulk concentration due to molecular interaction energies."

## 4.2 NRTL Equations for Binary Mixtures

### Excess Gibbs Energy

$$\frac{G^E}{RT} = x_1 x_2 \left( \frac{\tau_{21} G_{21}}{x_1 + x_2 G_{21}} + \frac{\tau_{12} G_{12}}{x_2 + x_1 G_{12}} \right)$$

### Activity Coefficients

$$\ln \gamma_1 = x_2^2 \left[ \tau_{21} \left( \frac{G_{21}}{x_1 + x_2 G_{21}} \right)^2 + \frac{\tau_{12} G_{12}}{(x_2 + x_1 G_{12})^2} \right]$$

$$\ln \gamma_2 = x_1^2 \left[ \tau_{12} \left( \frac{G_{12}}{x_2 + x_1 G_{12}} \right)^2 + \frac{\tau_{21} G_{21}}{(x_1 + x_2 G_{21})^2} \right]$$

### G and τ Parameter Relationships

$$G_{ij} = \exp(-\alpha_{ij} \tau_{ij})$$

$$\tau_{ij} = \frac{b_{ij}}{T} \quad \text{(ChemSep format, T in Kelvin)}$$

## 4.3 Physical Interpretation of Parameters

| Parameter | Symbol | Physical Meaning |
|-----------|--------|------------------|
| **Energy parameter** | $\tau_{ij}$ | Dimensionless energy difference between i-j and j-j interactions. Large positive $\tau$ indicates molecular repulsion. |
| **Non-randomness** | $\alpha_{ij}$ | Describes molecular distribution non-randomness. Typically 0.2-0.47 for polar systems. |

## 4.4 ChemSep Database Parameters (Validated)

**Source:** `thermo/Interaction Parameters/ChemSep/nrtl.json`  
**System:** Ethanol (CAS 64-17-5) / Water (CAS 7732-18-5)  
**Validation:** `tests/test_nrtl.py::test_NRTL_chemsep`

| Parameter | Value | Description |
|-----------|-------|-------------|
| α₁₂ = α₂₁ | **0.2937** | Symmetric non-randomness factor |
| b₁₂ | **-29.167 K** | Ethanol→Water interaction |
| b₂₁ | **+624.868 K** | Water→Ethanol interaction |

### Temperature-Dependent τ Calculation

At T = 343.15 K (70°C):

$$\tau_{12} = \frac{-29.167}{343.15} = -0.085$$

$$\tau_{21} = \frac{624.868}{343.15} = +1.821$$

The asymmetry in τ values indicates **stronger repulsion when water molecules are surrounded by ethanol** compared to the reverse.

---

# Part V: Phase Equilibrium Analysis

## 5.1 Bubble Point Calculation Algorithm

For a given liquid composition $x_1$ at pressure P:

1. **Guess** initial temperature T
2. **Calculate** $P_1^{sat}(T)$ and $P_2^{sat}(T)$ from Antoine equation
3. **Calculate** $\tau_{12}$, $\tau_{21}$, $G_{12}$, $G_{21}$ from NRTL parameters
4. **Calculate** $\gamma_1$ and $\gamma_2$ from NRTL equations
5. **Compute** $P_{calc} = x_1 \gamma_1 P_1^{sat} + x_2 \gamma_2 P_2^{sat}$
6. **Iterate** until $|P_{calc} - P| < \epsilon$
7. **Calculate** vapor composition: $y_1 = \frac{x_1 \gamma_1 P_1^{sat}}{P}$

## 5.2 T-x-y Diagram

![T-x-y Diagram](TxY_Ethanol_Water_1atm.png)

**Figure 1:** Temperature-composition diagram at 1 atm comparing Ideal (Raoult's Law), Van Laar, and NRTL models. Key observations:

- **NRTL** correctly predicts the minimum-boiling azeotrope
- **Ideal model** shows monotonic T-x relationship (no azeotrope)
- Maximum temperature deviation: **~8°C** in mid-composition range

## 5.3 y-x Diagram (McCabe-Thiele)

![y-x Diagram](Yx_Ethanol_Water_1atm.png)

**Figure 2:** Vapor-liquid equilibrium y-x diagram. The intersection with y=x diagonal indicates the **azeotropic composition** where no further separation is possible.

## 5.4 Azeotrope Characterization

| Property | Value | Engineering Implication |
|----------|-------|-------------------------|
| Type | Minimum-boiling | Lower T than either pure component |
| Temperature | **78.1°C** | Below pure ethanol (78.4°C) |
| Composition (x₁) | **0.874 mol** | 95.6 wt% ethanol |
| Relative volatility | **α₁₂ = 1.0** | Cannot separate by simple distillation |

---

# Part VI: Excess Property Analysis

## 6.1 Excess Gibbs Energy

![G^E Diagram](GE_vs_x_Ethanol_Water_1atm.png)

**Figure 3:** Excess Gibbs energy versus liquid composition along the bubble point curve.

### Key Results

| Property | Value | Composition |
|----------|-------|-------------|
| Maximum G^E | **924.55 J/mol** | x₁ = 0.425 |
| G^E at azeotrope | ~500 J/mol | x₁ = 0.874 |

**Interpretation:** Positive $G^E$ throughout confirms **positive deviation from ideality**, consistent with azeotrope formation.

## 6.2 Excess Enthalpy (Heat of Mixing)

![H^E Diagram](Hmix_vs_x_Ethanol_Water_1atm.png)

**Figure 4:** Excess enthalpy (heat of mixing) versus composition evaluated along the bubble curve.

### Key Results

| Property | Value | Composition |
|----------|-------|-------------|
| Maximum H^E | **663.22 J/mol** | x₁ = 0.342 |
| Thermal behavior | Endothermic | At bubble point T |

**Note:** The calculated H^E is positive (endothermic) at bubble point temperatures (78-100°C). At lower temperatures (~25°C), Ethanol-Water mixing is typically exothermic. This temperature dependence is captured through the τ = b/T relationship.

## 6.3 Engineering Implications

| Effect | Observation | Process Impact |
|--------|-------------|----------------|
| Heat of mixing | +663 J/mol max | Heat input required during mixing |
| Non-ideal VLE | γ₁ up to 5.7 | Larger column reflux ratio |
| Azeotrope limit | x₁ = 0.874 | Cannot exceed 95.6 wt% by distillation |

---

# Part VII: Engineering Recommendations

## 7.1 Azeotrope-Breaking Strategies

### Option A: Pressure-Swing Distillation

The azeotropic composition shifts with pressure:

| Pressure | Azeotrope x₁ (mol) | T_azeotrope |
|----------|-------------------|-------------|
| 1.0 atm | 0.874 | 78.1°C |
| 0.2 atm | ~0.91 | ~52°C |

**Strategy:** Two columns at different pressures can cross the azeotrope.

### Option B: Extractive Distillation

Add a high-boiling entrainer (e.g., ethylene glycol) that preferentially associates with water, breaking the azeotrope.

### Option C: Pervaporation Hybrid

1. Distill to ~90 mol% ethanol
2. Use membrane pervaporation for final dehydration

## 7.2 Heat Exchanger Redesign

Current design using ideal model **underestimates** heat duty. Corrections needed:

$$Q_{actual} = Q_{ideal} + Q_{mixing}$$

where $Q_{mixing}$ accounts for the 663 J/mol heat effect.

## 7.3 Recommended Action Summary

| Priority | Action | Expected Outcome |
|----------|--------|------------------|
| **1** | Update VLE model to NRTL | Accurate T-x-y predictions |
| **2** | Redesign for azeotrope | Enable >95.6 wt% ethanol |
| **3** | Recalculate heat duties | Proper exchanger sizing |

---

# Appendix

## A.1 Pure Component Properties

| Property | Ethanol | Water |
|----------|---------|-------|
| CAS Number | 64-17-5 | 7732-18-5 |
| Molecular Weight | 46.07 g/mol | 18.02 g/mol |
| Normal Boiling Point | 78.4°C | 100.0°C |
| Critical Temperature | 513.9 K | 647.1 K |
| Critical Pressure | 6.14 MPa | 22.06 MPa |
| Dipole Moment | 1.69 D | 1.85 D |

## A.2 Complete NRTL Parameter Set

**Source:** ChemSep Database (validated by thermo library)

| Parameter | Value | Definition |
|-----------|-------|------------|
| α₁₂ | 0.2937 | Non-randomness factor |
| α₂₁ | 0.2937 | (symmetric) |
| a₁₂ | 0.0 | Constant τ term (ChemSep: none) |
| a₂₁ | 0.0 | |
| b₁₂ | -29.167 K | Temperature coefficient |
| b₂₁ | +624.868 K | |

**τ Definition:** $\tau_{ij} = a_{ij} + \frac{b_{ij}}{T}$ (T in Kelvin)

## A.3 Calculation Validation

| Test Point | Expected γ₁ | Calculated γ₁ | Expected γ₂ | Calculated γ₂ |
|------------|-------------|---------------|-------------|---------------|
| T=343.15K, x₁=0.252 | 1.985 | 1.985 | 1.146 | 1.146 |
| T=343.15K, x₁=0.0 | 5.66 | 5.66 | 1.00 | 1.00 |
| T=343.15K, x₁=1.0 | 1.00 | 1.00 | 2.67 | 2.67 |

**Source:** `tests/test_nrtl.py::test_NRTL_chemsep`

## A.4 References

### Primary Sources

1. **Sandler, S.I.** (2017). *Chemical, Biochemical, and Engineering Thermodynamics* (5th ed.). Wiley.
   - Chapters 6-10: Activity coefficients, VLE theory, excess properties

2. **Renon, H., Prausnitz, J.M.** (1968). Local compositions in thermodynamic excess functions for liquid mixtures. *AIChE Journal*, 14(1), 135-144.
   - Original NRTL model derivation

### Data Sources

3. **ChemSep NRTL Parameter Database**
   - File: `thermo/Interaction Parameters/ChemSep/nrtl.json`
   - Format: τ = b/T (Kelvin)
   - Validated: `tests/test_nrtl.py::test_NRTL_chemsep`

4. **NIST Chemistry WebBook**
   - Antoine equation parameters
   - https://webbook.nist.gov/

### Computational Tools

5. **thermo Python Library** (Caleb Bell)
   - NRTL class implementation
   - https://github.com/CalebBell/thermo

6. **SciPy** - Brent's method for bubble point iteration

---

*Report generated using validated thermodynamic parameters from the ChemSep database. All calculations verified against thermo library test cases.*

**Document Control:**
- Version: 2.0
- Generated: January 2026
- Parameters validated: ✓
