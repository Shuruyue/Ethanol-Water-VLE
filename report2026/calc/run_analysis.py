"""Full analysis for the ethanol-water consulting report.

Outputs
  ../results/results.json   every number quoted in the report
  ../results/*.csv           tabulated curves
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq, least_squares

from vle_core import (DHVAP_NBP, NRTL, NRTL_LIT, R, EqCurve, Ideal, VanLaar, Virial,
                      azeotrope, bubble_T, mole_to_mass, nmin_total_reflux, psat,
                      rate_column, rmin_cmo, step_stages, tsat, txy)

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
RES.mkdir(exist_ok=True)
P_ATM = 101.325
VIR = Virial()
OUT: dict = {}

# ---------------------------------------------------------------------------
# 0. Data and pure-component checks
# ---------------------------------------------------------------------------
d = pd.read_csv(ROOT / "data" / "kamihama2012_ethanol_water_101p3kPa.csv", comment="#")
P_EXP = float(d.P_kPa.iloc[0])
OUT["pure"] = {
    "Tb_ethanol_K": tsat(0, P_ATM),
    "Tb_water_K": tsat(1, P_ATM),
}
# independent check of the ethanol Antoine set against the reference EOS of
# Schroeder et al. (2014) as implemented in CoolProp (Bell et al. 2014)
try:
    import CoolProp.CoolProp as _CP
    _Tg = np.linspace(293.15, 373.15, 81)
    _dev = [100 * (psat(0, t) / (_CP.PropsSI("P", "T", t, "Q", 0, "Ethanol") / 1000) - 1) for t in _Tg]
    OUT["pure"]["etoh_psat_dev_vs_eos_max_pct_293_373"] = float(np.max(np.abs(_dev)))
    OUT["pure"]["etoh_psat_dev_vs_eos_pct_393"] = float(100 * (psat(0, 393.15) / (_CP.PropsSI("P", "T", 393.15, "Q", 0, "Ethanol") / 1000) - 1))
    OUT["pure"]["etoh_Tb_eos_K"] = float(_CP.PropsSI("T", "P", P_ATM * 1000, "Q", 0, "Ethanol"))
except ImportError:  # CoolProp optional
    pass

# ---------------------------------------------------------------------------
# 1. Parameter regression (weighted least squares on T and y)
# ---------------------------------------------------------------------------
sT = d.U_T_K.values / 2.0
sy = d.U_y1.values / 2.0


def predict(model, P=P_EXP, vap=Virial(enabled=False)):
    r = np.array([bubble_T(x, P, model, vap) for x in d.x1.values])
    return r[:, 0], r[:, 1]


def residuals(model):
    T, y = predict(model)
    return np.concatenate([(T - d.T_K.values) / sT, (y - d.y1.values) / sy])


def fit(builder, p0):
    sol = least_squares(lambda p: residuals(builder(p)), p0, x_scale="jac")
    J = sol.jac
    n, p = J.shape
    s2 = 2 * sol.cost / (n - p)
    cov = np.linalg.inv(J.T @ J) * s2
    se = np.sqrt(np.diag(cov))
    corr = cov[0, 1] / (se[0] * se[1])
    return sol.x, se, corr, s2


vl_build = lambda p: VanLaar(A12=p[0], A21=p[1], name="Van Laar (this work)")
nr_build = lambda p: NRTL(a12=0.0, b12=p[0], a21=0.0, b21=p[1], alpha=0.3, name="NRTL (this work)")
p_vl, se_vl, c_vl, s2_vl = fit(vl_build, [1.7, 0.95])
p_nr, se_nr, c_nr, s2_nr = fit(nr_build, [-35.0, 630.0])
VL_FIT = vl_build(p_vl)
NR_FIT = nr_build(p_nr)
OUT["regression"] = {
    "van_laar": {"A12": p_vl[0], "A21": p_vl[1], "se": se_vl.tolist(), "corr": c_vl, "chi2_red": s2_vl},
    "nrtl_b_only": {"b12_K": p_nr[0], "b21_K": p_nr[1], "se": se_nr.tolist(), "corr": c_nr, "chi2_red": s2_nr},
    "n_points": int(len(d)),
}

# ---------------------------------------------------------------------------
# 2. Validation statistics
# ---------------------------------------------------------------------------
models = {
    "ideal": (Ideal(), Virial(enabled=False)),
    "van_laar_fit": (VL_FIT, Virial(enabled=False)),
    "nrtl_lit": (NRTL_LIT, Virial(enabled=False)),
    "nrtl_fit": (NR_FIT, Virial(enabled=False)),
    "nrtl_lit_virial": (NRTL_LIT, VIR),
}
val = {}
resid_table = {"x1": d.x1.values}
for key, (m, vap) in models.items():
    T, y = predict(m, vap=vap)
    dT, dy = T - d.T_K.values, y - d.y1.values
    val[key] = {
        "AAD_T_K": float(np.mean(np.abs(dT))),
        "max_abs_T_K": float(np.max(np.abs(dT))),
        "bias_T_K": float(np.mean(dT)),
        "AAD_y": float(np.mean(np.abs(dy))),
        "max_abs_y": float(np.max(np.abs(dy))),
        "frac_within_U_T": float(np.mean(np.abs(dT) <= d.U_T_K.values)),
        "frac_within_U_y": float(np.mean(np.abs(dy) <= d.U_y1.values)),
    }
    resid_table[f"dT_{key}"] = dT
    resid_table[f"dy_{key}"] = dy
OUT["validation"] = val
pd.DataFrame(resid_table).to_csv(RES / "residuals.csv", index=False)

# ---------------------------------------------------------------------------
# 3. Experimental activity coefficients and Herington test
# ---------------------------------------------------------------------------
g1e, g2e = [], []
for x, y, T in zip(d.x1, d.y1, d.T_K):
    P1, P2 = VIR.Phi(y, T, P_EXP)
    g1e.append(y * P1 * P_EXP / (x * psat(0, T)))
    g2e.append((1 - y) * P2 * P_EXP / ((1 - x) * psat(1, T)))
g1e, g2e = np.array(g1e), np.array(g2e)
pd.DataFrame({"x1": d.x1, "T_K": d.T_K, "gamma1_exp": g1e, "gamma2_exp": g2e}).to_csv(
    RES / "gamma_experimental.csv", index=False)

lnr = np.log(g1e / g2e)
coef = np.polyfit(d.x1.values, lnr, 3)
xx = np.linspace(0, 1, 20001)
f = np.polyval(coef, xx)
area_net = np.trapezoid(f, xx)
area_abs = np.trapezoid(np.abs(f), xx)
D_her = 100 * abs(area_net) / area_abs
Tmax, Tmin = tsat(1, P_EXP), d.T_K.min()
J_her = 150 * (Tmax - Tmin) / Tmin
OUT["herington"] = {"D": D_her, "J": J_her, "D_minus_J": D_her - J_her,
                    "consistent": bool(D_her - J_her < 10), "poly_coef": coef.tolist()}

# ---------------------------------------------------------------------------
# 4. Continuous curves at 1 atm: ideal vs NRTL, gamma, alpha
# ---------------------------------------------------------------------------
xg = np.concatenate([np.linspace(0, 0.05, 101), np.linspace(0.05, 1, 381)[1:]])
curves = {"x1": xg}
for key, m in {"ideal": Ideal(), "van_laar_fit": VL_FIT, "nrtl_lit": NRTL_LIT, "nrtl_fit": NR_FIT}.items():
    x_, y_, T_ = txy(m, P_ATM, xgrid=xg)
    curves[f"y_{key}"], curves[f"T_{key}"] = y_, T_
    l1, l2 = m.lngamma(xg, T_)
    curves[f"lng1_{key}"], curves[f"lng2_{key}"] = l1, l2
    with np.errstate(divide="ignore", invalid="ignore"):
        curves[f"alpha_{key}"] = np.exp(l1) * np.array([psat(0, t) for t in T_]) / (
            np.exp(l2) * np.array([psat(1, t) for t in T_]))
cv = pd.DataFrame(curves)
cv.to_csv(RES / "txy_curves_1atm.csv", index=False)

dT_id = cv.T_ideal - cv.T_nrtl_lit
k = int(np.argmax(dT_id))
az_lit = azeotrope(NRTL_LIT, P_ATM)
az_lit_v = azeotrope(NRTL_LIT, P_ATM, VIR)
az_fit = azeotrope(NR_FIT, P_ATM)
az_vl = azeotrope(VL_FIT, P_ATM)
inner = (cv.x1 > 0.0) & (cv.x1 < 1.0)
dd = (d.y1 - d.x1).values
ic = np.where(np.sign(dd[:-1]) != np.sign(dd[1:]))[0][0]
x_az_exp = d.x1.values[ic] + (d.x1.values[ic + 1] - d.x1.values[ic]) * dd[ic] / (dd[ic] - dd[ic + 1])
T_az_exp = d.T_K.values[ic] + (d.T_K.values[ic + 1] - d.T_K.values[ic]) * dd[ic] / (dd[ic] - dd[ic + 1])
OUT["azeotrope_exp_interp"] = {"x1": float(x_az_exp), "T_K": float(T_az_exp), "wt": mole_to_mass(x_az_exp)}
OUT["one_atm"] = {
    "max_dT_ideal_minus_nrtl_K": float(dT_id[k]), "at_x1": float(cv.x1[k]),
    "T_ideal_at": float(cv.T_ideal[k]), "T_nrtl_at": float(cv.T_nrtl_lit[k]),
    "dT_at_x_0p5": float(np.interp(0.5, cv.x1, dT_id)),
    "azeotrope_nrtl_lit": {"x1": az_lit[0], "T_K": az_lit[1], "wt": mole_to_mass(az_lit[0])},
    "azeotrope_nrtl_lit_virial": {"x1": az_lit_v[0], "T_K": az_lit_v[1]},
    "azeotrope_nrtl_fit": {"x1": az_fit[0], "T_K": az_fit[1]},
    "azeotrope_van_laar_fit": None if az_vl is None else {"x1": az_vl[0], "T_K": az_vl[1]},
    "gamma1_inf_nrtl": float(np.exp(cv.lng1_nrtl_lit.iloc[0])),
    "gamma2_inf_nrtl": float(np.exp(cv.lng2_nrtl_lit.iloc[-1])),
    "T_at_x0_K": float(cv.T_nrtl_lit.iloc[0]),
    "T_at_x1_K": float(cv.T_nrtl_lit.iloc[-1]),
    "alpha_nrtl_max": float(cv.alpha_nrtl_lit[inner].max()),
    "alpha_nrtl_at_x0p9": float(np.interp(0.9, cv.x1, cv.alpha_nrtl_lit)),
    "alpha_ideal_range": [float(cv.alpha_ideal[inner].min()), float(cv.alpha_ideal[inner].max())],
}

# ---------------------------------------------------------------------------
# 5. Vapour-phase non-ideality
# ---------------------------------------------------------------------------
vrows = []
for x in [0.02, 0.1, 0.3, 0.5, 0.7, 0.9]:
    T0, y0 = bubble_T(x, P_ATM, NRTL_LIT)
    T1, y1 = bubble_T(x, P_ATM, NRTL_LIT, VIR)
    P1, P2 = VIR.Phi(y1, T1, P_ATM)
    vrows.append({"x1": x, "T_ig_K": T0, "T_vir_K": T1, "dT_K": T1 - T0, "y_ig": y0, "y_vir": y1,
                  "dy": y1 - y0, "Phi1": P1, "Phi2": P2})
vdf = pd.DataFrame(vrows)
vdf.to_csv(RES / "vapour_correction.csv", index=False)
B = VIR.B(351.3)
OUT["virial"] = {"B11_B22_B12_cm3mol_at_351K": list(map(float, B)),
                 "max_abs_dT_K": float(vdf.dT_K.abs().max()), "max_abs_dy": float(vdf.dy.abs().max()),
                 "rows": vrows}

# ---------------------------------------------------------------------------
# 6. Excess functions
# ---------------------------------------------------------------------------
xe = np.linspace(0.0, 1.0, 501)
Tref = az_lit[1]
ex = {"x1": xe}
for key, m in {"nrtl_lit": NRTL_LIT, "nrtl_fit": NR_FIT}.items():
    gE = R * Tref * m.gE_RT(xe, Tref)
    hE = m.hE(xe, Tref)
    ex[f"gE_{key}"], ex[f"hE_{key}"], ex[f"TsE_{key}"] = gE, hE, hE - gE
    ex[f"cpE_{key}"] = m.cpE(xe, Tref)
    for Tk in (298.15, 323.15, 348.15, 373.15):
        ex[f"hE_{key}_{Tk:.2f}"] = m.hE(xe, Tk)
# bubble-line evaluation (what the previous report plotted)
Tb_line = np.interp(xe, cv.x1, cv.T_nrtl_lit)
ex["hE_nrtl_lit_bubbleline"] = NRTL_LIT.hE(xe, Tb_line)
exd = pd.DataFrame(ex)
exd.to_csv(RES / "excess_functions.csv", index=False)


def ext(series):
    i, j = int(np.argmin(series)), int(np.argmax(series))
    return {"min": float(series[i]), "x_min": float(xe[i]), "max": float(series[j]), "x_max": float(xe[j])}


hE298 = exd["hE_nrtl_lit_298.15"].values
hE298f = exd["hE_nrtl_fit_298.15"].values
zero = lambda s: [float(xe[i]) for i in np.where(np.diff(np.sign(s[1:-1])) != 0)[0] + 1]
OUT["excess"] = {
    "Tref_K": Tref,
    "gE_lit": ext(exd.gE_nrtl_lit.values),
    "hE_lit": ext(exd.hE_nrtl_lit.values), "hE_lit_zero": zero(exd.hE_nrtl_lit.values),
    "TsE_lit": ext(exd.TsE_nrtl_lit.values),
    "cpE_lit": ext(exd.cpE_nrtl_lit.values),
    "gE_fit": ext(exd.gE_nrtl_fit.values), "hE_fit": ext(exd.hE_nrtl_fit.values),
    "hE298_lit": ext(hE298), "hE298_fit": ext(hE298f),
    "hE298_exp_min": -775.0, "hE298_exp_xmin": 0.17,
    "hE298_lit_at_0p17": float(np.interp(0.17, xe, hE298)),
    "hE_bubbleline": ext(exd.hE_nrtl_lit_bubbleline.values),
    "hE_bubbleline_zero": zero(exd.hE_nrtl_lit_bubbleline.values),
    "max_abs_hE_over_dHvap_pct": float(100 * np.max(np.abs(exd.hE_nrtl_lit.values)) / DHVAP_NBP[1]),
    "dHvap_ratio_etoh_water": DHVAP_NBP[0] / DHVAP_NBP[1],
}
# adiabatic temperature change on isothermal-equivalent mixing at 25 C, x1 = 0.17
CP_L = (112.3, 75.3)  # J/mol/K, liquid ethanol and water at 298 K (NIST WebBook)
cp_mix = 0.17 * CP_L[0] + 0.83 * CP_L[1]
OUT["excess"]["cp_mix_x017"] = cp_mix
OUT["excess"]["dT_ad_exp_K"] = 775.0 / cp_mix
OUT["excess"]["dT_ad_nrtl_K"] = -float(np.interp(0.17, xe, hE298)) / cp_mix
for Tk in (298.15, 323.15, 348.15, 373.15):
    OUT["excess"][f"hE_lit_{Tk:.2f}"] = ext(exd[f"hE_nrtl_lit_{Tk:.2f}"].values)
    OUT["excess"][f"hE_lit_{Tk:.2f}_zero"] = zero(exd[f"hE_nrtl_lit_{Tk:.2f}"].values)

# ---------------------------------------------------------------------------
# 7. Azeotrope versus pressure and pressure-swing feasibility
# ---------------------------------------------------------------------------
Ps = np.array([5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 40, 50, 70, 101.325, 150, 200, 300])
arows = []
for P in Ps:
    a = azeotrope(NRTL_LIT, P, lo=0.5, hi=0.9999999)
    arows.append({"P_kPa": float(P), "x_az": None if a is None else a[0],
                  "T_az_K": None if a is None else a[1],
                  "wt_az": None if a is None else mole_to_mass(a[0])})
adf = pd.DataFrame(arows)
adf.to_csv(RES / "azeotrope_vs_pressure.csv", index=False)
OUT["azeotrope_vs_P"] = arows


def xaz(P):
    a = azeotrope(NRTL_LIT, P, lo=0.5, hi=0.9999999)
    return None if a is None else a[0]


psd = []
for PL, PH in [(20.0, 101.325), (10.0, 101.325), (101.325, 300.0), (20.0, 300.0)]:
    xl, xh = xaz(PL), xaz(PH)
    rr = (1 - xl) / (xl - xh) if (xl and xh and xl > xh) else None
    psd.append({"P_L_kPa": PL, "P_H_kPa": PH, "x_az_L": xl, "x_az_H": xh, "dx": None if rr is None else xl - xh,
                "recycle_per_product": rr})
OUT["psd"] = psd

# ---------------------------------------------------------------------------
# 8. Distillation case study (CMO McCabe-Thiele)
# ---------------------------------------------------------------------------
case = {"F_kmol_h": 100.0, "xF": 0.10, "q": 1.0, "xD": 0.85, "xB": 0.002, "P_kPa": P_ATM, "R_over_Rmin": 1.3}
eqs = {"ideal": EqCurve(Ideal(), P_ATM), "nrtl_lit": EqCurve(NRTL_LIT, P_ATM), "nrtl_fit": EqCurve(NR_FIT, P_ATM)}
Dfrac = (case["xF"] - case["xB"]) / (case["xD"] - case["xB"])
dhD = case["xD"] * DHVAP_NBP[0] + (1 - case["xD"]) * DHVAP_NBP[1]
design = {}
stage_dump = {}
for key, eq in eqs.items():
    Rm, xp = rmin_cmo(eq, case["xF"], case["xD"])
    Nm = nmin_total_reflux(eq, case["xD"], case["xB"])
    Rr = case["R_over_Rmin"] * Rm
    st, nf, ok = step_stages(eq, case["xF"], case["xD"], case["xB"], Rr)
    D = Dfrac * case["F_kmol_h"]
    V = (Rr + 1) * D
    QB = V * DHVAP_NBP[1] / 3600.0  # kW (saturated-liquid feed: Vbar = V)
    QC = V * dhD / 3600.0
    Tst = [float(eq.T_of_x(x)) for x, _ in st]
    design[key] = {"Rmin": Rm, "pinch_x": xp, "Nmin": int(Nm), "R": Rr, "N_stages_incl_reboiler": len(st),
                   "feed_stage": nf, "D_kmol_h": D, "V_kmol_h": V, "QB_kW": QB, "QC_kW": QC,
                   "QB_MJ_per_kmol_EtOH": QB * 3600 / (D * case["xD"]) / 1000.0,
                   "QB_MJ_per_kg_EtOH": QB * 3600 / (D * case["xD"] * 46.068) / 1000.0, "T_top_K": Tst[0], "T_bottom_K": Tst[-1]}
    stage_dump[key] = {"x": [s[0] for s in st], "y": [s[1] for s in st], "T": Tst}
OUT["column_case"] = case
OUT["column_design"] = design
json.dump(stage_dump, open(RES / "stages.json", "w"), indent=1)

# rating: ideal-designed hardware operated with the real (NRTL) phase behaviour
di = design["ideal"]
rate = rate_column(eqs["nrtl_lit"], di["N_stages_incl_reboiler"], di["feed_stage"], di["R"], Dfrac, case["xF"])
R_needed = brentq(lambda Rr: rate_column(eqs["nrtl_lit"], di["N_stages_incl_reboiler"], di["feed_stage"], Rr,
                                         Dfrac, case["xF"])[0] - case["xD"], di["R"], 200.0)
# tray temperatures of the ideal design as actually operated (true VLE)
def profile(eq_true, N, NF, Rr, xD, xF):
    B = 1 - Dfrac
    xB = (xF - Dfrac * xD) / B
    L, V = Rr * Dfrac, (Rr + 1) * Dfrac
    Lb, Vb = L + 1.0, V
    y = xD
    xs = []
    for n in range(1, N + 1):
        x = float(eq_true.x_of_y(y))
        xs.append(x)
        y = (L / V) * x + Dfrac * xD / V if n < NF else (Lb / Vb) * x - B * xB / Vb
    return xs, [float(eq_true.T_of_x(x)) for x in xs]


xs_act, T_act = profile(eqs["nrtl_lit"], di["N_stages_incl_reboiler"], di["feed_stage"], di["R"], rate[0], case["xF"])
stage_dump["ideal_actual"] = {"x": xs_act, "T": T_act}
json.dump(stage_dump, open(RES / "stages.json", "w"), indent=1)
OUT["column_rating_ideal_design"] = {
    "xD_achieved": rate[0], "xB_achieved": rate[1],
    "ethanol_loss_ratio": rate[1] / case["xB"],
    "R_needed_for_spec": R_needed,
    "QB_needed_kW": (R_needed + 1) * Dfrac * case["F_kmol_h"] * DHVAP_NBP[1] / 3600,
    "tray_T_shift_max_K": float(np.max(np.abs(np.array(stage_dump["ideal"]["T"]) - np.array(T_act)))),
}
# what the ideal model claims about high purity
Rm99, _ = rmin_cmo(eqs["ideal"], case["xF"], 0.995)
OUT["ideal_claims"] = {"xD_0p995_Rmin_ideal": Rm99, "nrtl_max_xD_limit": az_lit[0]}

# ---------------------------------------------------------------------------
# 9. Sensitivity of the NRTL design to the parameter set and to R/Rmin
# ---------------------------------------------------------------------------
sens = []
for frac in (1.1, 1.2, 1.3, 1.5, 2.0):
    eq = eqs["nrtl_lit"]
    Rm, _ = rmin_cmo(eq, case["xF"], case["xD"])
    st, nf, ok = step_stages(eq, case["xF"], case["xD"], case["xB"], frac * Rm)
    sens.append({"R_over_Rmin": frac, "R": frac * Rm, "N": len(st), "NF": nf,
                 "QB_kW": (frac * Rm + 1) * Dfrac * case["F_kmol_h"] * DHVAP_NBP[1] / 3600})
OUT["column_sensitivity"] = sens

# ---------------------------------------------------------------------------
# 10. Extractive agent: experimental ethanol/water volatility with ethylene glycol
#     alpha_12 = (y1/x1)/(y2/x2) from the ternary data of Kamihama et al. (2012),
#     compared with the binary (salt-free) NRTL value at the same solvent-free x1'.
# ---------------------------------------------------------------------------
tern = pd.read_csv(ROOT / "data" / "kamihama2012_ethanol_water_EG_ternary_101p3kPa.csv", comment="#")
tern["x_EG"] = 1.0 - tern.x_ethanol - tern.x_water
tern["xp"] = tern.x_ethanol / (tern.x_ethanol + tern.x_water)
tern["alpha"] = (tern.y_ethanol / tern.x_ethanol) / (tern.y_water / tern.x_water)
ab = []
for xp in tern.xp:
    Tb, yb = bubble_T(float(xp), P_EXP, NRTL_LIT)
    ab.append((yb / xp) / ((1 - yb) / (1 - xp)))
tern["alpha_binary_nrtl"] = ab
tern["enhancement"] = tern.alpha / tern.alpha_binary_nrtl
# points with x_water or y_water < 0.01 are dominated by analytical resolution (+/-0.001)
tern["reliable"] = (tern.x_water >= 0.01) & (tern.y_water >= 0.01)
tern.to_csv(RES / "eg_ternary_alpha.csv", index=False)
rel = tern[tern.reliable]
hi = rel[rel.xp >= 0.85]
OUT["extractive_eg"] = {
    "n_points": int(len(tern)),
    "n_reliable": int(len(rel)),
    "near_az_xp_range": [float(hi.xp.min()), float(hi.xp.max())],
    "near_az_xEG_range": [float(hi.x_EG.min()), float(hi.x_EG.max())],
    "near_az_alpha_range": [float(hi.alpha.min()), float(hi.alpha.max())],
    "near_az_alpha_binary_range": [float(hi.alpha_binary_nrtl.min()), float(hi.alpha_binary_nrtl.max())],
    "near_az_n": int(len(hi)),
    "low_x_enhancement_range": [float(rel[rel.xp < 0.25].enhancement.min()),
                                float(rel[rel.xp < 0.25].enhancement.max())],
}


# ---------------------------------------------------------------------------
# 11. Operating-variable sensitivity (T, P, composition) for Part V
# ---------------------------------------------------------------------------
_Tb = lambda x, P=P_ATM: bubble_T(x, P, NRTL_LIT)[0]
sensv = {
    "T_x070_C": _Tb(0.70) - 273.15, "T_x080_C": _Tb(0.80) - 273.15,
    "T_x085_C": _Tb(0.85) - 273.15, "T_x097_C": _Tb(0.97) - 273.15,
    "T_x002_C": _Tb(0.02) - 273.15, "T_x010_C": _Tb(0.10) - 273.15,
}
sensv["dT_top_070_to_az_K"] = sensv["T_x070_C"] - (OUT["one_atm"]["azeotrope_nrtl_lit"]["T_K"] - 273.15)
sensv["dTdx_dilute_K_per_molfrac"] = (sensv["T_x002_C"] - sensv["T_x010_C"]) / 0.08
sensv["dTdx_top_K_per_molfrac"] = (sensv["T_x070_C"] - sensv["T_x085_C"]) / 0.15
dP = 1.0
sensv["dTdP_x085_K_per_kPa"] = (_Tb(0.85, P_ATM + dP) - _Tb(0.85, P_ATM - dP)) / (2 * dP)
sensv["dTdP_x002_K_per_kPa"] = (_Tb(0.02, P_ATM + dP) - _Tb(0.02, P_ATM - dP)) / (2 * dP)
# composition change equivalent to a 1 kPa pressure disturbance in the rectifying section
sensv["dx_equiv_1kPa_top"] = sensv["dTdP_x085_K_per_kPa"] / sensv["dTdx_top_K_per_molfrac"]
_az = {round(a["P_kPa"], 3): a["x_az"] for a in OUT["azeotrope_vs_P"]}
sensv["dxaz_dlnP"] = (_az[200.0] - _az[50.0]) / np.log(200.0 / 50.0)
OUT["sensitivity_operating"] = sensv

def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    return o


json.dump(clean(OUT), open(RES / "results.json", "w"), indent=2, ensure_ascii=False)
print(json.dumps(clean({k: OUT[k] for k in ["regression", "validation", "herington", "one_atm", "virial"]}), indent=1)[:6000])
