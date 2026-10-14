"""Write every number quoted in the report as a LaTeX macro (tex/numbers.tex),
so that the text is generated from, and traceable to, results.json."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
r = json.load(open(ROOT / "results" / "results.json"))
K = 273.15
M = {}


def f(v, n):
    return f"{v:.{n}f}"


o, v, ex, reg = r["one_atm"], r["validation"], r["excess"], r["regression"]
cd, cr = r["column_design"], r["column_rating_ideal_design"]
M.update({
    "TbEtOH": f(r["pure"]["Tb_ethanol_K"] - K, 2), "TbWater": f(r["pure"]["Tb_water_K"] - K, 2),
    "dTmax": f(o["max_dT_ideal_minus_nrtl_K"], 1), "dTmaxX": f(o["at_x1"], 2),
    "dTmid": f(o["dT_at_x_0p5"], 1),
    "xAz": f(o["azeotrope_nrtl_lit"]["x1"], 4), "TAz": f(o["azeotrope_nrtl_lit"]["T_K"] - K, 2),
    "wAz": f(100 * o["azeotrope_nrtl_lit"]["wt"], 1),
    "xAzVir": f(o["azeotrope_nrtl_lit_virial"]["x1"], 4),
    "xAzFit": f(o["azeotrope_nrtl_fit"]["x1"], 4), "xAzVL": f(o["azeotrope_van_laar_fit"]["x1"], 4),
    "xAzExp": f(r["azeotrope_exp_interp"]["x1"], 3), "TAzExp": f(r["azeotrope_exp_interp"]["T_K"] - K, 2),
    "gOneInf": f(o["gamma1_inf_nrtl"], 2), "gTwoInf": f(o["gamma2_inf_nrtl"], 2),
    "alphaMax": f(o["alpha_nrtl_max"], 1), "alphaIdLo": f(o["alpha_ideal_range"][0], 2),
    "alphaIdHi": f(o["alpha_ideal_range"][1], 2),
    # regression
    "vlA": f(reg["van_laar"]["A12"], 3), "vlAse": f(reg["van_laar"]["se"][0], 3),
    "vlB": f(reg["van_laar"]["A21"], 3), "vlBse": f(reg["van_laar"]["se"][1], 3),
    "vlChi": f(reg["van_laar"]["chi2_red"], 2), "vlCorr": f(reg["van_laar"]["corr"], 2),
    "nrA": f(reg["nrtl_b_only"]["b12_K"], 1), "nrAse": f(reg["nrtl_b_only"]["se"][0], 1),
    "nrB": f(reg["nrtl_b_only"]["b21_K"], 1), "nrBse": f(reg["nrtl_b_only"]["se"][1], 1),
    "nrChi": f(reg["nrtl_b_only"]["chi2_red"], 2), "nrCorr": f(reg["nrtl_b_only"]["corr"], 3),
    "nPts": str(reg["n_points"]),
    # Herington
    "herD": f(r["herington"]["D"], 1), "herJ": f(r["herington"]["J"], 1),
    # virial
    "Bee": f(r["virial"]["B11_B22_B12_cm3mol_at_351K"][0], 0), "Bww": f(r["virial"]["B11_B22_B12_cm3mol_at_351K"][1], 0),
    "Bew": f(r["virial"]["B11_B22_B12_cm3mol_at_351K"][2], 0),
    "virdT": f(r["virial"]["max_abs_dT_K"], 2), "virdy": f(r["virial"]["max_abs_dy"], 4),
    # excess
    "gEmax": f(ex["gE_lit"]["max"], 0), "gEmaxX": f(ex["gE_lit"]["x_max"], 2),
    "hEmin": f(ex["hE_lit"]["min"], 0), "hEminX": f(ex["hE_lit"]["x_min"], 2),
    "hEmax": f(ex["hE_lit"]["max"], 0), "hEmaxX": f(ex["hE_lit"]["x_max"], 2),
    "hEzero": f(ex["hE_lit_zero"][0], 2),
    "TsEmin": f(ex["TsE_lit"]["min"], 0), "TsEminX": f(ex["TsE_lit"]["x_min"], 2),
    "cpEmax": f(ex["cpE_lit"]["max"], 1),
    "hEfitMax": f(ex["hE_fit"]["max"], 0), "hEfitMaxTwoNineEight": f(ex["hE298_fit"]["max"], 0),
    "hEtwoNineEightMin": f(ex["hE298_lit"]["min"], 0), "hEtwoNineEightMinX": f(ex["hE298_lit"]["x_min"], 2),
    "hEtwoNineEightAtExp": f(ex["hE298_lit_at_0p17"], 0),
    "hEratio": f(ex["max_abs_hE_over_dHvap_pct"], 2), "dHratio": f(ex["dHvap_ratio_etoh_water"], 3),
    "hEzeroTwoNineEight": f(ex["hE_lit_298.15_zero"][0], 2), "hEzeroThreeSevenThree": f(ex["hE_lit_373.15_zero"][0], 2),
    "dTadExp": f(ex["dT_ad_exp_K"], 1), "dTadNRTL": f(ex["dT_ad_nrtl_K"], 1), "cpMix": f(ex["cp_mix_x017"], 1),
    "TrefC": f(ex["Tref_K"] - K, 2),
})
for key, tag in [("ideal", "AAD"), ("van_laar_fit", "VL"), ("nrtl_lit", "NL"), ("nrtl_fit", "NF"),
                 ("nrtl_lit_virial", "NV")]:
    s = v[key]
    M[f"aadT{tag}"] = f(s["AAD_T_K"], 2 if key == "ideal" else 3)
    M[f"maxT{tag}"] = f(s["max_abs_T_K"], 2)
    M[f"aady{tag}"] = f(s["AAD_y"], 4)
    M[f"maxy{tag}"] = f(s["max_abs_y"], 3)
    M[f"biasT{tag}"] = f(s["bias_T_K"], 2 if key == "ideal" else 3)
# column
c = r["column_case"]
M.update({"xF": f(c["xF"], 2), "xD": f(c["xD"], 2), "xB": f(c["xB"], 3), "Fflow": f(c["F_kmol_h"], 0),
          "RR": f(c["R_over_Rmin"], 1)})
for key, tag in [("ideal", "Id"), ("nrtl_lit", "Nl"), ("nrtl_fit", "Nf")]:
    z = cd[key]
    M[f"Rmin{tag}"] = f(z["Rmin"], 2)
    M[f"R{tag}"] = f(z["R"], 2)
    M[f"Nmin{tag}"] = str(z["Nmin"])
    M[f"N{tag}"] = str(z["N_stages_incl_reboiler"])
    M[f"NF{tag}"] = str(z["feed_stage"])
    M[f"V{tag}"] = f(z["V_kmol_h"], 1)
    M[f"QB{tag}"] = f(z["QB_kW"], 0)
    M[f"QC{tag}"] = f(z["QC_kW"], 0)
    M[f"Esp{tag}"] = f(z["QB_MJ_per_kg_EtOH"], 2)
    M[f"pinch{tag}"] = f(z["pinch_x"], 3)
    M[f"Ttop{tag}"] = f(z["T_top_K"] - K, 1)
    M[f"Tbot{tag}"] = f(z["T_bottom_K"] - K, 1)
M["Dflow"] = f(cd["ideal"]["D_kmol_h"], 2)
M["QBratio"] = f(cd["ideal"]["QB_kW"] / cd["nrtl_lit"]["QB_kW"], 2)
M["xDach"] = f(cr["xD_achieved"], 3)
M["xBach"] = f(cr["xB_achieved"], 4)
M["lossRatio"] = f(cr["ethanol_loss_ratio"], 1)
M["Rneed"] = f(cr["R_needed_for_spec"], 0)
M["QBneed"] = f(cr["QB_needed_kW"], 0)
M["TshiftMax"] = f(cr["tray_T_shift_max_K"], 1)
M["RminIdNinetyNine"] = f(r["ideal_claims"]["xD_0p995_Rmin_ideal"], 2)
for s_ in r["column_sensitivity"]:
    t = {1.1: "a", 1.2: "b", 1.3: "c", 1.5: "d", 2.0: "e"}[s_["R_over_Rmin"]]
    M[f"sensR{t}"] = f(s_["R"], 2)
    M[f"sensN{t}"] = str(s_["N"])
    M[f"sensNF{t}"] = str(s_["NF"])
    M[f"sensQ{t}"] = f(s_["QB_kW"], 0)
# azeotrope vs P
azp = {round(a["P_kPa"], 3): a for a in r["azeotrope_vs_P"]}
for P, t in [(10.0, "Ten"), (20.0, "Twenty"), (50.0, "Fifty"), (101.325, "Atm"), (300.0, "ThreeHundred"),
             (5.0, "Five")]:
    M[f"xAzP{t}"] = f(azp[P]["x_az"], 3)
    M[f"TAzP{t}"] = f(azp[P]["T_az_K"] - K, 1)
for p in r["psd"]:
    t = {(20.0, 101.325): "a", (10.0, 101.325): "b", (101.325, 300.0): "c", (20.0, 300.0): "d"}[
        (p["P_L_kPa"], p["P_H_kPa"])]
    M[f"psdDx{t}"] = f(p["dx"], 3)
    M[f"psdRec{t}"] = f(p["recycle_per_product"], 1)

# vapour-pressure check
M["psatDevMax"] = f(r["pure"]["etoh_psat_dev_vs_eos_max_pct_293_373"], 1)
M["psatDevHi"] = f(r["pure"]["etoh_psat_dev_vs_eos_pct_393"], 1)
M["TbEtOHeos"] = f(r["pure"]["etoh_Tb_eos_K"] - K, 2)
# operating-variable sensitivity
so = r["sensitivity_operating"]
M["sTtopSeventy"] = f(so["T_x070_C"], 2)
M["sTtopEightyFive"] = f(so["T_x085_C"], 2)
M["sdTtop"] = f(so["dT_top_070_to_az_K"], 2)
M["sdTdxDil"] = f(so["dTdx_dilute_K_per_molfrac"], 0)
M["sdTdxTop"] = f(so["dTdx_top_K_per_molfrac"], 1)
M["sdTdP"] = f(so["dTdP_x085_K_per_kPa"], 2)
M["sdxEquiv"] = f(so["dx_equiv_1kPa_top"], 2)
M["sdxazdlnP"] = f(abs(so["dxaz_dlnP"]), 3)
M["sTdilTwo"] = f(so["T_x002_C"], 1)
M["sTdilTen"] = f(so["T_x010_C"], 1)
# extractive agent (ethylene glycol) from ternary data
eg = r["extractive_eg"]
M["egN"] = str(eg["n_points"])
M["egNrel"] = str(eg["n_reliable"])
M["egNhi"] = str(eg["near_az_n"])
M["egXpLo"], M["egXpHi"] = f(eg["near_az_xp_range"][0], 2), f(eg["near_az_xp_range"][1], 2)
M["egXegLo"], M["egXegHi"] = f(eg["near_az_xEG_range"][0], 2), f(eg["near_az_xEG_range"][1], 2)
M["egAlphaLo"], M["egAlphaHi"] = f(eg["near_az_alpha_range"][0], 2), f(eg["near_az_alpha_range"][1], 2)
M["egAlphaBinLo"], M["egAlphaBinHi"] = f(eg["near_az_alpha_binary_range"][0], 2), f(eg["near_az_alpha_binary_range"][1], 2)
M["egEnhLo"], M["egEnhHi"] = f(eg["low_x_enhancement_range"][0], 2), f(eg["low_x_enhancement_range"][1], 2)

out = ["% Auto-generated by calc/make_numbers.py from results/results.json -- do not edit by hand"]
for k_, val in M.items():
    assert k_.isalpha(), k_
    sval = val.replace("-", "\\ensuremath{-}", 1) if val.startswith("-") else val
    out.append(f"\\newcommand{{\\n{k_}}}{{{sval}}}")
(ROOT / "tex").mkdir(exist_ok=True)
(ROOT / "tex" / "numbers.tex").write_text("\n".join(out) + "\n", encoding="utf-8")
print(len(M), "macros")
