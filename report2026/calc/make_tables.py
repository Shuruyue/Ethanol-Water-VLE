"""LaTeX table bodies generated from results (appendix data table, vapour correction, azeotrope vs P)."""
import json
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parent.parent
K = 273.15
d = pd.read_csv(ROOT/"data"/"kamihama2012_ethanol_water_101p3kPa.csv", comment="#")
g = pd.read_csv(ROOT/"results"/"gamma_experimental.csv")
rs = pd.read_csv(ROOT/"results"/"residuals.csv")
rows = []
for i in range(len(d)):
    rows.append(f"{d.x1[i]:.3f} & {d.y1[i]:.3f} & {d.T_K[i]-K:.2f} & {d.U_T_K[i]:.2f} & {d.U_y1[i]:.3f} & "
                f"{g.gamma1_exp[i]:.3f} & {g.gamma2_exp[i]:.3f} & {rs.dT_nrtl_lit[i]:+.2f} & {rs.dy_nrtl_lit[i]:+.4f} \\\\")
(ROOT/"tex"/"tab_data.tex").write_text("\n".join(rows)+"\n", encoding="utf-8")
r = json.load(open(ROOT/"results"/"results.json"))
rows = []
for v in r["virial"]["rows"]:
    rows.append(f"{v['x1']:.2f} & {v['T_ig_K']-K:.2f} & {v['dT_K']:+.3f} & {v['y_ig']:.4f} & {v['dy']:+.4f} & "
                f"{v['Phi1']:.4f} & {v['Phi2']:.4f} \\\\")
(ROOT/"tex"/"tab_virial.tex").write_text("\n".join(rows)+"\n", encoding="utf-8")
rows = []
for a in r["azeotrope_vs_P"]:
    if a["P_kPa"] in (5, 10, 20, 50, 101.325, 200, 300):
        rows.append(f"{('101.3' if a['P_kPa']==101.325 else format(a['P_kPa'],'.0f'))} & {a['T_az_K']-K:.1f} & {a['x_az']:.3f} & {100*a['wt_az']:.1f} \\\\")
(ROOT/"tex"/"tab_azP.tex").write_text("\n".join(rows)+"\n", encoding="utf-8")
print("ok")
