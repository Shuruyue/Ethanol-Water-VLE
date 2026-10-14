"""Publication-style figures for the report (vector PDF)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

from vle_core import NRTL_LIT, EqCurve, Ideal, psat

ROOT = Path(__file__).resolve().parent.parent
RES, FIG = ROOT / "results", ROOT / "fig"
FIG.mkdir(exist_ok=True)

mpl.rcParams.update({
    "font.family": "serif", "font.serif": ["LM Roman 10", "DejaVu Serif"],
    "mathtext.fontset": "cm", "font.size": 9.5, "axes.labelsize": 10, "axes.titlesize": 10,
    "legend.fontsize": 8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.linewidth": 0.7, "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "lines.linewidth": 1.4, "legend.frameon": False, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "axes.grid": False, "pdf.fonttype": 42,
})
C = {"lit": "#1b4f8a", "fit": "#2a9d8f", "vl": "#d97706", "ideal": "#6b7280", "exp": "#111111",
     "vir": "#9d174d", "accent": "#b91c1c"}
W1, W2 = 3.4, 6.6  # inch

r = json.load(open(RES / "results.json"))
cv = pd.read_csv(RES / "txy_curves_1atm.csv")
d = pd.read_csv(ROOT / "data" / "kamihama2012_ethanol_water_101p3kPa.csv", comment="#")
ge = pd.read_csv(RES / "gamma_experimental.csv")
rs = pd.read_csv(RES / "residuals.csv")
ex = pd.read_csv(RES / "excess_functions.csv")
az = pd.DataFrame(r["azeotrope_vs_P"])
st = json.load(open(RES / "stages.json"))
K = 273.15
xa, Ta = r["one_atm"]["azeotrope_nrtl_lit"]["x1"], r["one_atm"]["azeotrope_nrtl_lit"]["T_K"] - K


def save(fig, name):
    fig.savefig(FIG / f"{name}.pdf")
    fig.savefig(FIG / f"{name}.png", dpi=220)
    plt.close(fig)


# ---------------------------------------------------------------- Fig. T-x-y
fig = plt.figure(figsize=(W2, 3.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.75, 1.0], wspace=0.28)
ax, ia = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
ax.plot(cv.x1, cv.T_ideal - K, "--", color=C["ideal"], lw=1.1, label="Ideal solution: bubble")
ax.plot(cv.y_ideal, cv.T_ideal - K, ":", color=C["ideal"], lw=1.1, label="Ideal solution: dew")
ax.plot(cv.x1, cv.T_nrtl_lit - K, "-", color=C["lit"], label="NRTL: bubble")
ax.plot(cv.y_nrtl_lit, cv.T_nrtl_lit - K, "-", color="#7aa7d9", lw=1.1, label="NRTL: dew")
ax.errorbar(d.x1, d.T_K - K, yerr=d.U_T_K, fmt="o", ms=3.4, mfc="white", mec=C["exp"], ecolor="0.55",
            elinewidth=0.6, capsize=1.5, label="Exp. $x_1$")
ax.plot(d.y1, d.T_K - K, "^", ms=3.4, mfc=C["exp"], mec=C["exp"], label="Exp. $y_1$")
ax.plot([xa], [Ta], "D", ms=4.6, mfc=C["accent"], mec="white", mew=0.6, zorder=6, label="NRTL azeotrope")
k = r["one_atm"]["at_x1"]
Ti, Tn = r["one_atm"]["T_ideal_at"] - K, r["one_atm"]["T_nrtl_at"] - K
ax.annotate("", xy=(k, Tn), xytext=(k, Ti), arrowprops=dict(arrowstyle="<->", lw=0.8, color=C["accent"]))
ax.text(k + 0.015, 0.5 * (Ti + Tn), f"$\\Delta T$ = {Ti - Tn:.1f} K\n($x_1$ = {k:.2f})", color=C["accent"],
        fontsize=8, va="center")
ax.add_patch(plt.Rectangle((0.65, 77.9), 0.35, 0.9, fill=False, ec="0.45", lw=0.6, ls="-"))
ax.set_xlim(0, 1)
ax.set_ylim(76, 101)
ax.set_xlabel("$x_1$ or $y_1$")
ax.set_ylabel("Temperature (°C)")
ax.text(0.02, 0.04, "(a)", transform=ax.transAxes)
m = cv.x1 > 0.6
ia.plot(cv.x1[m], cv.T_nrtl_lit[m] - K, "-", color=C["lit"])
ia.plot(cv.y_nrtl_lit[m], cv.T_nrtl_lit[m] - K, "-", color="#7aa7d9", lw=1.1)
mm = d.x1 > 0.6
ia.errorbar(d.x1[mm], d.T_K[mm] - K, yerr=d.U_T_K[mm], fmt="o", ms=3.2, mfc="white", mec=C["exp"],
            ecolor="0.55", elinewidth=0.6, capsize=1.2)
ia.plot(d.y1[mm], d.T_K[mm] - K, "^", ms=3.2, color=C["exp"])
ia.plot([xa], [Ta], "D", ms=4.6, mfc=C["accent"], mec="white", zorder=6)
xe_, Te_ = r["azeotrope_exp_interp"]["x1"], r["azeotrope_exp_interp"]["T_K"] - K
ia.annotate(f"NRTL: {xa:.4f}, {Ta:.2f} °C\nExp.: {xe_:.4f}, {Te_:.2f} °C", xy=(xa, Ta), xytext=(0.985, 78.70), ha="right",
            fontsize=7.2, arrowprops=dict(arrowstyle="-", lw=0.5, color="0.4"))
ia.set_xlim(0.65, 1.0)
ia.set_ylim(77.9, 78.8)
ia.set_xlabel("$x_1$ or $y_1$")
ia.text(0.04, 0.04, "(b) zoom near azeotrope", transform=ia.transAxes, fontsize=8)
h, l = ax.get_legend_handles_labels()
fig.legend(h, l, loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.02), fontsize=7.6)
save(fig, "fig_txy")

# ---------------------------------------------------------------- Fig. y-x + alpha
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 3.0))
a1.plot([0, 1], [0, 1], "-", color="0.75", lw=0.7)
a1.plot(cv.x1, cv.y_ideal, "--", color=C["ideal"], lw=1.1, label="Ideal solution")
a1.plot(cv.x1, cv.y_van_laar_fit, "-.", color=C["vl"], lw=1.1, label="Van Laar (this work)")
a1.plot(cv.x1, cv.y_nrtl_lit, "-", color=C["lit"], label="NRTL (Aspen/DDB)")
a1.errorbar(d.x1, d.y1, yerr=d.U_y1, fmt="o", ms=3.4, mfc="white", mec=C["exp"], ecolor="0.6", elinewidth=0.6,
            capsize=1.2, label="Exp. (Kamihama et al.)")
a1.plot([xa], [xa], "D", ms=4.5, mfc=C["accent"], mec="white")
a1.set_xlim(0, 1)
a1.set_ylim(0, 1)
a1.set_aspect("equal")
a1.set_xlabel("$x_1$")
a1.set_ylabel("$y_1$")
a1.legend(loc="lower right", fontsize=7.4)
a1.text(0.02, 0.95, "(a)", transform=a1.transAxes, va="top")
inner = (cv.x1 > 0.002) & (cv.x1 < 0.998)
a2.plot(cv.x1[inner], cv.alpha_ideal[inner], "--", color=C["ideal"], lw=1.1, label="Ideal solution")
a2.plot(cv.x1[inner], cv.alpha_nrtl_lit[inner], "-", color=C["lit"], label="NRTL")
ae = (d.y1 / d.x1) / ((1 - d.y1) / (1 - d.x1))
a2.plot(d.x1, ae, "o", ms=3.4, mfc="white", mec=C["exp"], label="Exp.")
a2.axhline(1, color=C["accent"], lw=0.7, ls=":")
a2.text(0.04, 0.88, r"$\alpha_{12}=1$ (azeotrope)", fontsize=7.5, color=C["accent"])
a2.set_yticks([1, 2, 5, 10, 20])
a2.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
a2.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
a2.set_yscale("log")
a2.set_xlim(0, 1)
a2.set_ylim(0.8, 20)
a2.set_xlabel("$x_1$")
a2.set_ylabel(r"$\alpha_{12}=(y_1/x_1)/(y_2/x_2)$")
a2.legend(loc="upper right")
a2.text(0.02, 0.95, "(b)", transform=a2.transAxes, va="top")
fig.tight_layout(w_pad=1.8)
save(fig, "fig_yx_alpha")

# ---------------------------------------------------------------- Fig. residuals
fig, (a1, a2) = plt.subplots(2, 1, figsize=(W2, 3.9), sharex=True)
sty = {"ideal": ("s", C["ideal"], "Ideal solution"), "van_laar_fit": ("v", C["vl"], "Van Laar (this work)"),
       "nrtl_fit": ("^", C["fit"], "NRTL, $b_{ij}$ only (this work)"), "nrtl_lit": ("o", C["lit"], "NRTL (Aspen/DDB)")}
a1.fill_between(d.x1, -d.U_T_K, d.U_T_K, color="0.88", step="mid", label="Exp. uncertainty $U(T)$")
a2.fill_between(d.x1, -d.U_y1, d.U_y1, color="0.88", step="mid", label="Exp. uncertainty $U(y_1)$")
for kkey, (mk, col, lab) in sty.items():
    if kkey == "ideal":
        continue
    a1.plot(rs.x1, rs[f"dT_{kkey}"], mk, ms=3.6, mfc="white" if kkey != "nrtl_lit" else col, mec=col, ls="-",
            lw=0.6, color=col, label=lab)
    a2.plot(rs.x1, rs[f"dy_{kkey}"], mk, ms=3.6, mfc="white" if kkey != "nrtl_lit" else col, mec=col, ls="-",
            lw=0.6, color=col)
a1.axhline(0, color="k", lw=0.5)
a2.axhline(0, color="k", lw=0.5)
a1.set_ylabel(r"$T_{\rm calc}-T_{\rm exp}$ (K)")
a2.set_ylabel(r"$y_{1,\rm calc}-y_{1,\rm exp}$")
a2.set_xlabel("$x_1$")
a1.set_ylim(-2, 2)
a2.set_ylim(-0.06, 0.06)
a1.legend(ncol=2, loc="upper right", fontsize=7.4)
a1.text(0.01, 0.92, "(a)", transform=a1.transAxes, va="top")
a2.text(0.01, 0.92, "(b)", transform=a2.transAxes, va="top")
fig.tight_layout(h_pad=0.4)
save(fig, "fig_residuals")

# ---------------------------------------------------------------- Fig. ln gamma
fig, ax = plt.subplots(figsize=(W1 + 0.4, 2.9))
ax.plot(ge.x1, np.log(ge.gamma1_exp), "o", ms=3.6, mfc=C["exp"], mec=C["exp"], label=r"$\ln\gamma_1$ exp.")
ax.plot(ge.x1, np.log(ge.gamma2_exp), "o", ms=3.6, mfc="white", mec=C["exp"], label=r"$\ln\gamma_2$ exp.")
ax.plot(cv.x1, cv.lng1_nrtl_lit, "-", color=C["lit"], label="NRTL (along bubble line)")
ax.plot(cv.x1, cv.lng2_nrtl_lit, "-", color=C["lit"])
ax.plot(cv.x1, cv.lng1_van_laar_fit, "-.", color=C["vl"], lw=1.1, label="Van Laar (this work)")
ax.plot(cv.x1, cv.lng2_van_laar_fit, "-.", color=C["vl"], lw=1.1)
ax.set_xlim(0, 1)
ax.set_ylim(-0.1, 1.9)
ax.set_xlabel("$x_1$")
ax.set_ylabel(r"$\ln\gamma_i$")
ax.legend(loc="upper center")
save(fig, "fig_gamma")

# ---------------------------------------------------------------- Fig. excess functions
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 3.0))
Tr = r["excess"]["Tref_K"]
a1.plot(ex.x1, ex.gE_nrtl_lit, "-", color=C["lit"], label="$G^E$")
a1.plot(ex.x1, ex.hE_nrtl_lit, "--", color=C["accent"], label="$H^E$")
a1.plot(ex.x1, ex.TsE_nrtl_lit, ":", color="k", label="$TS^E=H^E-G^E$")
a1.plot(ex.x1, ex.gE_nrtl_fit, "-", color=C["fit"], lw=0.9, alpha=0.9, label="$G^E$, NRTL $b_{ij}$-only")
a1.plot(ex.x1, ex.hE_nrtl_fit, "--", color=C["fit"], lw=0.9, alpha=0.9, label="$H^E$, NRTL $b_{ij}$-only")
a1.axhline(0, color="0.5", lw=0.5)
a1.set_xlim(0, 1)
a1.set_xlabel("$x_1$")
a1.set_ylabel("Excess property (J mol$^{-1}$)")
a1.set_ylim(-1150, 1950)
a1.legend(loc="upper center", ncol=2, fontsize=6.8, columnspacing=0.8, handlelength=1.8,
           title=f"isothermal, $T$ = {Tr - K:.2f} °C", title_fontsize=7.2)
a1.text(0.96, 0.04, "(a)", transform=a1.transAxes, ha="right")
cols = ["#93c5fd", "#3b82f6", "#1d4ed8", "#1e3a8a"]
for c, Tk in zip(cols, ["298.15", "323.15", "348.15", "373.15"]):
    a2.plot(ex.x1, ex[f"hE_nrtl_lit_{Tk}"], "-", color=c, label=f"NRTL, {float(Tk) - K:.0f} °C")
a2.plot(ex.x1, ex.hE_nrtl_lit_bubbleline, ":", color="k", lw=1.1, label="NRTL along bubble line")
a2.plot([0.17], [-775], "*", ms=9, color=C["accent"], mec="k", mew=0.4,
        label="Exp. minimum, 25 °C\n(Larkin 1975; Ott et al. 1986)")
a2.axhline(0, color="0.5", lw=0.5)
a2.set_xlim(0, 1)
a2.set_ylim(-850, 200)
a2.set_xlabel("$x_1$")
a2.set_ylabel("$H^E$ (J mol$^{-1}$)")
a2.legend(loc="lower right", fontsize=6.6, frameon=True, framealpha=1, edgecolor="0.8")
a2.text(0.02, 0.96, "(b)", transform=a2.transAxes, va="top")
fig.tight_layout(w_pad=1.6)
save(fig, "fig_excess")

# ---------------------------------------------------------------- Fig. azeotrope vs P
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W2, 2.7))
a1.semilogx(az.P_kPa, az.x_az, "o-", ms=3.2, color=C["lit"], lw=1.0, label="NRTL (extrapolated)")
a1.errorbar([101.3], [r["azeotrope_exp_interp"]["x1"]], yerr=[[0.01], [0.01]], fmt="s", ms=4, color=C["exp"],
            label="Exp. 101.3 kPa (interpolated)")
a1.set_xlim(4, 400)
a1.set_xlabel("$P$ (kPa)")
a1.set_ylabel(r"$x_{1}^{\rm az}$")
a1.legend(loc="upper right", fontsize=7)
a1.text(0.04, 0.06, "(a)", transform=a1.transAxes)
a2.semilogx(az.P_kPa, az.T_az_K - K, "o-", ms=3.2, color=C["lit"], lw=1.0)
a2.set_xlim(4, 400)
a2.set_xlabel("$P$ (kPa)")
a2.set_ylabel(r"$T^{\rm az}$ (°C)")
a2.text(0.02, 0.92, "(b)", transform=a2.transAxes)
fig.tight_layout(w_pad=1.6)
save(fig, "fig_azeotrope_P")

# ---------------------------------------------------------------- Fig. McCabe-Thiele
case = r["column_case"]
P = case["P_kPa"]
eqI, eqN = EqCurve(Ideal(), P, n=1500), EqCurve(NRTL_LIT, P, n=1500)
fig, axs = plt.subplots(1, 2, figsize=(W2, 3.35))
for ax, key, eq, title in [(axs[0], "ideal", eqI, "(a) Design with ideal-solution VLE"),
                           (axs[1], "nrtl_lit", eqN, "(b) Design with NRTL VLE")]:
    dz = r["column_design"][key]
    xF, xD, xB, Rr = case["xF"], case["xD"], case["xB"], dz["R"]
    ax.plot([0, 1], [0, 1], color="0.75", lw=0.6)
    ax.plot(eq.x, eq.y, "-", color=C["lit"] if key != "ideal" else C["ideal"], lw=1.2)
    xi = (xD / (Rr + 1)) / (1 - Rr / (Rr + 1))  # rectifying line meets q-line x = xF
    yF = Rr / (Rr + 1) * xF + xD / (Rr + 1)
    ax.plot([xD, xF], [xD, yF], color=C["accent"], lw=0.9)
    ax.plot([xF, xB], [yF, xB], color=C["accent"], lw=0.9)
    ax.plot([xF, xF], [xF, yF], color="0.4", lw=0.7, ls="--")
    # staircase
    xs, ys = st[key]["x"], st[key]["y"]
    px, py = [xD], [xD]
    for xx, yy in zip(xs, ys):
        px += [xx, xx]
        py += [yy, None]
    # rebuild proper stair: horizontal to (x_n, y_n), vertical to y_{n+1}
    sx, sy = [xD], [xD]
    for n, (xx, yy) in enumerate(zip(xs, ys)):
        sx.append(xx)
        sy.append(yy)
        if n + 1 < len(ys):
            sx.append(xx)
            sy.append(ys[n + 1])
    ax.plot(sx, sy, color="k", lw=0.55)
    # tangent-pinch marker
    if key != "ideal":
        xp = dz["pinch_x"]
        Rm = dz["Rmin"]
        ax.plot([xD, 0.55], [xD, xD + Rm / (Rm + 1) * (0.55 - xD)], ":", color=C["vir"], lw=0.9)
        ax.plot([xp], [eq.y_of_x(xp)], "o", ms=4, mfc="none", mec=C["vir"])
        ax.text(0.40, 0.87, "tangent pinch\n($R_{\\min}$ line)", fontsize=7, color=C["vir"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("$x_1$")
    ax.set_ylabel("$y_1$")
    ax.set_title(title, fontsize=8.8)
    ax.text(0.40, 0.08, f"$R_{{\\min}}$ = {dz['Rmin']:.2f}, $R$ = {dz['R']:.2f}\n"
            f"$N$ = {dz['N_stages_incl_reboiler']} (feed {dz['feed_stage']})", fontsize=7.6,
            bbox=dict(fc="white", ec="0.7", lw=0.5))
fig.tight_layout(w_pad=1.4)
save(fig, "fig_mccabe")

# ---------------------------------------------------------------- Fig. column temperature profiles
fig, ax = plt.subplots(figsize=(W1 + 0.5, 3.0))
ax.plot(np.array(st["ideal"]["T"]) - K, np.arange(1, len(st["ideal"]["T"]) + 1), "s--", ms=3, color=C["ideal"],
        lw=0.9, label="Ideal design: predicted")
ax.plot(np.array(st["ideal_actual"]["T"]) - K, np.arange(1, len(st["ideal_actual"]["T"]) + 1), "o-", ms=3,
        color=C["accent"], lw=0.9, label="Ideal design: actual (NRTL)")
ax.plot(np.array(st["nrtl_lit"]["T"]) - K, np.arange(1, len(st["nrtl_lit"]["T"]) + 1), "^-", ms=3,
        color=C["lit"], lw=0.9, label="NRTL design")
ax.invert_yaxis()
ax.set_xlabel("Stage temperature (°C)")
ax.set_ylabel("Stage number (top = 1)")
ax.legend(loc="center", bbox_to_anchor=(0.6, 0.45), fontsize=7.2, frameon=True, framealpha=1, edgecolor="0.85")
save(fig, "fig_temperature_profile")
print("figures written to", FIG)


# ---------------------------------------------------------------- Fig. EG extractive agent
eg = pd.read_csv(RES / "eg_ternary_alpha.csv")
egr = eg[eg.reliable]
xs_ = np.linspace(0.01, 0.995, 300)
ab_ = []
from vle_core import bubble_T as _bt
for xp_ in xs_:
    _T, _y = _bt(float(xp_), 101.3, NRTL_LIT)
    ab_.append((_y / xp_) / ((1 - _y) / (1 - xp_)))
fig, ax = plt.subplots(figsize=(W1 * 1.25, 2.9))
ax.semilogy(xs_, ab_, "-", color=C["lit"], label="Binary, NRTL (no solvent)")
sc = ax.scatter(egr.xp, egr.alpha, c=egr.x_EG, cmap="viridis_r", vmin=0.05, vmax=0.6, s=18,
                edgecolors="k", linewidths=0.4, zorder=3, label="With EG, exp. (Kamihama et al.)")
ax.axhline(1.0, color=C["accent"], lw=0.6, ls=":")
ax.text(0.02, 1.03, r"$\alpha_{12}=1$", color=C["accent"], fontsize=7)
ax.set_xlim(0, 1)
ax.set_ylim(0.8, 15)
ax.set_xlabel(r"solvent-free $x_1' = x_1/(x_1+x_2)$")
ax.set_ylabel(r"$\alpha_{12}=(y_1/x_1)/(y_2/x_2)$")
cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label(r"$x_{\rm EG}$", fontsize=8)
cb.ax.tick_params(labelsize=7)
ax.legend(loc="upper right", fontsize=6.8)
fig.tight_layout()
save(fig, "fig_eg_alpha")
