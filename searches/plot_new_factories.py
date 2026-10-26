"""Fig. 1: each new factory (orange) against the best prior construction for its target (gray), with the saving.
Run: python3 searches/plot_new_factories.py"""
import os, matplotlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.size": 8})
# (label, our N, best prior N, prior kind, level)   prior kind: "parts" = distill parts separately; "syn" = synthillation
# Read from outputs/error_coefficients.json (searches/error_coefficients.py): our n and output-error coefficient c
# (output error c p^2), and those of the best prior construction for the same target. A target is drawn only if
# our factory needs fewer inputs than that prior.
import json
LAB = {"CS12 T3": r"CS$_{12}$T$_3$", "CS12 CS34 CS56": r"CS$_{12}$CS$_{34}$CS$_{56}$", "CS12 CS34 T5": r"CS$_{12}$CS$_{34}$T$_5$",
       "CCZ123 CS45 T6": r"CCZ$_{123}$CS$_{45}$T$_6$", "CS12 T3 T4 T5": r"CS$_{12}$T$_3$T$_4$T$_5$", "CS12 T3...T7": r"CS$_{12}$T$_3\cdots$T$_7$",
       "CS12 CS34 CS56 T7": r"CS$_{12}$CS$_{34}$CS$_{56}$T$_7$", "CT12 sT3": r"CT$_{12}\sqrt{T}_3$", "CCS123 sT4": r"CCS$_{123}\sqrt{T}_4$",
       "CT12 sT3 sT4": r"CT$_{12}\sqrt{T}_3\sqrt{T}_4$", "CT12 CT34 sT5": r"CT$_{12}$CT$_{34}\sqrt{T}_5$",
       "CCS123 CT45 sT6": r"CCS$_{123}$CT$_{45}\sqrt{T}_6$", "CT12 CT34 sT5 sT6": r"CT$_{12}$CT$_{34}\sqrt{T}_5\sqrt{T}_6$"}
NQ = {"CS12 T3": 6, "CS12 CS34 CS56": 8, "CS12 CS34 T5": 8, "CCZ123 CS45 T6": 9, "CS12 T3 T4 T5": 8, "CS12 T3...T7": 10,
      "CS12 CS34 CS56 T7": 10, "CT12 sT3": 7, "CCS123 sT4": 8, "CT12 sT3 sT4": 8, "CT12 CT34 sT5": 9, "CCS123 CT45 sT6": 10, "CT12 CT34 sT5 sT6": 10}
PROVED = {"CS12 T3"}
E = [e for e in json.load(open(os.path.join(ROOT, "outputs", "error_coefficients.json"))) if e["n"] < e["prior_n"]]
# (label, our n, prior n, prior kind, level, proved, our c, prior c)
kind_of = lambda h: "syn" if "synth" in h and "+" not in h else ("mixed" if "synth" in h else "parts")
rows = [(LAB[e["target"]], e["n"], e["prior_n"], kind_of(e["prior"]),
         e["l"], e["target"] in PROVED, e["c"], e["prior_c"]) for e in E]
NS = [NQ[e["target"]] for e in E]
KS = [e["k"] for e in E]
NEW = "#eb6834"; PRIOR = "#8c8b87"
fig, ax = plt.subplots(figsize=(3.4, 3.9), constrained_layout=True)
ys = []
y = 0
for i, (lab, N, Np, kind, L, proved, c, cp) in enumerate(rows):
    if i == 0: y += 1.6
    if i > 0 and rows[i - 1][4] != L: y += 1.4
    ys.append(y)
    ax.plot([N, Np], [y, y], color="#cfcecb", lw=2.2, zorder=1, solid_capstyle="round")
    ax.scatter(Np, y, s=26, c="white", edgecolors=PRIOR, linewidths=1.2, zorder=2, marker={"syn": "s", "mixed": "^", "parts": "o"}[kind])
    ax.scatter(N, y, s=34, c=NEW, edgecolors="white", linewidths=0.6, zorder=3, marker="D")
    ax.text(Np + 1.2, y, f"$-{Np - N}$", va="center", fontsize=6.5, color="#52514e")
    ax.text(1.005, y, f"{c} / {cp}", va="center", ha="left", fontsize=6.0, color="#52514e", transform=ax.get_yaxis_transform())
    ax.text(N - (3.6 if proved else 2.0), y, lab, va="center", ha="right", fontsize=6.3, color="#0b0b0b")
    if proved: ax.scatter(N, y, s=95, facecolors="none", edgecolors="#0b0b0b", linewidths=0.8, zorder=4)
    y += 1
ax.set_yticks(ys); ax.set_yticklabels([f"[[{r[1]},{KS[j]},2]]  $N={NS[j]}$" for j, r in enumerate(rows)], fontsize=6.6)
ax.invert_yaxis()
ax.set_xlabel(r"$n$, number of input $T$ ($\ell=3$) or $\sqrt{T}$ ($\ell=4$) states", fontsize=7.5)
ax.set_xticks([20, 40, 60, 80])
ax.set_xlim(-22, 96); ax.tick_params(axis="x", labelsize=7)
ax.grid(True, axis="x", ls=":", lw=0.4, alpha=0.7)
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.tick_params(axis="y", length=0)
l3 = [y for y, r in zip(ys, rows) if r[4] == 3]; l4 = [y for y, r in zip(ys, rows) if r[4] == 4]
for yy, lab in ((l3[0] - 0.62, r"$\ell=3$"), (l4[0] - 0.85, r"$\ell=4$")):
    ax.text(-21, yy, lab, ha="left", va="center", fontsize=8.5, color="#0b0b0b", fontweight="bold")
ax.text(-23.5, ys[0] - 1.45, "factory", ha="right", va="center", fontsize=7.2, fontweight="bold", color="#0b0b0b")
ax.text(18, ys[0] - 1.45, "target", ha="right", va="center", fontsize=7.2, fontweight="bold", color="#0b0b0b")
ax.text(1.005, ys[0] - 1.45, "$c$ / $c_\\mathrm{prior}$", ha="left", va="center", fontsize=6.6, fontweight="bold", color="#0b0b0b", transform=ax.get_yaxis_transform())
h = [Line2D([0], [0], marker="D", color="w", markerfacecolor=NEW, markersize=6, label="this work"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor="white", markeredgecolor=PRIOR, markeredgewidth=1.2, markersize=6, label="distilling the parts separately"),
     Line2D([0], [0], marker="s", color="w", markerfacecolor="white", markeredgecolor=PRIOR, markeredgewidth=1.2, markersize=6, label="synthillation"),
     Line2D([0], [0], marker="^", color="w", markerfacecolor="white", markeredgecolor=PRIOR, markeredgewidth=1.2, markersize=6, label="synthillation + separate"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#0b0b0b", markersize=8, label="proved smallest for its circuit size")]
fig.legend(handles=h, loc="outside upper center", ncol=2, fontsize=6.3, frameon=False, handletextpad=0.3, columnspacing=0.9, labelspacing=0.3)
for out in (os.path.join(ROOT, "figures", "new_factories.png"),
            os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "new_factories.png")):
    if os.path.isdir(os.path.dirname(out)): fig.savefig(out, dpi=300); print("saved", out)
