"""Fig. 1: each new factory (orange) against the best prior construction for its target (gray), with the saving.
Run: python3 searches/plot_new_factories.py"""
import os, matplotlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.size": 8})
# (label, our N, best prior N, prior kind, level)   prior kind: "parts" = distill parts separately; "syn" = synthillation
rows = [
 (r"CS$_{12}$T$_3$",                     20, 26, "parts", 3, True),
 (r"CS$_{12}$CS$_{34}$CS$_{56}$",        24, 26, "syn",   3, False),
 (r"CS$_{12}$CS$_{34}$T$_5$",            26, 28, "syn",   3, False),
 (r"CCZ$_{123}$CS$_{45}$T$_6$",          26, 34, "parts", 3, False),
 (r"CS$_{12}$T$_3$T$_4$T$_5$",           30, 32, "parts", 3, False),
 (r"CS$_{12}$T$_3\cdots$T$_7$",          36, 38, "parts", 3, False),
 (r"CS$_{12}$CS$_{34}$CS$_{56}$T$_7$",   38, 50, "parts", 3, False),
 (r"CT$_{12}\sqrt{T}_3$",                44, 58, "parts", 4, False),
 (r"CCS$_{123}\sqrt{T}_4$",              48, 54, "parts", 4, False),
 (r"CT$_{12}\sqrt{T}_3\sqrt{T}_4$",      50, 58, "parts", 4, False),
 (r"CT$_{12}$CT$_{34}\sqrt{T}_5$",       68, 78, "parts", 4, False),
 (r"CCS$_{123}$CT$_{45}\sqrt{T}_6$",     68, 82, "parts", 4, False),
 (r"CT$_{12}$CT$_{34}\sqrt{T}_5\sqrt{T}_6$", 70, 78, "parts", 4, False),
]
NS = [6, 8, 8, 9, 8, 10, 10, 7, 8, 8, 9, 10, 10]
KS = [3, 6, 5, 6, 5, 7, 7, 3, 4, 4, 5, 6, 6]
NEW = "#eb6834"; PRIOR = "#8c8b87"
fig, ax = plt.subplots(figsize=(3.4, 4.0), constrained_layout=True)
ys = []
y = 0
for i, (lab, N, Np, kind, L, proved) in enumerate(rows):
    if i == 0: y += 1.6
    if i > 0 and rows[i - 1][4] != L: y += 1.4
    ys.append(y)
    ax.plot([N, Np], [y, y], color="#cfcecb", lw=2.2, zorder=1, solid_capstyle="round")
    ax.scatter(Np, y, s=26, c="white", edgecolors=PRIOR, linewidths=1.2, zorder=2, marker="s" if kind == "syn" else "o")
    ax.scatter(N, y, s=34, c=NEW, edgecolors="white", linewidths=0.6, zorder=3, marker="D")
    ax.text(Np + 1.2, y, f"$-{Np - N}$", va="center", fontsize=6.5, color="#52514e")
    ax.text(N - (3.6 if proved else 2.0), y, lab, va="center", ha="right", fontsize=6.3, color="#0b0b0b")
    if proved: ax.scatter(N, y, s=95, facecolors="none", edgecolors="#0b0b0b", linewidths=0.8, zorder=4)
    y += 1
ax.set_yticks(ys); ax.set_yticklabels([f"[[{r[1]},{KS[j]},2]]  $N={NS[j]}$" for j, r in enumerate(rows)], fontsize=6.6)
ax.invert_yaxis()
ax.set_xlabel(r"$n$, number of input $T$ ($\ell=3$) or $\sqrt{T}$ ($\ell=4$) states", fontsize=7.5)
ax.set_xticks([20, 40, 60, 80])
ax.set_xlim(-22, 90); ax.tick_params(axis="x", labelsize=7)
ax.grid(True, axis="x", ls=":", lw=0.4, alpha=0.7)
for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
ax.tick_params(axis="y", length=0)
l3 = [y for y, r in zip(ys, rows) if r[4] == 3]; l4 = [y for y, r in zip(ys, rows) if r[4] == 4]
for yy, lab in ((l3[0] - 0.62, r"$\ell=3$"), (l4[0] - 0.85, r"$\ell=4$")):
    ax.text(-21, yy, lab, ha="left", va="center", fontsize=8.5, color="#0b0b0b", fontweight="bold")
ax.text(-23.5, ys[0] - 1.45, "factory", ha="right", va="center", fontsize=7.2, fontweight="bold", color="#0b0b0b")
ax.text(18, ys[0] - 1.45, "target", ha="right", va="center", fontsize=7.2, fontweight="bold", color="#0b0b0b")
h = [Line2D([0], [0], marker="D", color="w", markerfacecolor=NEW, markersize=6, label="this work"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor="white", markeredgecolor=PRIOR, markeredgewidth=1.2, markersize=6, label="distilling the parts separately"),
     Line2D([0], [0], marker="s", color="w", markerfacecolor="white", markeredgecolor=PRIOR, markeredgewidth=1.2, markersize=6, label="synthillation"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#0b0b0b", markersize=8, label="proved smallest for its circuit size")]
fig.legend(handles=h, loc="outside upper center", ncol=2, fontsize=6.3, frameon=False, handletextpad=0.3, columnspacing=0.9, labelspacing=0.3)
for out in (os.path.join(ROOT, "figures", "new_factories.png"),
            os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "new_factories.png")):
    if os.path.isdir(os.path.dirname(out)): fig.savefig(out, dpi=300); print("saved", out)
