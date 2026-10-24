"""figures/landscape_l3.png (not in the paper): smallest factory for each l=3 and l=4 target state, marked by its relation to prior work.
Run: python3 searches/plot_landscape.py"""
import json, os, matplotlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
pts = json.load(open(os.path.join(ROOT, "outputs", "landscape_l3.json")))
plt.rcParams.update({"font.size": 8})
# three categories, the reference palette's first three slots (all-pairs validated for scatter)
COL = {"new": "#eb6834", "matched": "#2a78d6", "matched_worse": "#2a78d6", "known": "#8c8b87", "open": "#1baf7a", "pending": "#1baf7a"}
MK = {"new": "D", "matched": "o", "matched_worse": "o", "known": "o", "open": "s", "pending": "s"}
fig, (ax, ax4) = plt.subplots(2, 1, figsize=(3.4, 4.3), constrained_layout=True, gridspec_kw=dict(height_ratios=[1.35, 1]))
ax.axhspan(0.5, 8.5, xmin=0, xmax=1, color="none")
# proved-smallest region: rows with n <= 7 are certified; shade as a band behind their markers
for p in pts:
    if p["proved"]:
        y = p["k"] + {"known": -0.12, "matched": 0.0, "matched_worse": 0.0, "open": 0.12, "pending": 0.12, "new": 0.0}[p["status"]]
        ax.scatter(p["N"], y, s=70, facecolors="none", edgecolors="#0b0b0b", linewidths=0.9, zorder=4)
for st in ("known", "matched", "matched_worse", "open", "pending", "new"):
    P = [p for p in pts if p["status"] == st]
    for p in P:
        y = p["k"] + {"known": -0.12, "matched": 0.0, "matched_worse": 0.0, "open": 0.12, "pending": 0.12, "new": 0.0}[st]
        ax.scatter(p["N"], y, marker=MK[st], s=34 if st == "new" else 26, c=COL[st], edgecolors="white", linewidths=0.8, zorder=3)
LAB = {"CS T": (-17, 0), "CS CS CS": (-12, 9), "CS CS T": (0, -11), "CCZ CS T": (12, 9), "CS TTT": (14, 0),
       "CS CS CS T": (0, -11), "CS T^5": (0, 9)}
NICE = {"CS T": r"CS$\,$T", "CS CS CS": r"CS$^{3}$", "CS CS T": r"CS$^2$T", "CCZ CS T": r"CCZ$\,$CS$\,$T",
        "CS TTT": r"CS$\,$T$^3$", "CCZ T^4": r"CCZ$\,$T$^4$", "CS CS CS T": r"CS$^3$T", "CS T^5": r"CS$\,$T$^5$"}
for p in pts:
    if p["label"] in LAB:
        ax.annotate(NICE[p["label"]], (p["N"], p["k"]), xytext=LAB[p["label"]], textcoords="offset points", fontsize=6.5,
                    ha="center", va="center", color="#0b0b0b", zorder=5)
ax.set_xlim(6, 41); ax.set_ylim(0.4, 8.7); ax.set_yticks(range(1, 9))
ax.set_xlabel("$N$ (input $T$ states)", fontsize=8); ax.set_title(r"$\ell=3$", fontsize=9, pad=2); ax.set_ylabel("$k$ (output qubits)", fontsize=8)
ax.grid(True, ls=":", lw=0.4, alpha=0.7); ax.tick_params(labelsize=7.5)
for s in ("top", "right"): ax.spines[s].set_visible(False)
h = [Line2D([0], [0], marker="D", color="w", markerfacecolor=COL["new"], markersize=6, label="new: beats prior constructions"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["matched"], markersize=6, label="new here, matched by synthillation"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["known"], markersize=6, label="previously known"),
     Line2D([0], [0], marker="s", color="w", markerfacecolor=COL["open"], markersize=6, label="new, no prior construction to compare"),
     Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#0b0b0b", markeredgewidth=0.9, markersize=8, label="proved smallest on its $n$")]
ax.legend(handles=h, loc="lower right", bbox_to_anchor=(1.0, 0.0), fontsize=5.8, frameon=True, framealpha=0.95, edgecolor="#dddddd", handletextpad=0.3, labelspacing=0.35)

pts4 = json.load(open(os.path.join(ROOT, "outputs", "landscape_l4.json")))
SAV4 = {"CT12 sqrtT_{3}", "CCS123 sqrtT_{4}", "CT12 sqrtT_{3}sqrtT_{4}", "CT12 CT34 sqrtT_{5}", "CCS123 CT45 sqrtT_{6}", "CT12 CT34 sqrtT_{5}sqrtT_{6}"}
for p in pts4:
    st = "known" if p["status"] == "known" else ("new" if p["dec"] in SAV4 else "single4")
    c = {"known": COL["known"], "new": COL["new"], "single4": "#6250d6"}[st]
    mk = {"known": "o", "new": "D", "single4": "o"}[st]
    yj = {"new": 0.13, "single4": -0.13, "known": 0.0}[st]
    ax4.scatter(p["N"], p["k"] + yj, marker=mk, s=34 if st == "new" else 26, c=c, edgecolors="white", linewidths=0.8, zorder=3)
ax4.set_xlim(10, 102); ax4.set_ylim(0.4, 6.7); ax4.set_yticks(range(1, 7))
ax4.set_xlabel(r"$N$ (input $\sqrt{T}$ states)", fontsize=8); ax4.set_ylabel("$k$ (output qubits)", fontsize=8)
ax4.set_title(r"$\ell=4$", fontsize=9, pad=2)
ax4.grid(True, ls=":", lw=0.4, alpha=0.7); ax4.tick_params(labelsize=7.5)
for sp in ("top", "right"): ax4.spines[sp].set_visible(False)
h4 = [Line2D([0], [0], marker="D", color="w", markerfacecolor=COL["new"], markersize=6, label="new: beats distilling the parts"),
      Line2D([0], [0], marker="o", color="w", markerfacecolor="#6250d6", markersize=6, label="new: single gate type"),
      Line2D([0], [0], marker="o", color="w", markerfacecolor=COL["known"], markersize=6, label="previously known")]
ax4.legend(handles=h4, loc="lower right", fontsize=5.8, frameon=True, framealpha=0.95, edgecolor="#dddddd", handletextpad=0.3, labelspacing=0.35)
for out in (os.path.join(ROOT, "figures", "landscape_l3.png"),
            os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "landscape_l3.png")):
    if os.path.isdir(os.path.dirname(out)): fig.savefig(out, dpi=300)
