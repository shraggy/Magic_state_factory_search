"""Wall-clock time of each unified-search stage against circuit size N, one line per k_max (l=3), from
outputs/stage_scaling.json (searches/stage_scaling.py). Writes figures/stage_scaling.png.
Run: python3 searches/plot_stage_scaling.py"""
import json, os, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
plt.rcParams.update({"font.size": 8})
D = json.load(open(os.path.join(ROOT, "outputs", "stage_scaling.json")))
STAGES = [("two_group", "two-group"), ("gluing", "shared-check gluing"), ("padding", "Clifford padding"),
          ("classification", "classification (optional)")]
COL = {4: "#9ec5ef", 6: "#2a78d6", 8: "#0b3d7a"}         # one sequential hue, light to dark with k_max
fig, axes = plt.subplots(1, 4, figsize=(7.0, 2.2), sharey=True, constrained_layout=True)
for ax, (key, title) in zip(axes, STAGES):
    for k in (4, 6, 8):
        P = sorted((p["N_max"], max(p["seconds"][key], 0.05)) for p in D if p["k_max"] == k)
        ax.plot([a for a, _ in P], [b for _, b in P], marker="o", ms=3, lw=1.6, color=COL[k], zorder=3)
    ax.set_yscale("log"); ax.set_title(title, fontsize=8, pad=3)
    ax.set_xlabel("$N$, circuit qubits", fontsize=7.5); ax.set_xticks([6, 8, 10, 12, 14, 16]); ax.tick_params(labelsize=7)
    ax.grid(True, which="major", ls=":", lw=0.4, alpha=0.7)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
axes[0].set_ylabel("wall-clock time (s)", fontsize=7.5)
for yy, lab in ((60, "1 min"), (3600, "1 h")):
    for ax in axes: ax.axhline(yy, color="#bdbcb8", lw=0.6, zorder=1)
    axes[-1].text(16.3, yy, lab, fontsize=6.3, va="center", color="#6b6a66")
h = [Line2D([0], [0], color=COL[k], marker="o", ms=3, lw=1.6, label=f"$k\\leq{k}$") for k in (4, 6, 8)]
fig.legend(handles=h, loc="outside upper center", ncol=3, fontsize=7, frameon=False)
for out in (os.path.join(ROOT, "figures", "stage_scaling.png"),
            os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "stage_scaling.png")):
    if os.path.isdir(os.path.dirname(out)): fig.savefig(out, dpi=300); print("saved", out)
