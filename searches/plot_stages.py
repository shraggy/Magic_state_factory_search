"""Fig. S1: every factory the unified search finds with n <= 120, in one plot. Color is the Clifford level l,
marker shape the search stage; only factories at an (n, k) the two-group stage does not reach are labeled.
Run: python3 searches/plot_stages.py"""
import os, sys, csv, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.transforms as mt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.join(ROOT, "searches"))
from plot_supp_figures import two_group_points
_cwd = os.getcwd(); TG = two_group_points(); os.chdir(_cwd)
plt.rcParams.update({"font.size": 8})
rows = list(csv.DictReader(open(os.path.join(ROOT, "outputs", "smallest_per_class.csv"))))
COL = {"2": "#2a78d6", "3": "#eb6834", "4": "#1baf7a"}          # reference palette slots 1-3 (all-pairs validated)
MK = {"two-group": "o", "gluing": "s", "padding": "^", "sat": "D"}
NAME = {"two-group": "two-group", "gluing": "block gluing", "padding": "Clifford padding", "sat": "SAT"}
LJ = {"2": -0.2, "3": 0.0, "4": 0.2}
def stage(r):
    # search_stage column of unified_search.py: two-group, gluing, padding, SAT
    return r["search_stage"].lower()
fig, ax = plt.subplots(figsize=(3.4, 3.4), constrained_layout=True)
tg = set()
for N, k, l, o, d in TG:
    if N <= 120:
        L = str(l); tg.add((L, N, k))
        ax.scatter(N, k + LJ[L], s=9, marker="o", c=COL[L], alpha=0.22, edgecolors="none", zorder=1)
pts = []
for r in rows:
    st = stage(r)
    if st == "two-group" or int(r["n"]) > 120: continue
    L, N, k = r["l"], int(r["n"]), int(r["k"])
    pts.append((L, N, k, st))
    ax.scatter(N, k + LJ[L], s=30 if st != "sat" else 38, marker=MK[st], c=COL[L], edgecolors="white", linewidths=0.6, zorder=3)
ax.set_xscale("log"); ax.set_xlim(3, 125); ax.set_ylim(0.4, 8.6); ax.set_yticks(range(1, 9))
fig.canvas.draw(); rend = fig.canvas.get_renderer(); placed = []
for L, N, k, st in pts:
    x, y = ax.transData.transform((N, k + LJ[L])); placed.append(mt.Bbox([[x - 3, y - 3], [x + 3, y + 3]]))
axbb = ax.get_window_extent(rend)
lab = sorted({(L, N, k, st) for L, N, k, st in pts if st in ("gluing", "sat") and (L, N, k) not in tg and N <= 60})
for L, N, k, st in lab:
    y = k + LJ[L]
    for dx, dy in [(0, 8), (0, -8), (18, 6), (-18, 6), (18, -6), (-18, -6), (24, 0), (-24, 0), (0, 15), (0, -15), (26, 12), (-26, 12), (26, -12), (-26, -12)]:
        t = ax.annotate(f"[[{N},{k},2]]", (N, y), xytext=(dx, dy), textcoords="offset points", fontsize=4.9, ha="center", va="center")
        bb = t.get_window_extent(rend).expanded(1.0, 1.02); t.remove()
        if all(not bb.overlaps(b) for b in placed) and axbb.contains(bb.x0, bb.y0) and axbb.contains(bb.x1, bb.y1):
            placed.append(bb)
            ax.annotate(f"[[{N},{k},2]]", (N, y), xytext=(dx, dy), textcoords="offset points", fontsize=4.9, ha="center", va="center", color=COL[L],
                        arrowprops=dict(arrowstyle="-", color="#999999", lw=0.4, shrinkA=0, shrinkB=3))
            break
ax.set_xlabel("$n$ (input magic states)", fontsize=8); ax.set_ylabel("$k$ (output qubits)", fontsize=8)
ax.grid(True, ls=":", lw=0.4, alpha=0.6); ax.tick_params(labelsize=7.5)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
h = [Line2D([0], [0], marker="o", color="w", markerfacecolor=COL[L], markersize=6, label=rf"$\ell={L}$") for L in ("2", "3", "4")]
h += [Line2D([0], [0], marker=MK[s], color="w", markerfacecolor="#6b6a66", markersize=6, alpha=0.35 if s == "two-group" else 1, label=NAME[s]) for s in MK]
fig.legend(handles=h, loc="outside lower center", ncol=4, fontsize=6.3, frameon=False, handletextpad=0.2, columnspacing=0.8)
for out in (os.path.join(ROOT, "figures", "stages.png"), os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "stages.png")):
    if os.path.isdir(os.path.dirname(out)): fig.savefig(out, dpi=300); print("saved", out)
