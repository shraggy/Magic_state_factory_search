"""Fig. E1: every factory the two-group and unified searches find, all n (gray), with the literature factories they recover, by family.
Run: python3 searches/plot_literature.py"""
import json, csv, sys, os, contextlib, io, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
csv.field_size_limit(sys.maxsize)
plt.rcParams.update({"font.size": 8})
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "outputs") + os.sep
sys.path.insert(0, os.path.join(ROOT, "searches"))
from plot_supp_figures import two_group_points
# families: color slots 1..3 of the validated palette plus ink for the fourth, each also a distinct marker
FAM = {
 "triorthogonal / Reed–Muller": ("#2a78d6", "o"),   # Bravyi-Haah, Nezami-Haah, RM [[15,1,3]], [[7,1,3]], Litinski [[14,1,2]]
 "CCZ / CS / CCCZ codes":       ("#eb6834", "s"),   # cube [[8,3,2]], [[12,2,2]], iceberg/H-code, CCCZ [[16,4,2]]
 "synthillation":               ("#1baf7a", "^"),   # Campbell-Howard [[14,6,2]], [[18,4,2]]
 "catalytic (non-CSS)":         ("#0b0b0b", "D"),   # [[10,2,2]], [[11,1,2]]
 "malleable chains (this work)": ("#6250d6", "v"),   # sequential chains from one parent, e.g. l=4 [[16,4,2]]->[[24,3,2]]->[[28,2,2]]->[[30,1,2]]
}
known = {  # (l, N, k): family, label
 (2, 7, 1): ("triorthogonal / Reed–Muller", "[[7,1,3]]"), (3, 15, 1): ("triorthogonal / Reed–Muller", "[[15,1,3]]"),
 (3, 14, 1): ("triorthogonal / Reed–Muller", "[[14,1,2]]"), (3, 14, 2): ("triorthogonal / Reed–Muller", "[[14,2,2]]"),
 (3, 20, 4): ("triorthogonal / Reed–Muller", "[[20,4,2]]"), (3, 26, 6): ("triorthogonal / Reed–Muller", "[[26,6,2]]"),
 (3, 20, 3): ("triorthogonal / Reed–Muller", "[[20,3,2]]"), (3, 26, 5): ("triorthogonal / Reed–Muller", "[[26,5,2]]"),
 (3, 8, 3): ("CCZ / CS / CCCZ codes", "[[8,3,2]]"), (3, 12, 2): ("CCZ / CS / CCCZ codes", "[[12,2,2]]"),
 (2, 4, 2): ("CCZ / CS / CCCZ codes", "[[4,2,2]]"), (2, 6, 2): ("CCZ / CS / CCCZ codes", "[[6,2,2]]"),
 (4, 16, 4): ("CCZ / CS / CCCZ codes", "[[16,4,2]]"),
 (3, 14, 6): ("synthillation", "[[14,6,2]]"), (3, 18, 4): ("synthillation", "[[18,4,2]]"),
 (3, 10, 2): ("catalytic (non-CSS)", "[[10,2,2]]"), (3, 11, 1): ("catalytic (non-CSS)", "[[11,1,2]]"),
 (2, 6, 1): ("malleable chains (this work)", "[[6,1,2]]"),
 (3, 12, 3): ("malleable chains (this work)", "[[12,3,2]]"),
 (4, 24, 3): ("malleable chains (this work)", "[[24,3,2]]"), (4, 28, 2): ("malleable chains (this work)", "[[28,2,2]]"),
 (4, 30, 1): ("malleable chains (this work)", "[[30,1,2]]"), (4, 30, 2): ("malleable chains (this work)", "[[30,2,2]]"),
 (4, 31, 1): ("triorthogonal / Reed–Muller", "[[31,1,3]]"),
}
OFF = {(3, 15, 1): 0.22}
cwd = os.getcwd()
pts = two_group_points()
os.chdir(cwd)
allp = {(l, N, k) for N, k, l, o, d in pts}
csv.field_size_limit(sys.maxsize)
for r in csv.DictReader(open(R + "all_factories.csv")):          # n = input magic states
    allp.add((int(r["l"]), int(r["n"]), int(r["k"])))
fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.1), sharey=True, constrained_layout=True)
import matplotlib.transforms as mt
for ax, L in zip(axes, (2, 3, 4)):
    P = [p for p in allp if p[0] == L]
    ax.scatter([p[1] for p in P], [p[2] for p in P], s=14, c="#d6d5d1", edgecolors="none", zorder=1)
    ax.set_xscale("log"); ax.set_xlim(3, 2500); ax.set_ylim(0.4, 7.6)
    ax.axvspan(1, 100, color="#f1f0ed", zorder=0)
    fig.canvas.draw(); rend = fig.canvas.get_renderer(); placed = []
    for p in P:
        x, y = ax.transData.transform((p[1], p[2])); placed.append(mt.Bbox([[x - 2.5, y - 2.5], [x + 2.5, y + 2.5]]))
    items = sorted([(N, k, fam, lab) for (l, N, k), (fam, lab) in known.items() if l == L])
    skip = {"[[15,1,3]]"} | ({"[[28,2,2]]", "[[30,2,2]]", "[[30,1,2]]", "[[31,1,3]]"} if L == 4 else set())
    for N, k, fam, lab in items:
        c, m = FAM[fam]
        ax.scatter(N, k + OFF.get((L, N, k), 0.0), s=40, c=c, marker=m, edgecolors="white", linewidths=0.7, zorder=3)
    for N, k, fam, lab in items:
        x, y = ax.transData.transform((N, k + OFF.get((L, N, k), 0.0))); placed.append(mt.Bbox([[x - 7, y - 7], [x + 7, y + 7]]))
    axbb = ax.get_window_extent(rend)
    for N, k, fam, lab in [it for it in items if it[3] not in skip]:
        for dx, dy in [(0, 11), (0, -11), (24, 8), (-24, 8), (24, -8), (-24, -8), (32, 0), (-32, 0), (0, 19), (0, -19), (30, 16), (-30, 16), (30, -16), (-30, -16), (44, 6), (0, -27), (-44, 6), (-44, -6), (44, -6), (0, 27), (-36, 20), (36, 20), (-36, -20), (36, -20), (-50, 0), (50, 0), (0, -34), (0, 34)]:
            yk = k + OFF.get((L, N, k), 0.0)
            t = ax.annotate(lab, (N, yk), xytext=(dx, dy), textcoords="offset points", fontsize=5.6, ha="center", va="center", zorder=5)
            bb = t.get_window_extent(rend).expanded(1.02, 1.05)
            ok = all(not bb.overlaps(b) for b in placed) and axbb.contains(bb.x0, bb.y0) and axbb.contains(bb.x1, bb.y1)
            t.remove()
            if ok:
                placed.append(bb)
                ax.annotate(lab, (N, yk), xytext=(dx, dy), textcoords="offset points", fontsize=5.8, ha="center", va="center", zorder=5,
                            arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5, shrinkA=0, shrinkB=3.5))
                break
    ax.set_ylim(0.4, 7.6); ax.set_yticks(range(1, 8)); ax.set_xlabel("$n$, number of input magic states", fontsize=8)
    ax.text(0.015, 0.95, rf"$\ell={L}$", transform=ax.transAxes, fontsize=9, va="top")
    ax.grid(True, ls=":", lw=0.4, alpha=0.6); ax.tick_params(labelsize=7.5)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
axes[1].annotate("[[15,1,3]]", (15, 1 + OFF[(3, 15, 1)]), xytext=(40, -9), textcoords="offset points", fontsize=5.8, ha="center", va="center", zorder=5,
                 arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5, shrinkA=0, shrinkB=3.5))
for (N, k, lab, dx, dy) in ((28, 2, "[[28,2,2]]", -30, 10), (30, 2, "[[30,2,2]]", 34, 10), (30, 1, "[[30,1,2]]", -34, -2), (31, 1, "[[31,1,3]]", 36, -2)):
    axes[2].annotate(lab, (N, k), xytext=(dx, dy), textcoords="offset points", fontsize=5.8, ha="center", va="center", zorder=5,
                     arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5, shrinkA=0, shrinkB=3.5))
axes[0].set_ylabel("$k$, number of output qubits", fontsize=8)
h = [Line2D([0], [0], marker=m, color="w", markerfacecolor=c, markersize=6, label=f) for f, (c, m) in FAM.items()]
h.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#d6d5d1", markersize=5, label="all other factories found"))
h.append(plt.Rectangle((0, 0), 1, 1, fc="#f1f0ed", label=r"$n\leq100$"))
fig.legend(handles=h, loc="outside lower center", ncol=4, fontsize=6.5, frameon=False, handletextpad=0.3, columnspacing=1.0)
for out in (os.path.join(ROOT, "figures", "literature_recovered.png"),
            os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "literature_recovered.png")):
    if os.path.isdir(os.path.dirname(out)): fig.savefig(out, dpi=300); print("saved", out)
