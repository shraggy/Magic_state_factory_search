"""Regenerate figures/two_group.png and figures/two_group_ext.png (not in the paper) and Supplemental Fig. S2
(figures/sequential.png); two_group_points() also supplies the gray background of Figs. E1 and S1.

Both are drawn at one-column width with print-legible fonts. The data come from the search cells of
two_group.ipynb and sequential.ipynb, executed here unchanged. Run: python3 searches/plot_supp_figures.py
"""
import json, contextlib, io, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX_FIG = "Magic_State_Factories_tex/MSD_manuscript/FIGURES"


def two_group_points():
    """(N, k, l, output degree, d) for every nontrivial two-group factory of two_group.ipynb."""
    os.chdir(os.path.join(ROOT, "searches"))
    nb = json.load(open("two_group.ipynb"))
    src = "".join(nb["cells"][0]["source"]).replace('if __name__ == "__main__":', 'if True:')
    g = {}
    with contextlib.redirect_stdout(io.StringIO()):
        exec(src, g)
    res = g["all_results"]
    return sorted({(r["N"], r["k"], r["l"], r["output"], r["d"]) for r in res if r["output"] != "trivial"})


def fig_s1():
    pts = two_group_points()
    plt.rcParams.update({"font.size": 8})
    MK = {"deg1": "o", "deg2": "s", "deg3": "^", "deg4": "D"}
    CL = {"deg1": "#e41a1c", "deg2": "#377eb8", "deg3": "#2ca02c", "deg4": "#ff7f00"}
    LB = {"deg1": r"deg 1: $S/T/\sqrt{T}$", "deg2": "deg 2: CZ/CS/CT", "deg3": "deg 3: CCZ/CCS", "deg4": "deg 4: CCCZ"}
    JIT = {"deg1": -0.18, "deg2": -0.06, "deg3": 0.06, "deg4": 0.18}
    fig, axes = plt.subplots(3, 1, figsize=(3.4, 5.2), sharex=True, constrained_layout=True)
    for ax, l in zip(axes, (2, 3, 4)):
        ax.axvspan(1, 100, color="#eeeeee", zorder=0)
        for N, k, lv, out, d in pts:
            if lv != l: continue
            ax.scatter(N, k + JIT[out], marker=MK[out], s=11, c=CL[out], edgecolors="k" if d >= 3 else CL[out],
                       linewidths=0.8 if d >= 3 else 0.2, zorder=3)
        ax.set_xscale("log"); ax.set_xlim(3, 2500); ax.set_ylim(0.4, 7.6); ax.set_yticks(range(1, 8))
        ax.tick_params(labelsize=7.5); ax.grid(True, ls=":", lw=0.4, alpha=0.6)
        ax.set_ylabel("$k$", fontsize=8.5)
        ax.text(0.02, 0.92, rf"$\ell={l}$", transform=ax.transAxes, fontsize=9, va="top",
                bbox=dict(fc="white", ec="none", pad=1))
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    axes[-1].set_xlabel("$n$ (input magic states)", fontsize=8.5)
    h = [Line2D([0], [0], marker=MK[o], color="w", markerfacecolor=CL[o], markeredgecolor=CL[o], markersize=5, label=LB[o]) for o in MK]
    h.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="#e41a1c", markeredgecolor="k", markeredgewidth=1, markersize=5, label="distance 3"))
    h.append(plt.Rectangle((0, 0), 1, 1, fc="#eeeeee", label=r"$n\leq100$ (Fig. 1)"))
    fig.legend(handles=h, loc="outside lower center", ncol=2, fontsize=7, frameon=False, handletextpad=0.3, columnspacing=0.8)
    for out in ("figures/two_group_ext.png", os.path.join(TEX_FIG, "two_group_ext.png")):
        if os.path.isdir(os.path.join(ROOT, os.path.dirname(out))):
            fig.savefig(os.path.join(ROOT, out), dpi=300); print("saved", out)


def fig_s2():
    os.chdir(os.path.join(ROOT, "searches"))
    nb = json.load(open("sequential.ipynb"))
    g = {}
    with contextlib.redirect_stdout(io.StringIO()):
        exec("".join(nb["cells"][0]["source"]), g)
    pts = sorted({(r["N"], r["k"], r["l"]) for r in g["all_results"] if r["N"] <= 100})
    plt.rcParams.update({"font.size": 8})
    CL = {2: "#1fa354", 3: "#c0392b", 4: "#8e44ad"}; MK = {2: "o", 3: "s", 4: "^"}
    YJ = {2: -0.17, 3: 0.0, 4: 0.17}
    fig, ax = plt.subplots(figsize=(3.4, 2.6), constrained_layout=True)
    for N, k, l in pts:
        ax.scatter(N, k + YJ[l], marker=MK[l], s=24, c=CL[l], edgecolors="k", linewidths=0.4, zorder=3)
    LAB = [(8, 3, 3, (0.80, -0.45)), (8, 4, 3, (1.0, 0.45)), (12, 2, 3, (1.0, -0.5)), (14, 1, 3, (1.0, -0.45)),
           (12, 3, 2, (1.0, -0.45)), (4, 2, 2, (1.0, -0.5)), (16, 5, 4, (1.0, 0.45)), (30, 1, 4, (1.0, -0.4))]
    for N, k, l, (dx, dy) in LAB:
        ax.text(N * dx, k + YJ[l] + dy, f"[[{N},{k},2]]", ha="center", va="center", fontsize=6.3,
                color=CL[l], fontweight="bold", zorder=5)
    ax.set_xscale("log"); ax.set_xlim(3, 110); ax.set_ylim(0.4, 7.6); ax.set_yticks(range(1, 8))
    ax.set_xlabel("$n$ (input magic states)", fontsize=8); ax.set_ylabel("$k$ (output qubits)", fontsize=8)
    ax.grid(True, alpha=0.3); ax.tick_params(labelsize=7.5)
    h = [Line2D([0], [0], marker=MK[l], color="w", markerfacecolor=CL[l], markeredgecolor="k", markersize=6,
                label=rf"$\ell={l}$") for l in (2, 3, 4)]
    ax.legend(handles=h, loc="upper left", fontsize=6.5, framealpha=0.95, handletextpad=0.3, labelspacing=0.3)
    for out in ("figures/sequential.png", os.path.join(TEX_FIG, "sequential.png")):
        if os.path.isdir(os.path.join(ROOT, os.path.dirname(out))):
            fig.savefig(os.path.join(ROOT, out), dpi=300); print("saved", out)


def fig_1():
    """Main-text Fig. 1: two-group factories with N <= 100, smallest N per (k, distance) labeled."""
    pts = [p for p in two_group_points() if p[0] <= 100]
    plt.rcParams.update({"font.size": 9})
    MK = {"deg1": "o", "deg2": "s", "deg3": "^", "deg4": "D"}
    CL = {"deg1": "#e41a1c", "deg2": "#377eb8", "deg3": "#2ca02c", "deg4": "#ff7f00"}
    LB = {"deg1": r"deg 1: $S/T/\sqrt{T}$", "deg2": "deg 2: CZ/CS/CT", "deg3": "deg 3: CCZ/CCS", "deg4": "deg 4: CCCZ"}
    JIT = {"deg1": -0.16, "deg2": -0.05, "deg3": 0.05, "deg4": 0.16}
    CAT = [(10, 2, 3), (11, 1, 3)]
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.9), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.99, top=0.92, bottom=0.36, wspace=0.08)
    for ax, l in zip(axes, (2, 3, 4)):
        P = [p for p in pts if p[2] == l]
        for N, k, lv, out, d in P:
            ax.scatter(N, k + JIT[out], marker=MK[out], s=22, c=CL[out], edgecolors="k" if d >= 3 else CL[out],
                       linewidths=1.0 if d >= 3 else 0.3, zorder=3)
        for N, k, lc in CAT:
            if lc == l:
                ax.scatter(N, k + JIT["deg1"], marker="o", s=26, facecolors="white", edgecolors=CL["deg1"], linewidths=1.2, zorder=4)
        ax.set_xscale("log"); ax.set_xlim(3, 120); ax.set_ylim(0.4, 7.7)
        best = {}
        for N, k, lv, out, d in P:
            key = (k, d >= 3)
            if key not in best or N < best[key][0]: best[key] = (N, k, out, d)
        lab = [(N, k + JIT[o], f"[[{N},{k},{d}]]", CL[o], False) for N, k, o, d in best.values()]
        lab += [(N, k + JIT["deg1"], f"[[{N},{k},2]]", CL["deg1"], True) for N, k, lc in CAT if lc == l]
        fig.canvas.draw(); placed = []
        import matplotlib.transforms as mt
        for N, k, lv, out, d in P:
            x, y = ax.transData.transform((N, k + JIT[out])); placed.append(mt.Bbox([[x - 7, y - 7], [x + 7, y + 7]]))
        rend = fig.canvas.get_renderer()
        for N, y, txt, col, it in sorted(lab, key=lambda t: t[0]):
            for dx, dy in [(0, 11), (0, -11), (22, 9), (-22, 9), (22, -9), (-22, -9), (32, 0), (-32, 0), (0, 19), (0, -19), (26, 18), (-26, 18), (26, -18), (-26, -18)]:
                t = ax.annotate(txt, (N, y), xytext=(dx, dy), textcoords="offset points", fontsize=6.3, ha="center",
                                color=col, style="italic" if it else "normal", zorder=5, va="center")
                bb = t.get_window_extent(rend).expanded(1.04, 1.1)
                axbb = ax.get_window_extent(rend)
                ok = all(not bb.overlaps(b) for b in placed) and axbb.contains(bb.x0, bb.y0) and axbb.contains(bb.x1, bb.y1)
                if ok:
                    placed.append(bb); t.remove()
                    ax.annotate(txt, (N, y), xytext=(dx, dy), textcoords="offset points", fontsize=6.3, ha="center", va="center",
                                color=col, style="italic" if it else "normal", zorder=5,
                                arrowprops=dict(arrowstyle="-", color="#999999", lw=0.5, shrinkA=0, shrinkB=3.5))
                    break
                t.remove()
        ax.set_yticks(range(1, 8))
        ax.set_xlabel("$n$ (input magic states)", fontsize=9); ax.tick_params(labelsize=8.5)
        ax.grid(True, ls=":", lw=0.4, alpha=0.6)
        ax.set_title(rf"$\ell={l}$", fontsize=10, pad=3)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    axes[0].set_ylabel("$k$ (output qubits)", fontsize=9)
    h = [Line2D([0], [0], marker=MK[o], color="w", markerfacecolor=CL[o], markeredgecolor=CL[o], markersize=6, label=LB[o]) for o in MK]
    h.append(Line2D([0], [0], marker="o", color="w", markerfacecolor=CL["deg1"], markeredgecolor="k", markeredgewidth=1.2, markersize=6, label="distance 3 (ringed)"))
    h.append(Line2D([0], [0], marker="o", color="w", markerfacecolor="white", markeredgecolor=CL["deg1"], markeredgewidth=1.2, markersize=6, label="catalytic, non-CSS (open)"))
    fig.legend(handles=h, loc="lower center", ncol=3, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.0), columnspacing=1.2, handletextpad=0.3)
    for out in ("figures/two_group.png", os.path.join(TEX_FIG, "two_group.png")):
        if os.path.isdir(os.path.join(ROOT, os.path.dirname(out))):
            fig.savefig(os.path.join(ROOT, out), dpi=300); print("saved", out)


if __name__ == "__main__":
    fig_1()
    fig_s1()
    fig_s2()
