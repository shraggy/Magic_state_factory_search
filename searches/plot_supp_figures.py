"""Supplemental Fig. S2 (figures/sequential.png), from the search cell of sequential.ipynb, and two_group_points(),
the gray background of Figs. E1 and S1, from the search cell of two_group.ipynb; both cells are executed unchanged.
Run: python3 searches/plot_supp_figures.py
"""
import json, contextlib, io, os
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
    ax.set_xlabel("$n$, number of input magic states", fontsize=8); ax.set_ylabel("$k$, number of output qubits", fontsize=8)
    ax.grid(True, alpha=0.3); ax.tick_params(labelsize=7.5)
    h = [Line2D([0], [0], marker=MK[l], color="w", markerfacecolor=CL[l], markeredgecolor="k", markersize=6,
                label=rf"$\ell={l}$") for l in (2, 3, 4)]
    ax.legend(handles=h, loc="upper left", fontsize=6.5, framealpha=0.95, handletextpad=0.3, labelspacing=0.3)
    for out in ("figures/sequential.png", os.path.join(TEX_FIG, "sequential.png")):
        if os.path.isdir(os.path.join(ROOT, os.path.dirname(out))):
            fig.savefig(os.path.join(ROOT, out), dpi=300); print("saved", out)


if __name__ == "__main__":
    fig_s2()
