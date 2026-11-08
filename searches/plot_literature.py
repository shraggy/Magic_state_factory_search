"""Fig. 4 (End Matter, fig:search-ext): every factory found at l=2, 3, 4, with the literature factories among them marked by family.

Points (n, k), one per distinct factory size:
  main sweep   the smallest factory of each exact output class, every stage of the cascade (outputs/staged_best.json;
               l=2, 3 for k<=6, l=4 for k<=4), shaded region;
  larger sweep every verified factory (distance 2, and the distance-3 Reed--Muller members) of the two-group and disjoint-gluing stages beyond (and within) the
               class lists (outputs/large_unified_factories.csv), and at l=3, k=7..9 the stage-(i) and gluing outputs
               of the cascade beyond the class lists (outputs/staged_B, outputs/k78_pipeline/full).
At l=3 every class of the main sweep has a synthillation construction (App. S4.2): sizes (n, k) holding a class with
fewer inputs than synthillation are orange, the other main-sweep sizes, where the factory ties synthillation, green. Literature families, each member identified by its output, not only by (n, k):
  Reed--Muller [[2^{l+1}-1,1,3]] and [[2^{l+1}-2,1,2]]   one output, two-group;
  hypercube [[2^l,l,2]] (CCZ at l=3, CCCZ at l=4, CZ at l=2);
  Bravyi--Haah [[3k+8,k,2]], k even                      output T^{otimes k} (degree-1 terms only);
  H codes [[k+4,k,2]] (l=2), k even                      output S^{otimes k} (degree-1 terms only);
  Iceberg [[k+2,k,2]] (l=2), k even                      output CZ on k/2 disjoint pairs;
  CS layers [[6m+6,2m,2]] ((2m+2)CS->mCS of Ref. [45] with T inputs)  output CS on m disjoint pairs;
  CCZ layers [[6m+2,3m,2]] (synthillation's case 11 family)            output CCZ on m disjoint triples;
  Nezami--Haah T-only entries (n, k) = (20,3), (26,5)    output T^{otimes k}.
Each member is checked here against the files (two-group circuits rebuilt with merge_large_unified.regenerate_circuit
and their output polynomial read off); a member not found raises an error.
The 22 l=4 classes flagged third_level_target in staged_best.json (third-level targets, covered at l=3) are left out
of the main sweep unless --with-third-level-targets is given.
    python searches/plot_literature.py   -> figures/literature_recovered.png (+ manuscript copy), counts printed
"""
import os, glob, sys, re, csv, json, collections, matplotlib, matplotlib.patches
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
csv.field_size_limit(sys.maxsize)
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import argparse as _ap
_p = _ap.ArgumentParser(); _p.add_argument("--sb", default="staged_best.json", help="staged-best file under outputs/")
_p.add_argument("--out-dir", default=None, help="write every output here instead of outputs/ (and the manuscript)")
_p.add_argument("--with-third-level-targets", action="store_true",
                help="keep the 22 l=4 classes whose target is a level-3 gate (default: left out)")
_A = _p.parse_args(); SB_FILE = _A.sb; OUTD = _A.out_dir; WITH_T3 = _A.with_third_level_targets
_p2 = None
ALL_FOUND = os.path.join(OUTD, "all_found_l3.json") if OUTD and os.path.exists(os.path.join(OUTD, "all_found_l3.json")) else os.path.join(O, "checks", "all_found_l3.json")
FS = 8
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": FS, "axes.labelsize": FS,
                     "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": FS,
                     "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}"})
INK, MUTED, GRAY, ORANGE, BAND = "#0b0b0b", "#52514e", "#b9b8b3", "#eb6834", "#eef3fb"
FAM = {  # family -> (color, marker); slots of the repo palette, one marker each
    "Reed--Muller": ("#2a78d6", "o"), "Hypercube": ("#7a4fc9", "s"), "Bravyi--Haah": ("#1baf7a", "^"),
    r"H codes": ("#c2185b", "v"), r"Iceberg": ("#8d6e00", "D"), r"$\mathrm{CS}$ layers (AG codes)": ("#0097a7", "P"),
    r"$\mathrm{CCZ}$ layers": ("#37474f", "X"), "Nezami--Haah": ("#e0a800", "*")}
SYN = "#2e7d32"                                        # l=3 main-sweep sizes where the factory ties synthillation


def output_degrees(row):
    """Degrees of the output polynomial of a larger-sweep two-group row, from its rebuilt circuit."""
    from merge_large_unified import regenerate_circuit
    from classify import phase_polynomial
    l, k = int(row["l"]), int(row["k"])
    g = regenerate_circuit(l=l, k=k, N=int(row["N"]), params=row["params"])
    P = phase_polynomial([S for S, _ in g], [c for _, c in g], set(range(k)), l)
    return {len(T) for T, a in P.items() if a}


def layer(red, gate, k, size):
    toks = re.findall(r"\\mathrm\{([A-Z]+)\}_\{(\d+)\}", red)
    return bool(toks) and all(t == gate and len(s) == size for t, s in toks) and len(toks) == k // size \
        and len({q for _, s in toks for q in s}) == k


def literature(LU, SB):
    """{(l, n, k): family}, each member found in the files with the right output (else AssertionError)."""
    d2 = [r for r in LU if r["d"] == "2"]
    def find(l, n, k, test):
        rows = [r for r in d2 if int(r["l"]) == l and int(r["n"]) == n and int(r["k"]) == k]
        return next((r for r in rows if test(r)), None)
    def single(r):                                     # output T^{otimes k} / S^{otimes k}: degree-1 terms only
        if r["search_stage"] == "two-group" and "terms" in r["reduced"]:
            return output_degrees(r) == {1}
        toks = re.findall(r"\\mathrm\{([A-Z]+)\}_\{(\d+)\}", r["reduced"])
        return bool(toks) and all(len(s) == 1 for _, s in toks) and len(toks) == int(r["k"])
    out = collections.defaultdict(list)                # (l, n, k) -> families (one size can belong to several)
    def put(fam, key, ok):
        assert ok, f"{fam} member {key} not found with its output"
        out[key].append(fam)
    for l in (2, 3, 4):
        put("Reed--Muller", (l, 2 ** (l + 1) - 1, 1),                 # distance 3
            any(int(r["l"]) == l and int(r["n"]) == 2 ** (l + 1) - 1 and r["k"] == "1" and r["d"] != "2" for r in LU))
        put("Reed--Muller", (l, 2 ** (l + 1) - 2, 1), find(l, 2 ** (l + 1) - 2, 1, lambda _: True))
        put("Hypercube", (l, 2 ** l, l), find(l, 2 ** l, l, lambda _: True))
    for k in range(2, 13, 2):
        put("Bravyi--Haah", (3, 3 * k + 8, k), find(3, 3 * k + 8, k, single))
        put(r"H codes", (2, k + 4, k), find(2, k + 4, k, single))
    for k in range(2, 13, 2):
        r = find(2, k + 2, k, lambda r: layer(r["reduced"], "CZ", k, 2))
        if r: out[(2, k + 2, k)].append(r"Iceberg")
    for m in range(1, 7):
        r = find(3, 6 * m + 6, 2 * m, lambda r, m=m: layer(r["reduced"], "CS", 2 * m, 2))
        if r: out[(3, 6 * m + 6, 2 * m)].append(r"$\mathrm{CS}$ layers (AG codes)")
    for m in range(1, 5):
        r = find(3, 6 * m + 2, 3 * m, lambda r, m=m: layer(r["reduced"], "CCZ", 3 * m, 3))
        if r: out[(3, 6 * m + 2, 3 * m)].append(r"$\mathrm{CCZ}$ layers")
    # Nezami--Haah: every size (n, k) of their Fig. 1 at which the class-list search (stages (i)-(ii) of sweep A, stage (iii) of the final pipeline) recorded a T^{otimes k} factory, smallest or not
    # (outputs/checks/all_found_l3.json, from searches/all_found_l3.py)
    for n, k in json.load(open(ALL_FOUND))["hit"]:
        out[(3, n, k)].append("Nezami--Haah")
    return out



def proved_on_N(SB):
    """(l, n, k) of classes whose free-SAT query over all coefficients is unsatisfiable below n on every N up to the
    factory's own N (z3 coefficient form of the staged cells, free_coef* files and iii_f step records): no factory on N or
    fewer qubits has fewer inputs. Not affected by the history-dependent synthesis parts (no output parts are fixed)."""
    best = {(r["l"], r["k"], str(r["orbit"])): r for r in SB}
    mx = {}
    def note(key, steps):
        m = max([s["N"] for s in steps or [] if s.get("last") == "unsat"], default=0)
        mx[key] = max(mx.get(key, 0), m)
    for src in ("staged_A", "staged_pc"):
        for f in glob.glob(os.path.join(O, src, "l*_k*", "free_coef*", "*.json")):
            d = json.load(open(f)); L, k = int(d["L"]), int(d["k"])
            for o, v in d["per_class"].items(): note((L, k, o), v.get("steps"))
    for f in glob.glob(os.path.join(O, "staged_small", "l*_k*", "iii_f.json")) + glob.glob(os.path.join(O, "staged_l4fix", "l*_k*", "iii_f.json")):
        d = json.load(open(f))
        for o, v in (d.get("steps") or {}).items(): note((d["L"], d["k"], o), v if isinstance(v, list) else v.get("steps"))
    return {(b["l"], b["n"], b["k"]) for key, b in best.items() if mx.get(key, 0) >= b["N"]}

def main():
    SB = [r for r in json.load(open(os.path.join(O, SB_FILE)))["rows"] if WITH_T3 or not r.get("third_level_target")]
    LU = list(csv.DictReader(open(os.path.join(O, "large_unified_factories.csv"))))
    K78 = []                                                   # l=3, k=7..9 beyond the class lists: stage (i), both gluing forms
    for K in (7, 8, 9):                                        # (outputs/staged_B, outputs/k78_pipeline/full; every factory verified there)
        for f in glob.glob(os.path.join(O, "staged_B", f"l3_k{K}", "*.json")):
            K78 += [dict(n=v["n"], k=K) for v in json.load(open(f))["best"].values()]
        for p in ("ii_dj", "ii_eg"):
            f = os.path.join(O, "k78_pipeline", "full", f"l3_k{K}", f"{p}.json")
            if os.path.exists(f):
                K78 += [dict(n=v["n"], k=K) for v in json.load(open(f))["best"].values()]
    main_pts = {(r["l"], r["n"], r["k"]) for r in SB}
    fewer = {(r["l"], r["n"], r["k"]) for r in SB if r["l"] == 3 and r["n"] < r["syn_n"]}
    tie = {(r["l"], r["n"], r["k"]) for r in SB if r["l"] == 3 and r["n"] == r["syn_n"]} - fewer
    bound = {(r["l"], r["n"], r["k"]) for r in SB if r["tau"] is not None and r["n"] == r["tau"] + 1}
    onN = proved_on_N(SB) - bound                              # smallest on its own N (free SAT over all coefficients UNSAT)
    for K in (7, 8):                                           # k = 7, 8 CCZ layers on the bound (T-count decided)
        for p in ("ii_eg", "iii_f"):
            for v in json.load(open(os.path.join(O, "k78_pipeline", "full", f"l3_k{K}", f"{p}.json")))["best"].values():
                if v["tau"][0] == v["tau"][1] and v["n"] == v["tau"][0] + 1:
                    bound.add((3, v["n"], K))
    LU = [r for r in LU if int(r["essential_dim"]) == int(r["k"])]   # genuine k-output factories only (no Clifford padding)
    large = {(int(r["l"]), int(r["n"]), int(r["k"])) for r in LU} | {(3, r["n"], r["k"]) for r in K78}
    lit = literature(LU, SB)
    pts = main_pts | large
    kmax = max(p[2] for p in pts); nmax = max(p[1] for p in pts); NCUT = 50
    fig = plt.figure(figsize=(7.0, 3.25))
    # panels placed by hand: a narrow gap at each axis break (linear | log), a normal gap between levels
    L0, R0, B0, T0 = 0.06, 0.975, 0.36, 0.90                 # figure fractions; the legend sits below B0
    BREAK, LEVEL = 0.008, 0.045                              # gap at the axis break, gap between levels
    W = (R0 - L0 - 3 * BREAK - 2 * LEVEL) / 3                # width of one level (linear + log part)
    axes = []
    for c, L in enumerate((2, 3, 4)):
        x0 = L0 + c * (W + BREAK + LEVEL)
        wa = W * 0.75
        a = fig.add_axes([x0, B0, wa, T0 - B0], sharey=axes[0][0] if axes else None)
        b = fig.add_axes([x0 + wa + BREAK, B0, W - wa, T0 - B0], sharey=a)
        axes.append((a, b))
        kk = 6 if L < 4 else 4
        for ax in (a, b):
            ax.axhspan(0.4, kk + 0.5, color=BAND, zorder=0, lw=0)
            P = sorted(p for p in pts if p[0] == L and p not in fewer and p not in tie)
            ax.scatter([p[1] for p in P], [p[2] for p in P], s=7, c=GRAY, edgecolors="none", zorder=1)
            Tq = sorted(p for p in tie if p[0] == L)
            ax.scatter([p[1] for p in Tq], [p[2] for p in Tq], s=9, c=SYN, edgecolors="none", zorder=2)
            F = sorted(p for p in fewer if p[0] == L)
            ax.scatter([p[1] for p in F], [p[2] for p in F], s=24, c=ORANGE, marker="*", edgecolors="none", zorder=2)
            B = sorted(p for p in bound if p[0] == L and (p[1] <= NCUT) == (ax is a))
            ax.scatter([p[1] for p in B], [p[2] for p in B], s=70, facecolors="none", edgecolors=INK, linewidths=0.8, zorder=4)
            Q = sorted(p for p in onN if p[0] == L and (p[1] <= NCUT) == (ax is a))
            ax.scatter([p[1] for p in Q], [p[2] for p in Q], s=40, marker="s", facecolors="none", edgecolors=INK, linewidths=0.7, zorder=4)
            ax.grid(True, ls=":", lw=0.4, color="#d9d8d4"); ax.tick_params(colors=MUTED, length=2.5)
            for sp in ("top",): ax.spines[sp].set_visible(False)
        for (l, n, k), fams in lit.items():
            if l == L:
                for j, fam in enumerate(fams):
                    cc, m = FAM[fam]; y = k + 0.22 * (j - (len(fams) - 1) / 2)
                    a.scatter(n, y, s=30 if m != "*" else 55, c=cc, marker=m, edgecolors="white", linewidths=0.6, zorder=3)
        a.set_xlim(0, NCUT); a.set_xticks([0, 10, 20, 30, 40, 50])
        b.set_xscale("log"); b.set_xlim(NCUT, nmax * 2); b.set_xticks([1e3, 1e5]); b.set_xticklabels([r"$10^3$", r"$10^5$"])
        a.spines["right"].set_visible(False); b.spines["left"].set_visible(False); b.tick_params(axis="y", left=False, labelleft=False)
        for sp in ("left", "bottom"): a.spines[sp].set_color(MUTED)
        b.spines["bottom"].set_color(MUTED); b.spines["right"].set_visible(False)
        d = 0.012                                        # break marks on the x axis
        for ax, x in ((a, 1.0), (b, 0.0)):
            ax.plot([x - d * (3 if ax is b else 1), x + d * (3 if ax is b else 1)], [-d * 2, d * 2], transform=ax.transAxes,
                    color=MUTED, clip_on=False, lw=0.8)
        b.text(0.5, 1.01, r"no literature" "\n" r"factory", transform=b.transAxes, ha="center", va="bottom",
               color=MUTED, fontsize=FS - 1)
        a.text(0.03, 0.97, rf"$l={L}$", transform=a.transAxes, ha="left", va="top", fontsize=FS + 1, zorder=10, bbox=dict(fc="white", ec="none", pad=0.5))
        a.set_xlabel(r"$n$, input magic states", x=0.62)
        a.set_ylim(0.4, kmax + 0.6); a.set_yticks(range(1, kmax + 1))
        if c: a.tick_params(axis="y", labelleft=False)
    axes[0][0].set_ylabel(r"$k$, output qubits")
    h = [Line2D([0], [0], marker=m, color="w", markerfacecolor=c, markeredgecolor="white", markersize=6 if m != "*" else 9,
                label=f) for f, (c, m) in FAM.items() if any(f in v for v in lit.values())]
    h += [Line2D([0], [0], marker="*", color="w", markerfacecolor=ORANGE, markersize=8, label="fewer inputs than synthillation"),
          Line2D([0], [0], marker="o", color="w", markerfacecolor=SYN, markersize=4.5, label="same inputs as synthillation"),
          Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor=INK, markersize=7, label=r"on the bound $n=\tau+1$"),
          Line2D([0], [0], marker="s", color="w", markerfacecolor="none", markeredgecolor=INK, markersize=5.5, label=r"smallest on its $N$ (proved by free SAT)"),
          Line2D([0], [0], marker="o", color="w", markerfacecolor=GRAY, markersize=4.5, label="other factories found"),
          matplotlib.patches.Rectangle((0, 0), 1, 1, fc=BAND, label="classes listed")]
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=4, frameon=False, handletextpad=0.2, columnspacing=0.6, handlelength=1.2)
    for out in ((os.path.join(OUTD, "literature_recovered.png"),) if OUTD else (os.path.join(ROOT, "figures", "literature_recovered.png"),
                os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "literature_recovered.png"))):
        if os.path.isdir(os.path.dirname(out)):
            fig.savefig(out, dpi=300); print("saved", out)
    cnt = dict(main={L: sum(r["l"] == L for r in SB) for L in (2, 3, 4)},
               larger={L: sum(1 for r in LU if r["d"] == "2" and int(r["l"]) == L) for L in (2, 3, 4)},
               literature=dict(collections.Counter(f for v in lit.values() for f in v)), lit_nmax=max(n for (_, n, _) in lit),
               l3_sizes=dict(fewer=len(fewer), tie_only=len(tie)), on_bound=sorted(bound))
    json.dump(dict(counts=cnt, literature=[dict(l=l, n=n, k=k, families=f) for (l, n, k), f in sorted(lit.items())]),
              open(os.path.join(OUTD or os.path.join(O, "checks"), "literature_points.json"), "w"), indent=0)
    print(cnt)


if __name__ == "__main__":
    main()
