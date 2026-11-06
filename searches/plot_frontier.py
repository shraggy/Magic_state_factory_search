"""Fig. 2: distance-2 factories at l=3 against the bound n >= T-count + 1 and against the synthillation construction,
for every output class with k <= 6 and for every k = 7, 8 target that the search holds.

DATA.
  k <= 6: outputs/staged_best.json (rows with l = 3; 766 classes, the smallest factory of the staged cascade per class,
    exact T-count tau, syn_n = smallest valid synthillation construction, Campbell-Howard Thm 1 / Table I).
  k = 7, 8: outputs/checks/k78_frontier.json (searches/k78_frontier.py), one record per TARGET (reduced phase
    polynomial up to qubit permutation; no class list exists beyond k = 6), from both parts of the search: sweep B
    (stage (i) and the larger two-group / disjoint sweep) and the cascade on the CCZ-layer target list (gluing,
    restricted and free SAT). Per target: n of the smallest factory, T-count bounds lo <= tau <= hi, and a verdict
    against synthillation: fewer (n < syn_lo, below every Thm 1 / Table I construction), same (n = syn_hi, lo = hi),
    fewer_than_built (n < syn_hi only), more (n > syn_hi), open (n = syn_hi with lo < hi). 337 targets at k = 7,
    4,017 at k = 8.
(a) Each class / target at (n, T-count); x linear to n = 34 and compressed 3x beyond (break mark): every class with
    k <= 6 and 4,290 of the 4,354 targets at k = 7, 8 lie in the linear part. One glyph per (n, T-count, verdict,
    decided), standing for 1 to 60 classes / targets; the glyphs of one position sit side by side, centred on n (decided
    first, then upper-bound ones) and are moved apart along their row only as far as needed for no two to touch (the
    script prints the largest distance from n, 1.89). Same encoding in (a) and (b): orange fewer, gray same, dark more
    (triangle in (a)); filled / solid bar when the T-count is decided (lo = hi), hollow / dotted bar when the target is
    compared at its upper bound hi only (lo < hi: open, fewer_than_built, more). The 3,271 CCZ layers at k = 8 with
    n = hi + 1 and lo = 17 (asserted) are drawn as one hollow glyph per n at hi, 1.35 to the right of the line, with a
    bar down to lo = 17, so they are not read as on the bound. Dashed: Corollary 3, n = T-count + 1; shaded (neutral):
    excluded. Blue: on the bound with tau decided, one 8 pt label [[n, k, 2]] per point listing every k reaching it
    (N = k + 1 on all, asserted; every one 'same'); circles are in the synthillation paper (n = 8, 12, 14), diamonds
    are not.
(b) Per row (l = 3 unless marked), the share by verdict of the classes (k <= 6) or targets (k = 7, 8), with the total
    and the number with fewer inputs (inside the orange segment, or before the total where it is too narrow; with an
    asterisk when the comparison is with the construction built from the upper bound only, 'fewer_than_built').
    Segments under 2 % are drawn at 2 % (taken from the widest segment of the row). l = 4: classes with a factory, no
    synthillation benchmark, counting only the 4,145 classes that need a fourth-level rotation (the 22 flagged
    third_level_target in staged_best.json are third-level targets, covered at l = 3).
    python searches/plot_frontier.py                 # writes the paper's figure files (default, unchanged)
    python searches/plot_frontier.py --out PREFIX    # writes PREFIX.png and PREFIX.pdf only
"""
import os, json, collections, argparse, numpy as np, matplotlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
ap = argparse.ArgumentParser(); ap.add_argument("--out", default=None, help="write PREFIX.png and PREFIX.pdf only")
ARGS = ap.parse_args()
matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.patches, matplotlib.transforms
from matplotlib.lines import Line2D
FS = 8                                                     # caption size of the paper (revtex \small, 8pt)
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": FS, "axes.labelsize": FS,
                     "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": FS,
                     "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}"})
INK, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#d9d8d4", "#fcfcfb"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"          # validated categorical slots 1-3 (light)
GRAY, GRAYBAR, MORE, OPEN, LMORE = "#a8a7a2", "#c9c8c3", "#52514e", "#8f8e8a", "#9b9a96"   # same (dot, bar); more (x); open (hollow edge); more at upper bound (x)

# k <= 6: the staged cascade, one row per class
SB = json.load(open(os.path.join(O, "staged_best.json")))["rows"]
R = [dict(r, tcount=r["tau"], syn_strongest_n=r["syn_n"]) for r in SB if r["l"] == 3]
assert len(R) == 766
win = [r for r in R if r["n"] < r["syn_strongest_n"]]
same = [r for r in R if r["n"] == r["syn_strongest_n"]]
assert not [r for r in R if r["n"] > r["syn_strongest_n"]]
assert all(r["n"] >= r["tcount"] + 1 for r in R)
at = {r["n"]: r for r in R if r["n"] == r["tcount"] + 1}
KNOWN = {8, 12, 14}                                        # given in the synthillation paper (case 11)
FRONT = {n: (None, r["k"], r["N"], r["c"], n in KNOWN) for n, r in at.items()}
assert len(at) == len([r for r in R if r["n"] == r["tcount"] + 1]) == 5

# k = 7, 8: every target of both parts of the search (sweep B and the target-list cascade), outputs/checks/k78_frontier.json
KF = json.load(open(os.path.join(O, "checks", "k78_frontier.json")))
T78 = KF["targets"]
for K in (7, 8):                                           # the plotted records are exactly the file's totals
    TK = [r for r in T78 if r["k"] == K]; tot = KF["totals"][str(K)]
    assert len(TK) == tot["targets"], K
    assert collections.Counter(r["verdict"] for r in TK) == collections.Counter(tot["by_verdict"]), K
    assert sum(r["lo"] == r["hi"] for r in TK) == tot["lo_eq_hi"] and sum(r["on_bound"] for r in TK) == tot["on_bound"], K
assert all(r["n"] >= r["hi"] + 1 for r in T78)             # Corollary 3 against the upper bound, hence against tau
assert all(r["lo"] == r["hi"] for r in T78 if r["verdict"] in ("fewer", "same"))
assert all(r["lo"] < r["hi"] for r in T78 if r["verdict"] in ("fewer_than_built", "open"))

FIG_H = 4.45                                               # one column (246 pt); 0.15 in taller than before k = 7, 8 (two rows in (b))
fig = plt.figure(figsize=(246 / 72.27, FIG_H), constrained_layout=True, facecolor=SURF)
gs = fig.add_gridspec(2, 1, height_ratios=[1.75, 0.72])     # the key sits inside (a), in its empty lower right
ax = fig.add_subplot(gs[0]); bx = fig.add_subplot(gs[1])
for a in (ax, bx):
    a.set_facecolor(SURF)
    for sp in ("top", "right"): a.spines[sp].set_visible(False)
    for sp in ("left", "bottom"): a.spines[sp].set_color(MUTED)
    a.tick_params(colors=MUTED, labelsize=FS)

# (a) x axis linear to XB and compressed XC times beyond it (break mark on the axis): every class with k <= 6 and
# 4,290 of the 4,354 targets at k = 7, 8 lie in the linear part, which keeps room to set the glyphs of one position apart
XB, XC = 34.0, 3.0
fwd = lambda x: np.where(x <= XB, x, XB + (x - XB) / XC)
inv = lambda x: np.where(x <= XB, x, XB + (x - XB) * XC)
ax.set_xscale("function", functions=(fwd, inv))
XL = 2.5                                                   # left edge of (a): room for the labels in the excluded triangle
# Encoding, the same in (a) and (b): orange fewer, gray same, dark more (triangle in (a)); filled / solid when the
# T-count is decided (lo = hi), hollow / dotted when the target is compared only at its upper bound hi (lo < hi).
# Per position (n, T-count): one glyph per decided verdict, then one per upper-bound verdict, side by side and centred
# on n; a glyph stands for every class / target with that (n, T-count, verdict).
GAP = 0.70                                                 # in units of n, centre to centre within one position
G = collections.defaultdict(list)                         # (n, T) -> glyphs
for r in R:
    if r["n"] != r["tcount"] + 1:
        G[(r["n"], r["tcount"])].append(("fewer" if r["n"] < r["syn_strongest_n"] else "same", True))
CCZ_UB = [r for r in T78 if r["lo"] < r["hi"] and r["n"] == r["hi"] + 1]   # undecided CCZ layers next to the line
assert all(r["ccz_layer"] and r["lo"] == 17 and r["k"] == 8 for r in CCZ_UB) and len(CCZ_UB) == 3271
for r in T78:
    if r["on_bound"] or r in CCZ_UB: continue
    G[(r["n"], r["hi"])].append((r["verdict"], r["lo"] == r["hi"]))
# glyph styles: (scatter kwargs, outer diameter in pt)
STY = {("same", True): (dict(s=6, c=GRAY, alpha=0.9, linewidths=0, zorder=3), 2.45),
       ("fewer", True): (dict(s=6, c=ORANGE, linewidths=0, zorder=4), 2.45),
       ("more", True): (dict(s=9, marker="^", c=MORE, linewidths=0, zorder=4), 3.0),
       ("more", False): (dict(s=10, marker="^", facecolors=SURF, edgecolors=MORE, linewidths=0.7, zorder=6), 3.9),
       ("fewer_than_built", False): (dict(s=14, facecolors=SURF, edgecolors=ORANGE, linewidths=1.0, zorder=6), 4.75),
       ("open", False): (dict(s=10, facecolors=SURF, edgecolors=OPEN, linewidths=0.8, zorder=6), 4.0)}
VORD = ["same", "fewer", "more"]
GLYPHS = []                                                # (T, ideal x, verdict, decided) in drawing order along each row
for (n, t), gl in sorted(G.items()):
    vs = sorted({g[0] for g in gl if g[1]}, key=VORD.index)  # one glyph per decided verdict, then per upper-bound verdict
    items = [(v, True) for v in vs] + [(v, False) for v in sorted({g[0] for g in gl if not g[1]})]
    for i, it in enumerate(items):
        GLYPHS.append((t, n + GAP * (i - (len(items) - 1) / 2), it, n))
# the 3,271 CCZ layers at k = 8 whose factories meet the bound only if tau = hi: hollow at hi, a bar down to lo = 17
DXC = 1.35                                                 # set right of the line by DXC so they are not read as on it
for n in sorted({r["n"] for r in CCZ_UB}):
    x = n + DXC; hi = n - 1                                # every such n is <= 26, in the linear part
    ax.plot([x, x], [17, hi], color=OPEN, lw=0.8, zorder=5, solid_capstyle="butt")
    ax.plot([x - 0.1, x + 0.1], [17, 17], color=OPEN, lw=0.8, zorder=5)
    ax.scatter(x, hi, s=13, facecolors=SURF, edgecolors=OPEN, linewidths=0.8, zorder=6)
ax.text(27.9, 25.0, "{:,}".format(len(CCZ_UB)).replace(",", "{,}").join(["$", r"$ $\mathrm{CCZ}$ layers"]), ha="left", va="center", color=INK)
YTOP = 28.3
xs_ = [XL, 30]                                            # two points: a straight segment only in the linear part
ax.plot(xs_, [x - 1 for x in xs_], ls="--", lw=1.0, c=MUTED, zorder=1)
EXCL = "#ecebe7"                                          # neutral: blue marks only the factories on the bound
ax.fill_between([XL, YTOP + 1], [XL - 1, YTOP], YTOP, color=EXCL, lw=0, zorder=0)   # excluded by Corollary 3
# on the bound with tau decided: k <= 6 classes and k = 7, 8 targets; N = k + 1 on every one (asserted)
assert all(r["N"] == r["k"] + 1 for r in T78 if r["on_bound"]) and all(k + 1 == N for (_, k, N, _, _) in FRONT.values())
assert all(r["verdict"] == "same" for r in T78 if r["on_bound"])
labels = {}
for n, (lab, k, N, c, known) in FRONT.items():
    labels.setdefault(n, set()).add(k)
for r in T78:
    if r["on_bound"]:
        labels.setdefault(r["n"], set()).add(r["k"])
for n in labels:
    if n in FRONT and FRONT[n][4]:
        ax.scatter(n, n - 1, s=40, facecolors=SURF, edgecolors=BLUE, linewidths=1.4, zorder=7)
    else:
        ax.scatter(n, n - 1, s=44, marker="D", c=BLUE, edgecolors=SURF, linewidths=0.8, zorder=7)
for n, ks in labels.items():                                   # one 8 pt label per point: every k reaching it
    t = n - 1; lab = rf"$[[{n},{'/'.join(map(str, sorted(ks)))},2]]$"
    if n == 8:                                                 # right of and below the point
        ax.text(n + 0.8, t - 0.3, lab, ha="left", va="top", color=INK)
    else:                                                      # above-left, in the excluded triangle, stepping with the line
        ax.text(n - 0.55, t + 0.45, lab, ha="right", va="bottom", color=INK)
XMAX = 50.5
ax.set_xlim(XL, XMAX); ax.set_ylim(-3.6, YTOP)
ax.set_xlabel(r"$n$, input $T$ states", color=INK, labelpad=1); ax.set_ylabel(r"output $T$-count", color=INK)
ax.set_xticks([8, 16, 24, 32, 40, 48]); ax.set_yticks(range(0, 26, 5)); ax.grid(True, color=GRID, lw=0.4, zorder=0)
d_ = 0.012                                                 # break mark on the x axis at XB
xb = (float(fwd(XB)) - float(fwd(XL))) / (float(fwd(XMAX)) - float(fwd(XL)))
for dx in (-0.008, 0.008):
    ax.plot([xb + dx - d_ / 2, xb + dx + d_ / 2], [-0.018, 0.018], transform=ax.transAxes, color=MUTED, lw=0.8, clip_on=False, zorder=8)
ax.plot([xb - 0.006, xb + 0.006], [0, 0], transform=ax.transAxes, color=SURF, lw=1.6, clip_on=False, zorder=7.5)
h = [Line2D([0], [0], marker="o", ls="", mfc=SURF, mec=BLUE, mew=1.3, ms=5, label="synthillation paper"),
     Line2D([0], [0], marker="D", ls="", mfc=BLUE, mec=SURF, ms=5, label="not given there")]
lg = ax.legend(handles=h, loc="upper left", bbox_to_anchor=(-0.01, 1.02), ncol=1, frameon=False, handlelength=0.8,
               handletextpad=0.4, labelspacing=0.2, borderaxespad=0.1, labelcolor=INK)
ax.add_artist(lg)
ax.annotate("(a)", (0.005, 1.0), xycoords=("figure fraction", "axes fraction"), color=INK, va="top")

# (b) per row the share by verdict: k <= 6 classes, k = 7, 8 targets; l = 4 classes with a factory
def row_k6(lab, S):
    w = sum(r["n"] < r["syn_strongest_n"] for r in S)
    return (lab, len(S), dict(same=len(S) - w, fewer=w))
S4 = [r for r in R if r["k"] <= 4]
rows = [row_k6(r"$l=3$, $k\le4$", S4), row_k6(r"$l=3$, $k=5,6$", [r for r in R if r["k"] in (5, 6)])]   # rows merged (user)
for K in (7, 8):                                           # one row per k: k = 7 has every T-count exact, k = 8 mostly an upper bound
    TK = [r for r in T78 if r["k"] == K]
    rows.append((rf"$l=3$, $k={K}$", len(TK),
                 dict(collections.Counter(r["verdict"] + ("" if r["verdict"] != "more" or r["lo"] == r["hi"] else "_ub") for r in TK))))
rows.append((r"$l=4$, $k\le4$", sum(r["l"] == 4 and not r["third_level_target"] for r in SB), None))   # the 22 third-level targets left out
ORDER = [("same", GRAYBAR, None), ("open", GRAYBAR, "......"), ("more", MORE, None), ("more_ub", MORE, "......"),
         ("fewer_than_built", ORANGE, "......"), ("fewer", ORANGE, None)]
WMIN = 2.0                                                 # segments below WMIN % drawn at WMIN %, the excess taken from the widest
SHORT = dict(same="same", open="same, bound", more="more", more_ub="more, bound", fewer_than_built="fewer, bound", fewer="fewer")
PCT_PER_CHAR = 3.6                                         # width of one 7-8 pt digit in panel (b), in percent of the axis (about 3.4 in wide)
def textw_pct(t):
    return PCT_PER_CHAR * len(t.replace(",", "")) + 1.0
y = [5.3, 4.4, 3.1, 2.2, 0.9]                                # gaps: classes (k <= 6) | targets (k = 7, 8) | l = 4
H = 0.7
plt.rcParams["hatch.linewidth"] = 0.8
for yi, (lab, t, parts) in zip(y, rows):
    side = []
    if parts is None:
        bx.barh(yi, 100, color="#e7e6e2", height=H, edgecolor=SURF, linewidth=1.5, hatch="////")
    else:
        assert sum(parts.values()) == t, lab
        w = {v: 100 * c / t for v, c in parts.items() if c}
        small = [v for v in w if w[v] < WMIN]; big = max(w, key=w.get)
        for v in small: w[big] -= WMIN - w[v]; w[v] = WMIN
        left = 0.0
        for v, col, hat in ORDER:
            if v not in w: continue
            wdt = w[v]
            bx.barh(yi, wdt, left=left, color=col, height=H, edgecolor=SURF if hat is None else "white", linewidth=0.0 if hat else 1.0, hatch=hat)
            lab_ = rf"{parts[v]:,}"                         # every count inside its own segment; dark ink on light fills
            if wdt >= 0.75 * textw_pct(lab_):
                bx.text(left + wdt / 2, yi, lab_, va="center_baseline", ha="center",
                        color="white" if col in (ORANGE, MORE) else INK, fontsize=FS - 1 if wdt < 12 else FS, zorder=6,
                        bbox=None if hat is None else dict(fc=col, ec="none", pad=0.4))
            elif v in ("fewer", "fewer_than_built"):       # too narrow: only the fewer counts are named beside the bar
                side.append(rf"{parts[v]:,}" + ("" if v == "fewer" else r"$^{\ast}$"))
            left += wdt
        assert abs(left - 100) < 1e-9, lab
    bx.text(101.5, yi, " ".join(side + [f"({t:,})"]), va="center", ha="left", color=INK)
bx.set_yticks(y); bx.set_yticklabels([r[0] for r in rows], color=INK)
bx.set_xlim(0, 121); bx.set_ylim(0.35, 5.85)
bx.set_xticks([0, 50, 100]); bx.set_xticklabels([r"0", r"50", r"100"])
bx.set_xlabel(r"share (\%), with the total", color=INK, labelpad=1)
bx.annotate("(b)", (0.005, 1.0), xycoords=("figure fraction", "axes fraction"), color=INK, va="top")
# the key, inside (a) lower right: each entry the marker of (a) and the bar of (b); one column per kind of comparison
from matplotlib.patches import Patch
from matplotlib.legend_handler import HandlerTuple
hol = lambda c, m="o", ms=3.6: Line2D([0], [0], marker=m, ls="", mfc=SURF, mec=c, mew=0.9, ms=ms)
fil = lambda c, m="o", ms=3.4: Line2D([0], [0], marker=m, ls="", mfc=c, mec=c, mew=0, ms=ms)
dot = lambda c: Patch(fc=c, ec="white", lw=0, hatch="......")
HT = {tuple: HandlerTuple(ndivide=None, pad=0.3)}
KW = dict(frameon=False, handler_map=HT, handlelength=2.0, handletextpad=0.35, labelspacing=0.12, borderaxespad=0.0,
          labelcolor=INK, title_fontsize=FS, alignment="left")
k1 = ax.legend([(fil(ORANGE), Patch(fc=ORANGE)), (fil(GRAY), Patch(fc=GRAYBAR)), (fil(MORE, "^", 3.8), Patch(fc=MORE))],
               ["fewer", "same", "more"], loc="lower left", bbox_to_anchor=(0.36, 0.0), title=r"exact", **KW)
ax.add_artist(k1)
k2 = ax.legend([(hol(ORANGE, ms=4.0), dot(ORANGE)), (hol(OPEN), dot(GRAYBAR)), (hol(MORE, "^", 3.8), dot(MORE))],
               ["fewer", "same", "more"], loc="lower right", bbox_to_anchor=(1.0, 0.0), title=r"upper bound", **KW)

# place the glyphs: along each T-count row, the least-squares shift of the ideal x (group centred on n) such that
# neighbours are at least PAD pt apart edge to edge (pool-adjacent-violators on blocks); reports the largest shift
PAD = 0.3
fig.canvas.draw(); fig.set_layout_engine("none")   # freeze the layout so that the pt scale holds when saving
upt = float(fwd(ax.get_xlim()[1]) - fwd(ax.get_xlim()[0])) / (ax.get_window_extent().width * 72 / fig.dpi)   # display units per pt
SHIFT = 0.0; SHIFTS = []
for t in sorted({g[0] for g in GLYPHS}):
    row = sorted([(g[0], float(fwd(g[3])) + (g[1] - g[3]), g[2], float(fwd(g[3]))) for g in GLYPHS if g[0] == t], key=lambda g: g[1])
    rad = [STY[g[2]][1] / 2 * upt for g in row]
    off = [0.0]                                            # offsets of each glyph inside its block
    for k in range(1, len(row)): off.append(off[-1] + rad[k - 1] + rad[k] + PAD * upt)
    blocks = []                                            # [first, last, position of the block origin]
    for k in range(len(row)):
        blocks.append([k, k, row[k][1] - off[k] + off[k]])  # a block of one: origin so that glyph k sits at its ideal x
        blocks[-1][2] = row[k][1]
        def pos(b): return np.mean([row[q][1] - (off[q] - off[b[0]]) for q in range(b[0], b[1] + 1)])
        blocks[-1][2] = pos(blocks[-1])
        while len(blocks) > 1:
            b0, b1 = blocks[-2], blocks[-1]
            end0 = b0[2] + off[b0[1]] - off[b0[0]]
            if b1[2] - end0 >= rad[b0[1]] + rad[b1[0]] + PAD * upt - 1e-12: break
            blocks[-2:] = [[b0[0], b1[1], 0.0]]; blocks[-1][2] = pos(blocks[-1])
    for b in blocks:
        for q in range(b[0], b[1] + 1):
            x = b[2] + off[q] - off[b[0]]; SHIFT = max(SHIFT, abs(x - row[q][1])); SHIFTS.append(abs(x - row[q][3]))
            ax.scatter(float(inv(x)), t, **STY[row[q][2]][0])
print(f"panel (a): {len(GLYPHS)} glyphs; largest distance of a glyph from its n: {max(SHIFTS):.2f}")

if ARGS.out:
    for ext in ("png", "pdf"):
        fig.savefig(f"{ARGS.out}.{ext}", dpi=300, facecolor=SURF); print("saved", f"{ARGS.out}.{ext}")
else:
    for out in (os.path.join(ROOT, "figures", "frontier_l3.png"),
                os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "frontier_l3.png")):
        if os.path.isdir(os.path.dirname(out)):
            fig.savefig(out, dpi=300, facecolor=SURF); print("saved", out)
