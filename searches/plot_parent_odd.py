"""Supplement figure (one column, Fig. 1 style): the odd-weight parent of the first submission's Fig. 6, redrawn.

(a) l=3, five qubits, odd weights, every rotation pi/8: [[8,4,2]] -> [[12,3,2]] -> [[14,2,2]] -> [[15,1,3]];
(b) the six-qubit circuit of the first submission (13 rotations), with its two cuts [[8,5,2]] and [[12,3,2]]; it is not
    a borrowed identity, so the gates right of a cut do not synthesise its output: the labels give each EXACT output;
(c) conversions on the wires q1..q5 of (b): CNOT(q1->q2) CNOT(q3->q2) turn the [[8,5,2]] output into |+>_1 and the
    [[8,4,2]] output of (a) on q2..q5; sqrt(X)^dag_2, C_2(X_3X_4X_5) and an X measurement of q2 (outcome + with
    probability 1) give the [[12,3,2]] output on q3q4q5; one catalytic T on the parity q3q4q5 and S on each wire give |T>^3.
Every statement is checked by searches/lostv1/fig6_check.py -> outputs/checks/fig6_check.json (exact state vectors and
unified_search.is_factory on the full phase table). Coordinates are in points; text at the caption size (8 pt).
    python searches/plot_parent_odd.py   -> figures/parent_odd.{pdf,png}, MSD_manuscript/FIGURES/parent_odd.pdf
"""
import os, itertools, matplotlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
FS = 8
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": FS,
                     "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}"})
INK, WIRE = "#0b0b0b", "#c9c8c3"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#8a5cd6"
CUTC = [BLUE, ORANGE, AQUA, VIOLET]
LAYER = ["#dbe8f9", "#fbe3d8", "#d7f1e6", "#ebe3f8", "#eeeeec"]
WIDTH, DYW, GAPF, LINE, ROWGAP = 246.0, 8.6, 0.8, 8.6, 2.5

fig = plt.figure(figsize=(WIDTH / 72.27, 3.0))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, WIDTH); ax.set_ylim(-216, 0)
rend = fig.canvas.get_renderer()
def textw(s, **kw):
    t = ax.text(0, 0, s, **kw); bb = t.get_window_extent(renderer=rend); t.remove(); return bb.width * 72.0 / fig.dpi
def texth(s, **kw):
    t = ax.text(0, 0, s, **kw); bb = t.get_window_extent(renderer=rend); t.remove(); return bb.height * 72.0 / fig.dpi
Lw = textw(r"$q_0{:}\,|{+}\rangle$") + 2.5
Rw = textw(r"$|{+}\rangle$") + 2.0

def parent(ytop, x0, N, layers, cuts, letter, xright, p, plus_right=True, lgap=6.0, cutcol=CUTC, labels=True):
    """layers: list of lists of supports; cuts: list of (after_layer_index, first output wire, factory, output)."""
    g = GAPF * p
    y = {q: ytop - 12.0 - DYW * q for q in range(N)}
    L0 = x0 + Lw
    starts, xs = [], []
    for j, Lj in enumerate(layers):
        st = L0 + p / 2 + p * sum(len(M) for M in layers[:j]) + g * j; starts.append(st)
        xs += [st + p * i for i in range(len(Lj))]
    cols = [S for Lj in layers for S in Lj]
    top, bot = y[0] + 3.5, y[N - 1] - 3.5
    for j, Lj in enumerate(layers):
        ax.add_patch(FancyBboxPatch((starts[j] - p / 2, bot), p * len(Lj), top - bot,
                                    boxstyle="round,pad=0,rounding_size=2", fc=LAYER[j], ec="none", zorder=0))
    xend = xs[-1] + p / 2 + 1.0
    for q in range(N):
        ax.plot([L0 - 1, xend], [y[q], y[q]], color=WIRE, lw=0.5, zorder=1)
        ax.text(L0 - 1.5, y[q], rf"$q_{q}{{:}}\,|{{+}}\rangle$", ha="right", va="center", color=INK)
        if plus_right:
            ax.text(xend + 0.8, y[q], r"$|{+}\rangle$", ha="left", va="center", color=INK)
    for x, S in zip(xs, cols):
        if len(S) > 1:
            ax.plot([x, x], [y[max(S)], y[min(S)]], color=INK, lw=0.6, zorder=3)
        ax.scatter([x] * len(S), [y[q] for q in S], s=4.5, c=INK, linewidths=0, zorder=4)
    xcs = []
    for (j, q1, lab, outl), col in zip(cuts, cutcol):
        n = sum(len(M) for M in layers[:j])
        assert lab.startswith(rf"$[[{n},{N - q1},"), (lab, n, N - q1)
        xc = starts[j] - p / 2 - g / 2; xcs.append(xc)
        ax.plot([xc, xc], [bot - 1.0, y[q1] + 3.5], ls=(0, (2.0, 1.2)), lw=0.9, color=col, zorder=4)
    w = [max(textw(a), textw(b)) for _, _, a, b in cuts]
    cx = list(xcs)
    for i in range(1, len(cx)):
        cx[i] = max(cx[i], cx[i - 1] + (w[i - 1] + w[i]) / 2 + lgap)
    shift = sum(xcs) / len(xcs) - sum(cx) / len(cx)
    shift = max(shift, x0 + 2.0 - (cx[0] - w[0] / 2))
    cx = [c + shift for c in cx]
    over = cx[-1] + w[-1] / 2 - xright
    if over > 0:
        cx[-1] -= over
        for i in range(len(cx) - 2, -1, -1):
            cx[i] = min(cx[i], cx[i + 1] - (w[i] + w[i + 1]) / 2 - lgap)
    ylab = bot - 3.0
    for (_, _, lab, outl), xc, x, col in zip(cuts, xcs, cx, cutcol) if labels else ():
        if abs(x - xc) > 1:
            ax.plot([xc, x], [bot - 1.0, ylab + 0.3], color=col, lw=0.5, zorder=4)
        ax.text(x, ylab, outl, ha="center", va="top", color=col)
        ax.text(x, ylab - LINE, lab, ha="center", va="top", color=INK)
    ax.text(x0, ytop, letter, ha="left", va="top", color=INK)
    return dict(y=y, xcs=xcs, ylab=ylab, bottom=ylab - 2 * LINE, xend=xend, right=xend + (Rw if plus_right else 0), bot=bot)

def pitch(layers, avail, right=True):
    return (avail - Lw - (Rw if right else 1.5)) / (sum(map(len, layers)) + GAPF * (len(layers) - 1))

# ---- (a) five-qubit odd-weight parent, full width
five = [S for w in (5, 3, 1) for S in itertools.combinations(range(5), w)]
La = [[S for S in five if S[0] == j] for j in range(5)]
pa = pitch(La, WIDTH)
A = parent(0.0, 0.0, 5, La, [(1, 1, r"$[[8,4,2]]$", r"$\sim\!|\mathrm{CCZ}\rangle$"),
                             (2, 2, r"$[[12,3,2]]$", r"$\sim\!|\mathrm{CS}_{12}\mathrm{CS}_{13}\rangle$"),
                             (3, 3, r"$[[14,2,2]]$", r"$|T\rangle^{\otimes2}$"), (4, 4, r"$[[15,1,3]]$", r"$|T\rangle$")],
           "(a)", WIDTH, pa, lgap=7.0)

# ---- (b) six-qubit circuit of the first submission, full width; labels give the EXACT output of each cut
six = [(0, 1, 3, 4, 5), (0, 2, 3, 4, 5), (0, 1, 5), (0, 2, 5), (0, 1, 4), (0, 2, 4), (0, 1, 3), (0, 2, 3),
       (1, 2, 3, 4, 5), (1, 2, 5), (1, 2, 4), (1, 2, 3), (3, 4, 5)]
Lb = [six[:8], six[8:12], six[12:]]
pb = pa                                                    # same pitch as (a); the labels use the space on the right
ytop2 = A["bottom"] - ROWGAP
B = parent(ytop2, 0.0, 6, Lb, [(1, 1, r"$[[8,5,2]]$", ""), (2, 3, r"$[[12,3,2]]$", "")],
           "(b)", WIDTH, pb, plus_right=False, lgap=6.0, cutcol=[BLUE, ORANGE], labels=False)
# exact outputs (fig6_check.json exact_outputs_and_chain): first label from the wire start, second flush right
lb1 = r"$T_{12}T_{34}T_{35}T_{45}T_{1234}T_{1235}T_{1245}|{+}\rangle^{\otimes5}$"
lb2 = r"$T^{\dagger}_{3}T^{\dagger}_{4}T^{\dagger}_{5}T^{\dagger}_{345}|{+}\rangle^{\otimes3}$"
yb = B["bot"] - 5.0
xl1 = Lw - 3.0; wl1 = textw(lb1); wl2 = textw(lb2); xr2 = WIDTH - 0.5
assert xl1 + wl1 + 4.0 < xr2 - wl2, "(b) labels overlap"
ax.text(xl1, yb, lb1, ha="left", va="top", color=BLUE)
ax.text(xr2, yb, lb2, ha="right", va="top", color=ORANGE)
ax.text(xl1 + wl1 / 2, yb - LINE - 2.5, r"$[[8,5,2]]$", ha="center", va="top", color=INK)
ax.text(xr2 - wl2 / 2, yb - LINE - 2.5, r"$[[12,3,2]]$", ha="center", va="top", color=INK)
for xc, xt_, col in ((B["xcs"][0], xl1 + wl1 / 2, BLUE), (B["xcs"][1], xr2 - wl2 / 2, ORANGE)):
    ax.plot([xc, xt_], [B["bot"] - 1.0, yb + 0.2], color=col, lw=0.5, zorder=4)
B["bottom"] = yb - 2 * LINE - 2.5

# ---- (c) conversions, on the wires q1..q5 of (b)
ytop3 = B["bottom"] - ROWGAP
ry = {q: ytop3 - 12.0 - DYW * (q - 1) for q in range(1, 6)}
ax.text(0, ytop3, "(c)", ha="left", va="top", color=INK)
qlw = textw(r"$q_5$") + 2.0
x0 = qlw + 6.0
for q, yy in ry.items():
    ax.text(x0 - 1.5, yy, rf"$q_{q}$", ha="right", va="center", color=INK)
def box(x, w, yy, label, fc="white", ec=INK):
    w = max(w, textw(label, fontsize=6) + 2.4); h = max(7.2, texth(label, fontsize=6) + 2.0)
    ax.add_patch(Rectangle((x - w / 2, yy - h / 2), w, h, fc=fc, ec=ec, lw=0.5, zorder=3))
    ax.text(x, yy, label, ha="center", va="center", color=INK, zorder=4, fontsize=6)
    return x - w / 2, x + w / 2
def cnot_draw(x, c, ts, col=INK):
    ax.plot([x, x], [min(ry[t] for t in ts + [c]), max(ry[t] for t in ts + [c])], color=col, lw=0.5, zorder=2)
    ax.scatter([x], [ry[c]], s=4.5, c=col, zorder=4, linewidths=0)
    for t in ts:
        ax.add_patch(Circle((x, ry[t]), 2.1, fc="white", ec=col, lw=0.5, zorder=3))
        ax.plot([x - 2.1, x + 2.1], [ry[t], ry[t]], color=col, lw=0.5, zorder=4)
        ax.plot([x, x], [ry[t] - 2.1, ry[t] + 2.1], color=col, lw=0.5, zorder=4)
def marker(x, q1, col):
    ax.plot([x, x], [ry[5] - 3.5, ry[q1] + 3.5], ls=(0, (2.0, 1.2)), lw=0.9, color=col, zorder=4)
PW = textw(r"$|{+}\rangle$"); TW = textw(r"$|T\rangle^{\otimes3}$")
w0 = textw(r"$\sqrt{X}^{\dagger}$", fontsize=6) + 2.4
# x positions as offsets in units of sp (gate spacing) plus fixed widths; sp chosen to fill the column
def layout(gp):                                            # one free gap gp between consecutive elements, edges included
    X = {}; x = x0 + 3.0
    X["m1"] = x; x += gp + 2.1; X["g1"] = x; x += 4.2 + gp; X["g2"] = x; x += 2.1 + gp / 2; X["e1"] = x; x += PW + gp
    X["m2"] = x; x += gp + w0 / 2; X["g3"] = x; x += w0 / 2 + gp + 2.1; X["g4"] = x; x += 2.1 + gp / 2; X["e2"] = x
    x += PW + gp; X["m3"] = x; x += gp; X["t"] = x; x += gp + 3.0; X["s"] = x; x += 3.0 + gp / 2; X["e"] = x
    return X
lo, hi = 4.0, 60.0
for _ in range(50):
    mid = (lo + hi) / 2
    if layout(mid)["e"] + 2.3 + TW <= WIDTH - 1.0: lo = mid
    else: hi = mid
X = layout(lo)
xm1, xg1, xg2, xe1, xm2, xg3, xg4, xe2, xm3, xt, xs_, xe = (X[k] for k in ("m1", "g1", "g2", "e1", "m2", "g3", "g4", "e2", "m3", "t", "s", "e"))
for q, yy in ry.items():
    xz = xe1 if q == 1 else xe2 if q == 2 else xe
    ax.plot([x0, xz], [yy, yy], color=WIRE, lw=0.5, zorder=1)
for q, xz in ((1, xe1), (2, xe2)):
    ax.text(xz + 0.8, ry[q], r"$|{+}\rangle$", ha="left", va="center", color=INK)
cnot_draw(xg1, 1, [2]); cnot_draw(xg2, 3, [2])
box(xg3, 7.5, ry[2], r"$\sqrt{X}^{\dagger}$")
cnot_draw(xg4, 2, [3, 4, 5])
ax.plot([xt, xt], [ry[5], ry[3]], color=ORANGE, lw=0.8, zorder=3)
ax.scatter([xt] * 3, [ry[3], ry[4], ry[5]], s=6.5, c=ORANGE, linewidths=0, zorder=4)
ax.text(xt, ry[2], r"catalytic $T$", ha="center", va="center", color=ORANGE)
for t in (3, 4, 5):
    box(xs_, 6.0, ry[t], r"$S$")
ax.plot([xe + 0.5, xe + 1.5, xe + 1.5, xe + 0.5], [ry[3] + 2, ry[3] + 2, ry[5] - 2, ry[5] - 2], color=INK, lw=0.5)
ax.text(xe + 2.3, ry[4], r"$|T\rangle^{\otimes3}$", ha="left", va="center", color=INK)
marker(xm1, 1, BLUE); marker(xm2, 2, BLUE); marker(xm3, 3, ORANGE)
ylab3 = ry[5] - 4.0
labs = [(xm1, BLUE, r"(b) $[[8,5,2]]$"), (xm2, BLUE, r"(a) $[[8,4,2]]$"), (xm3, ORANGE, r"(a,b) $[[12,3,2]]$")]
for x, col, t in labs:
    ax.text(x - 1.5, ylab3, t, ha="left", va="top", color=col)
wl = [textw(t) for _, _, t in labs]
assert all(labs[i][0] + wl[i] + 3 < labs[i + 1][0] for i in range(2)) and labs[2][0] + wl[2] < WIDTH, "labels of (c) overlap"
assert xt - textw(r"catalytic $T$") / 2 > xe2 + PW + 1.0, "catalytic label hits |+> of q2"
assert xe + 2.3 + TW <= WIDTH, "(c) too wide"
C = dict(bottom=ylab3 - LINE - 2.0)
ytop, ybot = 0.5, C["bottom"] + 1.5
H = ytop - ybot
ax.set_xlim(0, WIDTH); ax.set_ylim(ybot, ytop)
fig.set_size_inches(WIDTH / 72.27, H / 72.27)
print(f"figure {WIDTH:.0f} x {H:.0f} pt; pitch (a) {pa:.2f}, (b) {pb:.2f} pt; (c) ends at {xe + 2.3 + TW:.0f} pt")
for out in (os.path.join(ROOT, "figures", "parent_odd.pdf"),
            os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "parent_odd.pdf"),
            os.path.join(ROOT, "figures", "parent_odd.png")):
    if os.path.isdir(os.path.dirname(out)):
        fig.savefig(out, dpi=300); print("saved", out)
