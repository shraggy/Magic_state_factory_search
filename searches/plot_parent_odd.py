"""Supplement figure (one column, Fig. 1 style): the odd-weight parent of the first submission's Fig. 6, redrawn.

(a) l=3, five qubits, odd weights, every rotation pi/8: [[8,4,2]] -> [[12,3,2]] -> [[14,2,2]] -> [[15,1,3]];
(b) conversions of its outputs: the [[8,4,2]] output on q1..q4 becomes the [[12,3,2]] output on q2q3q4 by
    sqrt(X)^dag_1, C_1(X_2X_3X_4) and an X measurement of q1 whose outcome is + with probability 1; one further T
    rotation on the parity q2q3q4 (the catalyst) and S on each wire then leave |T>^{x3};
(c) the six-qubit circuit of the first submission (13 rotations), drawn with its two cuts [[8,5,2]] and [[12,3,2]] only:
    it is not itself the identity, so its wires carry no |+> on the right.
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

def parent(ytop, x0, N, layers, cuts, letter, xright, p, plus_right=True, lgap=6.0, cutcol=CUTC):
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
    for (_, _, lab, outl), xc, x, col in zip(cuts, xcs, cx, cutcol):
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

# ---- second row: (b) conversion on the left, (c) six-qubit circuit on the right
six = [(0, 1, 3, 4, 5), (0, 2, 3, 4, 5), (0, 1, 5), (0, 2, 5), (0, 1, 4), (0, 2, 4), (0, 1, 3), (0, 2, 3),
       (1, 2, 3, 4, 5), (1, 2, 5), (1, 2, 4), (1, 2, 3), (3, 4, 5)]
Lc = [six[:8], six[8:12], six[12:]]
BW = 100.0                                                 # width given to (b); (c) takes the rest
SEP = 4.0
xc0 = BW + SEP
pc = pitch(Lc, WIDTH - xc0, right=False)
ytop2 = A["bottom"] - ROWGAP
C = parent(ytop2, xc0, 6, Lc, [(1, 1, r"$[[8,5,2]]$", r"$\sim\!|\mathrm{CCZ}\rangle$"),
                               (2, 3, r"$[[12,3,2]]$", r"$\sim\!|\mathrm{CS}_{12}\mathrm{CS}_{13}\rangle$")],
           "(c)", WIDTH, pc, plus_right=False, lgap=4.0, cutcol=[BLUE, ORANGE])

ry = {q: C["y"][q] for q in (1, 2, 3, 4)}
ax.text(0, ytop2, "(b)", ha="left", va="top", color=INK)
qlw = textw(r"$q_4$") + 2.0
x0 = qlw + 1.0
for q, yy in ry.items():
    ax.text(x0 - 1.5, yy, rf"$q_{q}$", ha="right", va="center", color=INK)
TW = textw(r"$|T\rangle^{\otimes3}$")
xe = BW - TW - 3.0
w0 = textw(r"$\sqrt{X}^{\dagger}$", fontsize=6) + 2.4; w2 = textw(r"$\langle{+}|$", fontsize=6) + 2.4
GX = [x0 + 1.0 + w0 / 2]; GX.append(GX[0] + w0 / 2 + 4.0); GX.append(GX[1] + 4.0 + w2 / 2)   # sqrtX^dag, fan-out, X meas.
xm = GX[2] + w2 / 2 + 5.0                                           # the [[12,3,2]] output: dashed marker
xs_ = xe - 5.5; xt = (xm + xs_ - 3.0) / 2                               # catalytic T parity rotation, S layer
def box(x, w, yy, label, fc="white", ec=INK):
    w = max(w, textw(label, fontsize=6) + 2.4); h = max(7.2, texth(label, fontsize=6) + 2.0)
    ax.add_patch(Rectangle((x - w / 2, yy - h / 2), w, h, fc=fc, ec=ec, lw=0.5, zorder=3))
    ax.text(x, yy, label, ha="center", va="center", color=INK, zorder=4, fontsize=6)
    return x - w / 2, x + w / 2
for q, yy in ry.items():
    xz = GX[2] if q == 1 else xe
    ax.plot([x0, xz], [yy, yy], color=WIRE, lw=0.5, zorder=1)
l0, r0 = box(GX[0], 7.5, ry[1], r"$\sqrt{X}^{\dagger}$")
ax.plot([GX[1], GX[1]], [ry[4], ry[1]], color=INK, lw=0.5, zorder=2)
ax.scatter([GX[1]], [ry[1]], s=4.5, c=INK, zorder=4, linewidths=0)
for t in (2, 3, 4):
    ax.add_patch(Circle((GX[1], ry[t]), 2.1, fc="white", ec=INK, lw=0.5, zorder=3))
    ax.plot([GX[1] - 2.1, GX[1] + 2.1], [ry[t], ry[t]], color=INK, lw=0.5, zorder=4)
    ax.plot([GX[1], GX[1]], [ry[t] - 2.1, ry[t] + 2.1], color=INK, lw=0.5, zorder=4)
l2, r2 = box(GX[2], 6.5, ry[1], r"$\langle{+}|$")
assert l2 - (GX[1] + 2.1) > 1.0 and GX[1] - 2.1 - r0 > 1.0, "boxes touch the fan-out"
ax.plot([xm, xm], [ry[4] - 3.5, ry[2] + 3.5], ls=(0, (2.0, 1.2)), lw=0.9, color=ORANGE, zorder=4)
ax.plot([xt, xt], [ry[4], ry[2]], color=ORANGE, lw=0.8, zorder=3)
ax.scatter([xt] * 3, [ry[2], ry[3], ry[4]], s=6.5, c=ORANGE, linewidths=0, zorder=4)
for t in (2, 3, 4):
    lS, _ = box(xs_, 6.0, ry[t], r"$S$")
assert lS - xt > 2.0 and xt - xm > 2.0, "catalytic T touches its neighbours"
ax.plot([xe + 0.5, xe + 1.5, xe + 1.5, xe + 0.5], [ry[2] + 2, ry[2] + 2, ry[4] - 2, ry[4] - 2], color=INK, lw=0.5)
ax.text(xe + 2.3, ry[3], r"$|T\rangle^{\otimes3}$", ha="left", va="center", color=INK)
# labels: the input above its wires (row of q0 in (c)), the [[12,3,2]] output and the catalyst under the wires
ax.text(x0, C["y"][0], r"$\sim\!|\mathrm{CCZ}\rangle$ of $[[8,4,2]]$", ha="left", va="center", color=BLUE)
xl = (xm + xt) / 2
ax.plot([xm, xl], [ry[4] - 3.5, C["ylab"] + 0.3], color=ORANGE, lw=0.5, zorder=4)
ax.text(xl, C["ylab"], r"$\sim\!|\mathrm{CS}_{12}\mathrm{CS}_{13}\rangle$", ha="center", va="top", color=ORANGE)
ax.text(xl, C["ylab"] - LINE, r"$+$ one catalytic $T$", ha="center", va="top", color=INK)
assert xe + 2.3 + TW <= BW + 0.5 and BW < C["xcs"][0], "(b) runs into (c)"
ytop, ybot = 0.5, C["bottom"] + 1.5
H = ytop - ybot
ax.set_xlim(0, WIDTH); ax.set_ylim(ybot, ytop)
fig.set_size_inches(WIDTH / 72.27, H / 72.27)
print(f"figure {WIDTH:.0f} x {H:.0f} pt; pitch (a) {pa:.2f}, (c) {pc:.2f} pt; (c) ends at {C['xend']:.0f} pt")
for out in (os.path.join(ROOT, "figures", "parent_odd.pdf"),
            os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "parent_odd.pdf"),
            os.path.join(ROOT, "figures", "parent_odd.png")):
    if os.path.isdir(os.path.dirname(out)):
        fig.savefig(out, dpi=300); print("saved", out)
