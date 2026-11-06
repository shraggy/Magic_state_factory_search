"""Fig. 1 (conceptual, one column): three symmetric borrowed identities drawn as support matrices and cut into
malleable factories, at one rotation pitch.

(a) l=3, four qubits, every weight: [[8,3,2]] -> [[12,2,2]] -> [[14,1,2]], with the [[10,2,2]] branch beside it, drawn as
    a small circuit on the |CCZ> output of [[8,3,2]];
(b) l=3, five qubits, odd weights: [[8,4,2]] -> [[12,3,2]] -> [[14,2,2]] -> [[15,1,3]];
(c) l=4, five qubits, every weight: [[16,4,2]] -> [[24,3,2]] -> [[28,2,2]] -> [[30,1,2]].

Each column is one pi/2^l parity rotation; dots mark the qubits it acts on. Layer j holds the rotations whose first qubit
is q_j, drawn widest first (they commute). Every wire starts in |+> and returns to |+>, asserted here from the symmetric
condition (Theorem 4). A cut after layer j measures q_0..q_{j-1}: the columns to its left are the inputs, the columns to
its right a synthesis of the output. Every cut is checked as a factory by searches/parent_cuts_check.py. The [[10,2,2]]
branch follows the gate list of App. S2.2: on the |CCZ> output (q1 q2 q3), sqrt(X)_1, C_1(X_2 X_3), the X-controlled
Z_2 Z_3 around T^dagger_1 (open circle: control on the X basis of q1), then T_1 S^dagger_1, and q1 is
accepted in |+>. Coordinates are in points, so the text is at the
caption size (8 pt) at print size.
    python searches/plot_parent.py     -> figures/parent_l3.{pdf,png}, MSD_manuscript/FIGURES/parent_l3.pdf
"""
import os, itertools, math, matplotlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle
FS = 8                                                     # caption size of the paper (revtex \small, 8pt)
plt.rcParams.update({"text.usetex": True, "font.family": "serif", "font.size": FS,
                     "text.latex.preamble": r"\usepackage{amsmath}\usepackage{amssymb}"})
INK, MUTED, WIRE = "#0b0b0b", "#52514e", "#c9c8c3"
BLUE, ORANGE, AQUA, VIOLET = "#2a78d6", "#eb6834", "#1baf7a", "#8a5cd6"
CUTC = [BLUE, ORANGE, AQUA, VIOLET]
LAYER = ["#dbe8f9", "#fbe3d8", "#d7f1e6", "#ebe3f8", "#eeeeec"]
WIDTH = 246.0                                              # column width of the paper, pt
DYW = 8.6                                                  # wire spacing, pt
GAPF = 0.8                                                 # layer gap in units of the rotation pitch
LINE = 8.6                                                 # text line height, pt
ROWGAP = 2.5                                               # between rows, pt (the 2.5 pt drop of each panel's wires is taken from here)

fig = plt.figure(figsize=(WIDTH / 72.27, 3.0))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, WIDTH); ax.set_ylim(-216, 0)                # 1 data unit = 1 pt (reset at the end)
rend = fig.canvas.get_renderer()
def textw(s, **kw):
    t = ax.text(0, 0, s, **kw); bb = t.get_window_extent(renderer=rend); t.remove()
    return bb.width * 72.0 / fig.dpi
def texth(s, **kw):
    t = ax.text(0, 0, s, **kw); bb = t.get_window_extent(renderer=rend); t.remove()
    return bb.height * 72.0 / fig.dpi

Lw = textw(r"$q_0{:}\,|{+}\rangle$") + 2.5                  # left wire labels
Rw = textw(r"$|{+}\rangle$") + 2.0                         # right wire labels

def layers_of(N, odd):
    weights = [w for w in range(1, N + 1) if not odd or w % 2]
    Ls = [[S for w in range(N, 0, -1) for S in itertools.combinations(range(N), w) if w in weights and S[0] == j]
          for j in range(N)]
    return weights, [x for x in Ls if x]

def symmetric_ok(N, l, weights):
    """Theorem 4: sum_{w>=i} C(N-i, w-i) in 2^(l-i+1) Z for i = 1..l."""
    return all(sum(math.comb(N - i, w - i) for w in weights if w >= i) % 2 ** (l - i + 1) == 0 for i in range(1, l + 1))

def span(layers, p):
    return Lw + Rw + p * (sum(map(len, layers)) + GAPF * (len(layers) - 1))

def parent(ytop, N, l, odd, p, cuts, letter, xright, lgap=2.5):
    """Draw one parent with its top at ytop (pt); cut labels stay left of xright. Return geometry."""
    weights, layers = layers_of(N, odd)
    assert symmetric_ok(N, l, weights)
    g = GAPF * p
    y = {q: ytop - 12.0 - DYW * q for q in range(N)}        # q0 sits below the panel letter
    xs, starts = [], []
    for j, Lj in enumerate(layers):
        start = Lw + p / 2 + p * sum(len(M) for M in layers[:j]) + g * j; starts.append(start)
        xs += [start + p * i for i in range(len(Lj))]
    cols = [S for Lj in layers for S in Lj]
    top, bot = y[0] + 3.5, y[N - 1] - 3.5
    for j, Lj in enumerate(layers):
        ax.add_patch(FancyBboxPatch((starts[j] - p / 2, bot), p * len(Lj), top - bot,
                                    boxstyle="round,pad=0,rounding_size=2", fc=LAYER[j], ec="none", zorder=0))
    xend = xs[-1] + p / 2 + 1.0
    for q in range(N):
        ax.plot([Lw - 1, xend], [y[q], y[q]], color=WIRE, lw=0.5, zorder=1)
        ax.text(Lw - 1.5, y[q], rf"$q_{q}{{:}}\,|{{+}}\rangle$", ha="right", va="center", color=INK)
        ax.text(xend + 0.8, y[q], r"$|{+}\rangle$", ha="left", va="center", color=INK)
    for x, S in zip(xs, cols):                             # one rotation: dots on its support, joined
        if len(S) > 1:
            ax.plot([x, x], [y[max(S)], y[min(S)]], color=INK, lw=0.6, zorder=3)
        ax.scatter([x] * len(S), [y[q] for q in S], s=4.5, c=INK, linewidths=0, zorder=4)
    xcs = []
    for j, (lab, out) in enumerate(cuts, start=1):
        n = sum(len(M) for M in layers[:j]); k = N - j
        assert lab.startswith(rf"$[[{n},{k},"), (lab, n, k)
        xc = starts[j] - p / 2 - g / 2; xcs.append(xc)
        # the cut after layer j leaves q_j..q_{N-1} as the outputs: the line spans only those wires
        ax.plot([xc, xc], [bot - 1.0, y[j] + 3.5], ls=(0, (2.0, 1.2)), lw=0.9, color=CUTC[j - 1], zorder=4)
    # labels under the cuts: centred on the cut, pushed apart, kept inside [Lw, xright]
    w = [max(textw(a), textw(b)) for a, b in cuts]
    cx = list(xcs)
    for i in range(1, len(cx)):
        cx[i] = max(cx[i], cx[i - 1] + (w[i - 1] + w[i]) / 2 + lgap)
    shift = sum(xcs) / len(xcs) - sum(cx) / len(cx)        # pull the packed row back left, centred on the cuts
    shift = max(shift, (cx[0] - w[0] / 2) and (2.0 - (cx[0] - w[0] / 2)))
    cx = [c + shift for c in cx]
    over = cx[-1] + w[-1] / 2 - xright
    if over > 0:
        cx[-1] -= over
        for i in range(len(cx) - 2, -1, -1):
            cx[i] = min(cx[i], cx[i + 1] - (w[i] + w[i + 1]) / 2 - lgap)
    ylab = bot - 3.0
    for (lab, out), xc, x, col in zip(cuts, xcs, cx, CUTC):
        if abs(x - xc) > 1:
            ax.plot([xc, x], [bot - 1.0, ylab + 0.3], color=col, lw=0.5, zorder=4)
        ax.text(x, ylab, out, ha="center", va="top", color=col)          # output state, under its cut
        ax.text(x, ylab - LINE, lab, ha="center", va="top", color=INK)   # factory: every block left of its cut
    ax.text(0, ytop, letter, ha="left", va="top", color=INK)
    return dict(y=y, xcs=xcs, ylab=ylab, bottom=ylab - 2 * LINE, xend=xend, right=xend + Rw)

# ---- the [[10,2,2]] branch, laid out first so (a) takes the rest of the row
BR_IN = textw(r"$|\mathrm{CCZ}\rangle$") + 5.0            # input label and brace
GX = [6.0, 16.0, 22.5, 30.0, 37.5, 44.5, 55.0]              # sqrtX, CX, CZ, Tdg, CZ, T, Sdg (offsets, pt)
BRL = 61.0                                                  # branch wire length, pt
BR_W = BR_IN + BRL + 1.0 + textw(r"$|T\rangle^{\otimes2}$")
SEP = 4.0

_, La = layers_of(4, False); _, Lb = layers_of(5, True); _, Lc = layers_of(5, False)
pa = (WIDTH - BR_W - SEP - Lw - Rw) / (sum(map(len, La)) + GAPF * (len(La) - 1))
pc = (WIDTH - Lw - Rw) / (sum(map(len, Lc)) + GAPF * (len(Lc) - 1))
P = pa                                                     # (a) fills the row beside the [[10,2,2]] branch
pb = (WIDTH - Lw - Rw) / (sum(map(len, Lb)) + GAPF * (len(Lb) - 1))   # (b) and (c) fill the column width
A = parent(0.0, 4, 3, False, P, [(r"$[[8,3,2]]$", r"$|\mathrm{CCZ}\rangle$"), (r"$[[12,2,2]]$", r"$|\mathrm{CS}\rangle$"),
                                 (r"$[[14,1,2]]$", r"$|T\rangle$")], "(a)",
           WIDTH - BR_W + BR_IN + BRL / 2 - textw(r"$8+2$ inputs") / 2 - 6.0, lgap=7.0)   # may run under the branch, left of its label

xb = WIDTH - BR_W                                          # branch on the right, its wires at the heights of q1..q3
ry = {1: A["y"][1], 2: A["y"][2], 3: A["y"][3]}
x0 = xb + BR_IN
ax.text(xb, ry[2], r"$|\mathrm{CCZ}\rangle$", ha="left", va="center", color=BLUE)
# arrow from the top of the [[8,3,2]] cut, over (a), down to the |CCZ> label of the branch
xc0, ytop_a = A["xcs"][0], A["y"][0] + 4.5
xl, yh = xb + textw(r"$|\mathrm{CCZ}\rangle$") / 2, ytop_a + 3.0
ax.plot([xc0, xc0, xl], [ytop_a, yh, yh], color=BLUE, lw=0.6, zorder=5, solid_joinstyle="round")
ax.annotate("", xy=(xl, ry[2] + 4.5), xytext=(xl, yh + 0.3),
            arrowprops=dict(arrowstyle="-|>", lw=0.6, color=BLUE, shrinkA=0, shrinkB=0, mutation_scale=5))
ax.plot([x0 - 1.5, x0 - 2.5, x0 - 2.5, x0 - 1.5], [ry[1] + 2, ry[1] + 2, ry[3] - 2, ry[3] - 2], color=INK, lw=0.5)
gx = [x0 + d for d in GX]
xe = x0 + BRL
for q, yy in ry.items():
    ax.plot([x0, xe], [yy, yy], color=WIRE, lw=0.5, zorder=1)
def box(x, w, yy, label, fc="white", ec=INK):
    w = max(w, textw(label, fontsize=6) + 2.4)             # at least 1.2 pt clear of the label on each side
    h = max(7.2, texth(label, fontsize=6) + 2.0)           # and 1 pt above and below it
    ax.add_patch(Rectangle((x - w / 2, yy - h / 2), w, h, fc=fc, ec=ec, lw=0.5, zorder=3))
    if "BOXES" in globals():
        BOXES[x] = (x - w / 2, x + w / 2)
    ax.text(x, yy, label, ha="center", va="center", color=INK, zorder=4, fontsize=6)
def ctrl(x, kind, open_=False):
    ax.plot([x, x], [ry[3], ry[1]], color=INK, lw=0.5, zorder=2)
    if open_:
        ax.add_patch(Circle((x, ry[1]), 1.3, fc="white", ec=INK, lw=0.5, zorder=4))
    else:
        ax.scatter([x], [ry[1]], s=4.5, c=INK, zorder=4, linewidths=0)
    for t in (2, 3):
        if kind == "X":
            ax.add_patch(Circle((x, ry[t]), 2.1, fc="white", ec=INK, lw=0.5, zorder=3))
            ax.plot([x - 2.1, x + 2.1], [ry[t], ry[t]], color=INK, lw=0.5, zorder=4)
            ax.plot([x, x], [ry[t] - 2.1, ry[t] + 2.1], color=INK, lw=0.5, zorder=4)
        else:
            ax.scatter([x], [ry[t]], s=4.5, c=INK, zorder=4, linewidths=0)
box(gx[0], 9.5, ry[1], r"$\sqrt{X}$")
ctrl(gx[1], "X")
ctrl(gx[2], "Z", open_=True)
box(gx[3], 6.5, ry[1], r"$T^\dagger$", fc="#fde7c7", ec=ORANGE)
ctrl(gx[4], "Z", open_=True)
BOXES = {}
box(gx[5], 5.5, ry[1], r"$T$", fc="#fde7c7", ec=ORANGE)
box(gx[6], 6.5, ry[1], r"$S^\dagger$")
assert BOXES[gx[6]][0] - BOXES[gx[5]][1] >= 1.0, f"T and Sdg boxes overlap or touch: gap {BOXES[gx[6]][0] - BOXES[gx[5]][1]:.2f} pt"
ax.text(xe + 0.8, ry[1], r"$|{+}\rangle$", ha="left", va="center", color=INK)
ax.text(xe + 0.8, (ry[2] + ry[3]) / 2, r"$|T\rangle^{\otimes2}$", ha="left", va="center", color=INK)
ax.text((x0 + xe) / 2, A["ylab"], r"$[[10,2,2]]$", ha="center", va="top", color=INK)
ax.text((x0 + xe) / 2, A["ylab"] - LINE, r"$8+2$ inputs", ha="center", va="top", color=INK)
assert span(La, P) <= xb - SEP + 0.01, "branch overlaps (a)"

B = parent(A["bottom"] - ROWGAP, 5, 3, True, pb, [(r"$[[8,4,2]]$", r"$\sim\!|\mathrm{CCZ}\rangle$"),
                                                 (r"$[[12,3,2]]$", r"$\sim\!|\mathrm{CS}_{12}\mathrm{CS}_{13}\rangle$"),
                                                 (r"$[[14,2,2]]$", r"$|T\rangle^{\otimes2}$"), (r"$[[15,1,3]]$", r"$|T\rangle$")],
           "(b)", WIDTH, lgap=7.0)
C = parent(B["bottom"] - ROWGAP, 5, 4, False, pc, [(r"$[[16,4,2]]$", r"$|\mathrm{CCCZ}\rangle$"),
                                                  (r"$[[24,3,2]]$", r"$|\mathrm{CCS}\rangle$"),
                                                  (r"$[[28,2,2]]$", r"$|\mathrm{CT}\rangle$"),
                                                  (r"$[[30,1,2]]$", r"$|\sqrt{T}\rangle$")], "(c)", WIDTH, lgap=7.0)
ytop, ybot = 0.5, C["bottom"] + 1.5          # the last label line has no descender below it
H = ytop - ybot
ax.set_xlim(0, WIDTH); ax.set_ylim(ybot, ytop)
fig.set_size_inches(WIDTH / 72.27, H / 72.27)
print(f"figure {WIDTH:.0f} x {H:.0f} pt; pitch (a) {P:.2f}, (b) {pb:.2f}, (c) {pc:.2f} pt; (b) ends at {B['right']:.0f}, (c) at {C['right']:.0f} pt; "
      f"APS estimate {150 * H / WIDTH + 20:.0f} words")

for out in (os.path.join(ROOT, "figures", "parent_l3.pdf"),
            os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "FIGURES", "parent_l3.pdf"),
            os.path.join(ROOT, "figures", "parent_l3.png")):
    if os.path.isdir(os.path.dirname(out)):
        fig.savefig(out, dpi=300); print("saved", out)
