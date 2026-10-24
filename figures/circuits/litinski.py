import re
"""Draw a parity-rotation circuit in the style of Litinski's figures, with the paper's colour scheme.

Each Z(pi/2^l) rotation on support S is one column: a box spanning the lowest to highest qubit of S, with "Z" on
every wire in S (wires inside the span but not in S pass through), and the angle written once below the box.
Clifford gates are plain boxes; post-selected X measurements end their wire. Pure TikZ, no quantikz."""
ROT_FILL, ROT_EDGE = "orange!18", "orange!70!black"
DEFAULT = r"\frac{\pi}{8}"
def tikz(nq, ops, labels=None, outs=None, cuts=None, dx=0.62, dy=0.5, title=None):
    """ops: ("H", qubits) | ("rot", support, angle_tex) | ("gate", q, tex) | ("cnot", c, targets) |
            ("xctrl", c, targets, target_tex) | ("meas", qubits) | ("cut", label)"""
    out = [r"\begin{tikzpicture}[x=1cm,y=1cm,font=\scriptsize]"]
    y = lambda q: -q * dy
    x = 0.6; alive = {q: 0.0 for q in range(nq)}; ends = {}
    cols = []
    for op in ops:
        k = op[0]
        if k == "H":
            for q in op[1]: out.append(r"\node[draw,fill=white,minimum size=0.36cm,inner sep=0pt] at (%.2f,%.2f) {$H$};" % (x, y(q)))
            x += dx
        elif k == "rot":
            S = sorted(op[1]); lo, hi = S[0], S[-1]
            w = 0.42; top = y(lo) + 0.2; bot = y(hi) - 0.2
            out.append(r"\draw[fill=%s,draw=%s,rounded corners=1pt] (%.2f,%.2f) rectangle (%.2f,%.2f);" % (ROT_FILL, ROT_EDGE, x - w/2, bot, x + w/2, top))
            for q in range(lo, hi + 1):
                if q in S: out.append(r"\node at (%.2f,%.2f) {$Z$};" % (x, y(q)))
                else: out.append(r"\draw (%.2f,%.2f) -- (%.2f,%.2f);" % (x - w/2, y(q), x + w/2, y(q)))
            # semicircular tab on top of the box carries the rotation angle (Litinski-style rotation, not a Pauli)
            r = w / 2 + (0.1 if op[2].startswith("-") or "16" in op[2] else 0.05)
            out.append(r"\draw[fill=%s,draw=%s] (%.2f,%.2f) arc (180:0:%.2f) -- cycle;" % (ROT_FILL, ROT_EDGE, x - r, top, r))
            lab = re.sub(r"\\frac\{\\pi\}\{(\d+)\}", r"\\pi/\1", op[2])
            out.append(r"\node[font=\fontsize{4.5}{5}\selectfont,inner sep=0pt] at (%.2f,%.2f) {$%s$};" % (x, top + r * 0.42, lab))
            x += dx * 0.9
        elif k == "gates":
            for q, t in op[1].items():
                out.append(r"\node[draw,fill=white,minimum size=0.36cm,inner sep=1pt] at (%.2f,%.2f) {$%s$};" % (x, y(q), t))
            x += dx
        elif k == "gate":
            out.append(r"\node[draw,fill=white,minimum size=0.36cm,inner sep=1pt] at (%.2f,%.2f) {$%s$};" % (x, y(op[1]), op[2]))
            x += dx
        elif k in ("cnot", "xctrl"):
            c, T = op[1], op[2]; lo, hi = min([c] + list(T)), max([c] + list(T))
            out.append(r"\draw (%.2f,%.2f) -- (%.2f,%.2f);" % (x, y(lo), x, y(hi)))
            if k == "cnot": out.append(r"\fill (%.2f,%.2f) circle (0.06);" % (x, y(c)))
            else: out.append(r"\draw[fill=white] (%.2f,%.2f) circle (0.06);" % (x, y(c)))
            for t in T:
                if k == "cnot":
                    out.append(r"\draw (%.2f,%.2f) circle (0.11); \draw (%.2f,%.2f) -- (%.2f,%.2f);" % (x, y(t), x - 0.11, y(t), x + 0.11, y(t)))
                    out.append(r"\draw (%.2f,%.2f) -- (%.2f,%.2f);" % (x, y(t) - 0.11, x, y(t) + 0.11))
                else:
                    out.append(r"\node[draw,fill=white,minimum size=0.34cm,inner sep=0pt] at (%.2f,%.2f) {$%s$};" % (x, y(t), op[3]))
            x += dx
        elif k == "meas":
            for q in op[1]:
                out.append(r"\node[draw,fill=white,minimum width=0.42cm,minimum height=0.32cm,inner sep=0pt,rounded corners=4pt] at (%.2f,%.2f) {$X$};" % (x, y(q)))
                ends[q] = x - 0.21
            x += dx
        elif k == "cut":
            out.append(r"\draw[blue!65,dashed,thick] (%.2f,%.2f) -- (%.2f,%.2f);" % (x - dx/2, 0.52, x - dx/2, y(nq - 1) - 0.3))
            out.append(r"\node[blue!65!black,font=\tiny] at (%.2f,%.2f) {$%s$};" % (x - dx/2, 0.64, op[1]))
    xend = x - dx * 0.4
    for q in range(nq):
        e = ends.get(q, xend)
        out.insert(1, r"\draw (0.25,%.2f) -- (%.2f,%.2f);" % (y(q), e, y(q)))
        lab = (labels or [f"q_{i}" for i in range(nq)])[q]
        out.append(r"\node[anchor=east] at (0.25,%.2f) {$%s$};" % (y(q), lab))
        if outs and q in outs: out.append(r"\node[anchor=west] at (%.2f,%.2f) {$%s$};" % (xend, y(q), outs[q]))
    out.append(r"\end{tikzpicture}")
    return "\n".join(out)
def doc(body):
    return r"\documentclass[border=3pt]{standalone}\usepackage{tikz}\usepackage{amsmath,amssymb}\newcommand{\ket}[1]{|#1\rangle}\begin{document}" + "\n" + body + "\n" + r"\end{document}"

# ---- Quirk gate list -> ops -------------------------------------------------------------
import re, json, urllib.parse
import os
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'Magic_State_Factories_tex', 'MSD_manuscript', 'circuits.tex')
SIGN = {"√i": r"\frac{\pi}{8}", "√-i": r"-\frac{\pi}{8}", "i": r"\frac{\pi}{4}", "-i": r"-\frac{\pi}{4}"}
TROT = {"Z^¼": r"\frac{\pi}{8}", "~d9b7": r"\frac{\pi}{8}", "~ntmg": r"\frac{\pi}{8}", "~fu1l": r"\frac{\pi}{8}",
        "Z^-¼": r"-\frac{\pi}{8}", "~703g": r"-\frac{\pi}{8}"}
CLIF = {"X^½": r"\sqrt{X}", "X^-½": r"\sqrt{X}^\dagger", "Z^½": "S", "Z^-½": r"S^\dagger", "~2gmh": "S", "~dvdi": "S", "~i2pe": r"S^\dagger", "Z": "Z"}
SKIP = ("Bloch", "Chance", "…")
def load(name):
    m = re.search(r'\\newcommand\{\\' + name + r'\}\{(.*?)\}\n', open(SRC).read())
    return json.loads(urllib.parse.unquote(m.group(1).replace('\\#', '#').replace('\\%', '%').split('circuit=')[1]))["cols"]
def from_quirk(name, cuts=None):
    cols = load(name); cuts = cuts or {}
    used = sorted({i for c in cols for i, x in enumerate(c) if x != 1 and x not in SKIP and x not in SIGN})
    off = used[0]; nq = used[-1] - off + 1
    ops = []
    for col in cols:
        g = {i - off: x for i, x in enumerate(col) if x != 1 and x not in SKIP and i >= off}
        if not g: continue
        sup = [q for q, x in g.items() if x == "zpar"]
        if sup:
            ph = [x for x in g.values() if x in SIGN]
            ops.append(("rot", sup, SIGN[ph[0]] if ph else r"\frac{\pi}{8}")); continue
        ctl = [q for q, x in g.items() if x in ("•", "⊖")]
        if ctl:
            T = [q for q, x in g.items() if x in ("X", "Z")]
            ops.append(("cnot", ctl[0], T) if g[ctl[0]] == "•" else ("xctrl", ctl[0], T, "Z")); continue
        meas = [q for q, x in g.items() if x == "|0⟩⟨0|"]
        if meas:
            # Quirk post-selects |0> after H (and possibly X): draw as an X-basis measurement
            for q in meas:
                for j in range(len(ops) - 1, -1, -1):
                    o = ops[j]
                    if o[0] == "gates" and q in o[1]:
                        if o[1][q] not in ("H", "X"): break
                        o[1].pop(q); continue
                    if o[0] == "rot" and q in o[1] or o[0] in ("cnot", "xctrl") and q in [o[1]] + list(o[2]): break
                ops[:] = [o for o in ops if not (o[0] == "gates" and not o[1])]
            ops.append(("meas", meas))
            if tuple(meas) in cuts: ops.append(("cut", cuts[tuple(meas)]))
            continue
        for q, x in g.items():
            if x in TROT: ops.append(("rot", [q], TROT[x]))
        d = {q: ("H" if x == "H" else "X" if x == "X" else CLIF[x]) for q, x in g.items() if x == "H" or x == "X" or x in CLIF}
        if d: ops.append(("gates", d))
    # merge separately-emitted H's into one leading column; collapse trailing empty
    return nq, ops

def malleable(nq, ops, last=None):
    """Reorder the (commuting) diagonal rotations of a parent circuit by anchor: rotations containing q0, its T, measure
    q0; then q1; ... Each cut is labelled [[n,k,2]] with n rotations so far and k unmeasured qubits."""
    head = [o for o in ops if o[0] == "gates"][:1]
    rots = [o for o in ops if o[0] == "rot"]
    out = list(head); n = 0
    for a in range(nq):
        mine = [r for r in rots if a in r[1]]
        rots = [r for r in rots if a not in r[1]]
        mine.sort(key=lambda r: (len(r[1]) == 1, -len(r[1])))
        out += mine; n += len(mine)
        out.append(("meas", [a]))
        if a < nq - 1: out.append(("cut", (last or {}).get(a, "[[%d,%d,2]]" % (n, nq - 1 - a))))
    return out
