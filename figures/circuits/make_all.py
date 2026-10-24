import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from litinski import tikz, doc, from_quirk, malleable
OUT = {"ten": {2: r"\ket{T}", 3: r"\ket{T}"}, "eleven": {0: r"\ket{T}"}, }
CUTS = {"six": {(0,): "[[8,5,2]]", (1, 2): "[[12,3,2]]"}}
for nm, q in {"four": "qparentfour", "eleven": "qelevenoneone", "ten": "qtentwotwo", "five": "qparentfive", "six": "qparentsix"}.items():
    nq, ops = from_quirk(q, CUTS.get(nm))
    if nm in ("four", "five"): ops = malleable(nq, ops, {3: "[[15,1,3]]"} if nm == "five" else None)
    if nm == "six":  # the Quirk link stops at [[12,3,2]]; complete the parent with T_3, T_4, T_5 (verified Clifford)
        ops += [("rot", [3], r"\frac{\pi}{8}"), ("meas", [3]), ("cut", "[[14,2,2]]"),
                ("rot", [4], r"\frac{\pi}{8}"), ("meas", [4]), ("cut", "[[15,1,3]]"), ("rot", [5], r"\frac{\pi}{8}"), ("meas", [5])]
    if ops[-1][0] == "cut": ops.pop()
    open(f"qk_{nm}.tex", "w").write(doc(tikz(nq, ops, outs=OUT.get(nm))))
# build:  python make_all.py && for f in qk_*.tex; do pdflatex $f; done
# l=4 analogue of the five-qubit parent (Table SII(c) [[16,4,2]], qVIfourh): all sixteen pi/16 rotations through the
# check qubit, drawn with the check on top; the later stages [[24,3,2]] -> [[28,2,2]] -> [[30,1,2]], [[31,1,3]] are elided.
import itertools
A16 = r"\frac{\pi}{16}"
rots = sorted((tuple([0] + [q + 1 for q in T]) for r in range(5) for T in itertools.combinations(range(4), r)),
              key=lambda S: (len(S) == 1, -len(S), S))
ops = [("gates", {q: "H" for q in range(5)})] + [("rot", list(S), A16) for S in rots] + [("meas", [0]), ("cut", "[[16,4,2]]")]
outs = {q: r"\cdots" for q in range(1, 5)}
open("qk_sixteen.tex", "w").write(doc(tikz(5, ops, outs=outs, dx=0.8)))
