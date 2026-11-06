"""Stage (iii) cascade with the multi-factor restricted stage (chain_iii2.py plus two stage kinds; same output schema).

New stages (XOR form, coefficient-1 rotations, output parts restricted to syntheses of the target, multi_parts.py):
  xsa:F:W:C   variant (a): output parts = union of the first F syntheses, one query sequence per N
  xsb:F:W:C   variant (b): one query sequence per synthesis f = 0..F-1 (its own parts), the cap lowered to the best so far
W = wmax (None, 5, 6), C = conflicts per query. Every query runs in a fresh spawned process (chain_iii2.xor_fresh), the
syntheses are enumerated in a fresh z3 context (multi_parts), so the result depends only on the inputs.
Existing kinds (xs:W:C, x:W:C) behave as in chain_iii2.py; z3 synthesis (syn) is still accepted but not used.
    python searches/chain_iii3.py k I M SEEDDIR OUT CHAIN [--classes o1,o2]   e.g. CHAIN xsb:4:None:2e6,xsb:4:6:2e6,x:5:2e6
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import chain_iii2 as C2
import multi_parts as MP
_base_stage = C2.stage


def _xor_seq(spec, k, cap, w, c, op):
    out = []
    for N in range(k + 1, k + 4):
        st, n, g = C2.xor_fresh(spec, k, N, cap, w, c, op)
        if g:
            out.append((n, N, g)); cap = n
    return out


def stage(name, k, spec, cap):
    kind = name.split(":")[0]
    if kind not in ("xsa", "xsb"):
        return _base_stage(name, k, spec, cap)
    _, F, w, c = name.split(":"); F = int(F); w = None if w == "None" else int(w); c = int(float(c))
    if kind == "xsa":
        return _xor_seq(spec, k, cap, w, c, MP.union_parts(spec, k, F))
    out = []
    for ps in MP.syntheses(spec, k, F):
        r = _xor_seq(spec, k, cap, w, c, ps)
        out += r
        if r:
            cap = min(cap, min(x[0] for x in r))
    return out


C2.stage = stage                                   # chain_iii2.run looks stage() up in its own module


if __name__ == "__main__":
    a = sys.argv[1:]; only = None
    if "--classes" in a:
        i = a.index("--classes"); only = a[i + 1].split(","); a = a[:i] + a[i + 2:]
    C2.run(int(a[0]), int(a[1]), int(a[2]), a[3], a[4], a[5].split(","), only)
