"""History-independent synthesis parts. synthillation_compare.opt_B solves the Lempel condition with z3 in the global
context, so the factor it returns (one of many valid ones) depends on every z3 query run earlier in the same process;
synthesis_sat.parts and pattern_sat.multiset call it, so the output parts of synthesis SAT, its XOR form and check-pattern
SAT could change with process history (measured: 4 successive SS.parts calls on one class give 4 different sets).
`install()` replaces synthillation_compare.opt_B, in this process only, by the same constraints in a fresh z3.Context
per call, so the result depends only on (k, target, mu). The module file synthillation_compare.py is not changed.
    python searches/det_parts.py        # self-test: same parts in any call order, repeated calls and after other z3 work
"""
import os, sys, hashlib
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import z3


def opt_B_fresh(k, target, mu):
    """synthillation_compare.opt_B with its own z3 context (same constraints, same variable order)."""
    ctx = z3.Context(); s = z3.Solver(ctx=ctx)
    b = [[z3.Bool(f"b{i}_{j}", ctx) for j in range(mu)] for i in range(k)]
    one, zero = z3.IntVal(1, ctx), z3.IntVal(0, ctx)
    x = lambda c: z3.If(c, one, zero)
    for i in range(k):
        s.add(z3.Sum([x(b[i][j]) for j in range(mu)]) % 2 == target.get(frozenset([i]), 0) % 2)
        for i2 in range(i + 1, k):
            q = (target.get(frozenset([i, i2]), 0) // 2) % 2
            s.add(z3.Sum([x(z3.And(b[i][j], b[i2][j])) for j in range(mu)]) % 2 == q)
    for j in range(mu):
        s.add(z3.Or([b[i][j] for i in range(k)]))
    if s.check() != z3.sat:
        return None
    m = s.model()
    return [sum(1 << i for i in range(k) if z3.is_true(m[b[i][j]])) for j in range(mu)]


def install():
    import synthillation_compare as sc
    sc.opt_B = opt_B_fresh


def _h(p):
    return hashlib.md5(repr(sorted(sorted(x) for x in p)).encode()).hexdigest()[:10]


if __name__ == "__main__":
    install()
    import json, synthesis_sat as SS, synthillation_compare as sc, class_table_l3 as CT
    sp = json.load(open(os.path.join(ROOT, "outputs", "staged_A", "l3_k6", "ii_eg.json")))["specs"]
    cls = ["2320:13064", "2352:1548", "2320:771", "41984:5"]
    ref = {o: _h(SS.parts(sp[o], 6)) for o in cls}
    ok = True
    for rep in range(3):
        for o in (cls if rep % 2 == 0 else cls[::-1]):
            z3.Solver().check(z3.Bool("noise") == True)                 # unrelated z3 work in the global context
            ok &= _h(SS.parts(sp[o], 6)) == ref[o]
    # every returned B satisfies the Lempel condition (the check opt_B's callers rely on)
    for o in cls:
        P = CT.parse(sp[o]); mu = 0
        while opt_B_fresh(6, P, mu) is None:
            mu += 1
        B = opt_B_fresh(6, P, mu)
        for i in range(6):
            ok &= sum(c >> i & 1 for c in B) % 2 == P.get(frozenset([i]), 0) % 2
    print("deterministic parts:", ok, ref, flush=True)
