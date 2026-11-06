"""Deterministic enumeration of several syntheses of a target, for the restricted (synthesis) SAT stage.

synthesis_sat.parts uses one optimal parity network A (synthillation_compare.opt_A) and one Lempel factor B for mu and
mu+1 (opt_B). Neither is unique. Here, all in a fresh z3.Context (no dependence on process history):
  lempel_factors(k, P, mu, F)  up to F distinct Lempel factors B with mu columns, columns strictly increasing as integers
                               (a canonical form: each factor is a set of distinct columns; repeated columns would
                               cancel in pairs and give a factor with fewer columns), blocking clause after each model,
                               so the list is fixed by (k, P, mu, F)
  opt_A_words(k, P, F)         up to F distinct minimum-weight coset words (optimal parity networks), in a fixed order
                               (ascending integer value of the word over the enumeration of coset elements)
  syntheses(spec, k, F)        list of part sets, factor f = (A_f, B_f(mu), B_f(mu+1)) for f = 0..F-1 (the i-th A
                               with the i-th B's; missing ones repeat the last), each set as synthesis_sat.parts does
  union_parts(spec, k, F)      the union of those part sets
F = 1 gives the first factor of each list.
    python searches/multi_parts.py --stats       # how often A and B are not unique at l=3, k<=6 (all 766 classes)
"""
import os, sys, itertools, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import z3


def lempel_factors(k, P, mu, F):
    if mu == 0:
        ok = all(P.get(frozenset([i]), 0) % 2 == 0 for i in range(k)) and \
             all((P.get(frozenset([i, j]), 0) // 2) % 2 == 0 for i in range(k) for j in range(i + 1, k))
        return [[]] if ok else []
    ctx = z3.Context(); s = z3.Solver(ctx=ctx)
    col = [z3.BitVec(f"c{j}", k, ctx) for j in range(mu)]
    bit = lambda j, i: z3.Extract(i, i, col[j]) == z3.BitVecVal(1, 1, ctx)
    one, zero = z3.IntVal(1, ctx), z3.IntVal(0, ctx)
    x = lambda c: z3.If(c, one, zero)
    for i in range(k):
        s.add(z3.Sum([x(bit(j, i)) for j in range(mu)]) % 2 == P.get(frozenset([i]), 0) % 2)
        for i2 in range(i + 1, k):
            q = (P.get(frozenset([i, i2]), 0) // 2) % 2
            s.add(z3.Sum([x(z3.And(bit(j, i), bit(j, i2))) for j in range(mu)]) % 2 == q)
    for j in range(mu):
        s.add(col[j] != z3.BitVecVal(0, k, ctx))
    for j in range(mu - 1):
        s.add(z3.ULT(col[j], col[j + 1]))                      # canonical: strictly increasing columns
    out = []
    while len(out) < F and s.check() == z3.sat:
        m = s.model(); B = [m.eval(c).as_long() for c in col]
        out.append(B)
        s.add(z3.Or([c != z3.BitVecVal(v, k, ctx) for c, v in zip(col, B)]))
    return out


def opt_A_words(k, P, F):
    """Up to F minimum-weight coset words (as column lists), as synthillation_compare.opt_A but all minima kept."""
    from tcount import snf_z8, _solve_z8, parity_matrix
    M = parity_matrix(k); U, W, D = snf_z8(M); nv = (1 << k) - 1
    f = []
    for xv in range(1, nv + 1):
        xs = {i for i in range(k) if xv >> i & 1}
        f.append(sum(c for T, c in P.items() if T <= xs) % 8)
    a0 = _solve_z8(U, W, D, f)
    basis = [sum(1 << r for r in range(nv) if W[r][i] & 1) for i in range(nv) if D[i] % 8 == 0]
    b = sum(1 << r for r in range(nv) if a0[r] & 1)
    best, words = None, []
    for c in itertools.product((0, 1), repeat=len(basis)):
        v = b
        for ci, bv in zip(c, basis):
            if ci:
                v ^= bv
        w = bin(v).count("1")
        if best is None or w < best:
            best, words = w, [v]
        elif w == best:
            words.append(v)
    words = sorted(set(words))[:F]
    return [[y for y in range(1, nv + 1) if v >> (y - 1) & 1] for v in words]


def _mu(k, P):
    mu = 0
    while not lempel_factors(k, P, mu, 1):
        mu += 1
    return mu


def syntheses(spec, k, F):
    import class_table_l3 as CT
    P = CT.parse(spec); mu = _mu(k, P)
    As = opt_A_words(k, P, F)
    Bs = [lempel_factors(k, P, m, F) for m in (mu, mu + 1)]
    c = sum(1 << i for i in range(k) if frozenset([i]) in P)
    out = []
    for f in range(F):
        cols = set(As[min(f, len(As) - 1)])
        for L in Bs:
            if L:
                cols |= set(L[min(f, len(L) - 1)])
        if c:
            cols.add(c)
        ps = frozenset(frozenset(q for q in range(k) if y >> q & 1) for y in cols) | {frozenset()}
        if ps not in out:
            out.append(ps)
    return out


def union_parts(spec, k, F):
    u = set()
    for ps in syntheses(spec, k, F):
        u |= ps
    return frozenset(u)


if __name__ == "__main__":
    if "--stats" in sys.argv:
        import class_table_l3 as CT
        rows = json.load(open(os.path.join(ROOT, "outputs", "solver_v2", "staged_A_v4", "best_per_class.json")))
        st = dict(classes=0, A_nonunique=0, B_nonunique=0, nA=[], nB=[], nparts={})
        for key, r in sorted(rows.items()):
            k = r["k"]; P = CT.parse(r["spec"]); mu = _mu(k, P)
            nA = len(opt_A_words(k, P, 64)); nB = len(lempel_factors(k, P, mu, 64))
            st["classes"] += 1; st["A_nonunique"] += nA > 1; st["B_nonunique"] += nB > 1
            st["nA"].append(nA); st["nB"].append(nB)
            st["nparts"][key] = [len(union_parts(r["spec"], k, F)) for F in (1, 2, 4, 8)]
        print({a: st[a] for a in ("classes", "A_nonunique", "B_nonunique")},
              "max A words", max(st["nA"]), "max B (cap 64)", max(st["nB"]), flush=True)
        json.dump(st, open(os.path.join(ROOT, "outputs", "solver_v2", "multi_parts_stats.json"), "w"))
