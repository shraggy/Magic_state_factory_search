"""Check-pattern search (App. S4 made constructive): fix the output parts of a factory's T rotations to a synthesis
multiset of the target (columns of an optimal network A once, a Lempel factor B twice, the linear-term column c four
times, as the fixed-check circuit of App. S4 must prepare the target) plus m rotations on checks alone, and solve
only for each rotation's check pattern p_j in {0,1}^R \\ {0} under the borrowed-identity condition (monomial form,
l=3, every coefficient 1). Boolean variables only. Each solution is verified with is_factory.
    python searches/pattern_sat.py --k 6 --spec T4,CS03,CS12,CCZ013 --R 3
"""
import os, sys, json, argparse, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import z3


def multiset(spec, k, mu):
    import class_table_l3 as CT, synthillation_compare as sc
    P = CT.parse(spec); A = sc.opt_A(k, P); B = sc.opt_B(k, P, mu)
    if B is None:
        return None, P
    c = sum(1 << i for i in range(k) if frozenset([i]) in P)
    return A + B + B + ([c] * 4 if c else []), P


def solve(spec, k, R, mu, m, rlimit=2e9):
    from unified_search import cliff_exp, is_factory, distance
    cols, P = multiset(spec, k, mu)
    if cols is None:
        return None
    J = len(cols) + m
    outs = cols + [0] * m                                   # output part of each rotation (bitmask on k outputs)
    ctx = z3.Context()                                      # a fresh context per query: no solver state carries
    p = [[z3.Bool(f"p{j}_{r}", ctx) for r in range(R)] for j in range(J)]   # over from earlier queries, so the
    s = z3.Solver(ctx=ctx); s.set("rlimit", int(rlimit))    # same query gives the same answer in any run order
    for j in range(J):
        s.add(z3.Or(p[j]))                                  # every rotation touches a check: distance 2
    bits = (lambda row: z3.Concat(*[z3.If(x, z3.BitVecVal(1, 1, ctx), z3.BitVecVal(0, 1, ctx)) for x in row]) if len(row) > 1
            else z3.If(row[0], z3.BitVecVal(1, 1, ctx), z3.BitVecVal(0, 1, ctx)))
    for j in range(len(cols), J - 1):                       # check-only rotations interchangeable: order them
        s.add(z3.ULE(bits(p[j]), bits(p[j + 1])))
    one = z3.BitVecVal(1, 3, ctx); zero = z3.BitVecVal(0, 3, ctx)
    for t in range(1, 4):
        for To in range(0, t + 1):
            for O in itertools.combinations(range(k), To):
                om = sum(1 << q for q in O)
                Js = [j for j in range(J) if outs[j] & om == om]
                for C in itertools.combinations(range(R), t - To):
                    if not C:
                        continue                            # output-only monomials: fixed by the synthesis
                    e = cliff_exp(t, 3); mod = 1 << e
                    terms = [z3.If(z3.And([p[j][r] for r in C]), one, zero) for j in Js]
                    if not terms:
                        continue
                    lhs = z3.Sum(terms) * z3.BitVecVal(((-2) ** (t - 1)) % 8, 3, ctx)
                    s.add(z3.Extract(e - 1, 0, lhs) == 0) if e > 0 else None
    res = s.check()
    if res != z3.sat:
        return str(res)
    mdl = s.model()
    gates = []
    for j in range(J):
        chk = {k + r for r in range(R) if z3.is_true(mdl[p[j][r]])}
        gates.append((frozenset({q for q in range(k) if outs[j] >> q & 1} | chk), 1))
    f = dict(stage="pattern-sat", k=k, r=R, gates=gates, params="")
    ok = is_factory(f, 3) and distance(f, 3) != 1
    return dict(n=J, N=k + R, ok=ok, gates=[[sorted(S), c] for S, c in gates])


def _main_single():
    ap = argparse.ArgumentParser(); ap.add_argument("--k", type=int); ap.add_argument("--spec")
    ap.add_argument("--R", type=int, default=3); ap.add_argument("--mu", type=int, default=None)
    a = ap.parse_args()
    import synthillation_compare as sc, class_table_l3 as CT
    P = CT.parse(a.spec); mu = a.mu
    if mu is None:
        mu = 0
        while sc.opt_B(a.k, P, mu) is None:
            mu += 1
    for m in range(0, 6):
        r = solve(a.spec, a.k, a.R, mu, m)
        print(a.spec, "mu", mu, "m", m, "->", r if not isinstance(r, dict) else (r["n"], r["N"], r["ok"]), flush=True)
        if isinstance(r, dict) and r["ok"]:
            break


def run_class(k, spec, cap, Rs=(1, 2, 3, 4), mus=None, mmax=12):
    """Smallest n over check counts R, Lempel sizes mu (smallest and next) and check-only rotations m, below cap."""
    import synthillation_compare as sc, class_table_l3 as CT
    P = CT.parse(spec)
    if mus is None:
        mu = 0
        while sc.opt_B(k, P, mu) is None:
            mu += 1
        mus = (mu, mu + 1)
    best = None; log = []
    for mu in mus:
        cols, _ = multiset(spec, k, mu)
        if cols is None:
            continue
        for R in Rs:
            for m in range(0, mmax + 1):
                n = len(cols) + m
                if n >= (best["n"] if best else cap):
                    break
                r = solve(spec, k, R, mu, m)
                log.append([mu, R, m, r if not isinstance(r, dict) else "sat"])
                if isinstance(r, dict) and r["ok"]:
                    if best is None or (r["n"], r["N"]) < (best["n"], best["N"]):
                        best = r
                    break
    return best, log


def _all(workers):
    """Every l=3 class with k <= 6 (outputs/class_table_l3.json): best n below our current best."""
    from multiprocessing import Pool
    C = json.load(open(os.path.join(ROOT, "outputs", "class_table_l3.json")))
    E = {(e["k"], e["spec"]): e for e in json.load(open(os.path.join(ROOT, "outputs", "entangled_gluing_l3.json")))}
    items = [(r["k"], r["spec"], min(x for x in (r["n"], E[(r["k"], r["spec"])]["eg_n"]) if x) + 1) for r in C]
    with Pool(workers) as p:
        res = p.starmap(run_class, items, chunksize=1)
    out = [dict(k=k, spec=s, best=b and dict(n=b["n"], N=b["N"], gates=b["gates"]), log=l) for (k, s, _), (b, l) in zip(items, res)]
    json.dump(out, open(os.path.join(ROOT, "outputs", "pattern_sat_l3.json"), "w"))
    print("done", len(out), flush=True)


if __name__ == "__main__":
    if "--index" in sys.argv:                               # cluster task: one class of outputs/class_table_l3.json
        i = int(sys.argv[sys.argv.index("--index") + 1])
        r = json.load(open(os.path.join(ROOT, "outputs", "class_table_l3.json")))[i]
        E = {(e["k"], e["spec"]): e for e in json.load(open(os.path.join(ROOT, "outputs", "entangled_gluing_l3.json")))}
        cap = min(x for x in (r["n"], E[(r["k"], r["spec"])]["eg_n"]) if x) + 1
        b, l = run_class(r["k"], r["spec"], cap)
        od = os.path.join(ROOT, "outputs", "pattern_sat"); os.makedirs(od, exist_ok=True)
        json.dump(dict(k=r["k"], spec=r["spec"], cap=cap, best=b and dict(n=b["n"], N=b["N"], gates=b["gates"]), log=l),
                  open(os.path.join(od, f"{i}.json"), "w"))
        print(r["spec"], cap - 1, "->", b and (b["n"], b["N"]), flush=True)
    elif "--all" in sys.argv:
        _all(int(sys.argv[sys.argv.index("--all") + 1]))
    else:
        _main_single()
