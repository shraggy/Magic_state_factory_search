"""SAT stage restricted by App. S4: fixing the checks of a distance-2 factory leaves a circuit preparing its target
from the rotations' output parts. Here every rotation's output part is drawn from a synthesis of the target (columns
of an optimal parity network A, of Lempel factors B for mu and mu+1, the linear-term column c, and the empty part),
with any check pattern; the condition solved is the full borrowed-identity condition of sat_certify. Each found
factory is re-verified with is_factory. Deterministic z3 step limit, fresh process per query.
    python searches/synthesis_sat.py --k 5 --spec T0,T1,CCZ014,CCZ023 --cap 24
"""
import os, sys, json, argparse, multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))


def parts(spec, k):
    import class_table_l3 as CT, synthillation_compare as sc
    P = CT.parse(spec); cols = set(sc.opt_A(k, P))
    mu = 0
    while sc.opt_B(k, P, mu) is None:
        mu += 1
    for m in (mu, mu + 1):
        B = sc.opt_B(k, P, m)
        if B: cols |= set(B)
    c = sum(1 << i for i in range(k) if frozenset([i]) in P)
    if c: cols.add(c)
    return {frozenset(q for q in range(k) if y >> q & 1) for y in cols} | {frozenset()}


def _q(spec, k, N, cap, rl, op, q):
    from coset_certify import parse_target
    from sat_certify import certify
    q.put(certify(parse_target(spec, 3), k, N, 3, cap, time_limit=None, rlimit=rl, output_parts=op))


def query(spec, k, N, cap, rl, op):
    ctx = mp.get_context("spawn"); q = ctx.Queue(); p = ctx.Process(target=_q, args=(spec, k, N, cap, rl, op, q))
    p.start(); r = q.get(); p.join(); return r


def solve(spec, k, cap, rl=3e9, span=(1, 2, 3)):
    import unified_search as us
    op = parts(spec, k); best = None; log = []
    for r in span:
        N = k + r
        while True:
            res, n, g = query(spec, k, N, cap, int(rl), op); log.append([N, cap, res])
            if res != "sat":
                break
            f = dict(stage="synthesis-sat", k=k, r=r, gates=[(frozenset(S), c) for S, c in g], params="")
            assert us.is_factory(f, 3) and us.distance(f, 3) != 1
            best = dict(n=n, N=N, gates=[[sorted(S), c] for S, c in g]); cap = n
    return best, log


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--k", type=int); ap.add_argument("--spec"); ap.add_argument("--cap", type=int)
    ap.add_argument("--rlimit", type=float, default=3e9); ap.add_argument("--index", type=int, default=None)
    ap.add_argument("--span", default="1,2,3", help="check-qubit counts r to try (N = k + r)")
    ap.add_argument("--retry", default=None, help="file of k,spec,cap lines: rerun these (index into it)")
    a = ap.parse_args()
    if a.retry is not None:                                # rerun of classes left 'unknown', from their best n
        a.k, a.spec, a.cap = open(a.retry).read().split()[a.index].split(";"); a.k, a.cap = int(a.k), int(a.cap)
    elif a.index is not None:                              # cluster task: the a.index-th class still above synthillation
        E = json.load(open(os.path.join(ROOT, "outputs", "entangled_gluing_l3.json")))
        lose = sorted([e for e in E if e["syn_n"] and min(x for x in (e["cascade_n"], e["eg_n"]) if x) > e["syn_n"]],
                      key=lambda e: (e["k"], e["spec"]))
        e = lose[a.index]; a.k, a.spec, a.cap = e["k"], e["spec"], min(x for x in (e["cascade_n"], e["eg_n"]) if x)
    best, log = solve(a.spec, a.k, a.cap, a.rlimit, tuple(int(x) for x in a.span.split(",")))
    out = os.path.join(ROOT, "outputs", "synthesis_sat"); os.makedirs(out, exist_ok=True)
    json.dump(dict(k=a.k, spec=a.spec, cap=a.cap, best=best, log=log), open(os.path.join(out, f"k{a.k}_{a.spec}" + ("_retry" if a.retry else "") + ".json"), "w"))
    print(a.spec, a.cap, "->", best and (best["n"], best["N"]), log, flush=True)
