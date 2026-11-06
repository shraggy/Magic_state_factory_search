"""Incremental XOR free SAT at l=3 (CryptoMiniSat): the problem of xor_sat.py / bool_sat.py, one solver per (class, N).

The factory condition is the same set of GF(2) rows (one XOR clause per set T of at most three qubits). The difference
from xor_sat.py: the cardinality bound is a sequential counter built once, up to the starting cap, and every later
bound n < c is passed as an assumption (the counter literal "at least c of the rotations"), so the solver keeps its
learned clauses while the cap is lowered. Each query has its own conflict limit `confl`, and the whole (class, N) run a
total limit `total` (sum of conflicts actually spent is not exposed by pycryptosat, so the total is counted as
queries x confl, an upper bound): deterministic, one thread, independent of machine speed.

Check symmetry: the N - k checks are interchangeable, so the first check is required to lie in some chosen rotation
(a valid break, as in xor_sat.py). With no starting cap, the first query runs without a counter and the counter is
built up to the first solution's n. With `op` (a set of output parts, e.g. synthesis_sat.parts), only rotations
whose output part is in `op` are allowed: the XOR form of synthesis SAT (coefficient-1 rotations only). Returned factories are verified with us.is_factory and us.distance != 1.
    python searches/xor_sat2.py --spec CS03,CS12,CCZ015,CCZ034 --k 6 --N 8 --cap 24 --wmax 5 --confl 2000000
"""
import os, sys, time, argparse
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
from bool_sat import digits


def counter(s, xs, K, nxt):
    """Sequential counter with outputs: returns (out, nxt) where out[j] (j=1..K+1) is implied true when at least j of
    xs are true (one-directional, enough for upper bounds). Asserting -out[c] forces at most c-1 true."""
    n = len(xs); K = min(K, n - 1)
    r = [[0] * (K + 2) for _ in range(n)]
    for i in range(n):
        for j in range(1, K + 2):
            r[i][j] = nxt; nxt += 1
    for i in range(n):
        x = xs[i]
        s.add_clause([-x, r[i][1]])
        if i > 0:
            for j in range(1, K + 2):
                s.add_clause([-r[i - 1][j], r[i][j]])
                if j > 1:
                    s.add_clause([-x, -r[i - 1][j - 1], r[i][j]])
    return r[n - 1], nxt


def build(spec, k, N, wmax, cap, op=None):
    import pycryptosat
    out = frozenset(range(k))
    sup = [frozenset(S) for w in range(1, N + 1) for S in combinations(range(N), w)
           if (wmax is None or w <= wmax) and not frozenset(S) <= out and (op is None or frozenset(S) & out in op)]
    var = {S: i + 1 for i, S in enumerate(sup)}; nxt = len(sup) + 1
    s = pycryptosat.Solver(threads=1)
    D = digits(spec, k)
    for t in range(1, 4):
        for T in combinations(range(N), t):
            T = frozenset(T)
            ids = [var[S] for S in sup if T <= S]
            want = bool(T <= out and T in D)
            if not ids:
                if want:
                    return None
                continue
            s.add_xor_clause(ids, want)
    if N - k >= 2:
        s.add_clause([var[S] for S in sup if k in S])
    cnt = None
    if cap is not None and cap - 1 < len(sup):
        cnt, nxt = counter(s, list(var.values()), cap - 1, nxt)
        s.add_clause([-cnt[cap]])                       # n <= cap-1 permanently
    return s, sup, var, cnt, nxt


def solve(spec, k, N, cap=None, wmax=None, confl=2_000_000, total=None, op=None):
    """Lower the cap until UNSAT or a query hits `confl`; returns (status, best n, gates, queries, trace)."""
    b = build(spec, k, N, wmax, cap, op)
    if b is None:
        return "UNSAT", None, None, 0, []
    s, sup, var, cnt, nxt = b
    best_n, best_g, q, trace = None, None, 0, []
    c = cap
    total = total if total is not None else 10 ** 18
    while True:
        if (q + 1) * confl > total:
            return "UNKNOWN", best_n, best_g, q, trace
        assum = [-cnt[c]] if (cnt is not None and c < len(cnt)) else []
        st, model = s.solve(assum, confl_limit=confl); q += 1
        if st is None:
            trace.append(("UNKNOWN", c)); return "UNKNOWN", best_n, best_g, q, trace
        if not st:
            trace.append(("UNSAT", c)); return "UNSAT", best_n, best_g, q, trace
        g = [(S, 1) for S in sup if model[var[S]]]
        f = dict(k=k, r=N - k, gates=g)
        assert us.is_factory(f, 3) and us.distance(f, 3) != 1, "solver returned a non-factory"
        best_n, best_g = len(g), g; trace.append(("SAT", len(g))); c = len(g)
        if cnt is None:                                  # first solution without a cap: build the counter to its n
            cnt, nxt = counter(s, list(var.values()), c - 1, nxt)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec"); ap.add_argument("--k", type=int); ap.add_argument("--N", type=int)
    ap.add_argument("--cap", type=int, default=None); ap.add_argument("--wmax", type=int, default=None)
    ap.add_argument("--confl", type=int, default=2_000_000); a = ap.parse_args()
    t = time.time(); st, n, g, q, tr = solve(a.spec, a.k, a.N, a.cap, a.wmax, a.confl)
    print(st, "best n", n, "queries", q, tr, round(time.time() - t, 1), "s", flush=True)
