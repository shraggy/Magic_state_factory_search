"""Free SAT at l=3 with native XOR constraints (CryptoMiniSat, pycryptosat), the same problem as bool_sat.py.

At l=3 a factory is a set of coefficient-1 rotations (bool_sat.py, docstring), and the factory condition is GF(2)
linear: for every set T of at most three qubits, the parity of #{chosen rotations containing T} is fixed (0 if T meets
a check, the target's digit if T lies in the outputs). CryptoMiniSat takes each such row as one XOR clause and runs
Gaussian elimination on them during search, instead of the integer encoding sum = 2z + b that CP-SAT needs. The bound
n < cap is a cardinality constraint, encoded as a sequential counter (Sinz 2005) in CNF. A query asks for a factory
with fewer inputs than `cap` and is repeated with the cap lowered until it is unsatisfiable or reaches its limit.

Deterministic: one thread, and the limit is a number of conflicts (`confl`), not a time; a query that reaches it
returns None (unknown). Every solution is expanded to rotations and checked with us.is_factory and us.distance != 1.
    python searches/xor_sat.py --spec CS03,CS12,CCZ015,CCZ024,CCZ123 --k 6 --N 8 --cap 26 --wmax 5 --confl 2000000
"""
import os, sys, time, argparse
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
from bool_sat import digits


def _atmost(s, xs, K, nxt):
    """Sequential counter: at most K of the literals xs are true. Returns the next free variable."""
    n = len(xs)
    if K >= n:
        return nxt
    if K == 0:
        for v in xs: s.add_clause([-v])
        return nxt
    r = [[0] * (K + 1) for _ in range(n)]                    # r[i][j]: at least j of xs[0..i] are true
    for i in range(n):
        for j in range(1, K + 1):
            r[i][j] = nxt; nxt += 1
    for i in range(n):
        x = xs[i]
        s.add_clause([-x, r[i][1]])
        if i > 0:
            for j in range(1, K + 1):
                s.add_clause([-r[i - 1][j], r[i][j]])
                if j > 1:
                    s.add_clause([-x, -r[i - 1][j - 1], r[i][j]])
            s.add_clause([-x, -r[i - 1][K]])                 # no (K+1)-th true literal
    return nxt


def query(spec, k, N, cap, wmax=None, confl=2_000_000, symmetry=True):
    """One query: a factory with n < cap. Returns ('SAT', n, gates), ('UNSAT', None, None) or ('UNKNOWN', None, None)."""
    import pycryptosat
    out = frozenset(range(k))
    sup = [frozenset(S) for w in range(1, N + 1) for S in combinations(range(N), w)
           if (wmax is None or w <= wmax) and not frozenset(S) <= out]
    var = {S: i + 1 for i, S in enumerate(sup)}; nxt = len(sup) + 1
    s = pycryptosat.Solver(threads=1, confl_limit=confl)
    D = digits(spec, k)
    for t in range(1, 4):
        for T in combinations(range(N), t):
            T = frozenset(T)
            ids = [var[S] for S in sup if T <= S]
            want = bool(T <= out and T in D)
            if not ids:
                if want:
                    return "UNSAT", None, None
                continue
            s.add_xor_clause(ids, want)
    if symmetry and N - k >= 2:                              # interchangeable checks: the first check is in at least
        q0 = [var[S] for S in sup if k in S]                 # one rotation (a valid break, cheap in CNF)
        s.add_clause(q0)
    if cap is not None:
        nxt = _atmost(s, list(var.values()), cap - 1, nxt)
    st, model = s.solve()
    if st is None:
        return "UNKNOWN", None, None
    if not st:
        return "UNSAT", None, None
    g = [(S, 1) for S in sup if model[var[S]]]
    return "SAT", len(g), g


def solve(spec, k, N, cap=None, wmax=None, confl=2_000_000):
    """Lower the cap until UNSAT or UNKNOWN; returns (last status, best n, best gates, queries)."""
    best_n, best_g, q = None, None, 0
    cap = cap if cap is not None else 10 ** 6
    while True:
        st, n, g = query(spec, k, N, cap, wmax, confl); q += 1
        if st != "SAT":
            return st, best_n, best_g, q
        f = dict(k=k, r=N - k, gates=g)
        assert us.is_factory(f, 3) and us.distance(f, 3) != 1, "solver returned a non-factory"
        best_n, best_g, cap = n, g, n


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec"); ap.add_argument("--k", type=int); ap.add_argument("--N", type=int)
    ap.add_argument("--cap", type=int, default=None); ap.add_argument("--wmax", type=int, default=None)
    ap.add_argument("--confl", type=int, default=2_000_000); a = ap.parse_args()
    t = time.time(); st, n, g, q = solve(a.spec, a.k, a.N, a.cap, a.wmax, a.confl)
    print(st, "best n", n, "queries", q, round(time.time() - t, 1), "s", flush=True)
