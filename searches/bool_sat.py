"""Free SAT at l=3 as a Boolean problem. At l=3 an even coefficient is a Clifford rotation, which costs no input and
changes neither the factory condition nor the output class (checked on all 17 catalogue factories whose best came from
free SAT: reducing every coefficient to its parity keeps the factory, its exact class and n). So a factory is a set of
rotations with coefficient 1, one Boolean x_Lambda per support Lambda that meets a check, and the condition is GF(2)
linear: every set T of at most three qubits that meets a check lies in an even number of chosen rotations, and for
T inside the outputs the parity of #{chosen rotations containing T} equals the target's digit on T (digits_of_spec,
the coordinates of searches/target_classes.py at l=3). The cost is the number of chosen rotations, n. So free SAT at
l=3 is: minimise |x| subject to linear equations over GF(2), a minimum-weight coset problem, solved here with CP-SAT
(one worker, fixed seed, deterministic time limit, so the outcome does not depend on the machine).

Every solution is expanded to rotations and checked with us.is_factory, us.distance != 1 and the exact class label.
    python searches/bool_sat.py --spec CS03,CS12,CCZ015,CCZ024,CCZ123 --k 6 --N 8 --cap 26 --wmax 5 --dlimit 120
"""
import os, sys, json, time, argparse
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us


def digits(spec, k):
    import orbit_group as OG
    return {frozenset(q for q in range(k) if m >> q & 1) for m in OG.digits_of_spec(spec, k)}


def solve(spec, k, N, cap=None, wmax=None, dlimit=60.0, symmetry=True):
    """Min-n factory for `spec` on N qubits with coefficient-1 rotations, n < cap. Returns (status, n, gates)."""
    from ortools.sat.python import cp_model
    out = frozenset(range(k)); chk = range(k, N)
    sup = [frozenset(S) for w in range(1, N + 1) for S in combinations(range(N), w)
           if (wmax is None or w <= wmax) and not frozenset(S) <= out]
    D = digits(spec, k)
    m = cp_model.CpModel(); x = {S: m.NewBoolVar("") for S in sup}
    for t in range(1, 4):
        for T in combinations(range(N), t):
            T = frozenset(T)
            ids = [x[S] for S in sup if T <= S]
            want = 1 if (T <= out and T in D) else 0
            if not ids:
                if want:
                    return "INFEASIBLE", None, None
                continue
            z = m.NewIntVar(0, len(ids), "")
            m.Add(sum(ids) == 2 * z + want)
    if symmetry and N - k >= 2:                      # checks interchangeable: order them by the count of rotations
        for q in range(k, N - 1):                    # through each (a valid symmetry break: some ordering exists)
            m.Add(sum(x[S] for S in sup if q in S) >= sum(x[S] for S in sup if q + 1 in S))
    n = sum(x.values())
    if cap is not None:
        m.Add(n <= cap - 1)
    m.Minimize(n)
    s = cp_model.CpSolver(); s.parameters.num_workers = 1; s.parameters.random_seed = 0
    s.parameters.max_deterministic_time = dlimit
    st = s.Solve(m); name = s.StatusName(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        g = [(S, 1) for S in sup if s.Value(x[S])]
        return name, len(g), g
    return name, None, None


def check(spec, k, N, g):
    import cascade_class as CC
    f = dict(k=k, r=N - k, gates=g)
    ok = us.is_factory(f, 3) and us.distance(f, 3) != 1
    return ok


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--spec"); ap.add_argument("--k", type=int); ap.add_argument("--N", type=int)
    ap.add_argument("--cap", type=int, default=None); ap.add_argument("--wmax", type=int, default=None)
    ap.add_argument("--dlimit", type=float, default=60.0); a = ap.parse_args()
    t = time.time(); st, n, g = solve(a.spec, a.k, a.N, a.cap, a.wmax, a.dlimit)
    print(st, n, "factory" if g and check(a.spec, a.k, a.N, g) else "", round(time.time() - t, 1), "s", flush=True)
