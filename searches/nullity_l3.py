"""Stabilizer nullity of every l=3 target class (k <= 6), and the lower bound it gives on any route that distils
known states and converts them with Clifford operations and Pauli measurements.

The stabilizer nullity nu(|psi>) = k - log2 |Stab(psi)|, with Stab the Pauli operators (up to sign) that fix |psi>,
cannot increase under stabilizer operations, including Pauli measurements and classical feed-forward (Beverland,
Campbell, Howard, Kliuchnikov, arXiv:1904.01124), and is additive on tensor products. A route that prepares the
target |U> from distilled states therefore needs states of total nullity >= nu(|U>). Resources and their costs in
input T states: a CCZ state from a batch factory [[6m+2,3m,2]] (nu = 3 each, 6m+2 inputs for m of them; [[8,3,2]]
for m = 1), a T state from catalysed CCZ -> 2T (nu = 1, 4 inputs each), a CS state (nu = 2, 8 inputs). LB(nu) is the
cheapest combination of total nullity >= nu. This is a lower bound on every such route from these resources: a factory with n < LB(nu)
beats all of them.
    python searches/nullity_l3.py      -> outputs/nullity_l3.json
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)


def parse(spec):
    P = {}
    for g in spec.split(","):
        name = g.rstrip("0123456789"); qs = frozenset(int(c) for c in g[len(name):])
        P[qs] = (P.get(qs, 0) + {"T": 1, "CS": 2, "CCZ": 4}[name]) % 8
    return P


def state(P, k):
    x = np.arange(1 << k); f = np.zeros(1 << k, dtype=np.int64)
    for T, c in P.items():
        m = sum(1 << q for q in T); f += c * ((x & m) == m)
    return np.exp(1j * np.pi / 4 * (f % 8)) / 2 ** (k / 2)


def nullity(psi, k):
    x = np.arange(1 << k); par = np.array([[bin(v & b).count("1") & 1 for v in x] for b in range(1 << k)])
    cnt = 0
    for a in range(1 << k):
        pa = psi[x ^ a]
        for b in range(1 << k):
            if abs(abs(np.vdot(psi, (1 - 2 * par[b]) * pa)) - 1) < 1e-9:
                cnt += 1
    return k - int(round(np.log2(cnt)))


def LB(nu):
    best = 10 ** 9
    for m in range(0, nu + 2):
        for t in range(0, nu + 1):
            for s in range(0, nu + 1):
                if 3 * m + t + 2 * s >= nu:
                    best = min(best, (6 * m + 2 if m else 0) + 4 * t + 8 * s)
    return best


if __name__ == "__main__":
    for spec, k, want in (("CCZ012", 3, 3), ("CS01", 2, 2), ("T0", 1, 1)):      # Beverland et al.: nu(CCZ)=3, nu(CS)=2, nu(T)=1
        assert nullity(state(parse(spec), k), k) == want, spec
    R = [dict(r, tcount=r["tau"]) for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]
         if r["l"] == 3]                                   # smallest factory per class, staged cascade (sweep A)
    out = []
    for r in R:
        nu = nullity(state(parse(r["spec"]), r["k"]), r["k"])
        out.append(dict(k=r["k"], spec=r["spec"], tcount=r["tcount"], nu=nu, n=r["n"], route_lower_bound=LB(nu),
                        beats_every_route=r["n"] < LB(nu)))
    json.dump(out, open(os.path.join(ROOT, "outputs", "nullity_l3.json"), "w"), indent=0)
    print("LB(nu) for nu = 1..6:", [LB(v) for v in range(1, 7)])
    print("nu == k for", sum(o["nu"] == o["k"] for o in out), "of", len(out))
    print("factory beats every distil-and-convert route:", [(o["spec"], o["n"], o["route_lower_bound"]) for o in out if o["beats_every_route"]])
