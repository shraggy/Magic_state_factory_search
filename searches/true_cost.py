"""True-cost comparison at l=3: our smallest factory for each target class versus the cheapest competitor route
built from [[8,3,2]] and synthesis, both fed by cultivated |T> states at error p.

Ours: one d=2 factory for the target (n inputs, accepted error c p^2, c = undetected weight-2 fault pairs).
Competitor: realise the target state from distilled |CCZ> states (8 inputs each, [[8,3,2]], error 28 p^2 each)
by Clifford circuits, measurements and gate teleportation, with any remaining T gates from the cheapest d=2 T
factory. Routes (the cheapest is taken):
  (a) T synthesis       T-count(target) T gates, each from Meier-Eastin-Knill (10 inputs -> 2 T, 5 per T)
  (b) CCZ/CS synthesis  write the target as a product of CCZ, CS and T terms in the cheapest CNOT frame found; each CCZ
                        costs one [[8,3,2]] (8), each CS one [[8,3,2]] + Y measurement (8), each T 5
The CS-from-CCZ conversion (Y measurement on one qubit of |CCZ>) is exact and keeps all 28 bad pairs (checked).
Error of a route = sum of its parts' errors (independent factories); T-gate errors for MEK are not included
9 p^2 per output T (their marginal error), so route (a) error = 9 * T-count.
Our c is recomputed from each factory's gates (pairs corrupting any output). Errors add across independent parts
(a union bound, as for any product of separately distilled states).

    python searches/true_cost.py --cells 3:3-5
"""
import os, sys, json, argparse
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search as us  # noqa: E402
from coset_certify import parse_target  # noqa: E402
from classify import essential_and_reduced  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
CCZ_N, CCZ_C = 8, 28                   # [[8,3,2]]: inputs, bad pairs (one CCZ)
CS_N, CS_C = 8, 28                     # [[8,3,2]] + Y measurement (one CS)
T_N = 5                                # Meier-Eastin-Knill (arXiv:1204.4221): 10 inputs -> 2 T
T_C = 9                                # its marginal error per output state, 9 p^2 (both outputs: 13 p^2)


def _terms_cost(poly):
    """Inputs and error of building a phase polynomial term by term: each CCZ, CS, T term is one state."""
    cnt = {1: 0, 2: 0, 3: 0}
    for S, c in poly.items():
        if c % (1 << len(S)):                         # genuine term (CCZ = 4, CS = 2 or 6, T = odd, mod 8)
            cnt[len(S)] += 1
    return cnt[3] * CCZ_N + cnt[2] * CS_N + cnt[1] * T_N, cnt[3] * CCZ_C + cnt[2] * CS_C + cnt[1] * T_C, cnt


def route_b(P, k, frames=4000, seed=7):
    """Cheapest term-by-term synthesis from CCZ, CS and T states over Clifford frames: the target as given, its
    classifier-reduced frame, and `frames` random CNOT frames (exact over all frames for k <= 4)."""
    import random
    from classify import _phase_table, _apply_linear, _tab_to_poly, _gl_matrices
    qs = list(range(k)); mod = 8
    tab = _phase_table(P, qs, mod)
    cands = [P, essential_and_reduced(P, qs, 3)[2]]
    if k <= 4:
        Ms = _gl_matrices(k)
    else:
        rnd = random.Random(seed); Ms = []
        while len(Ms) < frames:
            cols = [rnd.randrange(1, 1 << k) for _ in range(k)]
            basis = []
            for v in cols:
                w = v
                for b in basis:
                    w = min(w, w ^ b)
                if w:
                    basis.append(w)
            if len(basis) == k:
                Ms.append(tuple(cols))
    for M in Ms:
        cands.append(_tab_to_poly(_apply_linear(tab, qs, M), qs, mod))
    best = min((_terms_cost(p) for p in cands), key=lambda x: (x[0], x[1]))
    return best[0], best[1], best[2], k <= 4


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", nargs="+", default=["3:3-5"])
    ap.add_argument("--cascade-dir", default=os.path.join(ROOT, "outputs", "cascade"))
    a = ap.parse_args()
    out, wins = [], {"inputs": 0, "both": 0, "total": 0}
    for cell in a.cells:
        L, ks = cell.split(":"); L = int(L)
        lo, hi = map(int, ks.split("-"))
        for k in range(lo, hi + 1):
            d = json.load(open(os.path.join(a.cascade_dir, f"l{L}_k{k}.json")))
            rows = {str(r["orbit"]): r for r in d["rows"] if r.get("gates")}
            for orb, b in d["summary"]["best"].items():
                P = {T: c for T, c in parse_target(b["spec"], L).items() if c}
                tcount = us.zcount(P, k, L)
                nb, cb, cnt, exact = route_b(P, k)
                na = tcount * T_N if tcount is not None else None
                ca = tcount * T_C if tcount is not None else None
                comp = min([(na, ca, "T synth")] if na is not None else [] + [(nb, cb, "CCZ/CS synth")],
                           key=lambda x: (x[0], x[1]))
                comp = min([x for x in ((na, ca, "T synth"), (nb, cb, "CCZ/CS synth")) if x[0] is not None],
                           key=lambda x: (x[0], x[1]))
                comp_n, comp_c, comp_how = comp
                from benchmark_factories import fault_stats
                ours_c = None
                r_ = rows.get(orb)
                if r_ is not None:                    # SAT factory: circuit stored
                    ours_c = fault_stats([(frozenset(S), c) for S, c in r_["gates"]], k, L)["c"]
                else:                                 # stage (i)-(ii): rebuild the circuit from its parameters
                    import merge_large_unified as M
                    for x in d["rows"]:
                        if str(x["orbit"]) == orb and x["n"] == b["n"] and x["N"] == b["N"] and x["stage"] == b["stage"]:
                            g = M.regenerate_circuit(L, k, x["N"], x["params"], x["stage"], "")
                            if g:
                                ours_c = fault_stats(g, k, L)["c"]
                            break
                rec = dict(k=k, orbit=orb, spec=b["spec"], ours_n=b["n"], ours_N=b["N"], stage=b["stage"],
                           tcount=tcount, route_a_n=na, route_b_n=nb, route_b_c=cb, route_b_terms=cnt,
                           frame_exact=exact, competitor_n=comp_n, competitor_c=comp_c, competitor_route=comp_how,
                           ours_c=ours_c)
                out.append(rec); wins["total"] += 1
                if b["n"] < comp_n:
                    wins["inputs"] += 1
    json.dump(out, open(os.path.join(ROOT, "outputs", "true_cost_l3.json"), "w"), indent=1)
    both = sum(1 for r in out if r["ours_c"] is not None and r["ours_n"] < r["competitor_n"]
               and r["ours_c"] <= r["competitor_c"])
    print(f"{wins['total']} target classes; ours uses fewer inputs than the cheapest competitor route in "
          f"{wins['inputs']}, and fewer inputs with no larger error in {both}")
    for r in sorted(out, key=lambda r: r["ours_n"] - r["competitor_n"])[:15]:
        print(f"  k={r['k']} {r['spec'][:34]:34s} ours n={r['ours_n']:3d} c={r['ours_c']} | competitor "
              f"n={r['competitor_n']} c={r['competitor_c']} ({r['competitor_route']}) | T-count {r['tcount']}")


# ----------------------------------------------------------------- teleportation and its Clifford correction
def correction_cz(P, k, L=3):
    """Consuming |D> = D|+>^k teleports the diagonal D: CNOT data -> state, measure the state qubits in Z (outcome m),
    then correct with D(x) D(x+m)^{-1}. At l=3 that residual is a diagonal Clifford (checked for every m); return the
    worst-case and average number of CZ gates in it (its size-2 Newton terms with odd multiple of 4), over m != 0."""
    mod = 1 << L
    f = lambda x: sum(c for T, c in P.items() if all(x >> q & 1 for q in T)) % mod
    worst, tot, cnt = 0, 0, 0
    for m in range(1, 1 << k):
        a = us.poly_from_table([(f(x ^ m) - f(x)) % mod for x in range(1 << k)], k, mod)
        assert not any(c % (1 << us.cliff_exp(len(T), L)) for T, c in a.items()), "correction needs magic"
        cz = sum(1 for T, c in a.items() if len(T) == 2 and c % 8 == 4)
        worst, tot, cnt = max(worst, cz), tot + cz, cnt + 1
    return worst, tot / cnt


def teleport_cost(P, k):
    """(teleportations, CNOTs, worst-case correction CZs, mean correction CZs) for consuming one state |P>."""
    w, avg = correction_cz(P, k)
    return 1, k, w, avg
