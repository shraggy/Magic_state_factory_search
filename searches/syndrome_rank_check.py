"""Syndrome-rank bound n >= tau(U) + r on every factory with explicit gates and exact T-count, and the proof's map.

r is the GF(2) rank of the check supports of the input rotations (odd coefficient at l=3, unified_search.inputs > 0).
The bound: choose r input rotations with independent check supports b_i and a linear A with A b_i = a_i (their output
parts); CNOTs from outputs to checks, (x, y) -> (x, y + A^T x), send every rotation (a, b) -> (a + A b, b), fix
|+>^N and the X-basis record of the checks, so the factory is the same, and the r chosen rotations touch no output.
Fixing the checks to y = 0 then leaves a synthesis of the target with at most n - r rotations, so tau <= n - r.

Factories checked:
  outputs/l3_catalogue_k6.csv            766 classes at l=3, k<=6 (gates on qubits 0..N-1, outputs 0..k-1, exact tcount)
  outputs/k78_pipeline/all_isd/*.json    every CCZ layer at k=7, 8 of the cascade on a target list whose T-count is
                                         decided (exact at k=7; isotropy floor = synthesis at k=8), with its gates
For each: r, the bound, equality; at n = tau+1, r = 1. On every factory the map is applied and the result checked:
the r chosen rotations are output-free, is_factory holds (monomial form and brute force on 2^N states), the output
deposited is unchanged (entangled_gluing.deposits), and at y = 0 at most n - r input rotations act on the outputs.
    nice -n 19 python searches/syndrome_rank_check.py   -> outputs/checks/syndrome_rank_check.json
"""
import os, sys, csv, json
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
csv.field_size_limit(sys.maxsize)
import unified_search as us  # noqa: E402
from entangled_gluing import deposits  # noqa: E402
from classify import phase_polynomial  # noqa: E402
L = 3
O = os.path.join(ROOT, "outputs")
import argparse as _ap
_p = _ap.ArgumentParser(); _p.add_argument("--catalogue", default=os.path.join(O, "l3_catalogue_k6.csv"))
_p.add_argument("--out-dir", default=None, help="write syndrome_rank_check.json here instead of outputs/checks")
_A = _p.parse_args()


def mask(S):
    return sum(1 << q for q in S)


def independent(vecs):
    """Indices of a maximal linearly independent subset of vecs (GF(2) bitmasks), greedy in order."""
    basis, idx = {}, []          # pivot bit -> reduced vector
    for i, v in enumerate(vecs):
        w = v
        while w:
            p = w.bit_length() - 1
            if p in basis:
                w ^= basis[p]
            else:
                basis[p] = w; idx.append(i); break
    return idx


def linear_map(bs, as_, s):
    """A: GF(2)^s -> GF(2)^k (images of the unit vectors) with A b_i = a_i for independent b_i: reduce the pairs
    (b_i, a_i) to reduced row-echelon form in b, set A e_p = a' on each pivot p and A e_q = 0 elsewhere."""
    rows = [[b, a] for b, a in zip(bs, as_)]
    piv = []
    for t in range(len(rows)):
        p = rows[t][0].bit_length() - 1
        for u in range(len(rows)):
            if u != t and rows[u][0] >> p & 1:
                rows[u][0] ^= rows[t][0]; rows[u][1] ^= rows[t][1]
        piv.append(p)
        rows[t + 1:] = sorted(rows[t + 1:], key=lambda r: -r[0])
    A = [0] * s
    for (b, a), p in zip(rows, piv):
        A[p] = a
    for b, a in zip(bs, as_):
        assert apply(A, b) == a
    return A


def apply(A, x):
    out = 0
    j = 0
    while x:
        if x & 1:
            out ^= A[j]
        x >>= 1; j += 1
    return out


def check(k, N, gates, tau, name):
    km = (1 << k) - 1
    s = N - k
    n = sum(us.inputs(c, L) for _, c in gates)
    inp = [(mask(S), c) for S, c in gates if us.inputs(c, L)]
    bvec = [m >> k for m, _ in inp]
    sel = independent(bvec)
    r = len(sel)
    rec = dict(name=name, k=k, N=N, n=n, tau=tau, r=r, bound_ok=n >= tau + r, equality=n == tau + r)
    # the proof's map on this factory
    A = linear_map([bvec[i] for i in sel], [inp[i][0] & km for i in sel], s)
    new = []
    for S, c in gates:
        m = mask(S); a, b = m & km, m >> k
        a2 = a ^ apply(A, b)
        m2 = a2 | (b << k)
        new.append((frozenset(q for q in range(N) if m2 >> q & 1), c))
    newinp = [(mask(S), c) for S, c in new if us.inputs(c, L)]
    chosen_free = all(newinp[i][0] & km == 0 for i in sel)
    target = phase_polynomial([frozenset(S) for S, _ in gates], [c for _, c in gates], set(range(k)), L)
    f0 = dict(k=k, r=s, gates=[(frozenset(S), c) for S, c in gates])
    f2 = dict(k=k, r=s, gates=new)
    rec.update(orig_factory=us.is_factory(f0, L) and us.is_factory(f0, L, brute_force=True),
               map_factory=us.is_factory(f2, L) and us.is_factory(f2, L, brute_force=True),
               map_deposit_same=deposits(new, k, target, L) and deposits([(frozenset(S), c) for S, c in gates], k, target, L),
               map_chosen_output_free=chosen_free,
               synthesis_size=sum(1 for m, _ in newinp if m & km))
    rec["synthesis_le_n_minus_r"] = rec["synthesis_size"] <= n - r
    return rec


def main():
    recs = []
    for row in csv.DictReader(open(_A.catalogue)):
        g = [(tuple(S), c) for S, c in json.loads(row["gates"])]
        recs.append(dict(src="l3_catalogue_k6", spec=row["spec"],
                         **check(int(row["k"]), int(row["N"]), g, int(row["tcount"]), row["spec"])))
    import glob
    for f in sorted(glob.glob(os.path.join(O, "k78_pipeline", "all_isd", "*.json"))):   # k = 7, 8 CCZ layers
        for t in json.load(open(f))["results"]:
            if t["tcount_lo"] != t["tcount_hi"]:
                continue
            b = t["best"]; k = t["k"]
            g = [(tuple(S), c) for S, c in b["gates"]]
            rec = check(k, b["N"], g, t["tcount_lo"], str(t["ccz"]))
            rec["deposits_layer"] = deposits([(frozenset(S), c) for S, c in g], k, {frozenset(x): 4 for x in t["ccz"]}, L)
            recs.append(dict(src="k78_pipeline", **rec))
    flags = ["bound_ok", "orig_factory", "map_factory", "map_deposit_same", "map_chosen_output_free",
             "synthesis_le_n_minus_r"]
    summ = {}
    for src in ("l3_catalogue_k6", "k78_pipeline"):
        R = [x for x in recs if x["src"] == src]
        summ[src] = dict(checked=len(R), rebuild_failed=sum(1 for x in R if x.get("rebuild_failed")),
                         violations=sum(1 for x in R if not x.get("bound_ok", False)),
                         equality_n_eq_tau_plus_r=sum(1 for x in R if x.get("equality")),
                         at_tau_plus_1=sum(1 for x in R if x.get("n") == x.get("tau", -9) + 1),
                         r_at_tau_plus_1=dict(Counter(x["r"] for x in R if x.get("n") == x.get("tau", -9) + 1)),
                         r_distribution=dict(sorted(Counter(x.get("r") for x in R).items(), key=str)),
                         **{f + "_fail": sum(1 for x in R if not x.get(f, False)) for f in flags})
        if src == "k78_pipeline":
            summ[src]["deposits_layer_fail"] = sum(1 for x in R if not x.get("deposits_layer", False))
    tot = dict(checked=len(recs), violations=sum(1 for x in recs if not x.get("bound_ok", False)),
               equality_n_eq_tau_plus_r=sum(1 for x in recs if x.get("equality")),
               at_tau_plus_1=sum(1 for x in recs if x.get("n") == x.get("tau", -9) + 1),
               r_at_tau_plus_1=dict(Counter(x["r"] for x in recs if x.get("n") == x.get("tau", -9) + 1)),
               map_all_ok=all(all(x.get(f, False) for f in flags) for x in recs))
    out = dict(summary=dict(total=tot, **summ),
               equality_cases=[{k_: x[k_] for k_ in ("src", "name", "k", "N", "n", "tau", "r")} for x in recs
                               if x.get("equality")],
               rows=recs)
    OD = _A.out_dir or os.path.join(O, "checks"); os.makedirs(OD, exist_ok=True)
    json.dump(out, open(os.path.join(OD, "syndrome_rank_check.json"), "w"), indent=1)
    print(json.dumps(out["summary"], indent=1))


if __name__ == "__main__":
    main()
