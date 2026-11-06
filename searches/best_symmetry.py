"""Symmetry of the best factory per l=3 class (k <= 6): the group of qubit permutations that maps the factory's set of
rotations to itself and keeps outputs and checks in their roles (outputs among outputs, checks among checks).
For each factory: |Aut|, its action on the outputs (G_O) and on the checks (G_S), and whether Aut = G_O x G_S
(a product, the form the group-symmetric stage searches) or not (a joint symmetry it cannot express), and whether
the check part is fully symmetric (G_S = S_r, the form of stage (i)). Exact: backtracking over permutations with
pruning by each qubit's (role, degree, weight profile). Read-only on existing factories.
    python searches/best_symmetry.py --jobs 6   -> outputs/checks/best_symmetry_l3.json
"""
import os, sys, json, csv, argparse, collections
from itertools import permutations
from multiprocessing import Pool
from math import factorial
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")


def aut(k, N, gates):
    G = {frozenset(S) for S, c in gates if c % 2}            # rotations that consume inputs (odd c at l=3)
    prof = []
    for q in range(N):
        prof.append((q < k, tuple(sorted(len(S) for S in G if q in S))))
    perms = []

    def ext(p, used):
        q = len(p)
        if q == N:
            if {frozenset(p[x] for x in S) for S in G} == G:
                perms.append(tuple(p))
            return
        for t in range(N):
            if t not in used and prof[t] == prof[q] and (q < k) == (t < k):
                ok = True
                for S in G:                                   # partial check on rotations inside the mapped prefix
                    if max(S) == q:
                        if frozenset(p[x] if x < q else t for x in S) not in G:
                            ok = False; break
                if ok:
                    ext(p + [t], used | {t})
    ext([], frozenset())
    return perms


def one(row):
    k, N = int(row["k"]), int(row["N"]); gates = json.loads(row["gates"])
    P = aut(k, N, gates)
    GO = {p[:k] for p in P}; GS = {tuple(x - k for x in p[k:]) for p in P}
    r = N - k
    return dict(k=k, spec=row["spec"], n=int(row["n"]), N=N, source=row["source"], aut=len(P), G_O=len(GO), G_S=len(GS),
                product=len(P) == len(GO) * len(GS), checks_full=len(GS) == factorial(r))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=6); a = ap.parse_args()
    rows = list(csv.DictReader(open(os.path.join(O, "l3_catalogue_k6.csv"))))
    with Pool(a.jobs) as p:
        res = p.map(one, rows, chunksize=4)
    c = collections.Counter((x["checks_full"], x["G_O"] > 1, x["product"]) for x in res)
    by_src = collections.defaultdict(collections.Counter)
    for x in res:
        by_src[x["source"]][("checks S_r" if x["checks_full"] else "checks partial") + (", G_O>1" if x["G_O"] > 1 else ", G_O=1")
                            + ("" if x["product"] else ", joint")] += 1
    summ = dict(total=len(res), by_form={f"checks_full={a},G_O>1={b},product={c_}": n for (a, b, c_), n in c.items()},
                by_source={s: dict(v) for s, v in by_src.items()})
    json.dump(dict(summary=summ, rows=res), open(os.path.join(O, "checks", "best_symmetry_l3.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))
