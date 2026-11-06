"""Stage (ii) gluing with disjoint blocks (cascade.stage_gluing_multi: partitions of k into >= 2 symmetric blocks of
size <= l, r = 1..rmax checks, 2-adic solve, the smallest r per partition) at l = 3, k = 7, 8, 9, as its own pass of
sweep B. The same function and rmax = 4 as sweep A's ii_dj (staged.DEFAULT) and the larger sweep (large_unified.py
--rmax 4, which ran us.stage_gluing, the same loop including the one-block partition that stage (i) covers).
Every factory is verified: us.is_factory (monomial form; brute force on all 2^N states when N <= 14), us.distance,
and that it deposits its own phase polynomial (entangled_gluing.deposits against the reduced output); labelled by the
Clifford-reduced form (us.classify_factory, as sweep B and the larger sweep) and, when the output is a pure CCZ layer,
by the fingerprint of searches/k78_orbits.py (k <= 8). Cross-check: each row of outputs/large_unified_factories.csv
with search_stage 'gluing' at l=3, k=7..9 must be reproduced with the same (blocks, n, N).
    python searches/k789_ii_dj.py K OUTDIR      -> OUTDIR/l3_k{K}/ii_dj.json (staged.py schema, keyed by reduced form)
"""
import os, sys, csv, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
csv.field_size_limit(10 ** 9)
import unified_search as us
import cascade as CA
from classify import phase_polynomial
from entangled_gluing import deposits

if __name__ == "__main__":
    k, out = int(sys.argv[1]), sys.argv[2]; L = 3; RMAX = 4
    t0 = time.process_time(); best, rows = {}, []
    for f in CA.stage_gluing_multi(L, k, RMAX):
        g = [(frozenset(S), c % 8) for S, c in f["gates"]]; N = k + f["r"]; f = dict(f, gates=g)
        ok = us.is_factory(f, L) and us.distance(f, L) != 1
        bf = None
        if N <= 14:
            bf = bool(us.is_factory(f, L, brute_force=True)); ok = ok and bf
        P = phase_polynomial([S for S, _ in g], [c for _, c in g], set(range(k)), L)
        ok = ok and deposits(g, k, {T: v for T, v in P.items() if len(T) <= 3})
        cl = us.classify_factory(f, L) if ok else None
        n = sum(us.inputs(c, L) for _, c in g)
        tri = sorted(sorted(T) for T, v in P.items() if len(T) == 3 and v % 8 == 4)
        pure = all(v % 8 == 0 or (len(T) == 3 and v % 8 == 4) or (len(T) == 2 and v % 4 == 0) or (len(T) == 1 and v % 2 == 0) for T, v in P.items())
        row = dict(params=f["params"], n=n, N=N, verified=bool(ok), brute_force=bf, distance=us.distance(f, L),
                   reduced=cl and cl["reduced"], essential_dim=cl and cl["essential_dim"], ccz_layer=bool(pure and tri))
        rows.append(row)
        if ok and cl and cl["essential_dim"] == k:
            o = cl["reduced"]
            if o not in best or (n, N) < (best[o]["n"], best[o]["N"]):
                best[o] = dict(n=n, N=N, stage="ii_dj:" + f["params"], gates=[[sorted(S), c] for S, c in g])
    secs = round(time.process_time() - t0, 1)
    # cross-check against the larger sweep
    ref = [r for r in csv.DictReader(open(os.path.join(ROOT, "outputs", "large_unified_factories.csv")))
           if r["l"] == "3" and r["k"] == str(k) and r["search_stage"] == "gluing"]
    mine = {(r["params"], r["n"], r["N"]) for r in rows if r["verified"]}
    xc = [dict(params=r["params"], n=int(r["n"]), N=int(r["N"]), reproduced=(r["params"], int(r["n"]), int(r["N"])) in mine) for r in ref]
    rec = dict(L=L, k=k, pass_="ii_dj", settings=dict(rmax=RMAX, function="cascade.stage_gluing_multi"), seconds=secs,
               seconds_are_cpu=True, reached=len(best), best=best, rows=rows,
               crosscheck_large_unified=dict(rows=len(xc), reproduced=sum(x["reproduced"] for x in xc), detail=xc),
               ccz_layer_outputs=sum(r["ccz_layer"] for r in rows))
    d = os.path.join(out, f"l3_k{k}"); os.makedirs(d, exist_ok=True)
    json.dump(rec, open(os.path.join(d, "ii_dj.json"), "w"))
    print(f"k={k}: {len(rows)} partitions solved, {sum(r['verified'] for r in rows)} verified, {len(best)} reduced outputs, "
          f"{rec['ccz_layer_outputs']} CCZ layers; large_unified gluing rows reproduced {rec['crosscheck_large_unified']['reproduced']}"
          f"/{len(xc)}; {secs} CPU-s", flush=True)
