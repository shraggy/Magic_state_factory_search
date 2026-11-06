"""Stage (i) two-group and stage (ii) disjoint-gluing factories of the earlier larger sweep at k = 7, 8 (l = 3) whose
output is a CCZ layer: every row of outputs/large_unified_factories.csv at l=3, k=7,8, rebuilt with
merge_large_unified.regenerate_circuit, verified (us.is_factory, us.distance), its phase polynomial reduced mod
diagonal Cliffords and tested for being a pure CCZ layer with k essential outputs, then tagged with the
GL(k,2)-invariant fingerprint of searches/k78_orbits.py.  -> outputs/k78_pipeline/large_ccz.json
"""
import os, sys, csv, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
csv.field_size_limit(10 ** 9)
import unified_search as us
from classify import phase_polynomial
from merge_large_unified import regenerate_circuit
from k78_orbits import fingerprint, genuine
rows = [r for r in csv.DictReader(open(os.path.join(ROOT, "outputs", "large_unified_factories.csv"))) if r["l"] == "3" and r["k"] in ("7", "8")]
out, stats = [], dict(rows=len(rows), rebuilt=0, verified=0, ccz=0)
for r in rows:
    k = int(r["k"]); N = int(r["N"])
    if N > 14 and r["search_stage"] == "two-group":
        pass
    try:
        g = regenerate_circuit(3, k, N, r["params"], r["search_stage"], r.get("raw"))
    except Exception as e:
        g = None
    if not g:
        continue
    stats["rebuilt"] += 1
    f = dict(k=k, r=N - k, gates=g)
    if not (us.is_factory(f, 3) and us.distance(f, 3) != 1):
        continue
    stats["verified"] += 1
    P = phase_polynomial([S for S, _ in g], [c for _, c in g], set(range(k)), 3)
    pure = all((len(T) == 1 and v % 2 == 0) or (len(T) == 2 and v % 4 == 0) or (len(T) == 3 and v % 4 == 0) or len(T) != 3 and len(T) > 3
               for T, v in P.items() if v % 8)
    if not pure:
        continue
    tri = sorted(sorted(T) for T, v in P.items() if len(T) == 3 and v % 8 == 4)
    if not tri or not genuine(tri, k):
        continue
    stats["ccz"] += 1
    out.append(dict(k=k, n=int(r["n"]), N=N, stage=r["search_stage"], params=r["params"], ccz=tri, fp=fingerprint(tri, k)))
json.dump(dict(stats=stats, rows=out), open(os.path.join(ROOT, "outputs", "k78_pipeline", "large_ccz.json"), "w"))
print(stats, sorted((x["k"], x["n"], x["N"], x["stage"]) for x in out)[:30])
