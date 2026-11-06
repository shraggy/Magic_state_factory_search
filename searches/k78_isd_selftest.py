"""Known-answer check of the information-set decoder used to tighten T-count upper bounds in the k78 pass
(k78_tau_upper.isd, p = 3, seed 0, 200 iterations): on random CCZ layers at k = 5, 6 with exact T-count e from
classification/tcount.tcount_rm, every word it returns must lie in the coset (in_coset) and have weight >= e; we also
count how often it reaches e.   -> outputs/k78_pipeline/isd_selftest.json
"""
import os, sys, json, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import tcount7 as T, tcount as TC
from k78_tau_upper import isd, in_coset
rng = random.Random(9); rows = []
for _ in range(30):
    k = rng.randint(5, 6); tri = [t for t in itertools.combinations(range(k), 3) if rng.random() < 0.4]
    if not tri:
        continue
    e = TC.tcount_rm({frozenset(q + 1 for q in t): 4 for t in tri}, k); b, C0, nc = T.coset({frozenset(t): 4 for t in tri}, k)
    w, seen, its = isd(b, C0, nc, e, 200, 3, seed=0)
    rows.append(dict(k=k, ccz=tri, exact=e, found=seen if w is not None else None, in_coset=bool(w is not None and in_coset(w, b, C0)), iters=its))
bad = sum(r["found"] is not None and (r["found"] < r["exact"] or not r["in_coset"]) for r in rows)
reach = sum(r["found"] == r["exact"] for r in rows)
json.dump(dict(total=len(rows), reached_exact=reach, invalid=bad, rows=rows), open(os.path.join(ROOT, "outputs", "k78_pipeline", "isd_selftest.json"), "w"), indent=0)
print(f"{len(rows)} layers: reached the exact T-count on {reach} of {len(rows)}; invalid words {bad}")
