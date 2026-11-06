"""Known-answer check of the k=7 exact decision used by the k78 pass: tcount7.exists_below2 (one thread) against the
exact Reed-Muller decoder classification/tcount.tcount_rm on random CCZ layers at k <= 6: for each layer with exact
T-count e, exists_below2(cap = e) must report no word (and a proved bound >= e when complete) and exists_below2(cap =
e + 1) must find one.   -> outputs/k78_pipeline/tcount7_selftest.json
"""
import os, sys, json, random, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import tcount7 as T, tcount as TC
from k78_eg import _bin
_bin(T)
rng = random.Random(5); rows = []
for _ in range(40):
    k = rng.randint(4, 6); tri = [t for t in itertools.combinations(range(k), 3) if rng.random() < 0.35]
    if not tri:
        continue
    e = TC.tcount_rm({frozenset(q + 1 for q in t): 4 for t in tri}, k)
    b, C0, nc = T.coset({frozenset(t): 4 for t in tri}, k)
    f0, s0, p0 = T.exists_below2(b, C0, nc, e, 1); f1, s1, p1 = T.exists_below2(b, C0, nc, e + 1, 1)
    rows.append(dict(k=k, ccz=tri, exact=e, below_found=f0, below_bound=p0, at_found=f1, ok=(not f0) and f1))
ok = sum(r["ok"] for r in rows); complete = sum(r["below_bound"] is not None and r["below_bound"] >= r["exact"] for r in rows)
json.dump(dict(agree=ok, total=len(rows), complete_proofs=complete, rows=rows), open(os.path.join(ROOT, "outputs", "k78_pipeline", "tcount7_selftest.json"), "w"), indent=0)
print(f"{ok} of {len(rows)} agree with the exact decoder; {complete} of {len(rows)} with a complete proof below")
