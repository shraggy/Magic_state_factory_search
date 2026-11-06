"""Exact T-count decisions at k = 7 for the k78 pipeline pass: for each k = 7 target in a k78_eg.py run whose isotropy
floor is below its best synthesis w, decide with tcount7.exists_below2 (Brouwer-Zimmermann over two information sets,
compiled inner loop csrc/capped_minweight; complete when the returned bound reaches w) whether a coset word of weight
< w exists. 'exact' = no such word and the proved bound >= w: tau = w. Known-answer check first: the layers of
outputs/k7_candidates/exact2_20_22.json (no compilation below 17) and the k <= 6 self-test of tcount7.py.
    python searches/k78_exact7.py RUNDIR    -> outputs/k78_pipeline/exact7_<basename RUNDIR>.json
"""
import os, sys, json, glob, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import tcount7 as T
from k78_eg import _bin
_bin(T)
run = sys.argv[1]; out = []; seen = set()
for f in sorted(glob.glob(os.path.join(run, "c*.json"))):
    for r in json.load(open(f))["results"]:
        key = json.dumps(sorted(map(sorted, r["ccz"])))
        if r["k"] != 7 or key in seen:
            continue
        seen.add(key)
        from isotropy_floor import floor
        lo = floor(r["ccz"], 7); w = r["tcount_hi"]
        rec = dict(id=r.get("id"), fp=r.get("fp"), ccz=r["ccz"], floor=lo, synthesis=w, factory_n=r["best"] and r["best"]["n"])
        if lo < w:
            b, C0, nc = T.coset({frozenset(x): 4 for x in r["ccz"]}, 7); t = time.process_time()
            found, minseen, proved = T.exists_below2(b, C0, nc, w, 1)
            rec.update(found_below=found, min_seen=minseen, proved_bound=proved, cpu_s=round(time.process_time() - t, 1))
            rec["tau"] = w if (not found and proved is not None and proved >= w) else None
        else:
            rec["tau"] = w
        rec["on_bound"] = bool(rec["tau"] is not None and rec["factory_n"] == rec["tau"] + 1)
        out.append(rec); print(rec["id"], rec["fp"], "floor", lo, "synthesis", w, "tau", rec["tau"], "n", rec["factory_n"], flush=True)
json.dump(out, open(os.path.join(ROOT, "outputs", "k78_pipeline", f"exact7_{os.path.basename(run.rstrip('/'))}.json"), "w"), indent=1)
print(len(out), "k=7 layers;", sum(r["tau"] is not None for r in out), "with exact tau;", sum(r["on_bound"] for r in out), "on the bound")
