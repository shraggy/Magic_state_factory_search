"""Wall-clock time of each unified-search stage as the parameter space grows.

Runs unified_search.run() at one level for growing (k_max, N_max) and records the seconds spent in each stage
(two-group, gluing, padding, optional classification) and the number of factories found, so the trade-off between
stages can be plotted. Shared-check gluing is given at most N_max - k_max check qubits, as the two-group stage is.
Writes outputs/stage_scaling.json; each completed point is saved at once, so an interrupted run keeps its results.
Run: python3 searches/stage_scaling.py [l] [N_max_max]
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import unified_search as us
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "outputs", "stage_scaling.json")
L = int(sys.argv[1]) if len(sys.argv) > 1 else 3
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 16
res = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(p["l"], p["k_max"], p["N_max"]) for p in res}
for N in range(6, NMAX + 1):
    for k in (4, 6, 8):
        if k >= N or (L, k, N) in done:
            continue
        T = {}
        t = time.time()
        rows = us.run(L, k, N, min(4, N - k), 2, log=lambda *a, **kw: None, timing=T)
        T["total"] = time.time() - t
        stages = {}
        for r in rows:
            s = "padding" if r["stage"].endswith("+relabel") else r["stage"]
            stages[s] = stages.get(s, 0) + 1
        res.append(dict(l=L, k_max=k, N_max=N, seconds=T, factories=len(rows), per_stage=stages))
        json.dump(res, open(OUT, "w"), indent=1)
        print(f"l={L} k<={k} N<={N}: " + ", ".join(f"{a} {b:.1f}s" for a, b in T.items()) + f"; {len(rows)} factories", flush=True)
