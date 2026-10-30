"""Close the open steps of a quick stage (iv) run with longer, uncapped SAT proofs.

For every target in outputs/targets/quick_sat_l{L}_k{k}.json whose step at N ended 'unknown' (time limit), ask once,
with no gate-size cap and a longer limit, whether a factory with fewer inputs than the best found exists on N qubits.
UNSAT proves the reported n smallest on that many qubits; SAT gives a smaller factory. Only bounds n <= --nmax are
attempted (the exhaustive range of the known classifications).

    python searches/prove_quick.py --level 3 --k 4 --N 7 --limit 3600 --jobs 4
"""
import os, sys, json, time, argparse
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
TD = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs", "targets")


def prove(args):
    from coset_certify import parse_target
    from sat_certify import certify
    L, k, N, spec, bound, limit = args
    t = time.time()
    r, nf, g = certify(parse_target(spec, L), k, N, L, bound, limit)
    return dict(spec=spec, N=N, bound=bound, result=r, smaller_n=nf,
                gates=[[sorted(S), c] for S, c in g] if g else None, seconds=round(time.time() - t, 1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--N", type=int, required=True)
    ap.add_argument("--limit", type=float, default=3600.0)
    ap.add_argument("--nmax", type=int, default=54)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    path = os.path.join(TD, f"quick_sat_l{a.level}_k{a.k}.json")
    d = json.load(open(path))
    todo = []
    for t in d["targets"]:
        step = next((s for s in t["steps"] if s["N"] == a.N), None)
        if t["best"] and step and step["last"] == "unknown" and t["best"]["n"] <= a.nmax:
            todo.append((a.level, a.k, a.N, t["spec"], t["best"]["n"], a.limit))
    print(f"l={a.level} k={a.k} N={a.N}: {len(todo)} open steps, limit {a.limit:.0f}s each", flush=True)
    t0 = time.time()
    with ProcessPoolExecutor(a.jobs) as ex:
        res = []
        for r in ex.map(prove, todo):
            res.append(r)
            print(f"  {r['spec']}: n<{r['bound']} on N={r['N']}: {r['result']}"
                  + (f" -> n={r['smaller_n']}" if r["smaller_n"] else "") + f"  ({r['seconds']}s)", flush=True)
    out = os.path.join(TD, f"prove_l{a.level}_k{a.k}_N{a.N}.json")
    json.dump(dict(l=a.level, k=a.k, N=a.N, limit=a.limit, wall_seconds=round(time.time() - t0, 1), results=res),
              open(out, "w"), indent=1)
    c = {x: sum(r["result"] == x for r in res) for x in ("unsat", "sat", "unknown")}
    print(f"done: {c['unsat']} proved smallest, {c['sat']} smaller found, {c['unknown']} still open, "
          f"{time.time() - t0:.0f}s wall", flush=True)
