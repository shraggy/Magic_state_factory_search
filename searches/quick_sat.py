"""Stage (iv) in the quick regime of the entangled-output SAT search (distance_2/entangled_sat): for every exact
target class of outputs/targets, try N = k+1 .. k+SPAN qubits with a short limit per query and gate supports capped
at WMAX qubits, tightening n while queries keep succeeding. No minimality proofs: a timeout or UNSAT under the cap
only means 'not found in this regime'. Uses the Clifford-corrected borrowed-identity condition of sat_certify.

    python searches/quick_sat.py --cells 3:3-5 --jobs 4
writes outputs/targets/quick_sat_l{L}_k{k}.json (per target: best n, N, wall-clock per query) and prints a summary
against the catalogue (outputs/targets/catalogue_match.json). Each query is seeded with the best n stages (i)-(ii)
reached for that class, so it only asks for what they did not find; --unseeded reproduces the unseeded run.
"""
import os, sys, json, time, argparse
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
TD = os.path.join(ROOT, "outputs", "targets")


def solve_target(args):
    from coset_certify import parse_target
    from sat_certify import certify
    L, k, spec, budget, span, wmax, seed = args
    tgt = parse_target(spec, L)
    # seed = best n stages (i)-(ii) reached for this class: every query then asks only for something they did not
    # find (a smaller factory), so stage (iv) never re-finds a factory an earlier stage already has
    t0, cap, best, steps = time.time(), seed or 10 ** 6, None, []
    for NN in range(k + 1, k + span + 1):
        # gate-support caps, tightest first: a tighter cap is a much smaller SAT instance and often reaches a lower n
        # quickly (CS12 T3 T4 T5 on N=8: cap 5 gives n=26 in 0.4 s, cap 6 stalls after n=28). Any factory found under
        # a cap is a real factory. The uncapped pass (None) comes last, so at N <= wmax+1 an UNSAT there is a proof.
        caps = [w for w in (wmax - 1, wmax) if w < NN] + [None]
        q, t1, last = 0, time.time(), None
        for w in caps:
            while True:
                r, nf, g = certify(tgt, k, NN, L, cap, budget, max_weight=w); q += 1
                if r == "sat":
                    best, cap = dict(n=nf, N=NN, gates=[[sorted(S), c] for S, c in g]), nf
                else:
                    break
            last = r                                  # outcome of the last (least capped) pass
        steps.append(dict(N=NN, queries=q, last=last, seconds=round(time.time() - t1, 2)))
    return dict(spec=spec, best=best, steps=steps, seconds=round(time.time() - t0, 2))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", nargs="+", default=["3:3-5"])
    ap.add_argument("--budget", type=float, default=25.0, help="seconds per SAT query (entangled_sat used 25)")
    ap.add_argument("--span", type=int, default=3, help="try N = k+1 .. k+span (entangled_sat used 3)")
    ap.add_argument("--wmax", type=int, default=6, help="gate-support cap once N exceeds it (entangled_sat used 6)")
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--unseeded", action="store_true", help="ignore the stages' best n (re-finds their factories)")
    a = ap.parse_args()
    match = json.load(open(os.path.join(TD, "catalogue_match.json")))
    for c in a.cells:
        L, ks = c.split(":"); L = int(L)
        lo, hi = map(int, ks.split("-"))
        for k in range(lo, hi + 1):
            cell = match[f"l{L}_k{k}"]
            specs = [t["spec"] for t in cell["targets"]]
            stage_best = {t["spec"]: t["best_stage"] for t in cell["targets"]}
            t0 = time.time()
            with ProcessPoolExecutor(a.jobs) as ex:
                seeds = [stage_best[s][0] if (stage_best[s] and not a.unseeded) else None for s in specs]
                res = list(ex.map(solve_target, [(L, k, s, a.budget, a.span, a.wmax, sd) for s, sd in zip(specs, seeds)]))
            wall = time.time() - t0
            found = [r for r in res if r["best"]]
            new = [r for r in found if stage_best[r["spec"]] is None]
            better = [r for r in found if stage_best[r["spec"]] and r["best"]["n"] < stage_best[r["spec"]][0]]
            json.dump(dict(l=L, k=k, budget=a.budget, span=a.span, wmax=a.wmax, wall_seconds=round(wall, 1),
                           cpu_seconds=round(sum(r["seconds"] for r in res), 1), targets=res,
                           stage_best=stage_best), open(os.path.join(TD, f"quick_sat_l{L}_k{k}.json"), "w"), indent=1)
            print(f"l={L} k={k}: {len(specs)} classes | stages {sum(v is not None for v in stage_best.values())} | "
                  f"quick SAT found {len(found)} ({len(new)} the stages never reached, {len(better)} smaller than "
                  f"the stages) | {wall:.0f}s wall, {sum(r['seconds'] for r in res):.0f}s cpu", flush=True)
            for r in new + better:
                b, sb = r["best"], stage_best[r["spec"]]
                print(f"   {r['spec']}: n={b['n']} N={b['N']}" + (f"  (stages n={sb[0]} N={sb[1]})" if sb else ""), flush=True)
