"""Check-pattern SAT (pattern_sat, current code; z3 rlimit 2e9 per query, R <= 4, mmax 12, Lempel sizes mu and mu+1,
as staged.DEFAULT iii_r) with every z3 query in a fresh spawned process, seeded per class with the best n of
SEEDDIR/l3_k{k}/<SEEDPASS>.json (default ii_eg), or with --cap for a single class.
    python searches/pat_eg.py k I M SEEDDIR OUT [--seedpass ii_eg]   -> OUT/l3_k{k}/pat_<seedpass>/c{I}of{M}.json
    python searches/pat_eg.py --one k SPEC CAP                         (known-answer check, prints the result)
Classes sorted by label, index % M == I; classes at the bound n = tau+1 are skipped. Every factory is verified
(is_factory, distance != 1) and labelled with the exact class (staged.label_fix).
"""
import os, sys, json, time, multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))


def _q(args, q):
    import pattern_sat as PS
    q.put(PS.solve(*args))


def solve_fresh(spec, k, R, mu, m, rlimit=2e9):
    ctx = mp.get_context("spawn"); q = ctx.Queue()
    p = ctx.Process(target=_q, args=((spec, k, R, mu, m, rlimit), q)); p.start(); r = q.get(); p.join()
    return r


def run_class(k, spec, cap, Rs=(1, 2, 3, 4), mmax=12, rlimit=2e9):
    """pattern_sat.run_class with each query in its own process (same loop, same order)."""
    import pattern_sat as PS, synthillation_compare as sc, class_table_l3 as CT
    P = CT.parse(spec); mu = 0
    while sc.opt_B(k, P, mu) is None:
        mu += 1
    best = None; log = []
    for mu_ in (mu, mu + 1):
        cols, _ = PS.multiset(spec, k, mu_)
        if cols is None:
            continue
        for R in Rs:
            for m in range(0, mmax + 1):
                n = len(cols) + m
                if n >= (best["n"] if best else cap):
                    break
                r = solve_fresh(spec, k, R, mu_, m, rlimit)
                log.append([mu_, R, m, r if not isinstance(r, dict) else "sat"])
                if isinstance(r, dict) and r["ok"]:
                    if best is None or (r["n"], r["N"]) < (best["n"], best["N"]):
                        best = r
                    break
    return best, log


def run(k, I, M, seeddir, out, seedpass="ii_eg"):
    import staged as STG, staged_class as SC
    L = 3; cache, classes = SC.setup(L, k, seeddir)
    seed = json.load(open(os.path.join(seeddir, f"l{L}_k{k}", seedpass + ".json")))["best"]
    s = STG.DEFAULT["iii_r"]
    mine = [o for i, o in enumerate(sorted(classes)) if i % M == I]; res, per = [], {}
    for o in mine:
        spec = classes[o]["spec"]; cap = seed[o]["n"] if o in seed else 10 ** 6; tau = STG.tau_of(L, k, o, spec)
        if tau is not None and cap == tau + 1:
            per[o] = dict(seconds=0.0, skipped="at bound", seed=cap); continue
        t = time.time()
        b, log = run_class(k, spec, cap, Rs=tuple(range(1, s["R"] + 1)), mmax=s["mmax"], rlimit=s["pat_rlimit"])
        if b and b["ok"]:
            o2, gg = STG.label_fix(dict(k=k, r=b["N"] - k, gates=b["gates"]), L, k, cache)
            if o2 is not None:
                res.append(dict(orbit=o2, n=b["n"], N=b["N"], stage="iii_r:pattern:" + seedpass,
                                gates=[[sorted(S_), c] for S_, c in gg]))
        per[o] = dict(seconds=round(time.time() - t, 1), seed=cap, queries=len(log))
        print(f"pat k={k} {o} {spec}: {per[o]} -> {b and b['n']}", flush=True)
    SC.save(os.path.join(out, f"l{L}_k{k}", "pat_" + seedpass, f"c{I}of{M}.json"),
            dict(L=L, k=k, chunk=[I, M], classes=mine, host=os.uname().nodename, per_class=per, results=res))


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "--one":
        b, log = run_class(int(a[1]), a[2], int(a[3]))
        print("best", b and (b["n"], b["N"], b["ok"]), "queries", len(log), flush=True)
    else:
        sp = "ii_eg"
        if "--seedpass" in a:
            i = a.index("--seedpass"); sp = a[i + 1]; a = a[:i] + a[i + 2:]
        run(int(a[0]), int(a[1]), int(a[2]), a[3], a[4], sp)
