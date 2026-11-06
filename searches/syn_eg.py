"""Synthesis SAT (z3, synthesis_sat.solve, sweep A's settings: span N=k+1..k+3, step limit 3e9) seeded with the
gluing n (best after pass ii_eg) instead of the best after check-pattern SAT, per class at l=3.
    python searches/syn_eg.py k I M SEEDDIR OUT   -> OUT/l3_k{k}/syn_eg/c{I}of{M}.json
Classes sorted by label, index % M == I; classes at the bound n = tau+1 are skipped. Every factory is verified inside
synthesis_sat and labelled with the exact class (staged.label_fix).
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import staged as STG, staged_class as SC


def run(k, I, M, seeddir, out):
    L = 3; cache, classes = SC.setup(L, k, seeddir)
    seed = json.load(open(os.path.join(seeddir, f"l{L}_k{k}", "ii_eg.json")))["best"]
    s = STG.DEFAULT["iii_r"]; import synthesis_sat as SS
    mine = [o for i, o in enumerate(sorted(classes)) if i % M == I]; res, per = [], {}
    for o in mine:
        spec = classes[o]["spec"]; cap = seed[o]["n"] if o in seed else 10 ** 6; tau = STG.tau_of(L, k, o, spec)
        if tau is not None and cap == tau + 1:
            per[o] = dict(seconds=0.0, skipped="at bound"); continue
        t = time.time()
        b, log = SS.solve(spec, k, cap, rl=s["syn_rlimit"], span=tuple(range(1, s["syn_span"] + 1)))
        if b:
            o2, gg = STG.label_fix(dict(k=k, r=b["N"] - k, gates=b["gates"]), L, k, cache)
            if o2 is not None:
                res.append(dict(orbit=o2, n=b["n"], N=b["N"], stage="iii_r:synthesis:eg",
                                gates=[[sorted(S_), c] for S_, c in gg]))
        per[o] = dict(seconds=round(time.time() - t, 1), seed=cap, log=log)
        print(f"syn_eg k={k} {o} {spec}: {per[o]} -> {b and b['n']}", flush=True)
    SC.save(os.path.join(out, f"l{L}_k{k}", "syn_eg", f"c{I}of{M}.json"),
            dict(L=L, k=k, chunk=[I, M], classes=mine, host=os.uname().nodename, per_class=per, results=res))


if __name__ == "__main__":
    a = sys.argv[1:]
    run(int(a[0]), int(a[1]), int(a[2]), a[3], a[4])
