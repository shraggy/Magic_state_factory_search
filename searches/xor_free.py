"""Batch driver for incremental XOR free SAT (xor_sat2.py) over the l=3 classes of sweep A, as staged_class.py --free.

  python searches/xor_free.py L k I M VAR SEEDDIR OUT [--classes o1,o2,...]
VAR = x|xs:w<wmax|wNone>:c<confl>[:unseeded][:seed=<pass>]; xs restricts output parts to synthesis_sat.parts
(the XOR form of synthesis SAT); seed=<pass> reads the seed from SEEDDIR/l{L}_k{k}/<pass>.json (default iii_r). Classes sorted by label, index % M == I (or the explicit list). Seed cap per
class: best n of SEEDDIR/l{L}_k{k}/iii_r.json (sweep A), unless unseeded; classes at the bound n = tau+1 are skipped
(as staged_class.do_free). N = k+1..k+3; the cap is lowered to every n found for the class itself. Every factory is
verified (is_factory, distance != 1 inside xor_sat2) and labelled with the exact class (staged.label_fix): credited
to the class it deposits.  -> OUT/l{L}_k{k}/xor_<VAR>/c{I}of{M}.json with per-class seconds, steps and results.
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import staged as STG
import staged_class as SC
import xor_sat2 as X


def run(L, k, I, M, var, seeddir, out, only=None):
    cache, classes = SC.setup(L, k, seeddir)
    p = var.split(":")
    sp = [x[5:] for x in p if x.startswith("seed=")]; seedpass = sp[0] if sp else "iii_r"
    seedf = json.load(open(os.path.join(seeddir, f"l{L}_k{k}", seedpass + ".json")))["best"]
    restrict = p[0] == "xs"; w = None if p[1] == "wNone" else int(p[1][1:]); confl = int(float(p[2][1:]))
    unseeded = "unseeded" in p
    mine = only if only else [o for i, o in enumerate(sorted(classes)) if i % M == I]
    res, per = [], {}
    for o in mine:
        spec = classes[o]["spec"]; t = time.time(); c0 = time.process_time()
        cap = None if unseeded else (seedf[o]["n"] if o in seedf else None)
        tau = STG.tau_of(L, k, o, spec); steps = []
        if not unseeded and cap is not None and tau is not None and cap == tau + 1:
            per[o] = dict(seconds=0.0, skipped="at bound"); continue
        for N in range(k + 1, k + 4):
            tN = time.time()
            op = None
            if restrict:
                import synthesis_sat as SS
                op = SS.parts(spec, k)
            st, n, g, q, tr = X.solve(spec, k, N, cap, w, confl, op=op)
            steps.append(dict(N=N, status=st, n=n, queries=q, seconds=round(time.time() - tN, 2)))
            if g:
                o2, gg = STG.label_fix(dict(k=k, r=N - k, gates=g), L, k, cache)
                if o2 is not None:
                    res.append(dict(orbit=o2, n=n, N=N, stage=("iii_r:xorsyn:" if restrict else "iii_f:xor:") + var, gates=[[sorted(S_), c] for S_, c in gg]))
                    if o2 == o:
                        cap = n
        per[o] = dict(seconds=round(time.time() - t, 2), cpu=round(time.process_time() - c0, 2),
                      seed=None if unseeded else (seedf.get(o) or {}).get("n"), steps=steps)
        print(f"xor {var} l={L} k={k} {o} {spec}: {per[o]}", flush=True)
    tag = f"c{I}of{M}" if not only else "list_" + __import__("hashlib").md5(",".join(only).encode()).hexdigest()[:8]
    SC.save(os.path.join(out, f"l{L}_k{k}", "xor_" + var.replace(":", "_").replace("=", ""), tag + ".json"),
            dict(L=L, k=k, chunk=[I, M], var=var, classes=mine, host=os.uname().nodename, per_class=per, results=res))


if __name__ == "__main__":
    a = sys.argv[1:]
    only = None
    if "--classes" in a:
        i = a.index("--classes"); only = a[i + 1].split(","); a = a[:i] + a[i + 2:]
    run(int(a[0]), int(a[1]), int(a[2]), int(a[3]), a[4], a[5], a[6], only)
