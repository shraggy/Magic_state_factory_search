"""Stage (iii) of sweep A as one cascade per class, seeded from the gluing pass, with a chosen chain of SAT stages.

Per class (l=3): cap = best n after ii_eg (SEEDDIR/l3_k{k}/ii_eg.json). Each stage of CHAIN runs in order, seeded with the
running best (n of everything before it); a stage is skipped when the cap equals tau+1 (Corollary 3). Stages:
  pat          check-pattern SAT (pattern_sat, current code; staged.DEFAULT iii_r: R<=4, mmax 12, rlimit 2e9), every
               z3 query in a fresh spawned process (pat_eg.run_class)
  syn          synthesis SAT in z3 (synthesis_sat.solve: N=k+1..k+3, rlimit 3e9; each query a fresh spawned process)
  xs:W:C       synthesis SAT in XOR form (xor_sat2 with output parts synthesis_sat.parts), wmax W, C conflicts/query
  x:W:C        free XOR SAT (xor_sat2), wmax W (5, 6 or None), C conflicts per query
xs and x run N = k+1..k+3 with the cap lowered to each n found for the class itself. Every factory is verified
(is_factory, distance != 1) and labelled with the exact class (staged.label_fix); factories credited to another class are
kept in `results` (as in the other drivers), but only the class's own factories lower its running cap.
    python searches/chain_iii.py k I M SEEDDIR OUT CHAIN [--classes o1,o2]   CHAIN e.g. syn,xs:None:2e6,x:5:2e6
    -> OUT/l3_k{k}/chain/c{I}of{M}.json: per class, per stage: seed, n, seconds; all results with their stage
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import staged as STG, staged_class as SC


def stage(name, k, spec, cap):
    """Returns a list of (n, N, gates) found below cap, smallest last."""
    if name == "pat":
        import pat_eg as PE
        s = STG.DEFAULT["iii_r"]
        b, _ = PE.run_class(k, spec, cap, Rs=tuple(range(1, s["R"] + 1)), mmax=s["mmax"], rlimit=s["pat_rlimit"])
        return [(b["n"], b["N"], b["gates"])] if b and b["ok"] else []
    if name == "syn":
        import synthesis_sat as SS
        s = STG.DEFAULT["iii_r"]
        b, _ = SS.solve(spec, k, cap, rl=s["syn_rlimit"], span=tuple(range(1, s["syn_span"] + 1)))
        return [(b["n"], b["N"], b["gates"])] if b else []
    kind, w, c = name.split(":"); w = None if w == "None" else int(w); c = int(float(c))
    import xor_sat2 as X
    op = None
    if kind == "xs":
        import synthesis_sat as SS
        op = SS.parts(spec, k)
    out = []
    for N in range(k + 1, k + 4):
        st, n, g, q, tr = X.solve(spec, k, N, cap, w, c, op=op)
        if g:
            out.append((n, N, [[sorted(S), cc] for S, cc in g])); cap = n
    return out


def run(k, I, M, seeddir, out, chain, only=None):
    L = 3; cache, classes = SC.setup(L, k, seeddir)
    seed = json.load(open(os.path.join(seeddir, f"l{L}_k{k}", "ii_eg.json")))["best"]
    mine = only or [o for i, o in enumerate(sorted(classes)) if i % M == I]; res, per = [], {}
    for o in mine:
        spec = classes[o]["spec"]; tau = STG.tau_of(L, k, o, spec)
        cap = seed[o]["n"] if o in seed else 10 ** 6; per[o] = dict(seed=cap if o in seed else None, stages=[])
        for name in chain:
            if tau is not None and cap == tau + 1:
                per[o]["stages"].append(dict(stage=name, skipped="at bound")); continue
            t = time.time(); c0 = cap
            for n, N, g in stage(name, k, spec, cap):
                o2, gg = STG.label_fix(dict(k=k, r=N - k, gates=[(frozenset(S), cc) for S, cc in g]), L, k, cache)
                if o2 is None:
                    continue
                res.append(dict(orbit=o2, n=n, N=N, stage="iii:" + name, gates=[[sorted(S_), cc] for S_, cc in gg]))
                if o2 == o and n < cap:
                    cap = n
            per[o]["stages"].append(dict(stage=name, seed=c0, n=cap, seconds=round(time.time() - t, 1)))
        per[o]["final"] = cap
        print(f"chain k={k} {o} {spec}: {per[o]}", flush=True)
    tag = f"c{I}of{M}" if not only else "list_" + "_".join(o.replace(":", "-") for o in only)[:80]
    SC.save(os.path.join(out, f"l{L}_k{k}", "chain", tag + ".json"),
            dict(L=L, k=k, chunk=[I, M], chain=chain, classes=mine, host=os.uname().nodename, per_class=per, results=res))


if __name__ == "__main__":
    a = sys.argv[1:]; only = None
    if "--classes" in a:
        i = a.index("--classes"); only = a[i + 1].split(","); a = a[:i] + a[i + 2:]
    run(int(a[0]), int(a[1]), int(a[2]), a[3], a[4], a[5].split(","), only)
