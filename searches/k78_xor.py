"""Stage (iii) at k = 7, 8 on the CCZ-layer targets: the incremental XOR free SAT of searches/xor_sat2.py (one
solver per (target, N), per-query conflict limit and total limit: deterministic), seeded with the best n of the
target's fingerprint class after the gluing pass (outputs/k78_pipeline/full/l3_k{k}/ii_eg.json) and asked only for
fewer inputs, N = k+1 .. k+NSPAN. Classes already on the bound n = tau + 1 (tau fixed) are skipped, as in sweep A
(staged_class.do_free). A target is given to xor_sat2 as the spec "CCZabc,..." (qubits 0..k-1). Every factory is
re-verified here: us.is_factory (monomial form and brute force), distance 2, entangled_gluing.deposits.
    python searches/k78_xor.py TARGETS.json I M OUT --wmax W --confl C --total T --nspan 2
TARGETS.json: list of dict(id, k, ccz, fp) (one representative per class, or more).  -> OUT/c{I}of{M}.json
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import unified_search as us
import xor_sat2 as X
from entangled_gluing import deposits

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("targets"); ap.add_argument("I", type=int); ap.add_argument("M", type=int)
    ap.add_argument("out"); ap.add_argument("--wmax", type=int, default=None); ap.add_argument("--confl", type=int, default=2_000_000)
    ap.add_argument("--total", type=int, default=None); ap.add_argument("--nspan", type=int, default=2)
    ap.add_argument("--seeddir", default=os.path.join(ROOT, "outputs", "k78_pipeline", "full"))
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    dst = os.path.join(a.out, f"c{a.I}of{a.M}.json")
    if os.path.exists(dst):
        print("exists", dst); sys.exit()
    T = json.load(open(a.targets)); T = T["targets"] if isinstance(T, dict) else T
    seeds = {}
    for k in (7, 8):
        p = os.path.join(a.seeddir, f"l3_k{k}", "ii_eg.json")
        if os.path.exists(p):
            seeds[k] = json.load(open(p))["best"]
    res = []
    for i, t in enumerate(T):
        if i % a.M != a.I:
            continue
        k = t["k"]; spec = ",".join("CCZ" + "".join(map(str, sorted(x))) for x in t["ccz"])
        sb = seeds.get(k, {}).get(t["fp"]); cap = sb["n"] if sb else None
        rec = dict(id=t["id"], k=k, fp=t["fp"], ccz=t["ccz"], spec=spec, cap=cap, runs=[], best=None)
        if sb and sb.get("class_on_bound"):
            rec["skipped"] = "at bound"; res.append(rec); continue
        for N in range(k + 1, k + 1 + a.nspan):
            c0 = time.process_time()
            st, n, g, q, tr = X.solve(spec, k, N, cap, a.wmax, a.confl, a.total)
            run = dict(N=N, status=st, n=n, queries=q, trace=tr, cpu_s=round(time.process_time() - c0, 1))
            if g:
                f = dict(k=k, r=N - k, gates=g); tgt = {frozenset(x): 4 for x in t["ccz"]}
                ok = bool(us.is_factory(f, 3) and us.is_factory(f, 3, brute_force=True) and us.distance(f, 3) == 2 and deposits(g, k, tgt))
                run["verified"] = ok
                if ok and (rec["best"] is None or n < rec["best"]["n"]):
                    rec["best"] = dict(n=n, N=N, gates=[[sorted(S), c] for S, c in g]); cap = n
            rec["runs"].append(run)
            print(t["id"], k, t["fp"], "N", N, st, n, run["cpu_s"], "s", flush=True)
        res.append(rec)
    json.dump(dict(settings=vars(a), host=os.uname().nodename, results=res), open(dst + ".tmp", "w")); os.replace(dst + ".tmp", dst)
