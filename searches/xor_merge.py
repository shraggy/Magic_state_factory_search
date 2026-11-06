"""Drop-in pass iii_f from the XOR free-SAT runs (xor_free.py), the schema of staged_class.do_merge_free.

Seeds: SEEDDIR/l3_k{k}/iii_r.json (the published pass before iii_f). Adds every verified factory of every
OUT/l3_k{k}/xor_<VAR>/ directory whose name matches one of the --vars patterns (default: all), smallest (n, N) per
class, credited to the class it deposits (as labelled when found). Writes OUT/l3_k{k}/iii_f.json with the same keys
as sweep A's iii_f.json (best, lowered, new, variants, seconds as a CPU sum, specs).
    python searches/xor_merge.py k SEEDDIR OUT [--vars x_w5_c2e6,x_wNone_c2e6]
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import staged as STG, staged_class as SC


def merge(k, seeddir, out, vars_=None):
    L = 3; cache, classes = SC.setup(L, k, seeddir)
    seed = json.load(open(os.path.join(seeddir, f"l{L}_k{k}", "iii_r.json")))["best"]
    best = {o: dict(b) for o, b in seed.items()}; prevn = {o: b["n"] for o, b in best.items()}
    stats, secs = {}, 0.0
    for vd in sorted(glob.glob(os.path.join(out, f"l{L}_k{k}", "xor_*"))):
        v = os.path.basename(vd)[4:]
        if (vars_ and v not in vars_) or (not vars_ and "seed" in v):
            continue                                       # seeded-from-another-pass runs only when named
        done, low = set(), set()
        for f in glob.glob(os.path.join(vd, "*.json")):
            x = json.load(open(f)); done |= set(x["classes"])
            secs += sum(p.get("seconds", 0.0) for p in x["per_class"].values())
            for r in x["results"]:
                o = str(r["orbit"]); cur = best.get(o)
                if o in classes and (cur is None or (r["n"], r["N"]) < (cur["n"], cur["N"])):
                    best[o] = dict(n=r["n"], N=r["N"], stage=r["stage"], gates=r["gates"])
                if o in prevn and r["n"] < prevn[o]:
                    low.add(o)
        stats[v] = dict(classes_done=len(done), not_done=len(set(classes) - done), lowered_below_seed=len(low))
    lowd = {o: [a, best[o]["n"]] for o, a in prevn.items() if best[o]["n"] < a}
    new = sorted(o for o in best if o not in prevn)
    rec = dict(L=L, k=k, pass_="iii_f", settings=dict(solver="xor_sat2 (CryptoMiniSat 5.17.0, pycryptosat)",
               variants=sorted(stats), N=f"k+1..k+3"), variants=stats, seconds=round(secs, 1),
               seconds_are_cpu_sum=True, classes=len(classes), reached=len(best), new=new, lowered=lowd,
               inputs_saved=sum(a - b for a, b in lowd.values()), best=best, rows=[], steps=None,
               specs={o: c["spec"] for o, c in classes.items()})
    SC.save(os.path.join(out, f"l{L}_k{k}", "iii_f.json"), rec)
    print(f"xor merge l=3 k={k}: lowered {len(lowd)} below iii_r by {rec['inputs_saved']} inputs; {stats}", flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]; v = None
    if "--vars" in a:
        i = a.index("--vars"); v = set(a[i + 1].split(",")); a = a[:i] + a[i + 2:]
    merge(int(a[0]), a[1], a[2], v)
