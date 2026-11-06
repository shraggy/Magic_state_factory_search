"""Task list for the sweep A rerun of the classes a wall-clock limit may have touched (README, "Reproducing sweep A").

Sweep A ran synthesis SAT and coefficient free SAT with a 1 h wall-clock guard per z3 query in addition to the step
limit; the guard is now removed (sat_certify.certify(time_limit=None)). This writes one shell line per task to
tasks_A_rerun.txt (one single-core cluster task each), output under outputs/staged_A_rerun/:
  syn/A_*.json     synthesis SAT as in A (seed: best n after check-pattern SAT, step limit 3e9), no guard: the 9
                   classes that spent > 1 h in restricted SAT and the 9 that end above the earlier catalogue;
  syn/A6_*.json    the same 9 with A's seed and step limit 6e9 (the earlier run's), so the step limit is the only
                   change from A;
  syn/old_*.json   the 9 above the catalogue, with the earlier run's settings (seed: the gluing n, step limit 6e9),
                   to find the cause of the difference (AGENTS.md rule 8);
  freeold/*.json   the same 9, coefficient free SAT seeded at the gluing n (the earlier cascade's quick-sat form);
  l3_k6/free_coef_1.5e9_unseeded/  the 12 coefficient free-SAT steps that averaged > 3500 s per query.
    python searches/rerun_A_tasks.py           (from the repository root, after sweep A)
    python searches/rerun_A_tasks.py --merge   (after the rerun) -> outputs/staged_A_rerun/l3_k{5,6}/iii_f.json
"""
import os, json, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = os.path.join(ROOT, "outputs", "staged_A"); R = "outputs/staged_A_rerun"
NINE = {"1892", "2320:771", "2320:3224", "2320:9752", "2320:16992", "2320:82457", "2320:82505", "2320:82512",
        "41984:4228"}                                     # outputs/staged_A/compare_l3.json, worse_list


def per_class():
    info, specs = {}, {}
    for k in (5, 6):
        specs.update(json.load(open(os.path.join(A, f"l3_k{k}", "iii_r.json")))["specs"])
        for p in glob.glob(os.path.join(A, f"l3_k{k}", "cls", "*.json")):
            d = json.load(open(p)); I, M = d["chunk"]
            r = d["results"] if isinstance(d["results"], list) else eval(d["results"])
            for o in d["classes"]:
                n = lambda pre: [x["n"] for x in r if x["orbit"] == o and x["stage"].startswith(pre)]
                pre = [x["n"] for x in r if x["orbit"] == o and not x["stage"].startswith("iii_r")]
                info[o] = dict(k=k, I=I, M=M, sec=d["per_class_seconds"][o].get("iii_r", 0),
                               eg=min(n("ii_eg"), default=None), pat=min(n("iii_r:pattern"), default=None),
                               pre=min(pre, default=None))
    return info, specs


def main():
    info, specs = per_class()
    long_r = {o for o, i in info.items() if i["sec"] > 3600}
    long_f = set()
    for p in glob.glob(os.path.join(A, "l3_k*", "free_coef*", "*.json")):
        for o, v in json.load(open(p))["per_class"].items():
            if any(s["seconds"] / max(1, s["queries"]) >= 3500 for s in v.get("steps") or []):
                long_f.add(o)
    py = "cd searches && ../.venv/bin/python -u -c"
    L = []
    for o in sorted(NINE | long_r):
        i, sp, k, f = info[o], specs[o], info[o]["k"], o.replace(":", "_")
        caps = [("A", min(x for x in (i["pat"], i["pre"]) if x), 3e9)]
        if o in NINE:
            caps += [("A6", caps[0][1], 6e9), ("old", i["eg"] or i["pre"], 6e9)]
        for tag, cap, rl in caps:
            out = f"../{R}/syn/{tag}_k{k}_{f}.json"
            L.append(f"{py} \"import json, synthesis_sat as SS; b, l = SS.solve('{sp}', {k}, {cap}, rl={rl:.0f}); "
                     f"json.dump(dict(k={k}, spec='{sp}', orbit='{o}', cap={cap}, rlimit={rl:.0f}, best=b, log=l), "
                     f"open('{out}', 'w'))\"")
        if o in NINE:
            cap, out = i["eg"] or i["pre"], f"../{R}/freeold/k{k}_{f}.json"
            L.append(f"{py} \"import json; from quick_sat import solve_target as S; "
                     f"r = S((3, {k}, '{sp}', 25.0, 3, 6, {cap}, 600000000)); r.update(orbit='{o}', cap={cap}); "
                     f"json.dump(r, open('{out}', 'w'))\"")
    for o in sorted(long_f):
        i = info[o]
        L.append(f".venv/bin/python -u searches/staged_class.py --free 3 {i['k']} {i['I']} {i['M']} coef:1.5e9:unseeded {R}")
    for d in ("syn", "freeold", "l3_k5", "l3_k6"):
        os.makedirs(os.path.join(ROOT, R, d), exist_ok=True)
    for k in (5, 6):                                      # --free reads the iii_r seeds of its cell
        src, dst = os.path.join(A, f"l3_k{k}", "iii_r.json"), os.path.join(ROOT, R, f"l3_k{k}", "iii_r.json")
        if not os.path.exists(dst):
            open(dst, "w").write(open(src).read())
    open(os.path.join(ROOT, "tasks_A_rerun.txt"), "w").write("\n".join(L) + "\n")
    print(f"{len(L)} tasks: {len(NINE | long_r)} classes synthesis SAT (A settings), {len(NINE)} with step limit 6e9, "
          f"{len(NINE)} with earlier settings, "
          f"{len(NINE)} free SAT seeded at gluing n, {len(long_f)} coefficient free-SAT classes")


def merge():
    """After the rerun: every verified factory of outputs/staged_A_rerun/{syn,freeold,l3_k*/free_*}, credited to the
    class it deposits, written as outputs/staged_A_rerun/l3_k{k}/iii_f.json (the schema of a pass file, best per class
    only where the rerun found one), so searches/staged_best.py --extra outputs/staged_A_rerun takes the smaller n."""
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import staged as STG
    out = os.path.join(ROOT, R); cache = os.path.join(out, "_cache"); STG.CACHE[0] = cache
    found = {5: [], 6: []}
    for p in glob.glob(os.path.join(out, "syn", "*.json")) + glob.glob(os.path.join(out, "freeold", "*.json")):
        d = json.load(open(p)); b = d.get("best"); k = d.get("k") or int(os.path.basename(p).split("_k")[1][0])
        if b:
            tag = "iii_r:synthesis:" + os.path.basename(p).split("_k")[0] if "/syn/" in p else "iii_f:coef:gluing-seed"
            found[k].append((b, tag))
    for k in (5, 6):
        for p in glob.glob(os.path.join(out, f"l3_k{k}", "free_*", "*.json")):
            for r in json.load(open(p))["results"]:
                found[k].append((r, r["stage"]))
    for k, F in found.items():
        best = {}
        for b, tag in F:
            o, g = STG.label_fix(dict(k=k, r=b["N"] - k, gates=b["gates"]), 3, k, cache)
            n = STG.n_of(g, 3) if o is not None else None
            if o is not None and (o not in best or (n, b["N"]) < (best[o]["n"], best[o]["N"])):
                best[o] = dict(n=n, N=b["N"], stage=tag + " (rerun)", gates=[[sorted(S), c] for S, c in g])
        rec = dict(L=3, k=k, pass_="iii_f", note="rerun without the 1 h guard (rerun_A_tasks.py)", best=best)
        json.dump(rec, open(os.path.join(out, f"l3_k{k}", "iii_f.json"), "w"))
        print(f"l=3 k={k}: {len(best)} classes with a rerun factory")


if __name__ == "__main__":
    import sys
    merge() if sys.argv[1:] == ["--merge"] else main()
