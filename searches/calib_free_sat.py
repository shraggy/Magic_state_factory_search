"""Free-SAT step-limit calibration: for every l=3 class whose catalogue best came from a free-SAT run with a wall-clock
budget (no step limit), rerun free SAT (quick_sat.solve_target) with a deterministic z3 step limit, seeded with the best
n of the stages before free SAT (entangled gluing, check-pattern SAT, synthesis SAT), and record whether it reaches the
catalogue n. One fresh process per (class, rlimit).
    python searches/calib_free_sat.py --rlimits 6e8 3e9 --jobs 5   -> outputs/staged/calib_free_sat.json
"""
import os, sys, json, csv, glob, argparse, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))


def seeds():
    eg = {(e["k"], e["spec"]): e["eg_n"] for e in json.load(open(os.path.join(O, "entangled_gluing_l3.json")))}
    pre = dict(eg)
    for d in ("pattern_sat", "synthesis_sat"):
        for f in glob.glob(os.path.join(O, d, "*.json")):
            x = json.load(open(f))
            if x.get("best"):
                key = (x["k"], x["spec"]); n = x["best"]["n"]
                if pre.get(key) is None or n < pre[key]:
                    pre[key] = n
    out = []
    for r in csv.DictReader(open(os.path.join(O, "l3_catalogue_k6.csv"))):
        if r["source"] in ("cascade:quick-sat", "quick-sat"):
            key = (int(r["k"]), r["spec"]); out.append(dict(k=key[0], spec=r["spec"], cat_n=int(r["n"]), cat_N=int(r["N"]), seed=pre.get(key)))
    return out


def one(a):
    from quick_sat import solve_target
    c, rl = a; t = time.time()
    r = solve_target((3, c["k"], c["spec"], 25.0, 3, 6, c["seed"], int(rl)))
    b = r["best"]
    return dict(c, rlimit=rl, n=b and b["n"], N=b and b["N"], reached=bool((b and b["n"] <= c["cat_n"]) or (c["seed"] is not None and c["seed"] <= c["cat_n"])),
                steps=r["steps"], seconds=round(time.time() - t, 1))


if __name__ == "__main__":
    if sys.argv[1:2] == ["--one"]:                       # cluster task: class index, rlimit, out dir
        i, rl, od = int(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
        unseeded = "--unseeded" in sys.argv                # start free SAT with no cap (as the old unseeded runs)
        c = seeds()[i]; os.makedirs(od, exist_ok=True)
        if unseeded:
            c = dict(c, seed=None, seed_was=c["seed"])
        dst = os.path.join(od, f"c{i}_rl{rl:.0e}" + ("_unseeded" if unseeded else "") + ".json")
        if os.path.exists(dst):
            print("exists", dst); sys.exit(0)
        r = one((c, rl)); r["host"] = os.uname().nodename
        json.dump(r, open(dst + ".tmp", "w")); os.replace(dst + ".tmp", dst)
        print(r["k"], r["spec"], rl, "->", r["n"], "cat", r["cat_n"], r["reached"], r["seconds"], "s", flush=True)
        sys.exit(0)
    ap = argparse.ArgumentParser(); ap.add_argument("--rlimits", type=float, nargs="+", default=[6e8, 3e9])
    ap.add_argument("--jobs", type=int, default=5); a = ap.parse_args()
    C = seeds(); print(len(C), "classes:", [(c["k"], c["spec"], c["seed"], c["cat_n"]) for c in C], flush=True)
    tasks = [(c, rl) for rl in a.rlimits for c in C]
    res = []
    with Pool(a.jobs, maxtasksperchild=1) as p:
        for r in p.imap_unordered(one, tasks):
            res.append(r); print(r["k"], r["spec"], "rlimit", r["rlimit"], "seed", r["seed"], "->", r["n"], "cat", r["cat_n"], r["reached"], r["seconds"], "s", flush=True)
            json.dump(res, open(os.path.join(O, "staged", "calib_free_sat.json"), "w"), indent=1)
