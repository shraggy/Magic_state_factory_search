"""Boolean free SAT (bool_sat.py) on the 17 l=3 classes whose catalogue best came from wall-clock free SAT: for each
class, N = cat_N (8), cap = the seed of the stages before free SAT (as calib_free_sat), supports up to wmax, CP-SAT
deterministic time limit dlimit. Records whether it reaches the catalogue n and verifies every factory found.
    python searches/calib_bool_sat.py --one I WMAX DLIMIT OUTDIR
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import calib_free_sat as CF, bool_sat as B
if __name__ == "__main__":
    i, w, dl, od = int(sys.argv[2]), sys.argv[3], float(sys.argv[4]), sys.argv[5]
    w = None if w == "none" else int(w)
    c = CF.seeds()[i]; os.makedirs(od, exist_ok=True)
    dst = os.path.join(od, f"b{i}_w{w}_d{int(dl)}.json")
    t = time.time(); st, n, g = B.solve(c["spec"], c["k"], c["cat_N"], c["seed"], w, dl)
    ok = bool(g) and B.check(c["spec"], c["k"], c["cat_N"], g)
    r = dict(c, wmax=w, dlimit=dl, status=st, n=n, verified=ok, reached=bool(ok and n <= c["cat_n"]) or c["seed"] <= c["cat_n"],
             seconds=round(time.time() - t, 1), host=os.uname().nodename, gates=[[sorted(S), 1] for S, _ in g] if g else None)
    json.dump(r, open(dst + ".tmp", "w")); os.replace(dst + ".tmp", dst)
    print(c["k"], c["spec"], w, dl, st, n, "cat", c["cat_n"], ok, r["reached"], r["seconds"], "s", flush=True)
