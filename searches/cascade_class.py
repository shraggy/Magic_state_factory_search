"""The cascade for one output class, so a large (l, k) cell can be split across many cluster tasks.

Runs, for one target class of outputs/cascade_classes/l{L}_k{k}.json (written by --list), the per-class stages of
searches/cascade.py with the same functions and arguments: (ii') gluing on entangled support, (iv) check-pattern SAT
(l=3), (v) quick SAT over all rotations with the deterministic step limit, each seeded with the best n so far,
starting from the cell's stage (i)-(ii) best for that class. Every factory is verified (is_factory, distance) and
labelled with its exact class, as in cascade.run_cell. The merge (--merge) combines the per-class files with the
cell's stage (i)-(ii) file into the same l{L}_k{k}.json that cascade.py writes.

    python searches/cascade_class.py --list L k OUT            # class list + stages (i)-(ii) for the cell -> OUT/l{L}_k{k}_base.json
    python searches/cascade_class.py --one L k ORBIT OUT       # one class -> OUT/l{L}_k{k}/o{ORBIT}.json
    python searches/cascade_class.py --merge L k OUT           # -> OUT/l{L}_k{k}.json
"""
import os, sys, json, time
from concurrent.futures import ProcessPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cascade as CA
import unified_search as us
from classify import phase_polynomial
from enumerate_targets import enumerate_targets, index_of_poly
from quick_sat import solve_target

P = dict(nmax=11, smax=7, rmax=4, budget=25.0, span=3, wmax=6, rlimit=600000000)   # = cascade_main_v3


def classes_of(L, k, cache_dir=None):
    """Exact classes of cell (L, k) and the labelling context. With cache_dir, the result is pickled once to
    cache_dir/l{L}_k{k}_classes.pkl and reloaded by every later task (the table takes up to ~1 h to build)."""
    import pickle
    path = os.path.join(cache_dir, f"l{L}_k{k}_classes.pkl") if cache_dir else None
    if path and os.path.exists(path):
        return pickle.load(open(path, "rb"))
    out = _classes_of(L, k)
    if path:
        tmp = path + f".{os.getpid()}.tmp"
        pickle.dump(out, open(tmp, "wb"), protocol=pickle.HIGHEST_PROTOCOL); os.replace(tmp, path)
    return out


def _quiet(*a, **kw):
    """No-op logger; a module-level function so the labelling objects that keep it stay picklable."""
    return None


def _classes_of(L, k):
    try:
        res, (C, off, lab) = enumerate_targets(L, k, log=_quiet); fiber = None
    except ValueError:
        if L == 3:
            from enumerate_targets_fiber import enumerate_targets_fiber
            res, fiber = enumerate_targets_fiber(L, k, log=_quiet)
        else:
            from enumerate_targets_2adic import enumerate_targets_2adic
            res, fiber = enumerate_targets_2adic(L, k, log=_quiet)
        C = off = lab = None
    return {str(c["orbit"]): c for c in res["classes"] if c["genuine"]}, (C, off, lab, fiber)


def labeller(L, k, ctx, classes):
    C, off, lab, fiber = ctx
    def label(f):
        g = f["gates"]
        if not (us.is_factory(f, L) and us.distance(f, L) != 1):
            return None
        Pp = phase_polynomial([S for S, _ in g], [c for _, c in g], set(range(k)), L)
        if lab is not None:
            o = str(int(lab[index_of_poly(Pp, k, L, C, off)]))
        else:
            lb = fiber.label(fiber.digits_of_poly(Pp)); o = lb if isinstance(lb, str) else f"{lb[0]}:{lb[1]}"
        return o if o in classes else None
    return label


def do_list(L, k, out, jobs):
    t = time.time(); classes, ctx = classes_of(L, k, out); t_enum = time.time() - t
    label = labeller(L, k, ctx, classes); best = {}; T = dict(iii_enumerate=round(t_enum, 1))
    def rec(f, stage):
        o = label(f)
        if o is None:
            return
        n = sum(us.inputs(c, L) for _, c in f["gates"]); N = k + f["r"]
        if o not in best or (n, N) < best[o][:2]:
            best[o] = (n, N, stage, [[sorted(S), c] for S, c in f["gates"]])
    t = time.time()
    seen = set()
    for f in us.stage_two_group(L, k, P["nmax"], k, one_sign=True):
        key = (f["r"], frozenset(f["gates"]))
        if key not in seen:
            seen.add(key); rec(f, "two-group")
    T["i_two_group"] = round(time.time() - t, 1); t = time.time()
    us.S_MAX = P["smax"]
    for f in CA.stage_gluing_multi(L, k, P["rmax"]):
        rec(f, "gluing")
    T["ii_gluing"] = round(time.time() - t, 1)
    os.makedirs(out, exist_ok=True)
    json.dump(dict(L=L, k=k, params=P, seconds=T, classes={o: c["spec"] for o, c in classes.items()},
                   best={o: dict(n=v[0], N=v[1], stage=v[2], gates=v[3]) for o, v in best.items()}),
              open(os.path.join(out, f"l{L}_k{k}_base.json"), "w"))
    print(f"l={L} k={k}: {len(classes)} classes, (i)+(ii) reach {len(best)}; {T}", flush=True)


def do_one(L, k, orbit, out, jobs):
    dst = os.path.join(out, f"l{L}_k{k}", f"o{orbit.replace(':', '_')}.json")
    if os.path.exists(dst):
        print("exists", dst); return
    base = json.load(open(os.path.join(out, f"l{L}_k{k}_base.json")))
    spec = base["classes"][orbit]; b0 = base["best"].get(orbit)
    best = (b0["n"], b0["N"], b0["stage"], b0["gates"]) if b0 else None
    classes, ctx = classes_of(L, k, out); label = labeller(L, k, ctx, classes); T = {}
    def take(f, stage):
        nonlocal best
        if label(f) != orbit:
            return
        n = sum(us.inputs(c, L) for _, c in f["gates"]); N = k + f["r"]
        if best is None or (n, N) < best[:2]:
            best = (n, N, stage, [[sorted(S), c] for S, c in f["gates"]])
    t = time.time()
    f = CA._eglue((L, k, spec, P["rmax"]))
    if f:
        f["gates"] = [(frozenset(S), c) for S, c in f["gates"]]; take(f, "entangled-gluing")
    T["ii_entangled_gluing"] = round(time.time() - t, 1); t = time.time()
    if L == 3:
        r = CA._pattern((k, spec, best[0] if best else 10 ** 6))
        if r:
            take(dict(stage="pattern-sat", k=k, r=r[1] - k, gates=[(frozenset(S), c) for S, c in r[2]], params=spec),
                 "pattern-sat")
    T["iv_pattern_sat"] = round(time.time() - t, 1); t = time.time()
    r = solve_target((L, k, spec, P["budget"], P["span"], P["wmax"], best[0] if best else None, P["rlimit"]))
    if r["best"]:
        b = r["best"]
        if best is None or (b["n"], b["N"]) < best[:2]:
            best = (b["n"], b["N"], "quick-sat", b["gates"])
    T["v_quick_sat"] = round(time.time() - t, 1)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    tmp = dst + ".tmp"
    json.dump(dict(L=L, k=k, orbit=orbit, spec=spec, seconds=T, host=os.uname().nodename,
                   best=best and dict(n=best[0], N=best[1], stage=best[2], gates=best[3])), open(tmp, "w"))
    os.replace(tmp, dst)
    print(f"l={L} k={k} {orbit} {spec}: {best and best[:3]} {T}", flush=True)


def do_merge(L, k, out):
    base = json.load(open(os.path.join(out, f"l{L}_k{k}_base.json")))
    d = os.path.join(out, f"l{L}_k{k}"); best = {}; T = dict(base["seconds"]); missing = []
    for o, spec in base["classes"].items():
        p = os.path.join(d, f"o{o.replace(':', '_')}.json")
        if not os.path.exists(p):
            missing.append(o); b = base["best"].get(o)
        else:
            r = json.load(open(p)); b = r["best"]
            for s, v in r["seconds"].items():
                T[s] = round(T.get(s, 0) + v, 1)                  # summed over classes: CPU-seconds, not wall clock
        if b:
            best[o] = dict(n=b["n"], N=b["N"], stage=b["stage"], spec=spec)
    summary = dict(l=L, k=k, classes=len(base["classes"]), reached_i_ii=len(base["best"]), reached_total=len(best),
                   missing_class_files=missing, seconds=T, seconds_note="per-class stages summed over classes (CPU-s)",
                   params=base["params"], best=best)
    json.dump(dict(summary=summary), open(os.path.join(out, f"l{L}_k{k}.json"), "w"), indent=1)
    print(f"l={L} k={k}: {len(best)}/{len(base['classes'])} reached, {len(missing)} class files missing", flush=True)


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "--list":
        do_list(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(os.environ.get("LSB_DJOB_NUMPROC", "4")))
    elif m == "--one":
        do_one(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], 1)
    elif m == "--merge":
        do_merge(int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
