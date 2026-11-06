"""The staged cascade split for a cluster: cell-level sweeps once, then one task per chunk of classes, then a merge.

  --cell  L k OUT          cell-level passes, each from an empty book (what that pass alone finds, with gates):
                           i_tg two-group; i_gs sweep (orbit_group.run_k, every group in the menu); i_gh (G x H);
                           ii_dj disjoint gluing.  -> OUT/l{L}_k{k}/cell/<pass>.json
  --chunk L k I M OUT      classes sorted by label, those with index % M == I. Per class, seeded with the best n of
                           the cell files: i_gs target mode (G = Aut(U) and G = {e}, CP-SAT); ii_eg entangled
                           gluing; iii_r check-pattern then synthesis SAT (l=3; skipped at the bound n = tau+1). Each
                           factory is credited to the class it deposits (and, at l=4, its inverse to the inverse class).
                           -> OUT/l{L}_k{k}/cls/c{I}of{M}.json
  --merge L k OUT          cumulative best per class after each pass, in cascade order i_tg, i_gs, i_gh, ii_dj, ii_eg,
                           iii_r -> OUT/l{L}_k{k}/<pass>.json (the schema of staged.py, read by staged_compare.py)
  --free  L k I M VAR OUT  free SAT on the merged iii_r seeds, variant VAR: bool:w5:d600 (l=3 Boolean, supports<=5,
                           600 deterministic s per query), coef:1.5e9 (z3 over coefficients, step limit), add
                           :unseeded to start with no cap. N = k+1..k+3. -> OUT/l{L}_k{k}/free_<VAR>/c{I}of{M}.json
Per-class seeds come from the cell files only, so a class may start looser than in a one-process run; the merge takes
the minimum over everything found. Wall-clock and CPU seconds of every pass are recorded per class.
"""
import os, sys, json, time, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import staged as STG
import unified_search as us

CELL = ["i_tg", "i_gs", "i_gh", "ii_dj"]
ORDER = ["i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg", "iii_r"]
S = STG.DEFAULT


def setup(L, k, out):
    cache = os.path.join(out, "_cache"); os.makedirs(cache, exist_ok=True); STG.CACHE[0] = cache
    import orbit_group as OG, cascade_class as CC
    OG.load_classes = lambda k_, fc=None: CC.classes_of(3, k_, cache)      # one cached class table for every task
    classes, label = STG.ctx(L, k, cache)
    return cache, classes


def save(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    json.dump(obj, open(path + ".tmp", "w")); os.replace(path + ".tmp", path)


def p_gs_sweep(L, k, s, book, cache, jobs):
    if L != 3:
        return
    import orbit_group as OG
    labels = OG.Labels(k)
    sb, _ = OG.run_k(k, s["rmax"], labels, lambda *a: None)
    for o, b in sb.items():
        o2, g = STG.label_fix(dict(k=k, r=b["N"] - k, gates=b["gates"]), L, k, cache)
        book.add(o2, g, b["N"], "i_gs:sweep", keep_row=True)


def do_cell(L, k, out):
    cache, classes = setup(L, k, out)
    fn = dict(i_tg=STG.p_i_tg, i_gs=p_gs_sweep, i_gh=STG.p_i_gh, ii_dj=STG.p_ii_dj)
    for p in CELL:
        book = STG.Book(L, k, None); t = time.time()
        fn[p](L, k, S[p], book, cache, 1)
        save(os.path.join(out, f"l{L}_k{k}", "cell", f"{p}.json"),
             dict(L=L, k=k, pass_=p, settings=S[p], seconds=round(time.time() - t, 1), host=os.uname().nodename,
                  dropped=getattr(book, "dropped", 0), best=book.best, rows=book.rows))
        print(f"cell l={L} k={k} {p}: {len(book.best)} classes, {round(time.time() - t, 1)}s", flush=True)


def cell_seed(L, k, out):
    best = {}
    for p in CELL:
        f = os.path.join(out, f"l{L}_k{k}", "cell", f"{p}.json")
        for o, b in json.load(open(f))["best"].items():
            if o not in best or b["n"] < best[o]:
                best[o] = b["n"]
    return best


def gs_target(L, k, spec, s):
    """i_gs target mode for one class: G = Aut(U) and G = {e}, r = 1..rmax, CP-SAT deterministic limit."""
    import orbit_group as OG
    dig = OG.digits_of_spec(spec, k); A = OG.aut_of(dig, k)
    menu = [("Aut", A)] + ([("trivial", [tuple(range(k))])] if k <= s["trivial_kmax"] and len(A) > 1 else [])
    out = []
    for gname, G in menu:
        orbs, of = OG.subset_orbits(G, k)
        for r in range(1, s["rmax"] + 1):
            st, n, x, cells = OG.cpsat_target(k, r, orbs, of, dig, s["dlimit"])
            if n is not None:
                out.append((OG.gates_of(k, r, orbs, cells, x), k + r, "i_gs:target:" + gname))
    return out


def do_chunk(L, k, I, M, out):
    cache, classes = setup(L, k, out)
    seed = cell_seed(L, k, out)
    mine = [o for i, o in enumerate(sorted(classes)) if i % M == I]
    res, per = [], {}
    import cascade as CA

    class Rec(STG.Book):                       # records every credited factory (with the class it lands in)
        def add(self, o, g, N, pas, keep_row=False, partner=True):
            if o is not None:
                res.append(dict(orbit=o, n=STG.n_of(g, L), N=N, stage=pas, gates=[[sorted(S_), c] for S_, c in g]))
            super().add(o, g, N, pas, keep_row, partner)
    for o in mine:
        spec = classes[o]["spec"]; book = Rec(L, k, None); per[o] = {}
        cur = lambda: min([v for v in (seed.get(o), book.best.get(o, {}).get("n")) if v is not None], default=None)
        t = time.time()
        if L == 3:
            for g, N, st in gs_target(L, k, spec, S["i_gs"]):
                o2, gg = STG.label_fix(dict(k=k, r=N - k, gates=g), L, k, cache); book.add(o2, gg, N, st)
        per[o]["i_gs"] = round(time.time() - t, 1); t = time.time()
        f = CA._eglue((L, k, spec, S["ii_eg"]["rmax"]))
        if f:
            o2, gg = STG.label_fix(f, L, k, cache); book.add(o2, gg, k + f["r"], "ii_eg")
        per[o]["ii_eg"] = round(time.time() - t, 1); t = time.time()
        c = cur()
        tau = STG.tau_of(L, k, o, spec)
        if L == 3 and not (c is not None and tau is not None and c == tau + 1):
            for kind, b in STG._r_one((k, o, spec, c if c is not None else 10 ** 6, S["iii_r"]))[1]:
                o2, gg = STG.label_fix(dict(k=k, r=b["N"] - k, gates=b["gates"]), L, k, cache)
                book.add(o2, gg, b["N"], "iii_r:" + kind)
        per[o]["iii_r"] = round(time.time() - t, 1)
        print(f"l={L} k={k} {o} {spec}: seed {seed.get(o)} -> {cur()} {per[o]}", flush=True)
    save(os.path.join(out, f"l{L}_k{k}", "cls", f"c{I}of{M}.json"),
         dict(L=L, k=k, chunk=[I, M], classes=mine, host=os.uname().nodename, per_class_seconds=per, results=res))


def do_merge(L, k, out):
    cache, classes = setup(L, k, out)
    d = os.path.join(out, f"l{L}_k{k}")
    found = {p: [] for p in ORDER}                                # per pass: (orbit, n, N, stage, gates)
    secs = {p: 0.0 for p in ORDER}; cpu = {p: False for p in ORDER}
    for p in CELL:
        x = json.load(open(os.path.join(d, "cell", f"{p}.json"))); secs[p] += x["seconds"]
        found[p] += [dict(orbit=o, **b) for o, b in x["best"].items()]
    files = sorted(glob.glob(os.path.join(d, "cls", "*.json"))); covered = set()
    for f in files:
        x = json.load(open(f)); covered |= set(x["classes"])
        for r in x["results"]:
            found[r["stage"].split(":")[0]].append(r)
        for o, v in x["per_class_seconds"].items():
            for p, sec in v.items():
                secs[p] += sec; cpu[p] = True
    best, prevn = {}, {}
    for p in ORDER:
        prevn = {o: b["n"] for o, b in best.items()}
        for r in found[p]:
            o = str(r["orbit"]); cur = best.get(o)
            if o in classes and (cur is None or (r["n"], r["N"]) < (cur["n"], cur["N"])):
                best[o] = dict(n=r["n"], N=r["N"], stage=r["stage"], gates=r["gates"])
        new = sorted(o for o in best if o not in prevn); low = {o: [v, best[o]["n"]] for o, v in prevn.items() if best[o]["n"] < v}
        save(os.path.join(d, f"{p}.json"), dict(L=L, k=k, pass_=p, settings=S[p], seconds=round(secs[p], 1),
             seconds_are_cpu_sum=cpu[p], classes=len(classes), reached=len(best), new=new, lowered=low,
             inputs_saved=sum(a - b for a, b in low.values()), best=best, rows=[], steps=None,
             specs={o: c["spec"] for o, c in classes.items()}, missing_class_files=sorted(set(classes) - covered)))
        print(f"merge l={L} k={k} {p}: reached {len(best)}/{len(classes)} (+{len(new)}), lowered {len(low)}, "
              f"{round(secs[p], 1)}s{' cpu-sum' if cpu[p] else ''}", flush=True)
    if set(classes) - covered:
        print("MISSING class files for", len(set(classes) - covered), "classes", flush=True)


def do_free(L, k, I, M, var, out):
    cache, classes = setup(L, k, out)
    d = os.path.join(out, f"l{L}_k{k}")
    seedf = json.load(open(os.path.join(d, "iii_r.json")))["best"]
    parts = var.split(":"); kind = parts[0]; unseeded = "unseeded" in parts
    mine = [o for i, o in enumerate(sorted(classes)) if i % M == I]; res, per = [], {}
    for o in mine:
        spec = classes[o]["spec"]; t = time.time()
        cap = None if unseeded else (seedf[o]["n"] if o in seedf else None)
        tau = STG.tau_of(L, k, o, spec); steps = []
        if not unseeded and cap is not None and tau is not None and cap == tau + 1:
            per[o] = dict(seconds=0.0, skipped="at bound"); continue
        if kind == "bool" and L == 3:
            import bool_sat as B
            w = None if parts[1] == "wNone" else int(parts[1][1:]); dl = float(parts[2][1:])
            for N in range(k + 1, k + 4):
                st, n, g = B.solve(spec, k, N, cap, w, dl); steps.append(dict(N=N, status=st, n=n))
                if g:
                    o2, gg = STG.label_fix(dict(k=k, r=N - k, gates=g), L, k, cache)
                    if o2 is not None:
                        res.append(dict(orbit=o2, n=n, N=N, stage="iii_f:" + var, gates=[[sorted(S_), c] for S_, c in gg]))
                        if o2 == o:
                            cap = n
        else:
            from quick_sat import solve_target
            rl = float(parts[1])
            r = solve_target((L, k, spec, 25.0, 3, 6, cap, int(rl))); steps = r["steps"]
            if r["best"]:
                b = r["best"]; o2, gg = STG.label_fix(dict(k=k, r=b["N"] - k, gates=b["gates"]), L, k, cache)
                if o2 is not None:
                    res.append(dict(orbit=o2, n=b["n"], N=b["N"], stage="iii_f:" + var, gates=[[sorted(S_), c] for S_, c in gg]))
        per[o] = dict(seconds=round(time.time() - t, 1), seed=None if unseeded else (seedf.get(o) or {}).get("n"), steps=steps)
        print(f"free {var} l={L} k={k} {o} {spec}: {per[o]}", flush=True)
    save(os.path.join(d, "free_" + var.replace(":", "_"), f"c{I}of{M}.json"),
         dict(L=L, k=k, chunk=[I, M], var=var, classes=mine, host=os.uname().nodename, per_class=per, results=res))


def do_merge_free(L, k, out):
    """Final pass iii_f: the iii_r seeds plus every verified factory of every free_<VAR>/ variant, smallest per class
    (credited to the class it deposits). Writes iii_f.json with the schema of the other passes, the per-variant gain,
    and the classes each variant has not finished."""
    cache, classes = setup(L, k, out)
    d = os.path.join(out, f"l{L}_k{k}")
    seed = json.load(open(os.path.join(d, "iii_r.json")))["best"]
    best = {o: dict(b) for o, b in seed.items()}; prevn = {o: b["n"] for o, b in best.items()}
    var_stats, secs = {}, 0.0
    for vd in sorted(glob.glob(os.path.join(d, "free_*"))):
        v = os.path.basename(vd)[5:]; done = set(); lowered = set()
        for f in glob.glob(os.path.join(vd, "*.json")):
            x = json.load(open(f)); done |= set(x["classes"])
            secs += sum(p.get("seconds", 0.0) for p in x["per_class"].values())
            for r in x["results"]:
                o = str(r["orbit"]); cur = best.get(o)
                if o in classes and (cur is None or (r["n"], r["N"]) < (cur["n"], cur["N"])):
                    best[o] = dict(n=r["n"], N=r["N"], stage=r["stage"], gates=r["gates"])
                if o in prevn and r["n"] < prevn[o]:
                    lowered.add(o)
        var_stats[v] = dict(classes_done=len(done), not_done=len(set(classes) - done), lowered_below_seed=len(lowered))
    low = {o: [v, best[o]["n"]] for o, v in prevn.items() if best[o]["n"] < v}
    new = sorted(o for o in best if o not in prevn)
    save(os.path.join(d, "iii_f.json"), dict(L=L, k=k, pass_="iii_f", settings=S["iii_f"], variants=var_stats,
         seconds=round(secs, 1), seconds_are_cpu_sum=True, classes=len(classes), reached=len(best), new=new, lowered=low,
         inputs_saved=sum(a - b for a, b in low.values()), best=best, rows=[], steps=None,
         specs={o: c["spec"] for o, c in classes.items()}))
    print(f"merge-free l={L} k={k}: reached {len(best)}/{len(classes)} (+{len(new)}), lowered {len(low)} by "
          f"{sum(a - b for a, b in low.values())} inputs; variants {var_stats}", flush=True)


if __name__ == "__main__":
    m, a = sys.argv[1], sys.argv[2:]
    if m == "--cell":
        do_cell(int(a[0]), int(a[1]), a[2])
    elif m == "--chunk":
        do_chunk(int(a[0]), int(a[1]), int(a[2]), int(a[3]), a[4])
    elif m == "--merge":
        do_merge(int(a[0]), int(a[1]), a[2])
    elif m == "--free":
        do_free(int(a[0]), int(a[1]), int(a[2]), int(a[3]), a[4], a[5])
    elif m == "--merge-free":
        do_merge_free(int(a[0]), int(a[1]), a[2])
