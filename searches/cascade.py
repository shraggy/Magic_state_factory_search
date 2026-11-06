"""The cascaded search: four disjoint stages, each working only on what the stages before it did not cover.

  (i)   two-group circuits                      every symmetric-output circuit within the parameter bounds
  (ii)  shared-check gluing                     only partitions with >= 2 blocks (one block is a two-group circuit)
  (iii) exact target classes                    every Clifford class of output; each factory from (i)-(ii) is
                                                mapped to its class and the best n per class is recorded
  (ii') gluing on entangled support             for each class: an optimal Z(pi/2^l) network of the target as the
                                                output part, check couplings from the same 2-adic solve as (ii)
  (iv)  check-pattern SAT (l=3)                 output parts fixed to a synthesis of the target (App. S4: fixing the
                                                checks leaves a circuit that prepares it); solves only for which
                                                checks each rotation touches; asks for fewer inputs than (i)-(ii')
  (v)   quick SAT on the classes                over all rotations; asks only for fewer inputs than every stage
                                                before it, so it returns only factories outside (i)-(iv)

One (l, k) cell per call; wall-clock per stage is recorded. Writes outputs/cascade/l{L}_k{k}.json with every
factory (stage, n, N, exact class) and the per-stage times.

    python searches/cascade.py --cells 3:3-5 2:2-6 4:2-3 --jobs 4
"""
import os, sys, json, time, argparse
from concurrent.futures import ProcessPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search as us  # noqa: E402
from classify import phase_polynomial  # noqa: E402
from enumerate_targets import enumerate_targets, index_of_poly  # noqa: E402
from quick_sat import solve_target  # noqa: E402
from symfree_search import solve_binding, _partitions  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(ROOT, "outputs", "cascade")


def stage_gluing_multi(L, k, r_max):
    """Stage (ii) restricted to partitions with at least two blocks: a single block of size k is the fully symmetric
    output, which stage (i) already builds as a two-group circuit."""
    for parts in _partitions(k, L):
        if len(parts) < 2:
            continue
        for r in range(1, r_max + 1):
            sol = solve_binding(parts, r, L)
            if sol is not None:
                yield dict(stage="gluing", k=k, r=r, gates=[(S, c % (1 << L)) for S, c in sol["gates"]],
                           params="blocks=" + "+".join(map(str, parts)))
                break


SYN_RLIMIT = 3e9                                 # z3 step limit per synthesis-SAT query (as synthesis_sat.py)
APPROX_KMIN, APPROX_FRAMES = 11, 1              # as large_unified.py: at k >= 11 the class is the identity frame


def _eglue(args):
    """Stage (ii'): one class. Returns the factory dict (unverified) or None."""
    L, k, spec, rmax = args
    from coset_certify import parse_target
    import entangled_gluing as EG, entangled_gluing_l4 as EG4
    tgt = parse_target(spec, L)
    if L == 3:
        import synthillation_compare as sc
        oo = [(frozenset(q for q in range(k) if y >> q & 1), 1) for y in sc.opt_A(k, tgt)]
    else:
        net = EG4.network(tgt, k, L)
        if net is None:
            return None
        oo = [(frozenset(q for q in range(k) if y >> q & 1), c) for y, c in net[1].items()]
    best = None
    for nc in range(1, rmax + 1):
        # both forms of the output constraint: with free diagonal Cliffords absorbing output terms, and exact; the
        # 2-adic reduction of n is greedy, so either can give the smaller factory (e.g. T4,CS02,CS12,CS13,CCZ014,CCZ023:
        # 28 exact, 30 with free Cliffords); keep the smaller, both verified
        for fc in (True, False):
            g = EG.glue(oo, k, nc, L, free_clifford=fc)
            if g is None:
                continue
            f = dict(stage="entangled-gluing", k=k, r=nc, gates=g, params="")
            if not (us.is_factory(f, L) and us.distance(f, L) != 1):
                continue
            n = sum(us.inputs(c, L) for _, c in g)
            if best is None or n < best[0]:
                best = (n, f)
    return best and best[1]


def _pattern(args):
    """Restricted SAT, stage (iii) first form: one l=3 class, below cap. Two restrictions, both from the proof of
    Corollary 3 (fixing the checks leaves a synthesis of the target): check-pattern SAT (output parts fixed to the
    synthillation multiset, Boolean check patterns; pattern_sat.py) and synthesis SAT (output parts drawn from the
    columns of the target's syntheses, any check pattern; synthesis_sat.py, deterministic z3 step limit). Returns the
    smaller (n, N, gates) or None."""
    k, spec, cap = args
    import pattern_sat as PS, synthesis_sat as SS
    b, _ = PS.run_class(k, spec, cap)
    best = b and (b["n"], b["N"], b["gates"])
    c, _ = SS.solve(spec, k, best[0] if best else cap, rl=SYN_RLIMIT)
    if c and (best is None or (c["n"], c["N"]) < best[:2]):
        best = (c["n"], c["N"], c["gates"])
    return best


def run_cell(L, k, nmax, smax, rmax, budget, span, wmax, jobs, glue=True, sat=True, log=print, rlimit=None):
    us.S_MAX = smax
    if k >= APPROX_KMIN:                           # exact frame search costs ~1 h per factory at k >= 11
        import classify
        classify._GL_SAMPLE_CAP = APPROX_FRAMES
    T, rows = {}, []
    # (iii) first only to have the class map ready; its cost is reported as stage (iii)
    t = time.time()
    fiber = None
    if sat:
        try:
            res, (C, off, lab) = enumerate_targets(L, k, log=lambda *a: None)
        except ValueError:                               # orbit table too large: split the group in two
            if L == 3:                                   # degree filtration (l=3 k=6)
                from enumerate_targets_fiber import enumerate_targets_fiber
                res, fiber = enumerate_targets_fiber(L, k, log=lambda *a: None)
            else:                                        # 2-adic filtration (l=4 k=4)
                from enumerate_targets_2adic import enumerate_targets_2adic
                res, fiber = enumerate_targets_2adic(L, k, log=lambda *a: None)
            C = off = lab = None
        classes = {c["orbit"]: c for c in res["classes"] if c["genuine"]}
    else:                                                # beyond exact class lists: label by reduced output
        C = off = lab = None; classes = {}
    T["iii_enumerate"] = time.time() - t
    best = {}                                            # orbit -> (n, N, stage)
    T["verification"] = 0.0; T["classification"] = 0.0   # not search cost: reported separately

    def record(f, stage):
        t0 = time.time()
        ok = us.is_factory(f, L) and us.distance(f, L) != 1
        T["verification"] += time.time() - t0
        if not ok:
            return
        g = f["gates"]
        t0 = time.time()
        P = phase_polynomial([S for S, _ in g], [c for _, c in g], set(range(k)), L)
        if lab is not None:
            orb = int(lab[index_of_poly(P, k, L, C, off)]); genuine = orb in classes
            T["classification"] += time.time() - t0
        elif fiber is not None:                          # exact class via (quotient orbit, fiber orbit)
            lb = fiber.label(fiber.digits_of_poly(P))
            orb = lb if isinstance(lb, str) else f"{lb[0]}:{lb[1]}"; genuine = orb in classes
            T["classification"] += time.time() - t0
        else:                                            # no exact table: Clifford-reduced output (frame-sampled)
            info = us.classify_factory(f, L)
            if info is None:
                T["classification"] += time.time() - t0
                return
            orb = info["reduced"]; genuine = info["essential_dim"] == k
        T["classification"] += time.time() - t0
        n = sum(us.inputs(c, L) for _, c in g)
        rows.append(dict(stage=stage, orbit=orb, genuine=genuine, n=n, N=k + f["r"], params=f.get("params", "")))
        if genuine and (orb not in best or (n, k + f["r"]) < best[orb][:2]):
            best[orb] = (n, k + f["r"], stage)

    t = time.time(); seen = set()
    for f in us.stage_two_group(L, k, nmax, k, one_sign=True):
        key = (f["r"], frozenset(f["gates"]))
        if key not in seen:
            seen.add(key); record(f, "two-group")
    T["i_two_group"] = time.time() - t - T["verification"] - T["classification"]
    v0, c0 = T["verification"], T["classification"]
    t = time.time()
    if glue:
        for f in stage_gluing_multi(L, k, rmax):
            record(f, "gluing")
    T["ii_gluing"] = time.time() - t - (T["verification"] - v0) - (T["classification"] - c0)
    after_ii = dict(best)
    # (ii') gluing on entangled support, every class
    t = time.time(); v0, c0 = T["verification"], T["classification"]
    if sat and L in (2, 3, 4):
        with ProcessPoolExecutor(jobs) as ex:
            for f in ex.map(_eglue, [(L, k, c["spec"], rmax) for c in classes.values()]):
                if f:
                    f["gates"] = [(frozenset(S), cc) for S, cc in f["gates"]]
                    record(f, "entangled-gluing")
    T["ii_entangled_gluing"] = time.time() - t - (T["verification"] - v0) - (T["classification"] - c0)
    after_iip = dict(best)
    # (iv) check-pattern SAT (l=3), seeded with the best n so far
    t = time.time()
    if sat and L == 3:
        items = [(o, c) for o, c in classes.items()]
        with ProcessPoolExecutor(jobs) as ex:
            res = list(ex.map(_pattern, [(k, c["spec"], best[o][0] if o in best else 10 ** 6) for o, c in items]))
        for (o, c), r in zip(items, res):
            if r:
                f = dict(stage="pattern-sat", k=k, r=r[1] - k, gates=[(frozenset(S), cc) for S, cc in r[2]], params=c["spec"])
                record(f, "pattern-sat")
    T["iv_pattern_sat"] = time.time() - t
    after_iv = dict(best)
    # (v) quick SAT, seeded with the best n from every stage before it: asks only for what they did not find
    t = time.time()
    args = [(L, k, c["spec"], budget, span, wmax, best[o][0] if o in best else None, rlimit) for o, c in classes.items()]
    out = []
    if sat:
        with ProcessPoolExecutor(jobs) as ex:
            out = list(ex.map(solve_target, args))
    T["v_quick_sat"] = time.time() - t
    iv_new = iv_better = 0
    for (o, c), r in zip(list(classes.items())[:len(out)], out):
        if r["best"]:
            b = r["best"]
            rows.append(dict(stage="quick-sat", orbit=o, genuine=True, n=b["n"], N=b["N"], params=c["spec"],
                             gates=b["gates"]))
            if o in after_iv:
                iv_better += 1
            else:
                iv_new += 1
            if o not in best or (b["n"], b["N"]) < best[o][:2]:
                best[o] = (b["n"], b["N"], "quick-sat")
    summary = dict(l=L, k=k, classes=len(classes) if sat else None, reached_i_ii=len(after_ii), reached_iip=len(after_iip),
                   reached_iv=len(after_iv), v_new=iv_new, v_smaller=iv_better,
                   reached_total=len(best), seconds={x: round(v, 1) for x, v in T.items()},
                   best={str(o): dict(n=v[0], N=v[1], stage=v[2], spec=classes[o]["spec"] if o in classes else str(o))
                         for o, v in best.items()})
    log(f"l={L} k={k}: {len(classes)} classes | (i)+(ii) {len(after_ii)} | +(ii') {len(after_iip)} | +(iv) {len(after_iv)} | (v) +{iv_new} new, {iv_better} smaller | "
        f"total {len(best)} | " + ", ".join(f"{x} {v:.0f}s" for x, v in T.items()))
    return summary, rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", nargs="+", default=["3:3-5", "2:2-6", "4:2-3"])
    ap.add_argument("--nmax", type=int, default=11, help="max circuit qubits in two-group (main sweep 11, larger 18)")
    ap.add_argument("--nmax-l4", type=int, default=None, help="cap for l=4 (larger sweep 16)")
    ap.add_argument("--smax", type=int, default=7, help="two-group s_total, s_O (main sweep 7, larger 12)")
    ap.add_argument("--glue-kmax", type=int, default=99, help="run gluing only for k <= this (larger sweep 9)")
    ap.add_argument("--sat-kmax", type=int, default=99, help="run stages (iii)-(iv) only for k <= this")
    ap.add_argument("--out-dir", default=OUT)
    ap.add_argument("--rmax", type=int, default=4)
    ap.add_argument("--budget", type=float, default=25.0)
    ap.add_argument("--span", type=int, default=3)
    ap.add_argument("--wmax", type=int, default=6)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--rlimit", type=float, default=None, help="z3 step limit per stage-(v) query (deterministic)")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    for c in a.cells:
        L, ks = c.split(":"); L = int(L)
        lo, hi = map(int, ks.split("-"))
        for k in range(lo, hi + 1):
            nmax = a.nmax_l4 if (L == 4 and a.nmax_l4) else a.nmax
            summary, rows = run_cell(L, k, nmax, a.smax, a.rmax, a.budget, a.span, a.wmax, a.jobs,
                                     glue=k <= a.glue_kmax, sat=k <= a.sat_kmax,
                                     rlimit=int(a.rlimit) if a.rlimit else None)
            summary["params"] = dict(nmax=nmax, smax=a.smax, rmax=a.rmax, glue=k <= a.glue_kmax, sat=k <= a.sat_kmax,
                                     budget=a.budget, span=a.span, wmax=a.wmax, jobs=a.jobs, rlimit=a.rlimit)
            os.makedirs(a.out_dir, exist_ok=True)
            json.dump(dict(summary=summary, rows=rows), open(os.path.join(a.out_dir, f"l{L}_k{k}.json"), "w"), indent=1)
