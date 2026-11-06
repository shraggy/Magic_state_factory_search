"""Pass ii_eg (gluing on entangled support, staged.py / cascade._eglue) on a TARGET LIST of CCZ layers at k=7, 8,
where the pipeline has no class list (cascade_class.classes_of) and no exact synthesis (synthillation_compare.opt_A
enumerates the 2^dim Reed-Muller coset; dim 64 at k=7). The only changes from cascade._eglue:
  - the target comes from the list, not from a class table (no class label; the factory is checked against the layer
    itself with entangled_gluing.deposits, up to a diagonal Clifford and inverse, as _eglue's callers do);
  - the output synthesis is frontier_k7.bound (greedy descent over the coset, fixed seed, fixed restarts:
    deterministic), an upper bound w on tau(U), instead of the exact coset leader opt_A.
Added: the one-check form (all columns through one check), a particular solution for r = 1 that glue's greedy
reduction can miss (target front_k8_44 of outputs/k78_pipeline/targets_paper62.json: glue n = 20, one-check n = 18).
The check couplings are then solved by the same entangled_gluing.glue (2-adic linear system, both free-Clifford and
exact forms, r = 1..rmax checks), verified with us.is_factory, us.distance != 1 and deposits. The isotropy floor
(isotropy_floor.floor, a proved lower bound on tau) is recorded: n == floor + 1 means the factory is on the bound.
    python searches/k78_eg.py TARGETS.json I M OUTDIR [--rmax 4 --restarts 400 --seed 1 --isd ITERS]
--isd ITERS: also lower the compilation by information-set decoding (k78_tau_upper.isd, p = 3, fixed seeds, ITERS
iterations per step); the gluing and one-check forms are then run on both syntheses.
--exact7: at k = 7, decide tau exactly with tcount7.exists_below2 when the floor is below the synthesis (tcount_lo is
then the exact T-count); at k = 8 tcount_lo stays the isotropy floor.
TARGETS.json: list of dict(id, k, ccz, [ref_n, ref_tcount, src]). Targets with index % M == I -> OUTDIR/c{I}of{M}.json
"""
import os, sys, json, time, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import unified_search as us
import entangled_gluing as EG
from frontier_k7 import bound
from isotropy_floor import floor


def _bin(T):
    """tcount7.BIN is the binary compiled on the laptop (arm64 Mach-O); on another platform use a local build of the
    same source, searches/csrc_<machine>/capped_minweight (cc -O3 csrc/capped_minweight.c -lpthread), if present."""
    import platform
    alt = os.path.join(HERE, f"csrc_{platform.machine()}", "capped_minweight")
    if platform.system() != "Darwin" and os.path.exists(alt):
        T.BIN = alt


def isd_tighten(layer, k, w, cols, lo, iters, p=3):
    """Lower the compilation by information-set decoding (k78_tau_upper.isd, seeds 0, 1, ...: deterministic) until a
    fixed-budget search fails or w reaches the isotropy floor. Each word is checked to lie in the layer's coset."""
    import tcount7 as T
    from k78_tau_upper import isd, in_coset
    b, C0, nc = T.coset({frozenset(x): 4 for x in layer}, k); steps = []
    while w > lo:
        word, seen, its = isd(b, C0, nc, w - 1, iters, p, seed=len(steps))
        ok = bool(word is not None and seen < w and in_coset(word, b, C0))
        steps.append([w - 1, its, seen if ok else None])
        if not ok:
            break
        w, cols = seen, [y for y in range(1, nc + 1) if word >> (y - 1) & 1]
    return w, cols, steps


def eg_target(k, ccz, rmax=4, restarts=400, seed=1, isd_iters=0, exact7=False):
    layer = [tuple(x) for x in ccz]; tgt = {frozenset(x): 4 for x in layer}
    t0 = time.process_time()
    w, cols = bound([list(x) for x in layer], k, restarts=restarts, seed=seed)
    t_syn = time.process_time() - t0
    lo = floor([list(x) for x in layer], k)
    w0 = w; isd_steps = []
    syns = [cols]
    if isd_iters and w > lo:
        t2 = time.process_time(); w, cols2, isd_steps = isd_tighten(layer, k, w, cols, lo, isd_iters)
        t_syn += time.process_time() - t2
        if cols2 is not cols:
            syns.append(cols2)
    exact = lo == w
    if exact7 and k == 7 and w > lo:
        # exact decision at k = 7 (tcount7.exists_below2, two information sets, compiled search; one thread): no coset
        # word of weight < w, with the returned proved bound reaching w, fixes tau = w
        import tcount7 as T
        _bin(T)
        b, C0, nc = T.coset(tgt, k); t3 = time.process_time()
        found, seen, proved = T.exists_below2(b, C0, nc, w, 1); t_syn += time.process_time() - t3
        if not found and proved is not None and proved >= w:
            lo, exact = w, True
    best, tries = None, []
    t1 = time.process_time()

    def offer(g, nc, kind):
        nonlocal best
        f = dict(stage=kind, k=k, r=nc, gates=g, params="")
        ok = us.is_factory(f, 3) and us.distance(f, 3) != 1 and EG.deposits(g, k, tgt)
        ok = ok and us.is_factory(f, 3, brute_force=True)            # independent check on all 2^N basis states
        n = sum(us.inputs(c, 3) for _, c in g)
        tries.append([si, nc, kind, n if ok else -n])
        if ok and (best is None or (n, k + nc) < (best["n"], best["N"])):
            best = dict(n=n, N=k + nc, distance=us.distance(f, 3), kind=kind, synthesis=si,
                        gates=[[sorted(S), c] for S, c in g])

    for si, cols in enumerate(syns):
        oo = [(frozenset(q for q in range(k) if y >> q & 1), 1) for y in cols]
        for nc in range(1, rmax + 1):
            for fc in (True, False):
                g = EG.glue(oo, k, nc, 3, free_clifford=fc)
                if g is None:
                    tries.append([si, nc, "glue" + ("" if fc else ":exact"), None])
                else:
                    offer(g, nc, "glue" + ("" if fc else ":exact"))
        # the one-check form of Corollary 3: every column of the synthesis through one check (+ one check-only
        # rotation if the count is odd). A particular solution of glue's r = 1 system that glue's greedy kernel
        # reduction (symfree_search.reduce_support) does not always reach; offered as one more candidate.
        offer([(frozenset(oS | {k}), 1) for oS, _ in oo] + ([(frozenset({k}), 1)] if len(oo) % 2 else []), 1, "one-check")
    t_glue = time.process_time() - t1
    return dict(k=k, ccz=[list(x) for x in layer], tcount_hi=w, tcount_hi_greedy=w0, tcount_lo=lo, tcount_exact=exact, isd_steps=isd_steps,
                syn_cols=[len(c) for c in syns], best=best,
                tries=tries, on_bound=bool(best and lo == w and best["n"] == w + 1),
                cpu_s=dict(synthesis=round(t_syn, 2), glue_verify=round(t_glue, 2)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("targets"); ap.add_argument("I", type=int); ap.add_argument("M", type=int)
    ap.add_argument("out"); ap.add_argument("--rmax", type=int, default=4); ap.add_argument("--restarts", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1); ap.add_argument("--isd", type=int, default=0)
    ap.add_argument("--exact7", action="store_true"); a = ap.parse_args()
    T = json.load(open(a.targets)); T = T["targets"] if isinstance(T, dict) else T; os.makedirs(a.out, exist_ok=True)
    dst = os.path.join(a.out, f"c{a.I}of{a.M}.json")
    if os.path.exists(dst):
        print("exists", dst); sys.exit()
    res = []
    for i, t in enumerate(T):
        if i % a.M != a.I:
            continue
        r = eg_target(t["k"], t["ccz"], a.rmax, a.restarts, a.seed, a.isd, a.exact7); r.update({x: t[x] for x in t if x not in r})
        res.append(r)
        b = r["best"]
        print(f"{t.get('id')} k={t['k']}: tau in [{r['tcount_lo']},{r['tcount_hi']}] n={b and b['n']} N={b and b['N']}"
              f" ref_n={t.get('ref_n')} {r['cpu_s']}", flush=True)
    json.dump(dict(settings=vars(a), host=os.uname().nodename, results=res), open(dst + ".tmp", "w"))
    os.replace(dst + ".tmp", dst)
