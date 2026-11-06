"""The cascade run stage by stage: each stage is its own pass, reads the best factory per class from the pass before,
asks only for fewer inputs, and writes its own file, so the gain of every stage (and of every setting on its ladder)
is recorded separately and any pass can be rerun alone.

Passes, in cascade order (all verified with us.is_factory and us.distance != 1, and labelled with the exact class):
  i_tg   (i)   two-group, S_k x S_r symmetric                     setting: nmax (circuit qubits), smax
  i_gs   (i)   group-symmetric, outputs under G <= S_k (l=3)      setting: rmax; sweep over the group menu, then
                                                                  target mode G = Aut(U) and G = {e}, CP-SAT dlimit
  ii_dj  (ii)  gluing, disjoint blocks                            setting: rmax (checks)
  ii_eg  (ii)  gluing on entangled support, both output forms     setting: rmax
  iii_r  (iii) restricted SAT: check-pattern and synthesis SAT    settings: pattern R<=4, mmax; synthesis span, rlimit
  iii_f  (iii) free SAT over all rotations                        settings: span (N-k), wmax, rlimit
  Note: at l=3 the paper's stage (iii) is not iii_r/iii_f here but the XOR chain of searches/chain_iii2.py, assembled
  into outputs/staged_final by searches/build_staged_final.py; iii_r (check-pattern and z3 synthesis SAT) is the earlier
  form and not in the final pipeline. iii_f here (z3 quick_sat) is free SAT at l=2 and l=4, and the source of the
  all-coefficient minimality proofs at l=3.
A factory that deposits U^dagger instead of U is turned into one for U by negating every coefficient (the inverse
circuit is again a borrowed identity with the same n), and re-verified.

    python searches/staged.py --L 3 --k 4 --passes all --out outputs/staged --jobs 8
    -> outputs/staged/l3_k4/<pass>.json : settings, per-class best after the pass, gain over the previous pass, time
"""
import os, sys, json, time, argparse
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
import cascade as CA
import cascade_class as CC

PASSES = ["i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg", "iii_r", "iii_f"]
DEFAULT = dict(i_tg=dict(nmax=11, smax=7), i_gs=dict(rmax=4, dlimit=20.0, trivial_kmax=6), i_gh=dict(rmax=4, bits_max=24),
               ii_dj=dict(rmax=4, smax=7), ii_eg=dict(rmax=4),
               iii_r=dict(R=4, mmax=12, pat_rlimit=2e9, syn_span=3, syn_rlimit=3e9),
               iii_f=dict(span=3, wmax=6, rlimit=6e8))
_CTX = {}
CACHE = [None]                                   # class-cache dir of the current run (set in run)


def ctx(L, k, cache):
    if (L, k) not in _CTX:
        classes, c = CC.classes_of(L, k, cache)
        _CTX[(L, k)] = (classes, CC.labeller(L, k, c, classes))
    return _CTX[(L, k)]


def n_of(gates, L):
    return sum(us.inputs(c, L) for _, c in gates)


def label_fix(f, L, k, cache, want=None):
    """(class, gates) of a verified factory, or (None, None). As in cascade v4, a factory is credited to the class it
    actually deposits, whatever class it was searched for (`want` is accepted and ignored); Book.add also offers its
    inverse circuit to the inverse class at l>=4."""
    classes, label = ctx(L, k, cache)
    g = [(frozenset(S), c % (1 << L)) for S, c in f["gates"]]
    return label(dict(f, gates=g)), g


class Book:
    """Best factory per class, the pass that set it, and every factory seen (for the catalogue)."""
    @property
    def new(self):
        return {o for o in self.best if o not in self.prev}

    @property
    def lower(self):
        return {o: v for o, v in self.prev.items() if self.best[o]["n"] < v}

    def __init__(self, L, k, prev=None):
        self.L, self.k = L, k
        self.best = dict(prev["best"]) if prev else {}
        self.prev = {o: b["n"] for o, b in self.best.items()}   # gain is measured against this, not within the pass
        self.rows = []

    def add(self, o, g, N, pas, keep_row=False, partner=True):
        if o is None:
            return
        if partner and self.L >= 4:
            # the inverse circuit (every coefficient negated) is a factory with the same n for the inverse output,
            # whose class can differ from o at l=4 (at l<=3 a negated coefficient differs by a Clifford): offer it there too
            gi = [(S, (-c) % (1 << self.L)) for S, c in g]
            classes, label = ctx(self.L, self.k, CACHE[0])
            oi = label(dict(k=self.k, r=N - self.k, gates=gi))
            if oi is not None and oi != o:
                self.add(oi, gi, N, pas + ":inverse", keep_row, partner=False)
        n = n_of(g, self.L); cur = self.best.get(o)
        if keep_row:
            self.rows.append(dict(orbit=o, n=n, N=N, stage=pas))
        if cur is None or (n, N) < (cur["n"], cur["N"]):
            self.best[o] = dict(n=n, N=N, stage=pas, gates=[[sorted(S), c] for S, c in g])


# ----------------------------------------------------------------- passes
def p_i_tg(L, k, s, book, cache, jobs):
    us.S_MAX = s["smax"]; seen = set()
    for f in us.stage_two_group(L, k, s["nmax"], k, one_sign=True):
        key = (f["r"], frozenset(f["gates"]))
        if key in seen:
            continue
        seen.add(key); o, g = label_fix(f, L, k, cache)
        book.add(o, g, k + f["r"], "i_tg", keep_row=o is not None)


def p_ii_dj(L, k, s, book, cache, jobs):
    us.S_MAX = s["smax"]
    for f in CA.stage_gluing_multi(L, k, s["rmax"]):
        o, g = label_fix(f, L, k, cache)
        book.add(o, g, k + f["r"], "ii_dj", keep_row=o is not None)


def _eg_one(a):
    L, k, o, spec, rmax, cache = a
    f = CA._eglue((L, k, spec, rmax))
    if not f:
        return o, None
    oo, g = label_fix(f, L, k, cache)
    return oo, (g and dict(gates=[[sorted(S), c] for S, c in g], N=k + f["r"]))


def p_ii_eg(L, k, s, book, cache, jobs):
    classes, _ = ctx(L, k, cache)
    with ProcessPoolExecutor(jobs) as ex:
        for o, r in ex.map(_eg_one, [(L, k, o, c["spec"], s["rmax"], cache) for o, c in classes.items()]):
            if r:
                book.add(o, [(frozenset(S), c) for S, c in r["gates"]], r["N"], "ii_eg", keep_row=True)


def _r_one(a):
    k, o, spec, cap, s = a
    import pattern_sat as PS, synthesis_sat as SS
    out = []
    b, _ = PS.run_class(k, spec, cap, Rs=tuple(range(1, s["R"] + 1)), mmax=s["mmax"])
    if b:
        out.append(("pattern", b)); cap = min(cap, b["n"])
    c, _ = SS.solve(spec, k, cap, rl=s["syn_rlimit"], span=tuple(range(1, s["syn_span"] + 1)))
    if c:
        out.append(("synthesis", c))
    return o, out


def p_iii_r(L, k, s, book, cache, jobs):
    if L != 3:
        return
    classes, _ = ctx(L, k, cache)
    todo = [(k, o, c["spec"], book.best[o]["n"] if o in book.best else 10 ** 6, s) for o, c in classes.items()
            if not at_bound(book, o, c)]
    with ProcessPoolExecutor(jobs) as ex:
        for o, out in ex.map(_r_one, todo):
            for kind, b in out:
                f = dict(k=k, r=b["N"] - k, gates=b["gates"])
                oo, g = label_fix(f, L, k, cache)
                book.add(oo, g, b["N"], "iii_r:" + kind)


def _f_one(a):
    from quick_sat import solve_target
    L, k, o, spec, cap, s = a
    return o, solve_target((L, k, spec, 25.0, s["span"], s["wmax"], cap, int(s["rlimit"])))


def p_iii_f(L, k, s, book, cache, jobs):
    classes, _ = ctx(L, k, cache)
    todo = [(L, k, o, c["spec"], book.best[o]["n"] if o in book.best else None, s) for o, c in classes.items()
            if not at_bound(book, o, c)]
    with ProcessPoolExecutor(jobs) as ex:
        for o, r in ex.map(_f_one, todo):
            b = r["best"]
            if b:
                f = dict(k=k, r=b["N"] - k, gates=b["gates"])
                oo, g = label_fix(f, L, k, cache)
                book.add(oo, g, b["N"], "iii_f")
            book.steps = getattr(book, "steps", {}); book.steps[o] = r["steps"]


def p_i_gs(L, k, s, book, cache, jobs):
    if L != 3:
        return
    import orbit_group as OG
    labels = OG.Labels(k, os.path.join(cache, "l3_k6_fiber.pkl") if k == 6 else None)
    sb, _ = OG.run_k(k, s["rmax"], labels, lambda *a: None)
    tb, _ = OG.run_targets(k, s["rmax"], labels, lambda *a: None, s["dlimit"], trivial=k <= s["trivial_kmax"])
    for src in (sb, tb):
        for o, b in src.items():
            f = dict(k=k, r=b["N"] - k, gates=b["gates"])
            oo, g = label_fix(f, L, k, cache)
            book.add(oo, g, b["N"], "i_gs", keep_row=True)


def p_i_gh(L, k, s, book, cache, jobs):
    """(i) group-symmetric with a group on the checks too: outputs under G <= S_k, checks under H <= S_r (orbit_group2).
    Exhaustive DP sweep over each (G, H, r) whose syndrome space is <= 2^bits_max; every optimum is expanded, verified
    and credited to the class it deposits."""
    if L != 3:
        return
    import orbit_group as OG, orbit_group2 as OG2
    for name, G, orbs, of, _ in OG.group_menu(k):
        for r in range(1, s["rmax"] + 1):
            for hn, Hg, corbs, cof in OG2.check_groups(r):
                if hn == "H:S" + str(r):
                    continue                                  # H = S_r is pass i_gs
                res, info = OG2.with_cells(k, r, orbs, corbs, OG.dp_sweep, s["bits_max"])
                if res is None:
                    book.dropped = getattr(book, "dropped", 0) + 1
                    continue
                out, cells, cost, Dorb = res
                for dig, (n, x) in out.items():
                    if not x:
                        continue
                    g = OG2.gates2(k, r, orbs, corbs, cells, x)
                    f = dict(k=k, r=r, gates=g)
                    if not (us.is_factory(f, L) and us.distance(f, L) != 1):
                        continue
                    o, gg = label_fix(f, L, k, cache)
                    book.add(o, gg, k + r, "i_gh:" + name + "x" + hn[2:], keep_row=o is not None)


TAU = {}


def tau_of(L, k, o, spec):
    """Exact T-count (l=3, class_table_l3.json) or sqrt(T)-count (l=4, entangled_gluing_l4_k*.json 'tau'); else None."""
    if (L, k) not in TAU:
        t = {}
        if L == 3:
            for r in json.load(open(os.path.join(ROOT, "outputs", "class_table_l3.json"))):
                if r["k"] == k:
                    t[r["spec"]] = r["tcount"]
        elif L == 4:
            p = os.path.join(ROOT, "outputs", f"entangled_gluing_l4_k{k}.json")
            if os.path.exists(p):
                for r in json.load(open(p)):
                    t[r["spec"]] = r.get("tau")
        TAU[(L, k)] = t
    return TAU[(L, k)].get(spec)


def at_bound(book, o, c):
    """True if the best so far meets n = tau(U) + 1 (Corollary 3): no later stage can lower it."""
    t = tau_of(book.L, book.k, o, c["spec"]); b = book.best.get(o)
    return b is not None and t is not None and b["n"] == t + 1


FN = dict(i_tg=p_i_tg, i_gs=p_i_gs, i_gh=p_i_gh, ii_dj=p_ii_dj, ii_eg=p_ii_eg, iii_r=p_iii_r, iii_f=p_iii_f)


def run(L, k, passes, out, jobs, settings):
    d = os.path.join(out, f"l{L}_k{k}"); os.makedirs(d, exist_ok=True)
    cache = os.path.join(out, "_cache"); os.makedirs(cache, exist_ok=True)
    CACHE[0] = cache
    classes, _ = ctx(L, k, cache)
    prev = None
    for p in PASSES:
        path = os.path.join(d, f"{p}.json")
        if p not in passes:
            if os.path.exists(path):
                prev = json.load(open(path))
            continue
        s = dict(DEFAULT[p], **settings.get(p, {}))
        book = Book(L, k, prev); t = time.time()
        FN[p](L, k, s, book, cache, jobs)
        rec = dict(L=L, k=k, pass_=p, settings=s, seconds=round(time.time() - t, 1), classes=len(classes),
                   dropped=getattr(book, "dropped", 0),
                   reached=len(book.best), new=sorted(book.new), lowered={o: [v, book.best[o]["n"]] for o, v in book.lower.items()},
                   inputs_saved=sum(v - book.best[o]["n"] for o, v in book.lower.items()),
                   best=book.best, rows=book.rows, steps=getattr(book, "steps", None),
                   specs={o: c["spec"] for o, c in classes.items()})
        json.dump(rec, open(path + ".tmp", "w")); os.replace(path + ".tmp", path)
        print(f"l={L} k={k} {p}: reached {len(book.best)}/{len(classes)} (+{len(book.new)} new), lowered "
              f"{len(book.lower)} by {rec['inputs_saved']} inputs, {rec['seconds']}s", flush=True)
        prev = rec


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, required=True); ap.add_argument("--k", type=int, nargs="+", required=True)
    ap.add_argument("--passes", nargs="+", default=["all"]); ap.add_argument("--out", default=os.path.join(ROOT, "outputs", "staged"))
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--set", default="{}", help='JSON overrides, e.g. {"iii_f": {"rlimit": 2e9}}')
    a = ap.parse_args()
    ps = PASSES if a.passes == ["all"] else a.passes
    for k in a.k:
        run(a.L, k, ps, a.out, a.jobs, json.loads(a.set))
