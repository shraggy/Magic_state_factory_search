"""check_A.py for a sweep-A copy in another directory (e.g. outputs/solver_v2/staged_A_v2): the same checks (valid,
brute force for N <= 9, exact class, n recomputed, distance, leading coefficient c) on SRC/l{L}_k{k}/iii_f.json
(else iii_r.json), with the class cache CACHE (default outputs/staged_A/_cache). Results to
SRC/checks/l{L}_k{k}[_c{I}of{M}].json.
    python searches/check_A2.py SRC L k [I M] [--cache DIR]
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
import cascade_class as CC
from error_coefficients import coefficients

BRUTE_NMAX = 9


def main(srcdir, L, k, I=0, M=1, cache=None):
    d = os.path.join(srcdir, f"l{L}_k{k}")
    src = os.path.join(d, "iii_f.json") if os.path.exists(os.path.join(d, "iii_f.json")) else os.path.join(d, "iii_r.json")
    A = json.load(open(src))
    classes, ctx = CC.classes_of(L, k, cache or os.path.join(ROOT, "outputs", "staged_A", "_cache"))
    label = CC.labeller(L, k, ctx, classes)
    rows, t0 = {}, time.time()
    for i, o in enumerate(sorted(A["best"])):
        if i % M != I:
            continue
        b = A["best"][o]; g = [(frozenset(S), c) for S, c in b["gates"]]; N = b["N"]
        f = dict(k=k, r=N - k, gates=g)
        valid = us.is_factory(f, L)
        brute = us.is_factory(f, L, brute_force=True) if N <= BRUTE_NMAX else None
        n = sum(us.inputs(c, L) for _, c in g)
        dist = us.distance(f, L)
        try:
            n2, c = coefficients(g, k, L)
        except AssertionError:
            n2, c = n, None
        rows[o] = dict(spec=A["specs"][o], n=n, n_filed=b["n"], N=N, stage=b["stage"], valid=valid, brute=brute,
                       class_ok=label(f) == o, d=str(dist), c=c)
    summ = dict(source=os.path.relpath(src, ROOT), factories=len(rows),
                valid=sum(r["valid"] for r in rows.values()),
                brute_checked=sum(r["brute"] is not None for r in rows.values()),
                brute_ok=sum(bool(r["brute"]) for r in rows.values()),
                class_ok=sum(r["class_ok"] for r in rows.values()),
                n_ok=sum(r["n"] == r["n_filed"] for r in rows.values()),
                d_ge2=sum(r["d"] != "1" for r in rows.values()),
                d_dist={dd: sum(r["d"] == dd for r in rows.values()) for dd in sorted({r["d"] for r in rows.values()})},
                c_computed=sum(r["c"] is not None for r in rows.values()), seconds=round(time.time() - t0, 1))
    out = os.path.join(srcdir, "checks"); os.makedirs(out, exist_ok=True)
    name = f"l{L}_k{k}" + (f"_c{I}of{M}" if M > 1 else "") + ".json"
    json.dump(dict(summary=summ, rows=rows), open(os.path.join(out, name), "w"), indent=0)
    print(f"l={L} k={k} [{I}/{M}]: {json.dumps(summ)}", flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]; cache = None
    if "--cache" in a:
        i = a.index("--cache"); cache = a[i + 1]; a = a[:i] + a[i + 2:]
    main(a[0], *[int(x) for x in a[1:]], cache=cache)
