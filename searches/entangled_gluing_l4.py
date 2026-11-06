"""Gluing on entangled support at l=4: for every genuine l=4 output class (outputs/targets/l4_k{k}.json), the output
part is a sqrt(T)-count-optimal parity network of the target (Z(pi/16) rotations with coefficients mod 16, found by
the exact coset search of unified_search.zcount_exact), and the check couplings are solved by the 2-adic linear
system of entangled_gluing.glue at L=4. Every factory is verified with is_factory, the distance check, and by
comparing its deposited output with the target up to a diagonal Clifford (level 4). Compared with the cascade's best
n per class (outputs/cascade/l4_k{k}.json) and with the exact count tau (n >= tau + 1).
    python searches/entangled_gluing_l4.py --k 3 --workers 8      -> outputs/entangled_gluing_l4_k3.json
"""
import os, sys, json, argparse
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
L, MOD, BASE = 4, 16, 4


def network(P, k, L=L):
    """(count, {parity mask: coeff mod 2^L}) of a minimum-count network for P, or None (exact, coset <= cap)."""
    import unified_search as us
    MOD, BASE = 1 << L, 1 << max(L - 2, 1)
    n = (1 << k) - 1
    M = [[bin(x & y).count("1") & 1 for y in range(1, n + 1)] for x in range(1, n + 1)]
    f = [sum(c for T, c in P.items() if all(x >> q & 1 for q in T)) % MOD for x in range(1, n + 1)]
    U, W, D = us.snf_z8(M, MOD); cv = us._matvec_mod(U, f, MOD); ap = [0] * n; gens = []
    for i in range(n):
        d, ci = D[i] % MOD, cv[i] % MOD
        if d == 0:
            if ci:
                return None
            gens.append([W[r][i] % MOD for r in range(n)]); continue
        sols = [v for v in range(MOD) if (d * v) % MOD == ci]
        if not sols:
            return None
        ap[i] = sols[0]; st = MOD // (d & -d)
        if st < MOD:
            gens.append([(W[r][i] * st) % MOD for r in range(n)])
    a0 = tuple(x % MOD for x in us._matvec_mod(W, ap, MOD))
    gens = [g for g in gens if any(x % BASE for x in g)]
    # the count depends on a mod 4; enumerate the coset mod 4, keep a full mod-16 representative of the best
    cost = lambda v: sum(min(x % BASE, BASE - x % BASE) for x in v)
    span = {tuple(x % BASE for x in a0): a0}
    for g in gens:
        grown = dict(span)
        for key, v in span.items():
            w = list(v)
            for _ in range(BASE - 1):
                w = [(w[r] + g[r]) % MOD for r in range(n)]
                kk = tuple(x % BASE for x in w)
                if kk not in grown:
                    grown[kk] = tuple(w)
        span = grown
        if len(span) > 1 << 21:
            return None
    key = min(span, key=cost); v = span[key]
    return cost(key), {y + 1: c for y, c in enumerate(v) if c % MOD}


def deposits(gates, k, target):
    """Output-block phase polynomial equals target (or its inverse) up to a diagonal Clifford at level 4."""
    from classify import phase_polynomial
    import unified_search as us
    P = phase_polynomial([S for S, _ in gates], [c for _, c in gates], set(range(k)), L)
    keys = set(P) | set(target)
    for sgn in (1, -1):
        if all((P.get(T, 0) - sgn * target.get(T, 0)) % (1 << us.cliff_exp(len(T), L)) == 0 for T in keys if len(T) <= L):
            return True
    return False


def one(args):
    k, cls, rmax = args
    from coset_certify import parse_target
    import unified_search as us, entangled_gluing as EG
    tgt = parse_target(cls["spec"], L)
    net = network(tgt, k)
    if net is None:
        return dict(k=k, orbit=cls["orbit"], spec=cls["spec"], tau=None, eg_n=None)
    tau, par = net
    oo = [(frozenset(q for q in range(k) if y >> q & 1), c) for y, c in par.items()]
    best = None
    for nc in range(1, rmax + 1):
        g = EG.glue(oo, k, nc, L)
        if g is None:
            continue
        f = dict(stage="entangled-gluing", k=k, r=nc, gates=g, params="")
        if not (us.is_factory(f, L) and us.distance(f, L) != 1 and deposits(g, k, tgt)):
            continue
        n = sum(us.inputs(c, L) for _, c in g)
        if best is None or (n, k + nc) < (best["n"], best["N"]):
            best = dict(n=n, N=k + nc, gates=[[sorted(S), c] for S, c in g])
    return dict(k=k, orbit=cls["orbit"], spec=cls["spec"], tau=tau, eg_n=best and best["n"], eg_N=best and best["N"],
                gates=best and best["gates"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--workers", type=int, default=1); ap.add_argument("--rmax", type=int, default=3)
    ap.add_argument("--limit", type=int, default=0); a = ap.parse_args()
    cls = [c for c in json.load(open(os.path.join(ROOT, "outputs", "targets", f"l4_k{a.k}.json")))["classes"] if c["genuine"]]
    if a.limit:
        cls = cls[:a.limit]
    try:
        cas = json.load(open(os.path.join(ROOT, "outputs", "cascade", f"l4_k{a.k}.json")))["summary"]["best"]
    except FileNotFoundError:
        cas = {}
    with Pool(a.workers) as p:
        out = list(p.imap_unordered(one, [(a.k, c, a.rmax) for c in cls]))
    for o in out:
        b = cas.get(str(o["orbit"])); o["cascade_n"] = b and b["n"]
    json.dump(out, open(os.path.join(ROOT, "outputs", f"entangled_gluing_l4_k{a.k}.json"), "w"))
    found = [o for o in out if o["eg_n"]]
    print(f"l=4 k={a.k}: {len(out)} classes, factory for {len(found)}; at bound {sum(o['eg_n'] == o['tau'] + 1 for o in found)}; "
          f"cascade had {sum(1 for o in out if o['cascade_n'])}; EG < cascade {sum(1 for o in found if o['cascade_n'] and o['eg_n'] < o['cascade_n'])}; "
          f"EG new {sum(1 for o in found if not o['cascade_n'])}", flush=True)
