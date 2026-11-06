"""Gluing on entangled support: stage (ii) with the output-only gates set to a T-count-optimal parity network of the
target instead of symmetric all-weights blocks. The check couplings are then solved by the same 2-adic linear system
(symfree_search.solve_graded), so the factory comes from the borrowed-identity condition alone; no SAT.

For each class in outputs/class_table_l3.json: output-only gates = the optimal network A (tau rotations,
synthillation_compare.opt_A); for r = 1..R check qubits solve for the check-touching gates, keep the smallest n, and
verify with is_factory and by comparing the deposited output with the target up to a diagonal Clifford.
    python searches/entangled_gluing.py --workers 8 [--rmax 3]   -> outputs/entangled_gluing_l3.json
"""
import os, sys, json, argparse
from itertools import combinations
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))


def glue(oo, k, ncheck, L, free_clifford=True):
    """oo: the output to deposit, as parity rotations [(support, coeff)] on the output block; solve for the
    check-touching gates that deposit it (up to sign) and are the identity elsewhere. Returns those gates or None."""
    from symfree_search import solve_graded, reduce_support
    O = list(range(k)); C = list(range(k, k + ncheck)); allq = O + C
    cand = [frozenset(s) for r in range(1, len(allq) + 1) for s in combinations(allq, r) if set(s) & set(C)]
    cand.sort(key=len); nvars = len(cand)
    masks, rs, ms, varlists = [], [], [], []
    for t in range(1, L + 1):
        mT = L - t + 1
        for T in combinations(allq, t):
            Ts = frozenset(T)
            vl = [i for i, g in enumerate(cand) if Ts <= g]
            mask = 0
            for i in vl:
                mask |= 1 << i
            sig = sum(c for g, c in oo if Ts <= g)
            if free_clifford and Ts <= set(O):
                # output-only monomial: a free diagonal Clifford shifts its coefficient (-2)^(t-1) sig by multiples
                # of 2^cliff_exp(t, L), so only sig mod 2^(cliff_exp - t + 1) is fixed (mod 2 for every t < 3 at l=3)
                from unified_search import cliff_exp
                mT = min(mT, max(cliff_exp(t, L) - t + 1, 0))
            masks.append(mask); rs.append((-sig) % (1 << mT)); ms.append(mT); varlists.append(vl)
    sol = solve_graded(masks, rs, ms, nvars, L, varlists)
    if sol is None:
        return None
    x, nulls = sol
    x = reduce_support(rs, ms, x, nulls, nvars, L, varlists)
    return [(cand[i], x[i] % (1 << L)) for i in range(nvars) if x[i] % (1 << L)]   # only check-touching gates, as in
                                                                                     # solve_binding: oo fixes what they deposit


def deposits(gates, k, target, L=3):
    """Does the circuit's output-block phase polynomial equal target (or its inverse) up to a diagonal Clifford?"""
    from classify import phase_polynomial
    P = phase_polynomial([S for S, _ in gates], [c for _, c in gates], set(range(k)), L)
    keys = set(P) | set(target)
    for sgn in (1, -1):
        ok = True
        for T in keys:
            d = (P.get(T, 0) - sgn * target.get(T, 0)) % 8
            need = {1: 2, 2: 4, 3: 8}.get(len(T), 1)
            if len(T) <= 3 and d % need:
                ok = False; break
        if ok:
            return True
    return False


def one(args):
    r, rmax = args
    import class_table_l3 as CT, synthillation_compare as sc, unified_search as us
    k = r["k"]; tgt = CT.parse(r["spec"]); A = sc.opt_A(k, tgt)
    oo = [(frozenset(q for q in range(k) if y >> q & 1), 1) for y in A]
    best = None
    for nc in range(1, rmax + 1):
        g = glue(oo, k, nc, 3)
        if g is None:
            continue
        f = dict(stage="entangled-gluing", k=k, r=nc, gates=g, params="")
        if not (us.is_factory(f, 3) and us.distance(f, 3) != 1 and deposits(g, k, tgt)):
            continue
        n = sum(us.inputs(c, 3) for _, c in g)
        if best is None or (n, k + nc) < (best["n"], best["N"]):
            best = dict(n=n, N=k + nc, gates=[[sorted(S), c] for S, c in g])
    return dict(k=k, orbit=r["orbit"], spec=r["spec"], tcount=r["tcount"], cascade_n=r["n"], syn_n=r["syn_n"],
                eg_n=best and best["n"], eg_N=best and best["N"], gates=best and best["gates"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=1); ap.add_argument("--rmax", type=int, default=3)
    ap.add_argument("--only", default=""); a = ap.parse_args()
    R = json.load(open(os.path.join(ROOT, "outputs", "class_table_l3.json")))
    if a.only:
        R = [r for r in R if r["spec"] == a.only]
    with Pool(a.workers) as p:
        out = list(p.imap_unordered(one, [(r, a.rmax) for r in R]))
    if a.only:
        print(out); sys.exit()
    json.dump(out, open(os.path.join(ROOT, "outputs", "entangled_gluing_l3.json"), "w"))
    print("done", len(out))
