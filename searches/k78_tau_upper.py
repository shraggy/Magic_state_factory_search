"""T-count upper bounds read off the factories of searches/k78_eg.py, and a one-sided search for lower-weight
compilations.

(1) Corollary 3: fixing the check of a one-check factory to |0> leaves the output parts of its odd (input-consuming)
rotations as a synthesis of the target. For each factory in outputs/k78_pipeline/<run>/c*.json this script takes those
parities, checks independently that they are a synthesis of the CCZ layer, in two ways: (a) their phase polynomial
equals the layer up to a diagonal Clifford and inverse (entangled_gluing.deposits); (b) their indicator vector lies
in the layer's Reed-Muller coset b + span(C0) (tcount7.coset, GF(2) rank test). Its weight is then an upper bound on
tau(U). The isotropy floor (isotropy_floor.floor) is the proved lower bound.
(2) --isd ITER: Lee-Brickell information-set decoding (p <= 3, fixed seed, fixed iteration count: deterministic) for
a coset word of weight < the upper bound. One-sided: finding one lowers the bound (and is checked as in (1b));
not finding one proves nothing. Known-answer check: run on layers whose T-count is fixed (floor = bound).
    python searches/k78_tau_upper.py RUN [--isd 20000 --only ID,ID]   -> outputs/k78_pipeline/tau_upper_<RUN>.json
"""
import os, sys, json, glob, random, argparse, itertools, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import numpy as np
import tcount7 as T
from entangled_gluing import deposits
from isotropy_floor import floor


def rank(vs):
    B = []
    for v in vs:
        for b in B:
            v = min(v, v ^ b)
        if v:
            B.append(v)
    return len(B)


def in_coset(v, b, C0):
    return rank(C0 + [v ^ b]) == rank(C0)


def bits(v, n):
    return np.array([(v >> i) & 1 for i in range(n)], dtype=np.uint8)


def isd(b, C0, nc, target_w, iters, p=3, seed=0):
    """Lee-Brickell: search for a word of weight <= target_w in b + span(C0). Returns (int word or None, best seen)."""
    rng = np.random.default_rng(seed); G0 = np.array([bits(g, nc) for g in C0]); b0 = bits(b, nc); d = len(C0)
    best = (int(b0.sum()), b)                      # the coset base itself is a word
    for it in range(iters):
        perm = rng.permutation(nc); G = G0[:, perm].copy(); bb = b0[perm].copy(); piv = []; r = 0
        for c in range(nc):
            if r == d:
                break
            rows = np.nonzero(G[r:, c])[0]
            if len(rows) == 0:
                continue
            q = r + rows[0]
            if q != r:
                G[[r, q]] = G[[q, r]]
            m = G[:, c].astype(bool); m[r] = False; G[m] ^= G[r]
            if bb[c]:
                bb ^= G[r]
            piv.append(c); r += 1
        if r < d:
            continue
        rest = np.array([c for c in range(nc) if c not in set(piv)])
        R = G[:, rest]; br = bb[rest]
        # weight of b' + sum_{i in S} G_i = |S| + |br + sum R_i| (b' is zero on the pivots)
        Rp = np.packbits(R, axis=1); bp = np.packbits(br)
        pc = np.unpackbits(np.arange(256, dtype=np.uint8)[:, None], axis=1).sum(1)
        w0 = int(pc[bp].sum())
        cands = [((), w0)]
        x1 = Rp ^ bp; w1 = pc[x1].sum(1) + 1
        cands += [((i,), int(w1[i])) for i in np.nonzero(w1 <= target_w)[0]]
        if p >= 2:
            for i in range(d):
                x2 = x1[i] ^ Rp[i + 1:]; w2 = pc[x2].sum(1) + 2
                cands += [((i, i + 1 + j), int(w2[j])) for j in np.nonzero(w2 <= target_w)[0]]
                if p >= 3:
                    for j in range(i + 1, d):
                        x3 = x1[i] ^ Rp[j] ^ Rp[j + 1:]; w3 = pc[x3].sum(1) + 3
                        cands += [((i, j, j + 1 + l), int(w3[l])) for l in np.nonzero(w3 <= target_w)[0]]
        for S, w in cands:
            if w < best[0]:
                word = bb.copy()
                for i in S:
                    word ^= G[i]
                v = np.zeros(nc, dtype=np.uint8); v[perm] = word
                best = (int(v.sum()), sum(1 << i for i in range(nc) if v[i]))
        if best[0] <= target_w:
            return best[1], best[0], it + 1
    return None, best[0], iters


def one_check(cols, k, layer):
    """The one-check factory of a synthesis (every column through one check, + a check-only rotation if odd),
    verified with us.is_factory (monomial form and brute force on all 2^N states), distance, deposits."""
    import unified_search as us
    g = [(frozenset({q for q in range(k) if y >> q & 1} | {k}), 1) for y in cols]
    if len(g) % 2:
        g.append((frozenset({k}), 1))
    f = dict(stage="one-check", k=k, r=1, gates=g, params="")
    ok = bool(us.is_factory(f, 3) and us.is_factory(f, 3, brute_force=True) and us.distance(f, 3) == 2 and deposits(g, k, layer))
    return dict(n=sum(us.inputs(c, 3) for _, c in g), N=k + 1, verified=ok, gates=[[sorted(S), c] for S, c in g])


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("--isd", type=int, default=0)
    ap.add_argument("--p", type=int, default=2); ap.add_argument("--only", default=""); ap.add_argument("--tag", default="")
    a = ap.parse_args(); only = set(a.only.split(",")) - {""}
    R = []
    for f in sorted(glob.glob(os.path.join(ROOT, "outputs", "k78_pipeline", a.run, "c*.json"))):
        R += json.load(open(f))["results"]
    out = []
    for r in R:
        if only and r.get("id") not in only:
            continue
        k = r["k"]; layer = {frozenset(x): 4 for x in r["ccz"]}; g = r["best"]["gates"]
        syn = [(frozenset(q for q in S if q < k), 1) for S, c in g if c % 2 and any(q < k for q in S)]
        b, C0, nc = T.coset(layer, k)
        v = 0
        for S, _ in syn:
            v ^= 1 << (sum(1 << q for q in S) - 1)
        rec = dict(id=r.get("id"), k=k, ccz=r["ccz"], factory_n=r["best"]["n"], N=r["best"]["N"],
                   synthesis_weight=bin(v).count("1"), synthesis_ok_poly=bool(deposits(syn, k, layer)),
                   synthesis_ok_coset=in_coset(v, b, C0), compile_bound=r["tcount_hi"], floor=floor(r["ccz"], k),
                   ref_n=r.get("ref_n"), ref_tcount=r.get("ref_tcount"))
        rec["tau_upper"] = min(rec["synthesis_weight"], rec["compile_bound"]) if rec["synthesis_ok_poly"] and rec["synthesis_ok_coset"] else rec["compile_bound"]
        if a.isd and rec["tau_upper"] > rec["floor"]:
            t = time.time(); steps = []; word = None
            while rec["tau_upper"] > rec["floor"]:            # lower the target until a fixed-budget search fails
                w, seen, its = isd(b, C0, nc, rec["tau_upper"] - 1, a.isd, a.p, seed=len(steps))
                ok = bool(w is not None and in_coset(w, b, C0))
                steps.append(dict(target=rec["tau_upper"] - 1, iters=its, found=seen if ok else None))
                if not ok:
                    break
                rec["tau_upper"] = seen; word = w
            rec["isd"] = dict(p=a.p, budget=a.isd, steps=steps, seconds=round(time.time() - t, 1),
                              word=[y for y in range(1, nc + 1) if word >> (y - 1) & 1] if word else None)
            if word:
                rec["one_check"] = one_check(rec["isd"]["word"], k, layer)
        out.append(rec)
        print(rec["id"], "n", rec["factory_n"], "synthesis", rec["synthesis_weight"], rec["synthesis_ok_poly"], rec["synthesis_ok_coset"],
              "tau in", [rec["floor"], rec["tau_upper"]], rec.get("isd", ""), flush=True)
    dst = os.path.join(ROOT, "outputs", "k78_pipeline", f"tau_upper_{a.run}{'_isd' if a.isd else ''}{a.tag}.json")
    json.dump(dict(settings=vars(a), rows=out), open(dst, "w"), indent=1)
