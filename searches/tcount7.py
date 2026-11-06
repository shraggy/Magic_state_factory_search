"""Exact T-count of a diagonal level-3 gate by minimum-weight decoding of its Reed-Muller coset. The code handles up to
8 qubits (256 coordinates), but it has been run to completion only at k = 7; no k = 8 exact T-count was obtained.

The T-count of U (up to Clifford) is the minimum weight of a0 + RM(k-4, k)* over the 2^k - 1 nonzero parities, where
a0 is any parity solution of M a = f (mod 8) and the coset space is the punctured Reed-Muller code of the decoder in
classification/tcount.py (Amy-Mosca). That module materialises the whole coset (2^dim words), which stops at k = 6.
Here the minimum is found by a Brouwer-Zimmermann search over disjoint information sets with the compiled inner loop
of csrc/capped_minweight.c (vectors up to 256 coordinates): decide(w) asks whether a coset word of weight < w exists,
and a binary search gives the exact minimum.

    python searches/tcount7.py            # self-test against the exact k <= 6 decoder, then times k = 7 layers
"""
import os, sys, random, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
from coset_certify import gf2_rref  # noqa: E402
import tcount as TC  # noqa: E402
BIN = os.path.join(HERE, "csrc", "capped_minweight")


def coset(P, k):
    """(base, basis, ncoord) of the T-count coset of polynomial P (frozenset -> coeff mod 8) on qubits 0..k-1."""
    Pr = {frozenset(q + 1 for q in T): c for T, c in P.items()}       # tcount.py uses labels 1..k
    M = TC.parity_matrix(k); U, W, D = TC.snf_z8(M); fv = TC._fvec(Pr, k); a0 = TC._solve_z8(U, W, D, fv)
    if TC._matvec_mod(M, a0) != [v % 8 for v in fv]:
        raise RuntimeError("not a level-3 diagonal gate")
    n = (1 << k) - 1
    basis = [sum(1 << r for r in range(n) if W[r][i] & 1) for i in range(n) if D[i] % 8 == 0]
    base = sum(1 << r for r in range(n) if a0[r] & 1)
    return base, basis, n


def _info_sets(C0, nc, seed=1, tries=600):
    rnd = random.Random(seed); d = len(C0); best = []
    for _ in range(tries):
        free = list(range(nc)); rnd.shuffle(free); sets = []; used = set()
        while True:
            cols = [c for c in free if c not in used] + [c for c in free if c in used]
            perm = {c: i for i, c in enumerate(cols)}
            G = [sum(1 << perm[c] for c in range(nc) if v >> c & 1) for v in C0]
            R, piv = gf2_rref(G, nc)
            if len(R) < d:
                break
            I = {cols[p] for p in piv}
            if I & used:
                break
            sets.append(I); used |= I
        if len(sets) > len(best):
            best = sets
        if len(best) * d + d > nc:
            break
    return best


def exists_below(b, C0, nc, cap, thr):
    """Is there a word of weight < cap in b + span(C0)? Complete Brouwer-Zimmermann decision."""
    d = len(C0)
    if d == 0:
        return bin(b).count("1") < cap, bin(b).count("1")
    sets = _info_sets(C0, nc); m = len(sets); w = (cap - 1) // m
    h = lambda v: " ".join(f"{(v >> (64 * i)) & ((1 << 64) - 1):x}" for i in range(4))
    lb = 1 << 30
    for I in sets:
        order = sorted(I) + [c for c in range(nc) if c not in I]
        perm = {c: i for i, c in enumerate(order)}
        rm = lambda v: sum(1 << perm[c] for c in range(nc) if v >> c & 1)
        G, piv = gf2_rref([rm(v) for v in C0], nc)
        bb = rm(b)
        for g, p in zip(G, piv):
            if bb >> p & 1:
                bb ^= g
        inp = f"{d}\n{h(bb)}\n" + "\n".join(h(g) for g in G) + "\n"
        out = subprocess.run([BIN, str(thr), str(cap), str(w)], input=inp, capture_output=True, text=True,
                             check=True).stdout.split()
        lb = min(lb, int(out[1]))
        if out[0] == "1":
            return True, lb
    return False, lb


def tcount(P, k, thr=4):
    """Exact T-count: the smallest w with a coset word of weight w (binary search on exists_below)."""
    b, C0, nc = coset(P, k)
    ok, ub = exists_below(b, C0, nc, nc + 2, thr)             # any word: an upper bound on the minimum
    lo, hi = 0, ub                                             # minimum in [lo, hi]
    while lo < hi:
        mid = (lo + hi) // 2
        found, seen = exists_below(b, C0, nc, mid + 1, thr)    # weight <= mid exists?
        if found:
            hi = min(mid, seen)
        else:
            lo = mid + 1
    return lo


if __name__ == "__main__":
    import itertools
    if not os.path.exists(BIN):
        subprocess.run(["cc", "-O3", "-o", BIN, os.path.join(HERE, "csrc", "capped_minweight.c"), "-lpthread"], check=True)
    rng = random.Random(3); bad = 0; t0 = time.time()
    for _ in range(40):                                        # agree with the exact k <= 6 decoder
        k = rng.randint(3, 6)
        tri = [T for T in itertools.combinations(range(k), 3) if rng.random() < 0.35]
        if not tri:
            continue
        P = {frozenset(T): 4 for T in tri}
        a = tcount(P, k); e = TC.tcount_rm({frozenset(q + 1 for q in T): 4 for T in tri}, k)
        bad += a != e
        if a != e:
            print("MISMATCH", k, tri, a, e, flush=True)
    print(f"self-test vs exact decoder: {40 - bad}/40 agree, {time.time() - t0:.0f}s", flush=True)
    for tri in ([(0, 1, 2), (3, 4, 5), (0, 3, 6)], [(0, 1, 2), (2, 3, 4), (4, 5, 6), (0, 3, 6)]):
        t = time.time(); print("k=7", tri, "T-count", tcount({frozenset(T): 4 for T in tri}, 7),
                               f"{time.time() - t:.1f}s", flush=True)


def exists_below2(b, C0, nc, cap, thr):
    """Is there a word of weight < cap in b + span(C0)? Brouwer-Zimmermann with two disjoint information sets, the
    second allowed to be rank-deficient (rank r2 < d): a word not found by enumerating <= w generators in basis j has
    weight >= w + 1 - (d - r_j) on set j, so the two enumerations together rule out every word of weight
    < (w + 1) + (w + 1 - (d - r2)). Exact when that bound reaches cap. Returns (found, min weight seen, bound)."""
    d = len(C0); cols = list(range(nc)); I1 = None
    import random as _r
    rnd = _r.Random(1)
    for _ in range(200):                                   # first set: a full-rank information set
        rnd.shuffle(cols)
        perm = {c: i for i, c in enumerate(cols)}
        G, piv = gf2_rref([sum(1 << perm[c] for c in range(nc) if v >> c & 1) for v in C0], nc)
        I = {cols[p] for p in piv}
        rest = [c for c in range(nc) if c not in I]
        perm2 = {c: i for i, c in enumerate(rest + sorted(I))}
        G2, piv2 = gf2_rref([sum(1 << perm2[c] for c in range(nc) if v >> c & 1) for v in C0], nc)
        r2 = sum(1 for p in piv2 if p < len(rest))
        if I1 is None or r2 > best_r2:
            I1, best_r2 = I, r2
        if r2 == len(rest) or r2 == d:
            break
    r2 = best_r2; defect = d - r2
    w = 0
    while (w + 1) + max(0, w + 1 - defect) < cap:
        w += 1
    h = lambda v: " ".join(f"{(v >> (64 * i)) & ((1 << 64) - 1):x}" for i in range(4))
    lb = 1 << 30
    rest = [c for c in range(nc) if c not in I1]
    for first in (sorted(I1) + rest, rest + sorted(I1)):
        perm = {c: i for i, c in enumerate(first)}
        rm = lambda v: sum(1 << perm[c] for c in range(nc) if v >> c & 1)
        G, piv = gf2_rref([rm(v) for v in C0], nc)
        bb = rm(b)
        for g, p in zip(G, piv):
            if bb >> p & 1:
                bb ^= g
        inp = f"{d}\n{h(bb)}\n" + "\n".join(h(g) for g in G) + "\n"
        out = subprocess.run([BIN, str(thr), str(cap), str(w)], input=inp, capture_output=True, text=True,
                             check=True).stdout.split()
        lb = min(lb, int(out[1]))
        if out[0] == "1":
            return True, lb, None
    return False, lb, (w + 1) + max(0, w + 1 - defect)
