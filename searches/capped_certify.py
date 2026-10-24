"""Exact proof that no distance-2 factory with N < N* exists for a target on n circuit qubits (l=3).

Brouwer-Zimmermann decision on the affine space of coset_certify.affine_space: a vector of weight < N* has at
most floor((N*-1)/m) nonzeros on one of m disjoint information sets, so enumerating every generator subset of
that size, in systematic form on each set, is complete. The inner enumeration is compiled C
(csrc/capped_minweight.c, vectors up to 256 coordinates, multithreaded). This produced the n=7 minimality
results of App. S10.

Build once:  cc -O3 -o searches/csrc/capped_minweight searches/csrc/capped_minweight.c -lpthread
Run:         python3 searches/capped_certify.py
"""
import os, sys, subprocess, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
BIN = os.path.join(HERE, "csrc", "capped_minweight")
from coset_certify import affine_space, parse_target, gf2_rref
def info_sets(C0, nc, tries=4000, seed=1):
    rnd = random.Random(seed); d = len(C0); best = []
    for _ in range(tries):
        free = list(range(nc)); rnd.shuffle(free); sets = []
        used = set()
        while True:
            cols = [c for c in free if c not in used] + [c for c in free if c in used]
            perm = {c: i for i, c in enumerate(cols)}
            G = [sum(1 << perm[c] for c in range(nc) if v >> c & 1) for v in C0]
            R, piv = gf2_rref(G, nc)
            if len(R) < d: break
            I = {cols[p] for p in piv}
            if I & used: break
            sets.append(I); used |= I
        if len(sets) > len(best): best = sets
        if len(best) * d + d > nc: break
    return best
def decide(spec, k, n, Nstar, thr=16):
    b, C0, nc = affine_space(parse_target(spec), k, n); d = len(C0)
    sets = info_sets(C0, nc); m = len(sets); w = (Nstar - 1) // m
    h = lambda v: " ".join(f"{(v >> (64*i)) & ((1<<64)-1):x}" for i in range(4))
    assert nc <= 256
    lb = 1 << 30
    for I in sets:
        order = sorted(I) + [c for c in range(nc) if c not in I]
        perm = {c: i for i, c in enumerate(order)}
        rm = lambda v: sum(1 << perm[c] for c in range(nc) if v >> c & 1)
        G, piv = gf2_rref([rm(v) for v in C0], nc)
        assert sorted(piv) == list(range(d))
        bb = rm(b)
        for g, p in zip(G, piv):
            if bb >> p & 1: bb ^= g
        inp = f"{d}\n{h(bb)}\n" + "\n".join(h(g) for g in G) + "\n"
        out = subprocess.run([BIN, str(thr), str(Nstar), str(w)], input=inp, capture_output=True, text=True, check=True).stdout.split()
        lb = min(lb, int(out[1]))
        if out[0] == "1": return True, lb, d, m, w
    return False, lb, d, m, w
if __name__ == "__main__":
    # Sanity checks: a known factory exists at its own N (cap N+1 must find it) and not below it.
    for name, spec, k, n, N in [("CCZ CCZ", "CCZ012,CCZ345", 6, 7, 14), ("CCZ CCZ (cap 15)", "CCZ012,CCZ345", 6, 7, 15),
                                ("T T (cap 15)", "T0,T1", 2, 5, 15)]:
        t = time.time(); f, lb, d, m, w = decide(spec, k, n, N)
        print(f"check {name:18s} n={n} N*={N}: {'smaller exists' if f else 'none below'} (dim {d}, {m} sets, w<={w}, best seen {lb}) {time.time()-t:.1f}s", flush=True)
    # Table E1 targets on n = 7 circuit qubits (target qubits 0-based here).
    for name, spec, k, N in [("CCZ CCZ", "CCZ012,CCZ345", 6, 14), ("CCZ CS", "CCZ012,CS34", 5, 18), ("CCZ T", "CCZ012,T3", 4, 20),
                             ("CS T T", "CS01,T2,T3", 4, 22), ("T T T T", "T0,T1,T2,T3", 4, 20)]:
        t = time.time(); f, lb, d, m, w = decide(spec, k, 7, N)
        print(f"{name:8s} n=7 N*={N}: {'SMALLER FACTORY EXISTS' if f else 'proved: no factory with N<'+str(N)+' on 7 qubits'} (dim {d}, {m} sets, w<={w}, best seen {lb}) {time.time()-t:.0f}s", flush=True)
