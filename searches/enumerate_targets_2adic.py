"""Exact target classes at l=4 k=4 (|Q_4| = 2^29, too large for one orbit table) by the 2-adic filtration.

Write each output v in Q_k as v = v_lo + 2 v_hi digit by digit. The subgroup 2Q = {v : every digit even} is
invariant under the affine group AGL(k,2) (it is multiplication by 2 of a Z-module map), so AGL(k,2) acts on the
quotient Q/2Q (the low bits, 2^m elements, m = number of coordinates) and, for a fixed low part, on its fiber
{v_lo + 2w}. The low part of A v depends only on v_lo; the high part is A_hi(v) = (A v_lo div 2) + A w mod 2 per
coordinate of width 2, i.e. an affine map of w with offset from v_lo. So the same two-step method as
enumerate_targets_fiber.py applies: orbits of the low quotient with a BFS tree, then each root's stabilizer (random
Schreier generators, completed until the total matches the exact Burnside count) acting on its fiber.

    python searches/enumerate_targets_2adic.py --k 4          # writes outputs/targets/l4_k4.json
"""
import os, sys, json, time, random, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from target_classes import coords, action_matrix, orbits_burnside  # noqa: E402
from enumerate_targets import to_spec, poly_from_table  # noqa: E402
from enumerate_targets_fiber import gens_affine  # noqa: E402


class TwoAdic:
    def __init__(self, L, k, log=print, seed=1):
        self.L, self.k, self.log = L, k, log
        self.C = coords(L, k)
        self.E = np.array([e for _, _, e in self.C], dtype=np.int64)
        self.mods = 1 << self.E
        self.m = len(self.C)
        self.hi = [j for j in range(self.m) if self.E[j] == 2]          # coordinates with a high bit
        gens, self.inv = gens_affine(k)
        self.A = [np.array(action_matrix(g, L, k, self.C, b), dtype=np.int64) for g, b in gens]
        self.rng = random.Random(seed)

    # an element is (lo: m bits, hi: len(hi) bits); value per coordinate = lo_j + 2 hi_j
    def value(self, lo, hi):
        v = (lo >> np.arange(self.m)) & 1
        h = np.zeros(self.m, dtype=np.int64)
        for t, j in enumerate(self.hi):
            h[j] = (hi >> t) & 1
        return v + 2 * h

    def split(self, v):
        v = v % self.mods
        lo = int(sum(int(v[j] & 1) << j for j in range(self.m)))
        hi = int(sum(int((v[j] >> 1) & 1) << t for t, j in enumerate(self.hi)))
        return lo, hi

    def apply(self, M, lo, hi):
        return self.split(M @ self.value(lo, hi))

    def low_tables(self):
        t = time.time()
        n = 1 << self.m
        los = np.arange(n, dtype=np.int64)
        bits = ((los[None, :] >> np.arange(self.m)[:, None]) & 1)
        self.limg = [(((A @ bits) & 1) * (1 << np.arange(self.m))[:, None]).sum(axis=0).astype(np.int64) for A in self.A]
        lab = np.arange(n, dtype=np.int64)
        while True:
            new = lab.copy()
            for im in self.limg:
                np.minimum(new, lab[im], out=new)
            for _ in range(4):
                new = new[new]
            if np.array_equal(new, lab):
                break
            lab = new
        gates = bits.sum(axis=0)
        root = np.full(n, -1, dtype=np.int64); par = np.full(n, -1, dtype=np.int64); pg = np.full(n, -1, dtype=np.int8)
        self.roots = []
        for orb in np.unique(lab):
            mem = np.flatnonzero(lab == orb); r = int(mem[np.argmin(gates[mem])])
            self.roots.append(r); root[r] = r; fr = [r]
            while fr:
                nx = []
                for w in fr:
                    for s, im in enumerate(self.limg):
                        u = int(im[w])
                        if root[u] < 0:
                            root[u], par[u], pg[u] = r, w, s; nx.append(u)
                fr = nx
        self.lroot, self.lpar, self.lpg = root, par, pg
        self.log(f"low quotient 2^{self.m}: {len(self.roots)} orbits, {time.time() - t:.1f}s")

    def word_to_root(self, u):
        """List of generator indices g_1..g_t with g_t...g_1 (lift of u) in the fiber of root(u)."""
        w = []
        while self.lpg[u] >= 0:
            w.append(self.inv[int(self.lpg[u])]); u = int(self.lpar[u])
        return w

    def mat(self, word):
        M = np.eye(self.m, dtype=np.int64)
        for s in word:
            M = (self.A[s] @ M) % self.mods[:, None]
        return M

    def schreier(self, r, count):
        mem = np.flatnonzero(self.lroot == r); out = []
        for _ in range(count):
            u = int(mem[self.rng.randrange(len(mem))]); s = self.rng.randrange(len(self.A))
            su = int(self.limg[s][u])
            up = [int(self.lpg[x]) for x in self._path(u)]
            M = (self.mat(self.word_to_root(su)) @ self.A[s] @ self.mat(list(reversed(up)))) % self.mods[:, None]
            out.append(M)
        return out

    def _path(self, u):
        p = []
        while self.lpg[u] >= 0:
            p.append(u); u = int(self.lpar[u])
        return p

    def fiber_labels(self, r, Ms):
        nh = 1 << len(self.hi)
        imgs = []
        for M in Ms:
            lo0, _ = self.apply(M, r, 0); assert lo0 == r
            imgs.append(np.array([self.apply(M, r, h)[1] for h in range(nh)], dtype=np.int64))
        lab = np.arange(nh, dtype=np.int64)
        while True:
            new = lab.copy()
            for im in imgs:
                np.minimum(new, lab[im], out=new)
            for _ in range(4):
                new = new[new]
            if np.array_equal(new, lab):
                break
            lab = new
        return lab

    def run(self, ngen=24, rounds=8):
        self.low_tables()
        exact = orbits_burnside(self.L, self.k)
        self.hs = {r: self.schreier(r, ngen) for r in self.roots}
        for i in range(rounds):
            t = time.time()
            self.flab = {r: self.fiber_labels(r, self.hs[r]) for r in self.roots}
            tot = sum(len(np.unique(v)) for v in self.flab.values())
            self.log(f"round {i}: {tot} orbits (Burnside {exact}), {time.time() - t:.1f}s")
            if tot == exact:
                return
            assert tot > exact
            for r in self.roots:
                self.hs[r] += self.schreier(r, ngen)
        raise RuntimeError("stabilizers not generated")

    def digits_of_poly(self, P):
        mod = 1 << self.L
        tab = [sum(c for T, c in P.items() if all(x >> q & 1 for q in T)) % mod for x in range(1 << self.k)]
        a = poly_from_table(tab, self.k, mod)
        return np.array([(a.get(T, 0) >> sh) % (1 << e) for T, sh, e in self.C], dtype=np.int64)

    def label(self, dig):
        lo, hi = self.split(dig)
        r = int(self.lroot[lo])
        M = self.mat(self.word_to_root(lo))
        lo2, hi2 = self.apply(M, lo, hi)
        assert lo2 == r
        return f"{r}:{int(self.flab[r][hi2])}"


def enumerate_targets_2adic(L, k, log=print):
    t0 = time.time()
    F = TwoAdic(L, k, log); F.run()
    from enumerate_targets import enumerate_targets
    from coset_certify import parse_target
    sub, _ = enumerate_targets(L, k - 1, log=lambda *a: None)
    padded = {F.label(np.zeros(F.m, dtype=np.int64))}
    for c in sub["classes"]:
        padded.add(F.label(F.digits_of_poly(parse_target(c["spec"], L))))
    assert len(padded) == orbits_burnside(L, k - 1)
    out = []
    for r in F.roots:
        lab = F.flab[r]; seen = set()
        for h in np.argsort(lab, kind="stable"):
            key = f"{r}:{int(lab[h])}"
            if key in seen or key in padded:
                seen.add(key); continue
            seen.add(key)
            dig = F.value(r, int(h))
            out.append(dict(orbit=key, genuine=True, gates=int(dig.sum()), spec=to_spec([int(x) for x in dig], F.C, L)))
    log(f"l={L} k={k}: {len(out)} genuine classes, {time.time() - t0:.1f}s total")
    return dict(l=L, k=k, classes=out, seconds=round(time.time() - t0, 2)), F


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=4)
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--out-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs", "targets"))
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    res, _ = enumerate_targets_2adic(a.level, a.k)
    json.dump(res, open(os.path.join(a.out_dir, f"l{a.level}_k{a.k}.json"), "w"), indent=1)
