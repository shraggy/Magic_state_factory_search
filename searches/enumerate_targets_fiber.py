"""Exact target classes for cells whose output group Q_k is too large for one orbit table (l=3 k=6: 2^41).

The terms of degree >= s0 form a subgroup H of Q_k that the affine group AGL(k,2) maps to itself (lower-degree terms
can feed higher ones through the parity expansion x XOR y = x + y - 2xy, never the reverse), so AGL(k,2) also acts on
the quotient of the degree < s0 terms. (1) Label the orbits of that small 'top' quotient, with a BFS tree from one
root per orbit; (2) for each root c, generate its stabilizer from Schreier generators and label the orbits of the
'fiber' {c + terms of degree >= s0} under it. Orbits of Q_k are exactly the pairs (top orbit, fiber orbit). Random Schreier
generators can only miss stabilizer elements, which would split orbits, so the total is checked against the exact
Burnside count of target_classes.py and more generators are added until the two agree.

    python searches/enumerate_targets_fiber.py --level 3 --k 6     # writes outputs/targets/l3_k6.json
"""
import os, sys, json, time, random, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from target_classes import coords, action_matrix, orbits_burnside  # noqa: E402
from enumerate_targets import to_spec, poly_from_table  # noqa: E402


def gens_affine(k):
    """Generators of AGL(k,2) as (matrix, shift), and the index of each one's inverse."""
    I = np.eye(k, dtype=np.int64)
    cyc = np.roll(I, 1, axis=0)
    tv = I.copy(); tv[0, 1] = 1
    return [(cyc, 0), (cyc.T.copy(), 0), (tv, 0), (I, 1)], [1, 0, 2, 3]


class Fiber:
    def __init__(self, L, k, s0, log=print, seed=1):
        self.L, self.k, self.log = L, k, log
        C = coords(L, k)
        self.C = C
        self.mods = np.array([1 << e for _, _, e in C], dtype=np.int64)
        self.top = [j for j, (T, _, _) in enumerate(C) if len(T) < s0]      # quotient coordinates
        self.low = [j for j, (T, _, _) in enumerate(C) if len(T) >= s0]     # fiber coordinates (subgroup H)
        gens, self.inv = gens_affine(k)
        self.A = [np.array(action_matrix(g, L, k, C, b), dtype=np.int64) for g, b in gens]
        for A in self.A:                                  # H is invariant: its terms never feed lower degrees
            assert not A[np.ix_(self.top, self.low)].any()
        self.rng = random.Random(seed)
        self.tbits = [C[j][2] for j in self.top]
        self.lbits = [C[j][2] for j in self.low]
        self.ntop = 1 << sum(self.tbits)
        self.nlow = 1 << sum(self.lbits)

    # ---- packing helpers
    @staticmethod
    def _offsets(bits):
        off, o = [], 0
        for e in bits:
            off.append(o); o += e
        return off

    def unpack(self, idx, bits):
        off = self._offsets(bits)
        return [(idx >> o) & ((1 << e) - 1) for o, e in zip(off, bits)]

    @staticmethod
    def _pack_all(digs, bits):
        out = np.zeros(digs.shape[1], dtype=np.int64)
        o = 0
        for i, e in enumerate(bits):
            out |= (digs[i] & ((1 << e) - 1)) << o; o += e
        return out

    def _digits_all(self, n, bits):
        xs = np.arange(n, dtype=np.int64)
        off = self._offsets(bits)
        return np.array([(xs >> o) & ((1 << e) - 1) for o, e in zip(off, bits)])

    def mul(self, X, Y):
        return (X @ Y) % self.mods[:, None]

    # ---- top quotient: orbit labels and BFS trees
    def top_tables(self):
        t = time.time()
        D = self._digits_all(self.ntop, self.tbits)
        self.timg = []
        for A in self.A:
            M = A[np.ix_(self.top, self.top)]
            img = (M @ D) % self.mods[self.top][:, None]
            self.timg.append(self._pack_all(img, self.tbits).astype(np.int32))
        gates = D.sum(axis=0)
        # BFS from the fewest-gate element of each orbit (orbits found as connected components)
        par_gen = np.full(self.ntop, -1, dtype=np.int8)
        par = np.full(self.ntop, -1, dtype=np.int64)
        root = np.full(self.ntop, -1, dtype=np.int64)
        lab = np.arange(self.ntop, dtype=np.int64)
        while True:                                       # orbit labels by min-label propagation
            new = lab.copy()
            for im in self.timg:
                np.minimum(new, lab[im], out=new)
            for _ in range(4):
                new = new[new]
            if np.array_equal(new, lab):
                break
            lab = new
        self.roots = []
        for orb in np.unique(lab):
            members = np.flatnonzero(lab == orb)
            r = int(members[np.argmin(gates[members])])
            self.roots.append(r)
            root[r] = r
            frontier = [r]
            while frontier:
                nxt = []
                for w in frontier:
                    for s, im in enumerate(self.timg):
                        u = int(im[w])
                        if root[u] < 0:
                            root[u], par[u], par_gen[u] = r, w, s
                            nxt.append(u)
                frontier = nxt
        assert (root >= 0).all()
        self.troot, self.tpar, self.tpar_gen = root, par, par_gen
        self.log(f"top quotient 2^{sum(self.tbits)}: {len(self.roots)} orbits, {time.time() - t:.1f}s")

    def to_root(self, u):
        """Q-matrix B with B (u-lift) in the fiber of root(u): product of inverse generators up the BFS tree."""
        B = np.eye(len(self.C), dtype=np.int64)
        while self.tpar_gen[u] >= 0:
            B = self.mul(self.A[self.inv[self.tpar_gen[u]]], B)
            u = int(self.tpar[u])
        return B

    def from_root(self, u):
        F = np.eye(len(self.C), dtype=np.int64)
        path = []
        while self.tpar_gen[u] >= 0:
            path.append(int(self.tpar_gen[u])); u = int(self.tpar[u])
        for s in reversed(path):
            F = self.mul(self.A[s], F)
        return F

    def schreier(self, r, count):
        """`count` random Schreier generators of the stabilizer of root r in the top quotient."""
        members = np.flatnonzero(self.troot == r)
        out = []
        for _ in range(count):
            u = int(members[self.rng.randrange(len(members))])
            s = self.rng.randrange(len(self.A))
            su = int(self.timg[s][u])
            h = self.mul(self.to_root(su), self.mul(self.A[s], self.from_root(u)))
            out.append(h)
        return out

    # ---- fibers
    def fiber_labels(self, r, hs):
        """Orbit labels of {root r + lower terms} under the affine maps w -> off + H w induced by hs."""
        tdig = np.array(self.unpack(r, self.tbits), dtype=np.int64)
        D = self._digits_all(self.nlow, self.lbits)
        lm = self.mods[self.low][:, None]
        imgs = []
        for h in hs:
            assert np.array_equal((h[np.ix_(self.top, self.top)] @ tdig) % self.mods[self.top], tdig)
            off = h[np.ix_(self.low, self.top)] @ tdig
            img = (h[np.ix_(self.low, self.low)] @ D + off[:, None]) % lm
            imgs.append(self._pack_all(img, self.lbits).astype(np.int32))
        lab = np.arange(self.nlow, dtype=np.int32)
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

    def run(self, ngen=24, max_rounds=6):
        self.top_tables()
        exact = orbits_burnside(self.L, self.k)
        self.hs = {r: self.schreier(r, ngen) for r in self.roots}
        for rnd in range(max_rounds):
            t = time.time()
            self.flab = {r: self.fiber_labels(r, self.hs[r]) for r in self.roots}
            total = sum(len(np.unique(v)) for v in self.flab.values())
            self.log(f"round {rnd}: {total} orbits (Burnside {exact}), {time.time() - t:.1f}s")
            if total == exact:
                return
            assert total > exact
            for r in self.roots:
                self.hs[r] += self.schreier(r, ngen)
        raise RuntimeError("stabilizers not generated")

    # ---- canonical label of any output
    def digits_of_poly(self, P):
        mod = 1 << self.L
        tab = [sum(c for T, c in P.items() if all(x >> q & 1 for q in T)) % mod for x in range(1 << self.k)]
        a = poly_from_table(tab, self.k, mod)
        dig = []
        for T, sh, e in self.C:
            v = a.get(T, 0)
            assert v % (1 << sh) == 0
            dig.append((v >> sh) % (1 << e))
        return np.array(dig, dtype=np.int64)

    def label(self, dig):
        """(root, fiber orbit) of an output given by its coordinate digits."""
        tdig = dig[self.top]
        u = int(sum(int(d) << o for d, o in zip(tdig, self._offsets(self.tbits))))
        v = (self.to_root(u) @ dig) % self.mods
        r = int(self.troot[u])
        assert int(sum(int(d) << o for d, o in zip(v[self.top], self._offsets(self.tbits)))) == r
        w = int(sum(int(d) << o for d, o in zip(v[self.low], self._offsets(self.lbits))))
        return r, int(self.flab[r][w])


def enumerate_targets_fiber(L, k, s0=3, log=print):
    t0 = time.time()
    F = Fiber(L, k, s0, log)
    F.run()
    # non-genuine orbits: those containing an output that ignores the last qubit, i.e. the k-1 classes embedded.
    # Enumerate them as the orbits meeting that subgroup: label every element of the subgroup is too many, so label
    # the embedded k-1 representatives, whose count is the exact k-1 orbit count.
    from enumerate_targets import enumerate_targets
    sub, _ = enumerate_targets(L, k - 1, log=lambda *a: None)
    from coset_certify import parse_target
    padded = {F.label(np.zeros(len(F.C), dtype=np.int64))}
    for c in sub["classes"]:
        padded.add(F.label(F.digits_of_poly(parse_target(c["spec"], L))))
    assert len(padded) == orbits_burnside(L, k - 1), "embedded k-1 classes do not map to distinct orbits"
    out = []
    for r in F.roots:
        lab = F.flab[r]
        tdig = F.unpack(r, F.tbits)
        D = F._digits_all(F.nlow, F.lbits)
        gates = D.sum(axis=0) + sum(tdig)
        order = np.lexsort((np.arange(F.nlow), gates, lab))
        first = order[np.r_[True, lab[order][1:] != lab[order][:-1]]]
        for w in first:
            key = (r, int(lab[w]))
            if key in padded:
                continue
            dig = [0] * len(F.C)
            for j, d in zip(F.top, tdig):
                dig[j] = int(d)
            for j, d in zip(F.low, D[:, w]):
                dig[j] = int(d)
            out.append(dict(orbit=f"{r}:{int(lab[w])}", genuine=True, gates=int(gates[w]), spec=to_spec(dig, F.C, L)))
    log(f"l={L} k={k}: {len(out)} genuine classes, {time.time() - t0:.1f}s total")
    return dict(l=L, k=k, classes=out, seconds=round(time.time() - t0, 2)), F


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--s0", type=int, default=3, help="terms of degree >= s0 form the fiber subgroup")
    ap.add_argument("--out-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs", "targets"))
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    res, _ = enumerate_targets_fiber(a.level, a.k, a.s0)
    json.dump(res, open(os.path.join(a.out_dir, f"l{a.level}_k{a.k}.json"), "w"), indent=1)
