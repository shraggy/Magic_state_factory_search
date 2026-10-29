"""Exact number of SAT targets per (l, k): Clifford classes of k-qubit diagonal level-l outputs.

A target is a diagonal level-l state modulo the free corrections (cliff_exp), CNOT frames and Pauli X (the affine
group AGL(k,2), the classifier's 'GL x shift'). Modulo the corrections the outputs form a finite abelian group Q_k
on which AGL(k,2) acts linearly, so Burnside's lemma counts the orbits exactly. Targets that
genuinely need all k outputs number orb(k) - orb(k-1). Small cases are checked against brute-force orbit enumeration.

    python searches/target_classes.py --levels 2 3 4 --kmax 6
"""
import os, sys, json, time, argparse
from itertools import combinations, product
from fractions import Fraction
from math import prod

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_search import cliff_exp, poly_from_table  # noqa: E402


# ----------------------------------------------------------------- the group Q_k and the GL(k,2) action
def coords(L, k):
    """Coordinates of Q_k: supports T with Newton coefficient a_T in 2^{|T|-1} Z, taken modulo 2^cliff_exp."""
    out = []
    for s in range(1, k + 1):
        e = cliff_exp(s, L) - (s - 1)
        if e > 0:
            out += [(frozenset(T), s - 1, e) for T in combinations(range(k), s)]
    return out


def action_matrix(g, L, k, C, b=0):
    """Integer matrix of x -> g x + b on Q_k (CNOT frame g, then Pauli X on the qubits of b; the constant term is a
    global phase and is dropped): column j = coordinates of the j-th basis output after the map."""
    mod = 1 << L
    idx = {T: j for j, (T, _, _) in enumerate(C)}
    xs = np.arange(1 << k)
    gx = np.zeros_like(xs)                                  # bits of g x for every x
    for i in range(k):
        row = sum(int(g[i][j]) << j for j in range(k))
        gx |= (np.bitwise_count(xs & row) & 1).astype(np.int64) << i if hasattr(np, "bitwise_count") else \
              (np.array([bin(int(v) & row).count("1") & 1 for v in xs]) << i)
    gx ^= b
    A = np.zeros((len(C), len(C)), dtype=object)
    for j, (T, sh, _) in enumerate(C):
        m = sum(1 << q for q in T)
        tab = (((gx & m) == m).astype(np.int64) << sh) % mod
        for S, a in poly_from_table(tab.tolist(), k, mod).items():
            if S in idx:
                _, sh2, e2 = C[idx[S]]
                assert a % (1 << sh2) == 0
                A[idx[S], j] = (a >> sh2) % (1 << e2)
    return A


def fixed_points(A, C):
    """|{v in Q_k : (A - 1) v = 0}| = |Q_k / (A - 1) Q_k| = 2-adic index of the lattice spanned by the columns of
    A - 1 and diag(2^e). The lattice contains 2^E Z^m (E = max e), so Smith elimination over Z/2^E is exact."""
    m = len(C)
    E = max(e for _, _, e in C)
    mod = 1 << E
    cols = [[(int(A[i, j]) - (i == j)) % mod for i in range(m)] for j in range(m)]
    cols += [[(1 << C[i][2]) % mod if r == i else 0 for r in range(m)] for i in range(m)]
    M = [list(r) for r in zip(*cols)]                       # m x 2m over Z/2^E
    val = lambda x: E if x == 0 else (x & -x).bit_length() - 1
    logidx, rows, ncol = 0, list(range(m)), len(M[0])
    for _ in range(m):
        best = None
        for i in rows:
            for j in range(ncol):
                if M[i][j] and (best is None or val(M[i][j]) < best[0]):
                    best = (val(M[i][j]), i, j)
        if best is None:
            logidx += E * len(rows)
            break
        v, i, j = best
        logidx += v
        u = pow(M[i][j] >> v, -1, mod)                      # unit part of the pivot
        for r in rows:
            if r != i and M[r][j]:
                f = ((M[r][j] >> v) * u) % mod
                M[r] = [(a - f * b) % mod for a, b in zip(M[r], M[i])]
        for c in range(ncol):
            if c != j and M[i][c]:
                f = ((M[i][c] >> v) * u) % mod
                for r in range(m):
                    M[r][c] = (M[r][c] - f * M[r][j]) % mod
        rows.remove(i)
    return 1 << logidx


# ----------------------------------------------------------------- conjugacy classes of GL(k,2)
def irreducibles(dmax):
    """Irreducible polynomials over F_2 other than x, as bit masks (bit i = coefficient of x^i), by degree."""
    def mulmod2(a, b):
        r = 0
        while b:
            if b & 1:
                r ^= a
            a <<= 1; b >>= 1
        return r
    red = set()
    for da in range(1, dmax):
        for a in range(1 << da, 1 << (da + 1)):
            for db in range(da, dmax - da + 1):
                for b in range(1 << db, 1 << (db + 1)):
                    red.add(mulmod2(a, b))
    return {d: [p for p in range(1 << d, 1 << (d + 1)) if p & 1 and p not in red] for d in range(1, dmax + 1)}


def partitions(n, maxpart=None):
    maxpart = n if maxpart is None else maxpart
    if n == 0:
        yield (); return
    for p in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest


def poly_pow(p, e):
    r = 1
    for _ in range(e):
        a, b, out = r, p, 0
        while b:
            if b & 1:
                out ^= a
            a <<= 1; b >>= 1
        r = out
    return r


def companion(p):
    d = p.bit_length() - 1
    M = np.zeros((d, d), dtype=np.int64)
    for i in range(1, d):
        M[i, i - 1] = 1
    for i in range(d):
        M[i, d - 1] = (p >> i) & 1
    return M


def centralizer_order(blocks):
    """Order of the centralizer in GL(k,2) for elementary divisors {p: partition} (Macdonald's formula)."""
    out = Fraction(1)
    for p, lam in blocks:
        q = 2 ** (p.bit_length() - 1)
        conj = [sum(1 for x in lam if x > i) for i in range(lam[0])]
        out *= Fraction(q) ** sum(c * c for c in conj)
        for part in set(lam):
            mult = lam.count(part)
            for j in range(1, mult + 1):
                out *= 1 - Fraction(1, q ** j)
    assert out.denominator == 1
    return int(out)


def gl_classes(k):
    """(representative, class size) for every conjugacy class of GL(k,2), from elementary divisors."""
    order = prod((1 << k) - (1 << i) for i in range(k))
    irr = [p for d, ps in irreducibles(k).items() for p in ps]
    out = []

    def rec(i, left, chosen):
        if left == 0:
            blocks = []
            for p, lam in chosen:
                for part in lam:
                    blocks.append(companion(poly_pow(p, part)))
            g = np.zeros((k, k), dtype=np.int64); o = 0
            for B in blocks:
                d = B.shape[0]; g[o:o + d, o:o + d] = B; o += d
            out.append((g, order // centralizer_order(chosen)))
            return
        if i == len(irr):
            return
        p = irr[i]; d = p.bit_length() - 1
        rec(i + 1, left, chosen)
        for n in range(1, left // d + 1):
            for lam in partitions(n):
                rec(i + 1, left - n * d, chosen + [(p, lam)])
    rec(0, k, [])
    assert sum(s for _, s in out) == order, "class sizes do not sum to |GL(k,2)|"
    return out, order


# ----------------------------------------------------------------- counting
def orbits_burnside(L, k):
    """Orbits of the affine group AGL(k,2) (CNOT frames and Pauli X). sum over (g, b) of Fix(g, b) = sum over GL
    classes of |class| * sum_b Fix(g_rep, b), since conjugation by GL permutes the translations."""
    if k == 0:
        return 1
    C = coords(L, k)
    classes, order = gl_classes(k)
    total = sum(size * sum(fixed_points(action_matrix(g, L, k, C, b), C) for b in range(1 << k)) for g, size in classes)
    assert total % (order << k) == 0
    return total // (order << k)


def orbits_brute(L, k):
    """Orbit count by union-find over every element of Q_k under two generators of GL(k,2) (small k only)."""
    if k == 0:
        return 1
    C = coords(L, k)
    mods = [1 << e for _, _, e in C]
    size = prod(mods)
    if size > 1 << 22:
        return None
    gens = [np.eye(k, dtype=np.int64)]
    if k > 1:
        cyc = np.roll(np.eye(k, dtype=np.int64), 1, axis=0)
        tv = np.eye(k, dtype=np.int64); tv[0, 1] = 1
        gens = [cyc, tv]
    mats = [action_matrix(g, L, k, C) for g in gens] + [action_matrix(np.eye(k, dtype=np.int64), L, k, C, 1)]
    parent = list(range(size))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for v in product(*[range(m) for m in mods]):
        iv = 0
        for c, m in zip(v, mods):
            iv = iv * m + c
        for A in mats:
            w = [sum(int(A[i, j]) * v[j] for j in range(len(v))) % mods[i] for i in range(len(v))]
            iw = 0
            for c, m in zip(w, mods):
                iw = iw * m + c
            a, b = find(iv), find(iw)
            if a != b:
                parent[a] = b
    return sum(1 for x in range(size) if find(x) == x)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", type=int, nargs="+", default=[2, 3, 4])
    ap.add_argument("--kmax", type=int, default=6)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "outputs", "target_classes.json"))
    a = ap.parse_args()
    res = {}
    for L in a.levels:
        prev = 1
        for k in range(1, a.kmax + 1):
            t = time.time()
            orb = orbits_burnside(L, k)
            tb = time.time() - t
            brute = orbits_brute(L, k) if k <= 3 else None
            row = dict(l=L, k=k, group_size_log2=sum(e for _, _, e in coords(L, k)), orbits=orb,
                       genuine=orb - prev, brute_force=brute, seconds=round(tb, 2))
            prev = orb
            res[f"l{L}_k{k}"] = row
            print(json.dumps(row), flush=True)
            json.dump(res, open(a.out, "w"), indent=1)
