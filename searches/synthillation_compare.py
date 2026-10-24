"""Compare Table E1 with the synthillation construction (Campbell--Howard, arXiv:1606.01904, Theorem 1 and Table I).

For each target U the construction stacks an optimal gate-synthesis matrix A (tau[U] columns), two copies of a
Lempel factor B (mu[U] columns) and, for targets with T parts, four copies of the linear-term column c, onto the
constant rows of the Table I case fixed by the parities of tau and mu; it uses n = tau + 2 mu + Delta inputs. We
build that matrix and check it with the same borrowed-identity test as the rest of the repository: every monomial
touching a check row must be Clifford and the output monomials must equal the target up to a diagonal Clifford.
Run: python3 searches/synthillation_compare.py
"""
import itertools, os, sys

import z3
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "classification"))
from tcount import snf_z8, _solve_z8, parity_matrix
G = {"T": (1, {(0,): 1}), "CS": (2, {(0, 1): 2}), "CCZ": (3, {(0, 1, 2): 4})}
def tgt(names):
    t = {}; q = 0
    for nm in names:
        s, d = G[nm]
        for T, v in d.items(): t[frozenset(x + q for x in T)] = v
        q += s
    return q, t
def coeffs_ok(cols, n, k, target, upto):
    for t in range(1, upto + 1):
        mod = {1: 2, 2: 4, 3: 8}[t]
        for T in itertools.combinations(range(n), t):
            m = sum(1 << i for i in T)
            s = sum(1 for y in cols if y & m == m) * ((-2) ** (t - 1))
            want = target.get(frozenset(T), 0) if m < (1 << k) else 0
            if (s - want) % mod: return False
    return True
def opt_A(k, target):
    """A minimum-T-count parity network for the target (exact, coset leader)."""
    M = parity_matrix(k); U, W, D = snf_z8(M); nv = (1 << k) - 1
    f = []
    for x in range(1, nv + 1):
        xs = {i for i in range(k) if x >> i & 1}
        f.append(sum(c for T, c in target.items() if T <= xs) % 8)
    a0 = _solve_z8(U, W, D, f)
    basis = [sum(1 << r for r in range(nv) if W[r][i] & 1) for i in range(nv) if D[i] % 8 == 0]
    b = sum(1 << r for r in range(nv) if a0[r] & 1)
    best = None
    for c in itertools.product((0, 1), repeat=len(basis)):
        v = b
        for ci, bv in zip(c, basis):
            if ci: v ^= bv
        if best is None or bin(v).count("1") < bin(best).count("1"): best = v
    return [y for y in range(1, nv + 1) if best >> (y - 1) & 1]
def opt_B(k, target, mu):
    """mu columns with |B^T x| equal to the target up to cubic terms (Lempel condition Q = B B^T mod 2)."""
    s = z3.Solver(); b = [[z3.Bool(f"b{i}_{j}") for j in range(mu)] for i in range(k)]
    x = lambda c: z3.If(c, 1, 0)
    for i in range(k):
        s.add(z3.Sum([x(b[i][j]) for j in range(mu)]) % 2 == target.get(frozenset([i]), 0) % 2)
        for i2 in range(i + 1, k):
            q = (target.get(frozenset([i, i2]), 0) // 2) % 2
            s.add(z3.Sum([x(z3.And(b[i][j], b[i2][j])) for j in range(mu)]) % 2 == q)
    for j in range(mu): s.add(z3.Or([b[i][j] for i in range(k)]))
    if s.check() != z3.sat: return None
    m = s.model()
    return [sum(1 << i for i in range(k) if z3.is_true(m[b[i][j]])) for j in range(mu)]
T5 = {(0, 0): [[1,1,0],[0,1,1]], (0, 1): [[1,1,0,1,0],[0,1,1,1,1]], (1, 0): [[1,1,0,1],[1,0,1,1]], (1, 1): [[1,1,0,1,1,0],[1,0,1,1,0,1]]}
C14 = {1: [[1,1,0,1,0,0,1,1,0,0,1],[0,1,1,0,1,0,1,0,1,0,1],[0,0,0,1,1,1,1,1,1,1,1]],
       2: [[1,1,0,1,0,0,1,1,0,0,1,1],[0,1,1,0,1,0,1,0,1,0,1,1],[0,0,0,1,1,1,1,1,1,1,1,0]],
       3: [[1,1,0,1,0,0,1,1,0,0,1,1],[1,0,1,0,1,0,1,0,1,0,1,1],[0,0,0,1,1,1,1,1,1,1,1,0]],
       4: [[1,1,0,1,0,0,1,1,0,0,1,1,1,0],[1,0,1,0,1,0,1,0,1,0,1,1,0,1],[0,0,0,1,1,1,1,1,1,1,1,0,0,0]]}
def synth(names, mu):
    k, target = tgt(names); A = opt_A(k, target); tau = len(A)
    lin = any(len(T) == 1 for T in target)
    B = opt_B(k, target, mu)
    if B is None: return f"tau={tau}: no B with mu={mu}"
    if lin:
        tm = C14[{(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}[(tau % 2, mu % 2)]]
        c = sum(1 << i for i in range(k) if frozenset([i]) in target)
        blocks = [A, B, B] + [[c]] * 4
    else:
        tm = T5[(tau % 2, mu % 2)]; blocks = [A, B, B]
    blocks += [[0]] * (len(tm[0]) - len(blocks)); R = len(tm)
    cols = [v | sum(1 << (k + r) for r in range(R) if tm[r][bi]) for bi, blk in enumerate(blocks) for v in blk]
    ok = coeffs_ok(cols, k + R, k, target, 3) and all(y >> k for y in cols)
    return f"tau={tau} mu={mu} N={len(cols)} on {k + R} qubits: {'VALID' if ok else 'NOT a valid factory'}"
TARGETS = [(["CS"], 3), (["T"], 1), (["T", "T"], 2), (["CS", "CS"], 5), (["CCZ", "CS"], 3), (["CS", "T"], 3),
           (["CCZ", "T"], 1), (["T", "T", "T"], 3), (["CCZ", "T", "T"], 2), (["CS", "T", "T"], 4),
           (["CS", "CS", "CS"], 7), (["CS", "CS", "T"], 5), (["CCZ", "CS", "T"], 3), (["T"] * 4, 4),
           (["CCZ", "T", "T", "T"], 3), (["CS", "CS", "T", "T"], 6), (["CS", "T", "T", "T"], 5), (["T"] * 5, 5),
           (["CS"] + ["T"] * 4, 6), (["T"] * 6, 6)]

if __name__ == "__main__":
    # Pure-CCZ targets (mu = 0) use cases 9-11 with B empty: 1 CCZ -> 8, 2 CCZ -> 14, as in the synthillation paper.
    for names, mu in TARGETS:
        print(" ".join(names), "->", synth(names, mu), flush=True)
