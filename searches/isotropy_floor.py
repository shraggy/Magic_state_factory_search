"""Lower bound on the T-count of a CCZ layer (a pure-cubic diagonal level-3 gate): Mazumder's isotropy floor,
arXiv:2610.01024, Theorem 22 / Corollary 23. For a layer on k qubits with cubic form T_abc = 1 on its CCZ triples,
the radical R = {w : sum_a w_a T_abc = 0 for all b, c} has dimension d; then T-count >= 2(k - d) + 1 when the gate is
non-Clifford (d < k), over Hadamard-free {CNOT, T} circuits with clean ancillas. Computed by Gaussian elimination
over GF(2). Validated here against the exact Reed-Muller decoder for k <= 6.
"""
import itertools


def radical_dim(tri, k):
    """dim of {w in GF(2)^k : for all b<c, sum_a w_a [a,b,c in a triple] = 0}."""
    S = {tuple(sorted(T)) for T in tri}
    rows = []
    for b, c in itertools.combinations(range(k), 2):
        v = 0
        for a in range(k):
            if a != b and a != c and tuple(sorted((a, b, c))) in S:
                v |= 1 << a
        if v:
            rows.append(v)
    basis = []                                             # GF(2) rank of the constraint rows
    for v in rows:
        for b in basis:
            v = min(v, v ^ b)
        if v:
            basis.append(v)
    return k - len(basis)


def floor(tri, k):
    d = radical_dim(tri, k)
    return 0 if d == k else 2 * (k - d) + 1


if __name__ == "__main__":
    import random, sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "classification"))
    import tcount as TC
    rng = random.Random(7); worst = 0; tight = 0; n = 0
    for _ in range(300):
        k = rng.randint(3, 6); tri = [T for T in itertools.combinations(range(k), 3) if rng.random() < 0.3]
        if not tri:
            continue
        e = TC.tcount_rm({frozenset(q + 1 for q in T): 4 for T in tri}, k); f = floor(tri, k); n += 1
        assert f <= e, (k, tri, f, e)                       # a lower bound must never exceed the exact count
        tight += f == e; worst = max(worst, e - f)
    print(f"{n} random CCZ layers, k<=6: floor <= exact on all; floor == exact on {tight}; largest gap {worst}")
    print("checks: disjoint 2 CCZ ->", floor([(0, 1, 2), (3, 4, 5)], 6), "(13);",
          "cycle {012,034,135} ->", floor([(0, 1, 2), (0, 3, 4), (1, 3, 5)], 6), "(paper: floor 13, true 15)")
