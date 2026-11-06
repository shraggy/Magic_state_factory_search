"""Check the sequential cuts of the parents drawn in Figs. 1 and S1 against the factory condition, by direct simulation.

A parent on N qubits at level l is a set of supports (each a pi/2^l rotation). A cut after layer j measures q_0..q_{j-1}
(checks) and keeps as inputs every rotation touching a check; the factory is valid when, on |+>^N, the checks return
to |+> up to a diagonal Clifford, checked here on the full phase table: phi(x) restricted to the check-dependent part
must be a Clifford phase polynomial. Output: n inputs, k outputs, and the output T-count lower bound n - 1 is printed.
    python searches/parent_cuts_check.py
"""
import itertools, sys
def phase(cols, l, N, x):                                   # in units of 2*pi/2^l ... phase 2*theta per odd parity
    return sum(1 for S in cols if sum(x[q] for q in S) % 2) * 2 % (2 ** (l + 1))   # units of theta, mod 2pi = 2^(l+1) theta
def clifford_free(f, l, qs):
    """f: dict x -> phase in units of theta (mod 2^(l+1)). True if f is a diagonal Clifford phase polynomial on qs,
    i.e. its monomial coefficients are multiples of the Clifford moduli: degree 1 in units pi/2, degree 2 pi, degree 3 2pi."""
    import itertools as it
    M = 2 ** (l + 1)                                        # 2pi in units of theta
    # Moebius transform over the bits: coefficient a_T with f(x) = sum_T a_T prod x
    n = len(qs); a = {}
    for r in range(n + 1):
        for T in it.combinations(range(n), r):
            v = 0
            for rr in range(r + 1):
                for U in it.combinations(T, rr):
                    x = tuple(1 if i in U else 0 for i in range(n)); v += (-1) ** (r - rr) * f[x]
            a[T] = v % M
    mod = {1: M // 4, 2: M // 2}                            # pi/2, pi; degree >= 3: 2pi
    return all(a[T] % mod.get(len(T), M) == 0 for T in a if T)
def cut(cols, l, N, j):
    checks = list(range(j)); outs = list(range(j, N))
    keep = [S for S in cols if set(S) & set(checks)]
    # checks must return to |+> up to a diagonal Clifford for every output basis value
    ok = True
    for y in itertools.product((0, 1), repeat=len(outs)):
        f = {}
        for xc in itertools.product((0, 1), repeat=len(checks)):
            x = [0] * N
            for q, b in zip(checks, xc): x[q] = b
            for q, b in zip(outs, y): x[q] = b
            f[xc] = phase(keep, l, N, x)
        # subtract the xc = 0 value (a phase on the outputs), then the rest must be Clifford on the checks
        z = f[tuple([0] * len(checks))]; g = {k: (v - z) % 2 ** (l + 1) for k, v in f.items()}
        ok &= clifford_free(g, l, checks)
    return len(keep), len(outs), ok
def parent(N, l, odd):
    cols = [S for w in range(N, 0, -1) for S in itertools.combinations(range(N), w) if not odd or w % 2]
    return cols
for name, N, l, odd in (("l=3, 4 qubits, all weights (Fig. 1)", 4, 3, False),
                        ("l=3, 5 qubits, odd weights", 5, 3, True),
                        ("l=4, 5 qubits, all weights", 5, 4, False)):
    cols = parent(N, l, odd)
    print(name, len(cols), "rotations")
    for j in range(1, N):
        print("   cut after layer", j, "-> [[n,k]] =", cut(cols, l, N, j)[:2], "valid:", cut(cols, l, N, j)[2])
