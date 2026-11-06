"""Sequential condition (App. S2.3) in Phi_T form, checked against brute force.

A sequential parent on qubits q_1..q_N applies, at step j, every rotation on {q_j..q_N} that contains q_j and has
weight in W_j (signs sigma_j). After steps 1..j the anchor q_j returns to |+> up to a diagonal Clifford iff, for every
T that contains q_j, lies in {q_j..q_N} and has |T| <= l, Phi_T is a multiple of 2^{e_|T| - l} 2 pi (Eq. 4 with q_j as
the only new check). By symmetry Phi_T depends on |T| = i+1 only:
   Phi_T / theta = 2^{i+1} [ sigma_j sum_{w in W_j} C(N-j-i, w-1-i) + sum_{q<j} sigma_q sum_{w in W_q} C(N-q-1-i, w-2-i) ],
the first term counting step-j rotations that contain T, the second counting earlier rotations (anchor q_q, which must
also contain q_j and the i others). Checked against direct evaluation of the reduced state of q_j on all basis states.
    python searches/sequential_check.py     -> outputs/checks/sequential_check.json
"""
import os, json, random, itertools
from math import comb
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def e(t, l):
    return min(l, max(l - 3 + t, t))


def predicted(N, l, Ws, sig, j):
    """Does anchor q_j (1-based) return after steps 1..j, by the Phi_T formula?"""
    for i in range(0, l):                               # |T| = i + 1 <= l
        if i > N - j:
            break
        v = sig[j - 1] * sum(comb(N - j - i, w - 1 - i) for w in Ws[j - 1] if w - 1 - i >= 0)
        v += sum(sig[q - 1] * sum(comb(N - q - 1 - i, w - 2 - i) for w in Ws[q - 1] if w - 2 - i >= 0) for q in range(1, j))
        if (2 ** (i + 1) * v) % (2 ** (e(i + 1, l) + 1)):   # Phi_T/theta mod 2^{e+1}  (2pi = 2^{l+1} theta)
            return False
    return True


def direct(N, l, Ws, sig, j):
    """Brute force: after steps 1..j, is q_j in |+> up to a diagonal Clifford, for every value of q_{j+1}..q_N, given
    q_1..q_{j-1} already disentangled? Check: the phase difference between q_j=1 and q_j=0, as a function of the
    remaining qubits, is a Clifford phase function (degree-t coefficients in multiples of 2^{e_t}... in theta units)."""
    full = 2 ** (l + 1)
    rots = [(S, sig[t]) for t in range(j) for w in Ws[t] for S in itertools.combinations(range(t, N), w) if S[0] == t]
    rest = list(range(j, N))
    def ph(x):
        return sum(2 * c for S, c in rots if sum(x[q] for q in S) % 2) % full
    # Moebius coefficients of D(y) = ph(..., q_j=1, y) - ph(..., q_j=0, y), earlier anchors set to 0 (they are disentangled)
    f = {}
    for y in itertools.product((0, 1), repeat=len(rest)):
        x = [0] * N
        for q, b in zip(rest, y): x[q] = b
        x0 = list(x); x1 = list(x); x1[j - 1] = 1
        f[y] = (ph(x1) - ph(x0)) % full
    n = len(rest)
    for r in range(0, n + 1):
        for T in itertools.combinations(range(n), r):
            a = 0
            for rr in range(r + 1):
                for U in itertools.combinations(T, rr):
                    y = tuple(1 if t in U else 0 for t in range(n)); a += (-1) ** (r - rr) * f[y]
            a %= full
            t = r + 1                                    # the monomial is q_j * prod_{T} y
            if t > l:
                continue
            if a % (2 ** (e(t, l) + 1)):
                return False
    return True


rng = random.Random(5); agree = tot = pos = 0; bad = []
for _ in range(3000):
    l = rng.choice([3, 4]); N = rng.randint(3, 7); J = rng.randint(1, min(l + 1, N - 1))
    Ws = [sorted(rng.sample(range(1, N - t + 1), rng.randint(1, N - t))) for t in range(J)]
    sig = [rng.choice([1, -1]) for _ in range(J)]
    ok_prev = all(direct(N, l, Ws, sig, jj) for jj in range(1, J))
    if not ok_prev:
        continue
    a, b = predicted(N, l, Ws, sig, J), direct(N, l, Ws, sig, J)
    tot += 1; agree += a == b; pos += b
    if a != b and len(bad) < 5: bad.append(dict(l=l, N=N, W=Ws, sig=sig, pred=a, direct=b))
out = dict(tested=tot, agree=agree, anchor_returns=pos, mismatches=bad)
json.dump(out, open(os.path.join(ROOT, "outputs", "checks", "sequential_check.json"), "w"), indent=1)
print(out)
