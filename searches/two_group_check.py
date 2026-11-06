"""Check of Theorem 11 (two-group master condition) against Theorem 2 by direct evaluation of the phase table.

Random two-group circuits (l in {2,3,4}, N <= 7, k < N, random admissible pairs (w_O, w_S) and signs, plus structured
families with weight spacings s_total, s_O <= 4): a circuit is a borrowed identity by direct evaluation of
C|+>^N on all 2^N basis states, and by Eq. (asymmetric-master); the two must agree on every instance.
    python searches/two_group_check.py     -> outputs/checks/two_group_check.json
"""
import os, json, random, itertools
from math import comb
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def bi_direct(N, k, l, W, sig):
    full = 2 ** (l + 1)                                  # 2 pi in units of theta
    rots = [(So + Sc, sig[(a, b)]) for (a, b) in W for So in itertools.combinations(range(k), a)
            for Sc in itertools.combinations(range(k, N), b) if a + b]
    ph = {sum(2 * c for S, c in rots if sum(x[q] for q in S) % 2) % full for x in itertools.product((0, 1), repeat=N)}
    return len(ph) == 1


def master(N, k, l, W, sig):
    for io in range(k + 1):
        for js in range(N - k + 1):
            if (io, js) == (0, 0) or io + js > l:
                continue
            t = sum(sig[w] * comb(k - io, w[0] - io) * comb(N - k - js, w[1] - js) for w in W if w[0] >= io and w[1] >= js)
            if t % 2 ** (l - io - js + 1):
                return False
    return True


out = dict(random=[0, 0, 0], structured=[0, 0, 0])            # agree, total, borrowed identities
for mode, seed, reps in (("random", 1, 1500), ("structured", 7, 3000)):
    rng = random.Random(seed)
    for _ in range(reps):
        l = rng.choice([2, 3, 4]); N = rng.randint(3, 7); k = rng.randint(1, N - 1)
        if mode == "random":
            pairs = [(a, b) for a in range(k + 1) for b in range(N - k + 1) if a + b > 0]
            W = rng.sample(pairs, rng.randint(1, len(pairs))); sig = {w: rng.choice([1, -1]) for w in W}
        else:
            st, so = rng.randint(1, 4), rng.randint(1, 4)
            W = [(a, b) for a in range(k + 1) for b in range(N - k + 1)
                 if a + b > 0 and (a + b) % st == 1 % st and (a == 0 or a % so == 1 % so)]
            if not W:
                continue
            sig = {w: rng.choice([1, -1, 1, 1]) for w in W}
        a, b = bi_direct(N, k, l, W, sig), master(N, k, l, W, sig)
        out[mode][0] += a == b; out[mode][1] += 1; out[mode][2] += a
os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
json.dump(out, open(os.path.join(ROOT, "outputs", "checks", "two_group_check.json"), "w"), indent=1)
print(out)
