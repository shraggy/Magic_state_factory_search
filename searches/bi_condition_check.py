import os
"""Exhaustive check of the borrowed-identity condition used in the main text:
a circuit of parity rotations Z_S(c_S theta), theta = pi/2^l, is a borrowed identity iff
    n_T = sum_{S >= T} c_S  is in  2^{l-|T|+1} Z   for every support T with 1 <= |T| <= l,
compared with brute-force evaluation of C|+>^N on all 2^N basis states. Exhaustive over every coefficient vector for
(l, N) = (2, 3), (3, 2), (3, 3), (4, 2); random for larger cases.
    python searches/bi_condition_check.py
"""
import itertools, random
import numpy as np

def run(l, N, exhaustive, samples=0, seed=0):
    sup = [S for w in range(1, N + 1) for S in itertools.combinations(range(N), w)]
    par = np.array([[bin(x & sum(1 << q for q in S)).count("1") & 1 for S in sup] for x in range(1 << N)])
    Ts = [T for t in range(1, min(l, N) + 1) for T in itertools.combinations(range(N), t)]
    sup_has = np.array([[set(T) <= set(S) for S in sup] for T in Ts], dtype=np.int64)
    mods = np.array([1 << (l - len(T) + 1) for T in Ts])
    it = itertools.product(range(1 << l), repeat=len(sup)) if exhaustive else \
        (tuple(random.Random(seed + i).randrange(1 << l) for _ in sup) for i in range(samples))
    agree = tot = bis = 0
    for c in it:
        c = np.array(c)
        phase = (par @ c) % (1 << l)                # each odd overlap adds 2*theta*c = 2*pi*c/2^l
        brute = bool(np.all(phase == phase[0]))
        cond = bool(np.all((sup_has @ c) % mods == 0))
        agree += brute == cond; tot += 1; bis += brute
    return tot, bis, agree

if __name__ == "__main__":
    import json
    res = []
    for l, N in ((2, 3), (3, 2), (3, 3), (4, 2)):
        tot, bis, agree = run(l, N, True)
        res.append(dict(l=l, N=N, circuits=tot, borrowed_identities=bis, agree=agree))
        print(f"l={l} N={N} exhaustive: {tot} circuits, {bis} borrowed identities, condition agrees on {agree}", flush=True)
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs", "checks", "bi_condition_check.json")
    json.dump(dict(cells=res, circuits=sum(r["circuits"] for r in res), agree=sum(r["agree"] for r in res),
                   borrowed_identities=sum(r["borrowed_identities"] for r in res)), open(out, "w"), indent=1)
