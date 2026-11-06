"""Check the factory form of the borrowed-identity condition against brute force.

A circuit of parity rotations Z_S(c_S theta) on k outputs and N-k checks, with no odd rotation on the outputs alone,
is a distance-2 factory when its checks return to |+> up to a diagonal Clifford. The claimed condition: for every
support T with 1 <= |T| <= l that touches a check,
    n_T = sum_{S >= T} c_S  lies in  2^{e_|T| - |T| + 1} Z,   e_t = min(l, max(l - 3 + t, t)),
where 2^{e_t} is the coefficient modulus at which a degree-t monomial is a free diagonal Clifford. At l = 3 every
modulus is 2: every set of at most three qubits touching a check lies in an even number of rotations (weighted by c).
Compared with unified_search.is_factory(brute_force=True), which evaluates the phase table on all 2^N basis states.
    python searches/factory_condition_check.py
"""
import sys, os, random, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search as us

def cond(gates, k, N, l):
    for t in range(1, min(l, N) + 1):
        e = min(l, max(l - 3 + t, t)); m = e - t + 1
        for T in itertools.combinations(range(N), t):
            if set(T) <= set(range(k)): continue
            if m > 0 and sum(c for S, c in gates if set(T) <= S) % (1 << m): return False
    return True

if __name__ == "__main__":
    rnd = random.Random(11); res = {}
    for l in (2, 3, 4):
        agree = tot = fac = 0
        for _ in range(4000 if l < 4 else 2500):
            N = rnd.randint(2, 6); k = rnd.randint(0, N - 1)
            sup = [frozenset(S) for w in range(1, N + 1) for S in itertools.combinations(range(N), w) if not set(S) <= set(range(k))]
            g = [(S, rnd.randrange(1, 1 << l)) for S in rnd.sample(sup, rnd.randint(1, min(len(sup), 12)))]
            if rnd.random() < 0.5:   # bias toward factories: add the complement needed at random supports
                g += [(S, (-c) % (1 << l)) for S, c in rnd.sample(g, min(3, len(g)))]
            f = dict(stage="x", k=k, r=N - k, gates=g, params="")
            b = us.is_factory(f, l, brute_force=True); a = cond(g, k, N, l)
            agree += a == b; tot += 1; fac += b
        res[l] = (tot, fac, agree)
        print(f"l={l}: {tot} random circuits, {fac} factories (brute force), condition agrees on {agree}", flush=True)
    import json
    json.dump({str(l): dict(circuits=t, factories=f, agree=a) for l, (t, f, a) in res.items()},
              open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs", "checks",
                                "factory_condition_check.json"), "w"), indent=1)
