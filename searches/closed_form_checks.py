"""Checks behind the closed forms quoted in the main text (writes outputs/checks/closed_forms.json).
- binomial sums: all weights on N qubits give sum_{w>=i} C(N-i,w-i) = 2^{N-i}; odd weights give 2^{N-i-1}
- Theorem 4 holds for all weights at N = l+1 and odd weights at N = l+2, l = 2,3,4
- families built and verified with is_factory and distance: [[2^{l+1}-2,1,2]], [[2^{l+1}-1,1,3]], [[2^l,l,2]]
"""
import sys, os, json, itertools
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search as us
out = {"binomial": all(sum(comb(N - i, w - i) for w in range(i, N + 1)) == 2 ** (N - i) and
                       (N == i or sum(comb(N - i, w - i) for w in range(i, N + 1) if w % 2) == 2 ** (N - i - 1))
                       for N in range(2, 12) for i in range(1, N + 1))}
def sym(N, W): return [(frozenset(S), 1) for w in W for S in itertools.combinations(range(N), w)]
fam = {}
for l in (2, 3, 4):
    N = l + 1; g = sym(N, range(1, N + 1)); f = [(S, c) for S, c in g if not S <= {0}]
    a = dict(stage="x", k=1, r=N - 1, gates=f, params="")
    N2 = l + 2; g2 = sym(N2, range(1, N2 + 1, 2)); f2 = [(S, c) for S, c in g2 if not S <= {0}]
    b = dict(stage="x", k=1, r=N2 - 1, gates=f2, params="")
    h = [(frozenset(q - 1 if q else N - 1 for q in S), c) for S, c in g if 0 in S]
    hc = dict(stage="x", k=N - 1, r=1, gates=h, params="")
    fam[l] = {"[[2^{l+1}-2,1,2]]": (len(f), us.is_factory(a, l), str(us.distance(a, l))),
              "[[2^{l+1}-1,1,3]]": (len(f2), us.is_factory(b, l), str(us.distance(b, l))),
              "[[2^l,l,2]]": (len(h), us.is_factory(hc, l), str(us.distance(hc, l)))}
out["families"] = fam
json.dump(out, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs", "checks", "closed_forms.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
