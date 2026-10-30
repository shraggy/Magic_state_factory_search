"""Exact minimal cost of synthesising a k-qubit level-3 diagonal target from CCZ, CS and T gates with arbitrary CNOTs
and Pauli X in between (unitary circuits, no ancillas or measurements), up to a final diagonal Clifford.

With CNOTs anywhere, each gate acts on affine parities of the inputs: CCZ(l1,l2,l3), CS(l1,l2) or T(l1). Modulo free
diagonal Cliffords the outputs form Q_k = F_2^m (at l=3 every coordinate is one bit) and composing gates adds their
elements. The gate set G is invariant under CNOT frames and X, so the minimal cost is constant on Clifford classes:
Dijkstra on the class graph (edges a -> class(rep(a) + g), weight of g) gives the exact minimum for every class.
Costs: CCZ and CS 8 inputs each ([[8,3,2]], CS by a Y measurement), T 5 (Meier-Eastin-Knill); error 28, 28, 9.

    python searches/min_synthesis.py --k 5        # writes outputs/min_synthesis_l3_k5.json
"""
import os, sys, json, heapq, argparse
from itertools import combinations, product
import numpy as np
P_TIME = None                                  # set to p to minimise time = inputs / acceptance instead of inputs

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_search import poly_from_table  # noqa: E402
from enumerate_targets import orbit_table  # noqa: E402

COST = {3: (8, 28), 2: (8, 28), 1: (5, 9)}     # gate size -> (inputs, error coefficient)
# time per gate = injected states / acceptance: [[8,3,2]] accepts with 1 - 8p; MEK (10 inputs, 2 T) with 1 - 10p
def time_cost(size, p):
    return 8 / (1 - 8 * p) if size in (2, 3) else 5 / (1 - 10 * p)
COEF = {3: 4, 2: 2, 1: 1}                      # phase pi/4 * COEF * product of the parities


def gate_elements(k, C, off):
    """Q-index of every CCZ / CS / T on affine parities, with its (inputs, error)."""
    xs = range(1 << k)
    forms = [(v, b) for v in range(1, 1 << k) for b in (0, 1)]
    val = {f: [((bin(x & f[0]).count("1") & 1) ^ f[1]) for x in xs] for f in forms}
    idx_of = {}
    pos = {T: (o, sh) for o, (T, sh, e) in zip(off, C)}
    for size in (1, 2, 3):
        for fs in combinations(forms, size):
            tab = [(COEF[size] * int(all(val[f][x] for f in fs))) % 8 for x in xs]
            a = poly_from_table(tab, k, 8)
            idx = 0
            for T, c in a.items():
                if T in pos:
                    o, sh = pos[T]
                    idx |= ((c >> sh) & 1) << o
            w = (time_cost(size, P_TIME), COST[size][1]) if P_TIME else COST[size]
            if idx and (idx not in idx_of or w < idx_of[idx]):
                idx_of[idx] = w
    return idx_of


def min_costs(k):
    C, off, bits, lab = orbit_table(3, k, log=lambda *a: None)
    G = gate_elements(k, C, off)
    reps = {}
    for i in np.unique(lab):
        reps[int(i)] = int(i)                       # the label is the smallest index in the class: a member
    dist = {int(lab[0]): (0, 0)}
    heap = [(0, 0, int(lab[0]))]
    gi = np.array(list(G.keys()), dtype=np.int64); gw = list(G.values())
    while heap:
        n, c, a = heapq.heappop(heap)
        if dist.get(a, (1 << 30, 0)) < (n, c):
            continue
        nb = lab[reps[a] ^ gi]
        for b, (wn, wc) in zip(nb.tolist(), gw):
            cand = (n + wn, c + wc)                 # wn: inputs, or time if P_TIME is set
            if cand < dist.get(b, (1 << 30, 1 << 30)):
                dist[b] = cand
                heapq.heappush(heap, (cand[0], cand[1], b))
    return dist, len(G)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--p", type=float, default=None, help="minimise time = inputs / acceptance at this input error")
    a = ap.parse_args()
    globals()["P_TIME"] = a.p
    dist, ng = min_costs(a.k)
    tag = "" if a.p is None else f"_time_p{a.p:g}"
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs", f"min_synthesis_l3_k{a.k}{tag}.json")
    json.dump({str(o): dict(n=v[0], c=v[1]) for o, v in dist.items()}, open(out, "w"), indent=1)
    print(f"k={a.k}: {ng} distinct gates on affine parities; exact minimal synthesis cost for {len(dist)} classes")
