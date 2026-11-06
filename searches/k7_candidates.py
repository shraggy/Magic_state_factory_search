"""k=7 frontier candidates: for random CCZ layers on 7 qubits, a heuristic upper bound on the T-count (repeated
greedy descent, many restarts so the bound is as low as we can make it) and the isotropy floor. A layer is a candidate
for n = bound + 1 (+1 if odd, one check) when its bound exceeds the frontier below n. Keeps, for each n, the layers
with the largest bound; these are then decided exactly by tcount7.exists_below.
    python searches/k7_candidates.py --seed 0 --layers 20000
"""
import os, sys, json, random, itertools, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import frontier_k7 as F, isotropy_floor as I
ap = argparse.ArgumentParser(); ap.add_argument("--seed", type=int, default=0); ap.add_argument("--layers", type=int, default=20000)
ap.add_argument("--k", type=int, default=7)
ap.add_argument("--out", default="outputs/k7_candidates"); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
K = a.k; rng = random.Random(a.seed); tri = list(itertools.combinations(range(K), 3)); keep = {}
for _ in range(a.layers):
    m = rng.randint(4, 4 + 2 * K); L = sorted(rng.sample(tri, m))
    if len({q for x in L for q in x}) != K: continue
    ub, _ = F.bound(L, K, restarts=200, seed=rng.randrange(1 << 30)); lo = I.floor(L, K)
    n = ub + 1 + (ub + 1) % 2
    keep.setdefault(n, [])
    if len(keep[n]) < 5: keep[n].append(dict(n=n, upper=ub, floor=lo, ccz=L))
json.dump(keep, open(os.path.join(a.out, f"cand_k{K}_s{a.seed}.json"), "w"), indent=1)
print("seed", a.seed, "| n with candidates:", sorted(keep), flush=True)
