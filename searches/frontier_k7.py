"""k>=7 frontier sweep, certified. Every CCZ layer on k qubits with m CCZ gates, split into shards: for each layer, a
T-count upper bound by greedy descent with restarts over its Reed-Muller coset (tcount7.coset), and the one-check
factory n = bound + 1 (rounded up to even) built from that compilation and verified with is_factory. Exact T-counts are
not computable at k=7 (the coset has dimension 64), so a result is a verified factory whose target's T-count is at most
n - 1; it is on the frontier only if that T-count is certified separately.
    python searches/frontier_k7.py --k 7 --m 3 --shard 0 --nshards 64
"""
import os, sys, json, time, random, argparse, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tcount7 as T  # noqa: E402
from isotropy_floor import floor  # noqa: E402
import unified_search as us  # noqa: E402
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def bound(tri, k, restarts=40, seed=0):
    b, C0, nc = T.coset({frozenset(x): 4 for x in tri}, k)
    rng = random.Random(seed); best, bv = bin(b).count("1"), b
    for _ in range(restarts):
        v = b
        for g in C0:
            if rng.random() < 0.5:
                v ^= g
        w = bin(v).count("1"); improved = True
        while improved:
            improved = False
            for g in C0:
                u = v ^ g; wu = bin(u).count("1")
                if wu < w:
                    v, w, improved = u, wu, True
        if w < best:
            best, bv = w, v
    return best, [y for y in range(1, nc + 1) if bv >> (y - 1) & 1]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=7); ap.add_argument("--m", type=int, required=True)
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--out-dir", default=os.path.join(ROOT, "outputs", "frontier_k7"))
    a = ap.parse_args(); os.makedirs(a.out_dir, exist_ok=True)
    k = a.k; triples = list(itertools.combinations(range(k), 3)); t0 = time.time()
    best = {}; done = 0
    for i, tri in enumerate(itertools.combinations(triples, a.m)):
        if i % a.nshards != a.shard:
            continue
        if len({q for x in tri for q in x}) != k:
            continue                                      # fewer than k qubits used: covered at smaller k
        lo = floor(tri, k)                                # proved lower bound (Mazumder isotropy floor)
        w, cols = bound(tri, k); done += 1
        if w != lo:
            continue                                      # T-count not certified: skip (could be a poor compilation)
        n = w + 1 + (w + 1) % 2
        if n not in best or w > best[n]["tbound"]:
            cols2 = cols + ([0] if len(cols) % 2 else [])
            g = [(frozenset({q for q in range(k) if c >> q & 1} | {k}), 1) for c in cols2]
            if us.is_factory(dict(stage="frontier", k=k, r=1, gates=g, params=""), 3):
                best[n] = dict(n=n, tbound=w, tcount_certified=True, k=k, N=k + 1, ccz=[list(x) for x in tri],
                               gates=[[sorted(S), c] for S, c in g])
    out = os.path.join(a.out_dir, f"k{k}_m{a.m}_s{a.shard}.json")
    json.dump(dict(k=k, m=a.m, shard=a.shard, nshards=a.nshards, layers=done, note="only layers whose heuristic upper bound equals the isotropy floor (T-count exact) are kept", seconds=round(time.time() - t0, 1),
                   best=best), open(out, "w"), indent=1)
    print(f"k={k} m={a.m} shard {a.shard}: {done} layers, {time.time() - t0:.0f}s, best (n, T-count bound):",
          sorted((v["n"], v["tbound"]) for v in best.values()), flush=True)
