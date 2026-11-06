"""GL(k,2)-invariant fingerprints of CCZ layers (pure-cubic level-3 targets; up to CNOTs on the outputs and a diagonal
Clifford, a CCZ layer is a GL(k,2)-orbit of cubic forms in RM(3,k)/RM(2,k)). For a cubic form f with triples T and
a in GF(2)^k, the derivative D_a f is a quadratic form whose alternating matrix is B_a[j][l] = sum_{(i,j,l) in T} a_i.
Fingerprint: the multiset of rank(B_a) over a != 0 and the multiset of sorted (rank B_a, rank B_b, rank B_a+b) over
pairs a < b. Equal layers in one orbit have equal fingerprints; different fingerprints prove different orbits (the
converse is not proved in general: distinct fingerprints give a lower bound on the number of orbits; at k = 7, 8 the
Burnside count of searches/k78_class_count.py gives exactly 6 and 20 classes, which the fingerprints of the target list
reach, so there they separate every class). Known-answer check: the
genuine pure-CCZ classes of outputs/class_table_l3.json (k <= 6).
    python searches/k78_orbits.py --sample 20000 --seed 0   -> outputs/k78_pipeline/orbits.json
"""
import os, sys, json, random, itertools, argparse
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)


def rank2(rows):
    B = []
    for v in rows:
        for b in B:
            v = min(v, v ^ b)
        if v:
            B.append(v)
    return len(B)


def ranks(tri, k):
    out = [0] * (1 << k)
    for a in range(1, 1 << k):
        M = [0] * k
        for T in tri:
            for i in T:
                if a >> i & 1:
                    j, l = [q for q in T if q != i]
                    M[j] ^= 1 << l; M[l] ^= 1 << j
        out[a] = rank2(M)
    return out


def fingerprint(tri, k):
    r = ranks(tri, k)
    c1 = Counter(r[1:])
    c2 = Counter(tuple(sorted((r[a], r[b], r[a ^ b]))) for a in range(1, 1 << k) for b in range(a + 1, 1 << k))
    return json.dumps([sorted(c1.items()), sorted(c2.items())])


def genuine(tri, k):
    """Essential dimension k (no CNOT circuit reduces the layer to fewer outputs): f mod RM(2) depends on fewer than
    k variables iff some a != 0 has deg D_a f <= 1, i.e. B_a = 0 (the radical of the cubic form, isotropy_floor)."""
    r = ranks(tri, k)
    return all(r[a] for a in range(1, 1 << k))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--sample", type=int, default=20000); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--kmax", type=int, default=8); ap.add_argument("--kmin", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(ROOT, "outputs", "k78_pipeline", "orbits.json")); a = ap.parse_args()
    rng = random.Random(a.seed); res = {}
    # known answer: genuine pure-CCZ classes at k <= 6 in the exact class table
    import re
    tab = json.load(open(os.path.join(ROOT, "outputs", "class_table_l3.json")))
    known = {}
    for r in tab:
        s = r["spec"].split(",")
        if all(re.fullmatch(r"CCZ\d{3}", x) for x in s):
            known.setdefault(r["k"], []).append(r["spec"])
    for k in range(a.kmin, a.kmax + 1):
        tri = list(itertools.combinations(range(k), 3)); fps = {}
        for _ in range(a.sample):
            L = [T for T in tri if rng.random() < rng.random()]
            if not L or not genuine(L, k):
                continue
            fp = fingerprint(L, k)
            fps.setdefault(fp, L)
        res[k] = dict(distinct_fingerprints=len(fps), reps=list(fps.values()),
                      class_table_genuine=len(known.get(k, [])) if k <= 6 else None)
        print(k, "distinct fingerprints (lower bound on genuine orbits):", len(fps),
              "| class table:", res[k]["class_table_genuine"], flush=True)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(dict(settings=vars(a), result=res), open(a.out, "w"), indent=0)
