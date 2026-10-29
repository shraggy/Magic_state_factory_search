"""Stage (iii) of the unified search: one target per Clifford class of k-qubit diagonal level-l outputs.

Every element of Q_k (outputs modulo free diagonal Cliffords, see target_classes.py) is labelled by its orbit under
CNOT frames and Pauli X (the affine group AGL(k,2)), by min-label propagation over the Schreier graph of its
generators and their inverses. The orbit
count is checked against the exact Burnside count. Each orbit keeps the representative with the fewest gates, and
orbits that some CNOT frame reduces to fewer than k outputs are marked non-genuine (they are Clifford padding of a
smaller-k target). `canonical_id` maps any output polynomial to its orbit, so catalogue factories can be matched to
targets exactly (also an exact Clifford classification wherever the orbit table fits in memory).

    python searches/enumerate_targets.py --level 3 --k 5      # writes outputs/targets/l3_k5.json
"""
import os, sys, json, time, argparse
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from unified_search import cliff_exp, poly_from_table  # noqa: E402
from target_classes import coords, action_matrix, orbits_burnside  # noqa: E402

GEN = {2: ("S", "CZ"), 3: ("T", "CS", "CCZ"), 4: ("sT", "CT", "CCS", "CCCZ")}
MAX_BITS = 27                                            # orbit table of 2^27 int32 entries = 512 MB


def layout(C):
    """Bit offset of each coordinate in the packed index (coordinate j holds e_j bits)."""
    off, o = [], 0
    for _, _, e in C:
        off.append(o); o += e
    return off, o


def image_index(A, C, off, bits, chunk=1 << 22):
    """Packed index of A v for every packed v, as an int32 array of length 2^bits (computed in chunks)."""
    out = np.empty(1 << bits, dtype=np.int32)
    for s0 in range(0, 1 << bits, chunk):
        xs = np.arange(s0, min(s0 + chunk, 1 << bits), dtype=np.int64)
        res = np.zeros(len(xs), dtype=np.int64)
        for i, (_, _, ei) in enumerate(C):
            acc = np.zeros(len(xs), dtype=np.int64)
            for j, (_, _, ej) in enumerate(C):
                a = int(A[i, j])
                if a:
                    acc += a * ((xs >> off[j]) & ((1 << ej) - 1))
            res |= (acc & ((1 << ei) - 1)) << off[i]
        out[s0:s0 + len(xs)] = res
    return out


def generators(k):
    """Generators of AGL(k,2) as (matrix, shift): CNOT frames and a Pauli X (X|+> = |+>, so X is free)."""
    I = np.eye(k, dtype=np.int64)
    if k == 1:
        return [(I, 1)]
    cyc = np.roll(np.eye(k, dtype=np.int64), 1, axis=0)
    tv = np.eye(k, dtype=np.int64); tv[0, 1] = 1        # a transvection is its own inverse
    return [(cyc, 0), (cyc.T.copy(), 0), (tv, 0), (I, 1)]       # cyc.T = cyc^{-1}; X and tv are involutions


def orbit_table(L, k, log=print):
    C = coords(L, k)
    off, bits = layout(C)
    if bits > MAX_BITS:
        raise ValueError(f"|Q_{k}| = 2^{bits} at l={L}: orbit table too large")
    t = time.time()
    imgs = [image_index(action_matrix(g, L, k, C, b), C, off, bits) for g, b in generators(k)]
    lab = np.arange(1 << bits, dtype=np.int32)
    it = 0
    while True:
        new = lab.copy()
        for im in imgs:
            np.minimum(new, lab[im], out=new)
        for _ in range(4):                               # pointer jumping
            new = new[new]
        it += 1
        if np.array_equal(new, lab):
            break
        lab = new
    log(f"l={L} k={k}: |Q| = 2^{bits}, orbits labelled in {it} rounds, {time.time() - t:.1f}s")
    return C, off, bits, lab


def digits_of(idx, C, off):
    return [(idx >> o) & ((1 << e) - 1) for o, (_, _, e) in zip(off, C)]


def to_spec(dig, C, L):
    """Target string for coset_certify.parse_target: a_T = v * 2^{|T|-1} is v copies of the size-|T| gate."""
    toks = []
    for v, (T, sh, _) in zip(dig, C):
        name = GEN[L][len(T) - 1] + "".join(map(str, sorted(T)))
        toks += [name] * v
    return ",".join(toks)


def index_of_poly(P, k, L, C, off):
    """Packed index of an output polynomial {support: coefficient mod 2^L} (Newton/AND-monomial form)."""
    mod = 1 << L
    tab = [sum(c for T, c in P.items() if all(x >> q & 1 for q in T)) % mod for x in range(1 << k)]
    a = poly_from_table(tab, k, mod)
    idx = 0
    for o, (T, sh, e) in zip(off, C):
        v = a.get(T, 0)
        assert v % (1 << sh) == 0, "not a level-l diagonal output"
        idx |= ((v >> sh) % (1 << e)) << o
    return idx


def enumerate_targets(L, k, log=print):
    t0 = time.time()
    C, off, bits, lab = orbit_table(L, k, log)
    reps = np.unique(lab)
    exact = orbits_burnside(L, k)
    assert len(reps) == exact, f"orbit table gives {len(reps)} orbits, Burnside {exact}"
    xs = np.arange(1 << bits, dtype=np.int64)
    # cost of each element: number of gates, then number of nonzero coordinates
    gates = np.zeros(1 << bits, dtype=np.int64)
    for o, (_, _, e) in zip(off, C):
        gates += (xs >> o) & ((1 << e) - 1)
    # non-genuine: the orbit meets the subgroup of outputs that do not involve the last output qubit
    last = sum(((1 << e) - 1) << o for o, (T, _, e) in zip(off, C) if (k - 1) in T)
    sub = (xs & last) == 0
    padded = np.zeros(1 << bits, dtype=bool)
    padded[np.unique(lab[sub])] = True
    order = np.lexsort((xs, gates, lab))                 # by orbit, then fewest gates, then index
    first = order[np.r_[True, lab[order][1:] != lab[order][:-1]]]
    size = np.bincount(lab, minlength=1 << bits)
    out = []
    for r in first:
        orb = int(lab[r])
        if orb == 0:
            continue                                     # the trivial output
        dig = digits_of(int(r), C, off)
        out.append(dict(orbit=orb, genuine=not bool(padded[orb]), gates=int(gates[r]), orbit_size=int(size[orb]),
                        spec=to_spec(dig, C, L)))
    log(f"l={L} k={k}: {len(out)} nontrivial classes ({sum(o['genuine'] for o in out)} genuine), "
        f"{time.time() - t0:.1f}s total")
    return dict(l=L, k=k, classes=out, seconds=round(time.time() - t0, 2)), (C, off, lab)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--out-dir", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "outputs", "targets"))
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    res, _ = enumerate_targets(a.level, a.k)
    json.dump(res, open(os.path.join(a.out_dir, f"l{a.level}_k{a.k}.json"), "w"), indent=1)
