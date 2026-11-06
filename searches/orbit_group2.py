"""Group-symmetric stage with a group on the checks as well: outputs symmetric under G <= S_k, checks under H <= S_r,
a circuit symmetric under G x H. A cell is (G-orbit of an output subset A, H-orbit of a nonempty check subset W); it
holds every rotation A' u W' with A' in the first orbit and W' in the second, coefficient 1 (l=3). H = S_r gives the
cells of orbit_group.py (check orbits = check weights). The factory condition at l=3 (every set T of <=3 qubits that
meets a check lies in an even number of rotations) is one GF(2) row per (G-orbit of T∩O, H-orbit of T∩S), since the
count is G x H invariant; the deposited output digit on T ⊆ O is the parity of #{rotations ⊇ T}. Same DP sweep and
CP-SAT target mode as orbit_group.py, on these cells.
    python searches/orbit_group2.py --selftest
"""
import os, sys, json, time, argparse
from itertools import combinations, permutations
from math import comb
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import orbit_group as OG
import unified_search as us
L = 3


def check_groups(r):
    """Menu of check groups H <= S_r: S_r (as orbit_group), every Young subgroup, C_r, D_r, trivial (deduplicated by
    their orbits on check subsets)."""
    out = []
    for p in OG.partitions(r):
        out.append(("H:" + "x".join(f"S{q}" for q in p), OG.young(p)))
    if r >= 3:
        out.append((f"H:C{r}", OG.cyclic(r))); out.append((f"H:D{r}", OG.dihedral(r)))
    seen, menu = set(), []
    for name, H in out:
        orbs, of = OG.subset_orbits(H, r)
        key = tuple(of)
        if key in seen:
            continue
        seen.add(key); menu.append((name, H, orbs, of))
    return menu


def build2(k, r, orbs, corbs):
    """cells (oi, ci) over output orbits orbs (masks on k) and nonempty check orbits corbs (masks on r)."""
    cnz = [ci for ci, o in enumerate(corbs) if o[0] != 0]
    cells = [(oi, ci) for oi in range(len(orbs)) for ci in cnz]
    cost = [len(orbs[oi]) * len(corbs[ci]) for oi, ci in cells]
    cntO = [[sum(1 for A in o if A & B == B) for B in range(1 << k)] for o in orbs]
    cntS = [[sum(1 for W in o if W & J == J) for J in range(1 << r)] for o in corbs]
    H, D, Dorb = [], [], []
    for bi, o in enumerate(orbs):
        B = o[0]; i = bin(B).count("1")
        if i <= 2:
            for cj, co in enumerate(corbs):
                J = co[0]; j = bin(J).count("1")
                if 1 <= j <= 3 - i:
                    row = 0
                    for x, (oi, ci) in enumerate(cells):
                        if (cntO[oi][B] * cntS[ci][J]) & 1:
                            row |= 1 << x
                    H.append(row)
        if 1 <= i <= 3:
            row = 0
            for x, (oi, ci) in enumerate(cells):
                if (cntO[oi][B] * len(corbs[ci])) & 1:
                    row |= 1 << x
            D.append(row); Dorb.append(bi)
    return cells, cost, H, D, Dorb


def gates2(k, r, orbs, corbs, cells, x):
    g = []
    for c in x:
        oi, ci = cells[c]
        for A in orbs[oi]:
            As = frozenset(i for i in range(k) if A >> i & 1)
            for W in corbs[ci]:
                g.append((As | frozenset(k + j for j in range(r) if W >> j & 1), 1))
    return g


_orig_build = OG.build


def with_cells(k, r, orbs, corbs, fn, *a):
    """Run an orbit_group routine (dp_sweep / cpsat_target) on G x H cells by swapping OG.build for build2."""
    OG.build = lambda k_, r_, orbs_: build2(k_, r_, orbs_, corbs)
    try:
        return fn(k, r, orbs, *a)
    finally:
        OG.build = _orig_build


def selftest(trials=400, seed=11):
    import random
    rnd = random.Random(seed); agree = tot = fac = 0; same = 0
    for _ in range(trials):
        k = rnd.randint(1, 4); r = rnd.randint(1, 3)
        name, G, orbs, of, _ = rnd.choice(OG.group_menu(k)) if k >= 1 else None
        hn, Hg, corbs, cof = rnd.choice(check_groups(r))
        cells, cost, H, D, Dorb = build2(k, r, orbs, corbs)
        x = [i for i in range(len(cells)) if rnd.random() < 0.3]
        g = gates2(k, r, orbs, corbs, cells, x)
        if not g:
            continue
        pred = all(bin(row & sum(1 << i for i in x)).count("1") % 2 == 0 for row in H)
        f = dict(k=k, r=r, gates=g, stage="t", params="")
        truth = us.is_factory(f, L) and us.is_factory(f, L, brute_force=True)
        tot += 1; agree += pred == truth; fac += truth
    # H = S_r must reproduce orbit_group's cells exactly
    for k in (2, 3, 4):
        for r in (1, 2, 3):
            for name, G, orbs, of, _ in OG.group_menu(k):
                full = [h for h in check_groups(r) if h[0] == "H:S" + str(r)][0]
                c1 = _orig_build(k, r, orbs); c2 = build2(k, r, orbs, full[2])
                same += (sorted(c1[1]) == sorted(c2[1]) and len(c1[2]) == len(c2[2]))
    print(f"condition: {agree}/{tot} random G x H circuits agree with is_factory (monomial and brute force); {fac} factories;"
          f" H=S_r reproduces orbit_group cell costs and row counts on {same} (k,r,G) cases", flush=True)
    return agree == tot


def dp_check(kmax=3, rmax=3):
    """Every DP optimum over G x H cells, expanded to rotations, must pass is_factory (monomial and brute force)."""
    nf = nv = 0
    for k in range(1, kmax + 1):
        for r in range(1, rmax + 1):
            for name, G, orbs, of, _ in OG.group_menu(k):
                for hn, Hg, corbs, cof in check_groups(r):
                    res, info = with_cells(k, r, orbs, corbs, OG.dp_sweep)
                    if res is None:
                        continue
                    out, cells, cost, Dorb = res
                    for dig, (n, x) in out.items():
                        if not x:
                            continue
                        g = gates2(k, r, orbs, corbs, cells, x)
                        f = dict(k=k, r=r, gates=g, stage="t", params="")
                        nf += 1; nv += us.is_factory(f, L) and us.is_factory(f, L, brute_force=True)
    print(f"DP optima over G x H cells, k<={kmax}, r<={rmax}: {nv}/{nf} pass is_factory (monomial and brute force)", flush=True)
    return nv == nf


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--selftest", action="store_true"); a = ap.parse_args()
    if a.selftest:
        selftest(); dp_check()

