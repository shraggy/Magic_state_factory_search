"""Orbit-group stage (l=3, d=2): checks interchangeable (S_r), outputs symmetric under a chosen group G <= S_k.

Generalizes stage (i), two-group (S_k x S_r, searches/Two_group.py), toward gluing. Qubits 0..k-1 are outputs,
k..k+r-1 checks. A CELL is (omega, w): omega a G-orbit of output subsets A (the empty set included: check-only
rotations), w = 1..r a check weight; the cell holds every rotation Z_{A u W}(theta) with A in omega and W any w-subset
of the checks. A circuit is a set of cells, each all-in (coefficient 1) or out. At l=3 coefficient 1 is general: an
odd c differs from 1 by c-1 even, a Clifford rotation (inputs(2,3)=0) that changes neither the condition below nor the
output class; output-only rotations (w=0) are not part of a factory (as in unified_search.two_group_gates).

Condition (Eq. factory-condition / unified_search.is_factory at l=3): for every T, 1<=|T|<=3, that meets a check,
n_T = #{included rotations S >= T} is even. Write T = B u J (B outputs, J checks, |J| = j >= 1). A cell (omega, w)
contributes #{A in omega : A >= B} * C(r-j, w-j) rotations, which depends only on (G-orbit of B, j) because omega is
G-invariant and the checks are interchangeable. So the condition is one GF(2) row per (orbit beta of output subsets
with |B| <= 2, j = 1..min(3-|B|, r)):  sum_cells x_cell * cnt(omega, B) * C(r-j, w-j) = 0 (mod 2).
Output: the deposited Newton coefficient on T <= O is a_T = (-1)^{|T|+1} 2^{|T|-1} n_T (mod 8), so its digit in
Q_k = (Z_2)^{bits} (coords of searches/target_classes.py, e_T = 1 at l=3) is n_T mod 2
= sum_cells x_cell * cnt(omega, T) * C(r, w) (mod 2), constant on G-orbits of T. Both maps are GF(2)-linear in x.
Cost n = sum_cells x_cell |omega| C(r, w) (Z(pi/8) inputs, unified_search.inputs).

Two modes, both exact within their space:
  sweep  : for one (k, G, r), min n for EVERY G-invariant output, by a dynamic program over the GF(2) syndrome
           (condition part, output part): dp[s] = min cost of a cell subset with syndrome s; exhaustive over all
           2^cells subsets. Run when the syndrome space is <= 2^BITS_MAX. Every output is labelled with its exact
           class (cascade_class.labeller / orbit tables of enumerate_targets(_fiber)).
  target : for one class representative P (the spec of outputs/targets/l3_k*.json), deposit exactly P:
           G = Aut(P) (output permutations fixing P) or G = {e}; min n by CP-SAT (one worker, deterministic step
           limit, OPTIMAL or not is recorded). G = {e} contains every G' for that representative.
Every reported factory is expanded to explicit rotations and re-verified with unified_search.is_factory (monomial
form and brute-force phase table), unified_search.distance != 1, and the exact class label.

The orbit-folded condition is the d=2 counterpart (roles flipped: checks fully symmetric, outputs under G) of
Distillation_factories/distance_3/check_symmetry_search/orbit_bi.py (is_BI_orbit: sigma(T) = sum_omega sigma_omega
N(omega, T)) and orbit_search.py (slots = (output part, check orbit)); the code here is written for this repo.

    python searches/orbit_group.py --selftest
    python searches/orbit_group.py --k 3 4 5 6 --rmax 4 --out outputs/checks/orbit_group_l3.json
"""
import os, sys, json, time, argparse, pickle, random
from itertools import combinations, permutations
from math import comb

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us  # noqa: E402
import cascade_class as CC  # noqa: E402
from target_classes import coords  # noqa: E402
from enumerate_targets import layout  # noqa: E402

L = 3
BITS_MAX = 24
INF = 1 << 30


# ----------------------------------------------------------------- groups on the k outputs
def closure(gens, k):
    e = tuple(range(k))
    G, fr = {e}, [e]
    while fr:
        g = fr.pop()
        for h in gens:
            gh = tuple(h[g[i]] for i in range(k))
            if gh not in G:
                G.add(gh); fr.append(gh)
    return sorted(G)


def partitions(k, m=None):
    m = k if m is None else m
    if k == 0:
        yield ()
        return
    for p in range(min(k, m), 0, -1):
        for rest in partitions(k - p, p):
            yield (p,) + rest


def young(parts):
    """S_{p1} x S_{p2} x ... on consecutive blocks of outputs."""
    k, gens, s = sum(parts), [], 0
    for p in parts:
        if p >= 2:
            t = list(range(k)); t[s], t[s + 1] = s + 1, s; gens.append(tuple(t))
            c = list(range(k))
            for i in range(p):
                c[s + i] = s + (i + 1) % p
            gens.append(tuple(c))
        s += p
    return closure(gens, k)


def cyclic(k):
    return closure([tuple((i + 1) % k for i in range(k))], k)


def dihedral(k):
    return closure([tuple((i + 1) % k for i in range(k)), tuple((-i) % k for i in range(k))], k)


def img(mask, g, k):
    return sum(1 << g[i] for i in range(k) if mask >> i & 1)


def subset_orbits(G, k):
    """G-orbits on all subsets of the k outputs (bitmasks), each sorted; orbit index of every mask."""
    of, orbs = [-1] * (1 << k), []
    for m in range(1 << k):
        if of[m] >= 0:
            continue
        o = sorted({img(m, g, k) for g in G})
        for a in o:
            of[a] = len(orbs)
        orbs.append(o)
    return orbs, of


def digits_of_spec(spec, k):
    """Q_k digits (support mask -> 1) of a spec 'CCZ012,CS34,T0' at l=3: digit of T = number of tokens on T, mod 2."""
    from coset_certify import parse_target
    P = parse_target(spec, L)
    d = {}
    for T, c in P.items():
        m = sum(1 << q for q in T)
        sh = len(T) - 1
        assert c % (1 << sh) == 0
        if (c >> sh) & 1:
            d[m] = 1
    return d


def aut_of(dig, k):
    """Output permutations fixing the digit vector (support masks with digit 1)."""
    S = set(dig)
    return [g for g in permutations(range(k)) if {img(m, g, k) for m in S} == S]


# ----------------------------------------------------------------- the GF(2) system on cells
def build(k, r, orbs):
    cells = [(oi, w) for oi in range(len(orbs)) for w in range(1, r + 1)]
    cost = [len(orbs[oi]) * comb(r, w) for oi, w in cells]
    # cnt[oi][B] = #{A in orbit oi : A >= B}
    cnt = [[sum(1 for A in o if A & B == B) for B in range(1 << k)] for o in orbs]
    H, D, Dorb = [], [], []
    for bi, o in enumerate(orbs):
        B = o[0]; i = bin(B).count("1")
        if i <= 2:
            for j in range(1, min(3 - i, r) + 1):
                row = 0
                for ci, (oi, w) in enumerate(cells):
                    if w >= j and (cnt[oi][B] * comb(r - j, w - j)) & 1:
                        row |= 1 << ci
                H.append(row)
        if 1 <= i <= 3:
            row = 0
            for ci, (oi, w) in enumerate(cells):
                if (cnt[oi][B] * comb(r, w)) & 1:
                    row |= 1 << ci
            D.append(row); Dorb.append(bi)
    return cells, cost, H, D, Dorb


def rowbasis(H, D):
    """Row basis of [H; D] with H rows first. Returns (basis rows, h = #basis rows from H, combo of every D row over
    the basis indices) so that, for x with H x = 0, D_t x = parity(state & combo_t), state bit b = basis_b . x."""
    piv, basis, combos = {}, [], []          # piv: pivot bit -> (vector with that pivot only, its combo)

    def reduce(row):
        v, cmb = row, 0
        changed = True
        while v and changed:
            changed = False
            for p, (pv, pc) in piv.items():
                if v & p:
                    v ^= pv; cmb ^= pc; changed = True
        return v, cmb

    def add2(row):
        v, cmb = reduce(row)
        if v:
            b = len(basis); basis.append(row); cmb ^= 1 << b
            p = v & -v
            # eliminate p from existing pivot vectors so every pivot bit appears in exactly one vector
            for q in list(piv):
                pv, pc = piv[q]
                if pv & p:
                    piv[q] = (pv ^ v, pc ^ cmb)
            piv[p] = (v, cmb)
            return None
        return cmb
    for row in H:
        add2(row)
    h = len(basis)
    for row in D:
        add2(row)
    for row in D:
        v, cmb = reduce(row)
        assert v == 0
        combos.append(cmb)
    return basis, h, combos


def kernel(H, nv):
    """GF(2) null space basis (bitmasks over nv variables) of the rows H."""
    rows, piv = [], []
    for row in H:
        for pv, pc in zip(rows, piv):
            if row >> pc & 1:
                row ^= pv
        if row:
            pc = (row & -row).bit_length() - 1
            for i in range(len(rows)):
                if rows[i] >> pc & 1:
                    rows[i] ^= row
            rows.append(row); piv.append(pc)
    free = [c for c in range(nv) if c not in piv]
    out = []
    for f in free:
        v = 1 << f
        for pv, pc in zip(rows, piv):
            if pv >> f & 1:
                v |= 1 << pc
        out.append(v)
    return out


def parity_arr(a):
    a = a.copy()
    for s in (32, 16, 8, 4, 2, 1):
        a ^= a >> s
    return a & 1


def dp_sweep(k, r, orbs, bits_max=BITS_MAX):
    """Exact min cost for every reachable G-invariant output. Returns (outputs, info) where outputs maps a digit
    tuple (over Dorb) to (cost, x as list of cell indices)."""
    cells, cost, H, D, Dorb = build(k, r, orbs)
    basis, h, combos = rowbasis(H, D)
    bits = len(basis)
    info = dict(cells=len(cells), cond_rows=len(H), out_rows=len(D), h=h, d=bits - h, bits=bits)
    if bits > bits_max:
        return None, info
    cols = [sum(((basis[b] >> ci) & 1) << b for b in range(bits)) for ci in range(len(cells))]
    S = 1 << bits
    idx = np.arange(S, dtype=np.int64)
    dp = np.full(S, INF, dtype=np.int64); dp[0] = 0
    takes = []
    cheapest = {}                              # equal syndrome columns: only the cheapest can be in an optimum
    for ci in range(len(cells)):               # (two equal ones cancel at positive cost)
        if cols[ci] and (cols[ci] not in cheapest or cost[ci] < cost[cheapest[cols[ci]]]):
            cheapest[cols[ci]] = ci
    use = set(cheapest.values())
    info["cells_used"] = len(use)
    for ci in range(len(cells)):
        c = cols[ci]
        if c == 0 or ci not in use:
            takes.append(None); continue
        cand = dp[idx ^ c] + cost[ci]
        t = cand < dp
        dp = np.where(t, cand, dp)
        takes.append(np.packbits(t))
    finals = (np.arange(1 << (bits - h), dtype=np.int64) << h)
    ok = finals[dp[finals] < INF]
    out = {}
    for s in ok.tolist():
        dig = tuple(int(bin(s & cm).count("1") & 1) for cm in combos)
        x, ss = [], s
        for ci in range(len(cells) - 1, -1, -1):
            t = takes[ci]
            if t is not None and (t[ss >> 3] >> (7 - (ss & 7))) & 1:
                x.append(ci); ss ^= cols[ci]
        assert ss == 0
        n = sum(cost[ci] for ci in x)
        assert n == int(dp[s])
        out[dig] = (n, sorted(x))
    info["reachable_outputs"] = len(out)
    return (out, cells, cost, Dorb), info


def cpsat_target(k, r, orbs, of, dig, dlimit=20.0):
    """Min cost with deposit exactly `dig` (support mask -> 1, G-invariant). Returns (status, n, x, cells)."""
    from ortools.sat.python import cp_model
    cells, cost, H, D, Dorb = build(k, r, orbs)
    rhs = [dig.get(orbs[bi][0], 0) for bi in Dorb]
    m = cp_model.CpModel()
    x = [m.NewBoolVar(f"x{i}") for i in range(len(cells))]
    for row, b in [(hr, 0) for hr in H] + list(zip(D, rhs)):
        ids = [c for c in range(len(cells)) if row >> c & 1]
        if not ids:
            if b:
                return "INFEASIBLE", None, None, cells
            continue
        z = m.NewIntVar(0, len(ids), "")
        m.Add(sum(x[c] for c in ids) == 2 * z + b)
    m.Minimize(sum(cost[c] * x[c] for c in range(len(cells))))
    s = cp_model.CpSolver()
    s.parameters.num_workers = 1
    s.parameters.random_seed = 0
    s.parameters.max_deterministic_time = dlimit
    st = s.Solve(m)
    name = s.StatusName(st)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        xs = [c for c in range(len(cells)) if s.Value(x[c])]
        return name, sum(cost[c] for c in xs), xs, cells
    return name, None, None, cells


def gates_of(k, r, orbs, cells, x):
    g = []
    chk = range(k, k + r)
    for ci in x:
        oi, w = cells[ci]
        for A in orbs[oi]:
            As = frozenset(i for i in range(k) if A >> i & 1)
            for W in combinations(chk, w):
                g.append((As | frozenset(W), 1))
    return g


# ----------------------------------------------------------------- exact classes
def load_classes(k, fiber_cache=None):
    if k == 6 and fiber_cache and os.path.exists(fiber_cache):
        res, F = pickle.load(open(fiber_cache, "rb"))
        F.log = lambda *a: None
        return {str(c["orbit"]): c for c in res["classes"] if c["genuine"]}, (None, None, None, F)
    return CC.classes_of(L, k)


class Labels:
    def __init__(self, k, fiber_cache=None):
        self.k = k
        self.classes, self.ctx = load_classes(k, fiber_cache)
        self.C = coords(L, k)
        self.off, _ = layout(self.C)
        self.masks = [sum(1 << q for q in T) for T, _, _ in self.C]
        self.verify_label = CC.labeller(L, k, self.ctx, self.classes)

    def of_digits(self, dmask):
        """Class label (str) of the output with digit 1 exactly on the masks in dmask; None if not genuine."""
        C_, off_, lab, F = self.ctx
        dig = [1 if m in dmask else 0 for m in self.masks]
        if lab is not None:
            o = str(int(lab[sum(d << o for d, o in zip(dig, self.off))]))
        else:
            lb = F.label(np.array(dig, dtype=np.int64)); o = lb if isinstance(lb, str) else f"{lb[0]}:{lb[1]}"
        return o if o in self.classes else None


def statevector_check(k, r, gates):
    """Independent of the monomial algebra: build the 2^N state H^N then the rotations exp(i pi/8 c parity_S) (global
    phase dropped), undo every Newton monomial on a support touching a check (a diagonal Clifford when the factory
    condition holds; skipped and False returned otherwise), and test that the result is |psi_out> (x) |+>^r with
    |psi_out> = the deposited output. Returns (checks_returned, clifford_correction)."""
    N = k + r
    xs = np.arange(1 << N)
    # Newton coefficients a_T of the phase table (units of 2 pi / 8) on every T touching a check
    tab = [int(v) for v in us.phases(gates, N, L)]
    P = us.poly_from_table(tab, N, 1 << L)
    out = set(range(k)); cliff = True
    corr = np.zeros(1 << N)
    for T, aT in P.items():
        if T <= out:
            continue
        if aT % (1 << us.cliff_exp(len(T), L)):
            cliff = False
        m = sum(1 << q for q in T)
        corr -= aT * ((xs & m) == m)
    psi = np.exp(2j * np.pi * (np.array(tab) + corr) / (1 << L)) / np.sqrt(1 << N)
    M = psi.reshape(1 << r, 1 << k)                  # rows: check bits (high), columns: output bits
    ok = np.allclose(M, M[0][None, :])               # every check value carries the same output state => |+>^r
    return bool(ok), cliff


def verify(k, r, gates, orbit, n_expected, labels):
    f = dict(stage="orbit-group", k=k, r=r, gates=gates, params="")
    n = sum(us.inputs(c, L) for _, c in gates)
    res = dict(is_factory=us.is_factory(f, L), is_factory_bruteforce=us.is_factory(f, L, brute_force=True),
               distance=str(us.distance(f, L)), n_ok=n == n_expected, label=labels.verify_label(f))
    sv, cl = statevector_check(k, r, gates)
    res.update(statevector_checks_plus=sv, correction_clifford=cl)
    res["verified"] = bool(res["is_factory"] and res["is_factory_bruteforce"] and res["distance"] != "1"
                           and res["n_ok"] and res["label"] == orbit and sv and cl)
    return res


# ----------------------------------------------------------------- groups to sweep
def group_menu(k):
    out = []
    for p in partitions(k):
        out.append(("x".join(f"S{q}" for q in p), young(p)))
    if k >= 3:
        out.append((f"C{k}", cyclic(k)))
        out.append((f"D{k}", dihedral(k)))
    seen, menu = {}, []
    for name, G in out:
        orbs, of = subset_orbits(G, k)
        key = tuple(of)
        if key in seen:
            seen[key].append(name); continue
        seen[key] = [name]
        menu.append((name, G, orbs, of, seen[key]))
    return menu


# ----------------------------------------------------------------- self tests
def selftest(trials=1000, seed=7):
    """(1) the orbit-folded condition and output digits against is_factory (both forms) and the deposited polynomial
    on random cell subsets; (2) DP optimum against CP-SAT on the same cells and output."""
    from classify import phase_polynomial
    rnd = random.Random(seed); agree = tot = fac = 0; dig_ok = 0
    for _ in range(trials):
        k = rnd.randint(1, 4); r = rnd.randint(1, 3)
        menu = group_menu(k); name, G, orbs, of, _ = rnd.choice(menu)
        cells, cost, H, D, Dorb = build(k, r, orbs)
        if rnd.random() < 0.5:
            x = [ci for ci in range(len(cells)) if rnd.random() < 0.3]
        else:                                            # a random element of the kernel of H (a factory)
            ker = kernel(H, len(cells))
            xm = 0
            for v in ker:
                if rnd.random() < 0.5:
                    xm ^= v
            x = [ci for ci in range(len(cells)) if xm >> ci & 1]
        xm = sum(1 << ci for ci in x)
        cond = all(bin(row & xm).count("1") % 2 == 0 for row in H)
        g = gates_of(k, r, orbs, cells, x)
        if not g:
            continue
        f = dict(k=k, r=r, gates=g, stage="t", params="")
        b1 = us.is_factory(f, L); b2 = us.is_factory(f, L, brute_force=True)
        tot += 1; agree += (cond == b1 == b2); fac += b1
        P = phase_polynomial([S for S, _ in g], [c for _, c in g], set(range(k)), L)
        want = {sum(1 << q for q in T): 1 for T, c in P.items() if (c >> (len(T) - 1)) & 1}
        got = {}
        for bi, row in zip(Dorb, D):
            if bin(row & xm).count("1") & 1:
                for m in orbs[bi]:
                    got[m] = 1
        dig_ok += want == got
    print(f"condition: {agree}/{tot} random orbit circuits agree with is_factory (monomial and brute force); "
          f"{fac} are factories; output digits agree on {dig_ok}/{tot}", flush=True)
    # every DP optimum at k=2,3 must be a factory (the random subsets above are rarely factories)
    nfac = nver = 0
    for k in (2, 3):
        labels = Labels(k)
        for name, G, orbs, of, _ in group_menu(k):
            for r in (1, 2, 3):
                res, info = dp_sweep(k, r, orbs)
                out, cells, cost, Dorb = res
                for dig, (n, x) in out.items():
                    g = gates_of(k, r, orbs, cells, x)
                    if not g:
                        continue
                    f = dict(k=k, r=r, gates=g, stage="t", params="")
                    nfac += 1; nver += us.is_factory(f, L) and us.is_factory(f, L, brute_force=True)
    print(f"DP solutions k=2,3: {nver}/{nfac} pass is_factory (monomial and brute force)", flush=True)
    json.dump(dict(trials=trials, seed=seed, circuits=tot, agree=agree, factories=fac, digits_agree=dig_ok,
                   dp_solutions=nfac, dp_verified=nver),
              open(os.path.join(ROOT, "outputs", "checks", "orbit_group_selftest.json"), "w"), indent=1)
    return agree == tot and dig_ok == tot and nver == nfac


# ----------------------------------------------------------------- main run
def run_k(k, rmax, labels, log, bits_max=BITS_MAX, cpsat_check=3):
    """Sweep mode for every group in the menu and r = 1..rmax. Returns per-class best and per-group stats."""
    best, groups = {}, []
    for name, G, orbs, of, aliases in group_menu(k):
        for r in range(1, rmax + 1):
            t = time.time()
            res, info = dp_sweep(k, r, orbs, bits_max)
            rec = dict(k=k, group=name, aliases=aliases, order=len(G), r=r, N=k + r, **info)
            if res is None:
                rec.update(status=f"dropped: syndrome space 2^{info['bits']} > 2^{bits_max}", seconds=round(time.time() - t, 2))
                groups.append(rec); log(f"k={k} {name} r={r}: {rec['status']}"); continue
            out, cells, cost, Dorb = res
            per = {}
            for dig, (n, x) in out.items():
                dmask = set()
                for bi, dbit in zip(Dorb, dig):
                    if dbit:
                        dmask.update(orbs[bi])
                o = labels.of_digits(dmask)
                if o is None:
                    continue
                if o not in per or n < per[o][0]:
                    per[o] = (n, x, dig)
            nver = 0; checked = 0
            for o, (n, x, dig) in per.items():
                cur = best.get(o)
                if cur is not None and (cur["n"], cur["N"]) <= (n, k + r):
                    continue
                g = gates_of(k, r, orbs, cells, x)
                v = verify(k, r, g, o, n, labels)
                if not v["verified"]:
                    log(f"  VERIFY FAILED k={k} {name} r={r} class {o}: {v}")
                    continue
                best[o] = dict(n=n, N=k + r, group=name, order=len(G), mode="sweep", verified=True,
                               checks=v, gates=[[sorted(S), c] for S, c in g])
            # second method: CP-SAT on the same cells must reproduce the DP optimum for a few outputs
            for o, (n, x, dig) in list(per.items())[:cpsat_check]:
                dmask = {}
                for bi, dbit in zip(Dorb, dig):
                    if dbit:
                        for m in orbs[bi]:
                            dmask[m] = 1
                st, n2, _, _ = cpsat_target(k, r, orbs, of, dmask)
                checked += 1; nver += (st == "OPTIMAL" and n2 == n)
            rec.update(status="ok", classes_reached=len(per), cpsat_agrees=f"{nver}/{checked}",
                       seconds=round(time.time() - t, 2))
            groups.append(rec)
            log(f"k={k} {name} r={r}: bits {info['bits']}, {len(per)} classes, CP-SAT agrees {nver}/{checked}, "
                f"{rec['seconds']}s")
    return best, groups


def run_targets(k, rmax, labels, log, dlimit=20.0, trivial=True, rmin=1, shard=(0, 1)):
    """Target mode: G = Aut(rep) and (if trivial) G = {e}, deposit exactly the spec representative."""
    best, recs = {}, []
    triv = [tuple(range(k))]
    for i, (o, c) in enumerate(sorted(labels.classes.items())):
        if i % shard[1] != shard[0]:
            continue
        dig = digits_of_spec(c["spec"], k)
        A = aut_of(dig, k)
        menu = [("Aut", A)] + ([("trivial", triv)] if trivial and len(A) > 1 else [])
        for gname, G in menu:
            orbs, of = subset_orbits(G, k)
            for r in range(rmin, rmax + 1):
                t = time.time()
                st, n, x, cells = cpsat_target(k, r, orbs, of, dig, dlimit)
                rec = dict(k=k, orbit=o, spec=c["spec"], group=gname, order=len(G), r=r, N=k + r, status=st, n=n,
                           seconds=round(time.time() - t, 2))
                if n is not None:
                    g = gates_of(k, r, orbs, cells, x)
                    v = verify(k, r, g, o, n, labels)
                    rec["verified"] = v["verified"]
                    if v["verified"]:
                        cur = best.get(o)
                        if cur is None or (n, k + r) < (cur["n"], cur["N"]):
                            best[o] = dict(n=n, N=k + r, group=f"{gname}(order {len(G)})", order=len(G),
                                           mode="target", status=st, verified=True, checks=v,
                                           gates=[[sorted(S), cc] for S, cc in g])
                    else:
                        log(f"  VERIFY FAILED target k={k} {c['spec']} {gname} r={r}: {v}")
                recs.append(rec)
        b = best.get(o)
        log(f"k={k} target {c['spec']}: |Aut|={len(A)} best {b and (b['n'], b['N'], b['group'])}")
    return best, recs


# ----------------------------------------------------------------- comparison with the current best per class
def stage_first(O):
    """First cascade stage reaching each class's best n: the logic of searches/stage_contributions.py, recomputed
    here (that script writes outputs/stage_contributions_l3.json on import)."""
    import glob
    C = json.load(open(os.path.join(O, "class_table_l3.json")))
    key_of = {(r["k"], str(r["orbit"])): (r["k"], r["spec"]) for r in C}
    ORDER = ["two-group", "gluing", "entangled-gluing", "restricted-SAT", "unrestricted-SAT"]
    n = {s: {} for s in ORDER}

    def put(s, key, v):
        if v and (key not in n[s] or v < n[s][key]):
            n[s][key] = v
    for k in range(1, 7):
        for r in json.load(open(os.path.join(O, "cascade", f"l3_k{k}.json")))["rows"]:
            key = key_of.get((k, str(r["orbit"])))
            if key and r.get("genuine"):
                put({"two-group": "two-group", "gluing": "gluing", "quick-sat": "unrestricted-SAT"}[r["stage"]], key, r["n"])
    for e in json.load(open(os.path.join(O, "entangled_gluing_l3.json"))):
        put("entangled-gluing", (e["k"], e["spec"]), e["eg_n"])
    for f in glob.glob(os.path.join(O, "pattern_sat", "*.json")) + glob.glob(os.path.join(O, "pattern_sat_local", "*.json")) \
            + glob.glob(os.path.join(O, "synthesis_sat", "*.json")):
        d = json.load(open(f))
        if d.get("best"):
            put("restricted-SAT", (d["k"], d["spec"]), d["best"]["n"])
    keys = [(r["k"], r["spec"]) for r in C]
    best = {key: min(n[s][key] for s in ORDER if key in n[s]) for key in keys}
    first = {key: next(s for s in ORDER if n[s].get(key) == best[key]) for key in keys}
    return first, best, n, ORDER


def report(run_files, out_json, out_md, fiber_cache=None):
    O = os.path.join(ROOT, "outputs")
    first, sbest, n_stage, ORDER = stage_first(O)
    cur = {(r["k"], r["spec"]): r for r in json.load(open(os.path.join(O, "l3_best_per_class.json")))}
    ours, meta = {}, {}
    for fpath in run_files:
        d = json.load(open(fpath))
        for k, v in d["per_k"].items():
            k = int(k)
            meta.setdefault(k, []).append(dict(file=os.path.relpath(fpath, ROOT), args=d.get("args"),
                                               seconds=v["seconds"], groups=v["groups"],
                                               target_runs=len(v["target_runs"]),
                                               target_status={st: sum(r["status"] == st for r in v["target_runs"])
                                                              for st in sorted({r["status"] for r in v["target_runs"]})}))
            for src in ("sweep_best", "target_best"):
                for o, b in v[src].items():
                    key = (k, v["classes"][o])
                    if key not in ours or (b["n"], b["N"]) < (ours[key]["n"], ours[key]["N"]):
                        ours[key] = dict(b, orbit=o)
    # re-verify every factory kept (both is_factory forms, distance, statevector, exact label), whatever run made it
    labs = {}
    for (k, spec), b in ours.items():
        if k not in labs:
            labs[k] = Labels(k, fiber_cache)
        g = [(frozenset(S), c) for S, c in b["gates"]]
        b["checks"] = verify(k, b["N"] - k, g, b["orbit"], b["n"], labs[k])
        b["verified"] = b["checks"]["verified"]
    recs, cnt = [], dict(classes=len(cur), reached=0, equal=0, smaller=0, larger=0)
    by_stage = {s: dict(classes=0, reached=0, equal=0, smaller=0, larger=0, examples_equal=[], examples_smaller=[])
                for s in ORDER}
    tg = set(n_stage["two-group"])
    for key, c in sorted(cur.items()):
        b = ours.get(key)
        st = first[key]
        by_stage[st]["classes"] += 1
        rec = dict(k=key[0], spec=key[1], best_n=c["n"], best_N=c["N"], best_source=c["source"], best_first_stage=st,
                   stage_contrib_best_n=sbest[key], tcount=c["tcount"])
        if b:
            rel = "smaller" if b["n"] < c["n"] else "equal" if b["n"] == c["n"] else "larger"
            rec.update(og_n=b["n"], og_N=b["N"], og_group=b["group"], og_mode=b["mode"], og_status=b.get("status", "exhaustive"),
                       verified=b["verified"], og_checks=b["checks"], relation=rel, og_gates=b["gates"])
            cnt["reached"] += 1; cnt[rel] += 1
            cnt["unverified"] = cnt.get("unverified", 0) + (not b["verified"])
            S = by_stage[st]; S["reached"] += 1; S[rel] += 1
            ex = dict(spec=key[1], k=key[0], n=b["n"], N=b["N"], G=b["group"], best_n=c["n"], verified=b["verified"])
            if rel in ("equal", "smaller") and len(S["examples_" + rel]) < 5:
                S["examples_" + rel].append(ex)
        recs.append(rec)
    for r in recs:
        if "og_n" in r:
            r["og_group_order"] = ours[(r["k"], r["spec"])]["order"]
    # the same counts restricted to factories whose group G is nontrivial (|G| > 1): the orbit-group part proper.
    # Factories with G = {e} (target mode, trivial Aut or the trivial-group run) only have the checks symmetric.
    cnt["nontrivial_G"] = {rel_: sum(1 for r in recs if r.get("og_group_order", 1) > 1 and r.get("relation") == rel_)
                           for rel_ in ("equal", "smaller", "larger")}
    cnt["nontrivial_G"]["reached"] = sum(cnt["nontrivial_G"].values())
    for st_, S in by_stage.items():
        S["nontrivial_G"] = {rel_: sum(1 for r in recs if r["best_first_stage"] == st_ and r.get("og_group_order", 1) > 1
                                       and r.get("relation") == rel_) for rel_ in ("equal", "smaller", "larger")}
    cnt["by_mode_status"] = {}
    for r in recs:
        if "og_n" in r:
            key_ = f"{r['og_mode']}:{r['og_status']}"
            cnt["by_mode_status"][key_] = cnt["by_mode_status"].get(key_, 0) + 1
    reached = {(r["k"], r["spec"]) for r in recs if "og_n" in r}
    cnt["reached_not_two_group"] = len(reached - tg)
    cnt["two_group_reached"] = len(tg)
    cnt["per_k"] = {k: dict(classes=sum(1 for x in cur if x[0] == k), reached=sum(1 for x in reached if x[0] == k),
                            equal=sum(1 for r in recs if r["k"] == k and r.get("relation") == "equal"),
                            smaller=sum(1 for r in recs if r["k"] == k and r.get("relation") == "smaller"))
                    for k in sorted({x[0] for x in cur})}
    json.dump(dict(summary=cnt, by_first_stage=by_stage, runs=meta, classes=recs), open(out_json, "w"), indent=1)
    if out_md:
        write_md(out_md, out_json, cnt, by_stage, recs, meta)
    return cnt, by_stage, recs, meta


METHOD_NOTE = """Cells (omega, w): omega a G-orbit of output subsets (empty set = check-only), w = 1..r check weight; each cell
all-in with coefficient 1 or out (at l=3 general: odd c = 1 + even Clifford). Factory condition folded onto (orbit of
T cap O, |T cap S|) rows over GF(2); deposited Q_k digit of T <= O = n_T mod 2 per output orbit. Two modes:
sweep (exhaustive over all 2^cells subsets by a DP on the GF(2) syndrome, exact min n for every G-invariant output;
groups: Young subgroups S_p1 x ... for every partition of k, C_k, D_k, r = 1..4) and target (deposit the class
representative of outputs/targets/l3_k*.json exactly; G = Aut(representative), and G = {e} for k <= 5; CP-SAT, one
worker, deterministic time limit; status OPTIMAL is exact within the space, FEASIBLE is an upper bound only).
Every factory reported is re-verified here: unified_search.is_factory (monomial and brute-force phase table),
unified_search.distance != 1, a 2^N statevector check that the checks return to |+> after the diagonal Clifford
correction, n recomputed with unified_search.inputs, and the exact class label of cascade_class.labeller equal to the
class it is credited to. Self-test (python searches/orbit_group.py --selftest): the folded condition and digits agree
with is_factory on random orbit circuits; the DP optimum is cross-checked against CP-SAT on 3 outputs per sweep cell.
Not covered: l=4 (signs), r > 4 (N > k+4), sweep cells dropped for size (listed above), Aut of other representatives
of a class (Aut depends on the representative), FEASIBLE-only target runs (not proved optimal in their space)."""


def write_md(path, out_json, cnt, by_stage, recs, meta):
    rel = lambda p: os.path.relpath(p, ROOT)
    L_ = ["# Orbit-group stage at l=3 (d=2): checks S_r, outputs under G <= S_k", "",
          f"Code: `searches/orbit_group.py`. Per-class records: `{rel(out_json)}`. Current best per class: "
          "`outputs/l3_best_per_class.json`; first stage reaching it: logic of `searches/stage_contributions.py` "
          "(reproduces 9/15/425/301/16).", "",
          "## Counts (N of M classes, k = 1..6)", "",
          f"- reached: {cnt['reached']} of {cnt['classes']} (unverified among them: {cnt.get('unverified', 0)})",
          f"- equal to current best n: {cnt['equal']} of {cnt['reached']} reached",
          f"- smaller than current best n: {cnt['smaller']} of {cnt['reached']} reached",
          f"- larger: {cnt['larger']} of {cnt['reached']} reached",
          f"- reached and not reached by two-group: {cnt['reached_not_two_group']} (two-group reaches {cnt['two_group_reached']})",
          f"- with a nontrivial output group (|G| > 1) only: reached {cnt['nontrivial_G']['reached']}, equal "
          f"{cnt['nontrivial_G']['equal']}, smaller {cnt['nontrivial_G']['smaller']}, larger {cnt['nontrivial_G']['larger']}; "
          "the rest have G = {e} (checks symmetric, outputs unstructured)",
          f"- best factory per class by mode and CP-SAT status: {cnt['by_mode_status']} (FEASIBLE = upper bound within "
          "the space, not proved minimal there; sweep = exhaustive)",
          "", "| k | classes | reached | equal | smaller |", "|---|---|---|---|---|"]
    for k, v in cnt["per_k"].items():
        L_.append(f"| {k} | {v['classes']} | {v['reached']} | {v['equal']} | {v['smaller']} |")
    L_ += ["", "## By the stage that currently gives the best n first", "",
           "| first stage | classes | reached | equal n | smaller n | larger n | equal / smaller with |G|>1 |",
           "|---|---|---|---|---|---|---|"]
    for s_, v in by_stage.items():
        L_.append(f"| {s_} | {v['classes']} | {v['reached']} | {v['equal']} | {v['smaller']} | {v['larger']} | "
                  f"{v['nontrivial_G']['equal']} / {v['nontrivial_G']['smaller']} |")
    for s_ in ("entangled-gluing", "restricted-SAT"):
        v = by_stage[s_]
        for kind in ("smaller", "equal"):
            if v["examples_" + kind]:
                L_ += ["", f"{s_}-first classes, orbit-group n {kind} (up to 5):", ""]
                for e in v["examples_" + kind]:
                    L_.append(f"- k={e['k']} `{e['spec']}`: n={e['n']}, N={e['N']}, G={e['G']}, current best n={e['best_n']}, "
                              f"verified={e['verified']}")
    L_ += ["", "## Runs", ""]
    for k, ms in sorted(meta.items()):
        for m in ms:
            gs = m["groups"]
            dropped = [f"{g['group']} r={g['r']}" for g in gs if g["status"].startswith("dropped")]
            line = f"- k={k}: `{m['file']}`, seconds {m['seconds']}"
            if gs:
                line += (f"; sweep (group, r) cells: {sum(g['status'] == 'ok' for g in gs)} exhaustive, {len(dropped)} "
                         f"dropped (syndrome space > 2^{BITS_MAX}): {', '.join(dropped) or 'none'}")
            if m["target_runs"]:
                line += f"; target CP-SAT queries: {m['target_runs']} ({m['target_status']})"
            L_.append(line)
    L_ += ["", "## Method and checks", "", METHOD_NOTE]
    sm = [r for r in recs if r.get("relation") == "smaller"]
    L_ += ["", "## Every factory smaller than the current best", ""]
    for r in sm:
        L_.append(f"- k={r['k']} `{r['spec']}`: n={r['og_n']}, N={r['og_N']}, G={r['og_group']} ({r['og_mode']}), "
                  f"current best n={r['best_n']} N={r['best_N']} ({r['best_source']}), T-count {r['tcount']}; checks "
                  + ", ".join(f"{a}={b}" for a, b in r["og_checks"].items()))
    if not sm:
        L_.append("- none")
    open(path, "w").write("\n".join(L_) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--k", type=int, nargs="+", default=[3, 4, 5, 6])
    ap.add_argument("--rmax", type=int, default=4)
    ap.add_argument("--bits-max", type=int, default=BITS_MAX)
    ap.add_argument("--targets", action="store_true", help="also run target mode (Aut(rep) and trivial G, CP-SAT)")
    ap.add_argument("--trivial-kmax", type=int, default=5, help="target mode with G = {e} only for k <= this")
    ap.add_argument("--no-sweep", action="store_true")
    ap.add_argument("--rmin", type=int, default=1, help="smallest r in target mode")
    ap.add_argument("--shard", default="0/1", help="target mode: classes i with i %% m == s, given as s/m")
    ap.add_argument("--dlimit", type=float, default=20.0, help="CP-SAT deterministic time limit per query")
    ap.add_argument("--fiber-cache", default=None, help="pickled (res, Fiber) of enumerate_targets_fiber(3, 6); without it the "
                                                        "k=6 class table is rebuilt (~14 min); the cache used was built by "
                                                        "pickle.dump(enumerate_targets_fiber(3, 6)) and its labels match "
                                                        "outputs/targets/l3_k6.json on all 679 classes")
    ap.add_argument("--out", required=False)
    ap.add_argument("--md", default=None, help="with --report: also write this markdown summary")
    ap.add_argument("--report", nargs="+", default=None, help="run JSON files to compare with the current best")
    a = ap.parse_args()
    if a.selftest:
        sys.exit(0 if selftest() else 1)
    if a.report:
        cnt, by_stage, recs, meta = report(a.report, a.out, a.md, a.fiber_cache)
        print(json.dumps(cnt, indent=1))
        for s_, v in by_stage.items():
            print(s_, {x: y for x, y in v.items() if not x.startswith("examples")})
        sys.exit(0)
    log = lambda *s: print(*s, flush=True)
    result = dict(l=L, args=vars(a), per_k={})
    for k in a.k:
        t0 = time.time()
        labels = Labels(k, a.fiber_cache)
        t_lab = time.time() - t0
        t = time.time()
        sb, groups = ({}, []) if a.no_sweep else run_k(k, a.rmax, labels, log, a.bits_max)
        t_sweep = time.time() - t
        tb, trecs = {}, []
        t = time.time()
        if a.targets:
            tb, trecs = run_targets(k, a.rmax, labels, log, a.dlimit, trivial=k <= a.trivial_kmax, rmin=a.rmin,
                                     shard=tuple(map(int, a.shard.split("/"))))
        t_tgt = time.time() - t
        result["per_k"][k] = dict(classes={o: c["spec"] for o, c in labels.classes.items()}, sweep_best=sb,
                                  target_best=tb, groups=groups, target_runs=trecs,
                                  seconds=dict(labels=round(t_lab, 1), sweep=round(t_sweep, 1), target=round(t_tgt, 1)))
        log(f"k={k}: sweep reaches {len(sb)}, target reaches {len(tb)} of {len(labels.classes)}; "
            f"labels {t_lab:.0f}s sweep {t_sweep:.0f}s target {t_tgt:.0f}s")
        if a.out:
            json.dump(result, open(a.out, "w"))
