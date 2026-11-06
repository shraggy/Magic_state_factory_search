"""k = 7, 8 (l = 3) factories of sweep B against the synthillation construction, per target, for Fig. 2.

WHAT. Every verified distance-2 factory at l=3, k=7,8 that sweep B holds, grouped by TARGET: the output phase
polynomial reduced modulo diagonal Cliffords (T coefficients mod 2, CS mod 4, CCZ mod 8; classify.phase_polynomial, as
searches/k78_large_ccz.py). Beyond k=6 there is no class list, so targets equivalent under GL(k,2) are separate keys
("targets", never "classes"). A target is the key canonicalised under qubit permutations (exact: the minimum over
every relabelling compatible with a permutation-invariant colouring of the qubits); each factory is relabelled into
that frame (re-verified there), so every frame of a target shares one best n. Raw-frame and canonical counts are both
reported. GL(k,2) equivalence is used only for CCZ layers (fingerprint classes, k78_orbits.py, which at k=7,8 number
the Burnside count of k78_class_count.py and so are the classes), as a separate per-class summary. Targets with
fewer than k essential outputs (some CNOT frame leaves an output with no T, CS or CCZ term: tested exactly over all
2^k - 1 directions) are dropped and counted. Per target the smallest n is kept (ties: smallest N).

SOURCES (factories). (1) larger sweep: every l=3, k=7,8 row of outputs/large_unified_factories.csv (two-group and
disjoint gluing), rebuilt with merge_large_unified.regenerate_circuit as k78_large_ccz.py does. (2) sweep B stage (i):
outputs/staged_B/l3_k{7,8}/c*.json (two-group, group-symmetric and check-group circuits; added because the user asked
for every sweep B factory, and the target-list pass only keeps its CCZ-layer class bests). (3) the cascade on the
CCZ-layer target list: outputs/k78_pipeline/all_isd/*.json (best, tcount_lo, tcount_hi per layer) and
outputs/k78_pipeline/full/l3_k{7,8}/{ii_dj,ii_eg,iii_f}.json (best). Every factory used is re-verified here:
us.is_factory (monomial form), us.distance == 2, n = sum of us.inputs; for target-list records also that the factory
deposits the record's layer (same reduced key).

T-COUNT BOUNDS (tau = least number of T rotations, Hadamard-free {CNOT, T} with clean ancillas, diagonal Cliffords free).
  hi = min of: the target list's tcount_hi (k78_eg synthesis, a coset word); the synthesis read off each factory (odd
  rotations' output parts, k78_tau_upper.one_check logic); the same after the CNOT map of the proof of Corollary 3,
  which clears r output parts, so weight <= n - r (k78_on_bound.rank_of_checks); Lee-Brickell ISD from the best word
  (k78_tau_upper.isd, p=3, fixed seeds and iterations); at k=7 the weight a complete Brouwer-Zimmermann search found.
  Every explicit word is checked to lie in the target's Reed-Muller coset (tcount7.coset logic, GF(2) rank test).
  lo = max of PROVEN lower bounds: isotropy floor 2(k-d)+1 for CCZ layers (Mazumder arXiv:2610.01024 Cor. 23,
  isotropy_floor.py); Thm 25 of the same paper, tau >= 2(k-d*) - r (d* = dim of the radical of the full moment tensor,
  r = rank of the quadratic moment matrix; any diagonal level-3 gate; same circuit model, Cor. 13 for ancillas), on the
  target and on its restrictions to every hyperplane; restriction monotonicity (Lemmas 10, 11 there: tau(U) >= tau(U|W)
  for every subspace W after a CNOT relabelling) with the exact k<=6 decoder (tcount.tcount_rm) on every hyperplane at
  k=7 and on the 28 coordinate 6-subsets at k=8; and at k=7 the exact T-count by tcount7.exists_below2 (two
  information sets, complete when its bound reaches the cap, deterministic), run whenever lo < hi <= CAPMAX = 21
  (a fixed rule on the cap, not a time limit: the search cost grows as C(64, cap/2)). No exact k=8 T-count is used.
  These bounds do not cover circuits with Hadamards or measurements (Mazumder: CCZ from 4 T with measurement).

SYNTHILLATION (Campbell-Howard arXiv:1606.01904, Thm 1 and Table I, text read: n = col(A) + 2 col(B) + Delta; Delta =
  8, 9, 9, 11 (cases 1-4, linear terms) or 0, 1, 1, 3 (cases 5-8, no linear terms) by the parities of (col A, col B);
  B empty and F = 4C: case 9 (col A even) 0, case 10 2, case 11 (col A odd) 1).
  syn_hi = smallest VALID construction built here: synthillation_strongest.build on A in {A, A+0}, B in {B*, B*+0,
  B(mu*+1)}, A = the lightest explicit synthesis (weight hi when one was found), B from det_parts.opt_B_fresh at the
  smallest feasible mu* (searched upward), each checked with synthillation_compare.coeffs_ok (borrowed identity).
  syn_lo = min_{a in {lo, lo+1}, b in {mu*, mu*+1}} a + 2b + Delta_min(lin, a mod 2, b mod 2), and lo + (lo mod 2) for
  CCZ layers (cases 9/11). Delta_min (delta_min()) is the least Delta of any matrix with the printed Table I head of
  that case (its check rows on A, B, B, and c c c c with linear terms) plus extra columns with zero data part: every
  monomial on check rows alone must have an even count (a necessary condition of the borrowed identity), checked
  exhaustively. It gives 8, 10, 9, 11 (linear terms) and 0, 2, 1, 3 (none). Two printed values are not attainable:
  case 2 prints Delta = 9 (its 5 trailing columns: 0 of 7^5 completions pass), case 6 prints Delta = 1 but its own
  matrix 'A B B 0 0' has 2 padding columns and Example IV.3 counts 18 = 6 + 2*5 + 2 (the v5 text notes errors fixed in
  case 6 and Example IV.3). A proven lower bound on every Thm 1 / Table I construction (not on Thm 4's general
  quasitransversal G matrices): col(A) >= tau >= lo (A is a gate synthesis of U), col(B) >= mu* (B must match the
  order-1 and order-2 moments of U, the constraints of opt_B, so mu* is exact).
VERDICT: fewer (n < syn_lo), same (n == syn_hi and lo == hi; same_as_every_table1 when also n == syn_lo),
  fewer_than_built (n < syn_hi only), more (n > syn_hi), open otherwise; on_bound = (lo == hi == n - 1).

HOW VERIFIED (known answers, run first; the main run stops on any failure): T-count code on the 142 layers of
outputs/k78_pipeline/exact7_all_isd.json (tau there) and on 20 l=3 classes of outputs/staged_best.json (exact tau);
lower bounds <= exact tau, the essential-output test and permutation invariance of the canonical form on all 766
classes of outputs/class_table_l3.json; the k>=7-only restriction bounds on NRESTRICT k=6 classes with T or CS terms,
each embedded at k=7 and at k=8 in a random CNOT frame (fixed seeds): every bound <= tau; synthillation: syn_hi with the
exact A (synthillation_compare.opt_A) equals outputs/synthillation_strongest_l3.json best_n, syn_lo <= best_n, on all
766 (and the number with syn_lo = best_n is reported). Every known answer is recomputed in the same run (no reuse).
    python searches/k78_frontier.py [--workers W]  -> outputs/checks/k78_frontier.json, k78_frontier.md,
                                                      k78_frontier_known_answers.json
"""
import os, sys, json, glob, csv, time, argparse, itertools, math, random
from collections import Counter, defaultdict
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
csv.field_size_limit(10 ** 9)
import numpy as np

MODC = {1: 2, 2: 4, 3: 8}; CANON = {1: 1, 2: 2, 3: 4}
NRESTRICT = 60                  # k=6 classes embedded at k=7 and at k=8 for the restriction-bound known answer
CAPMAX = 21                     # exists_below2 is run only for caps <= this (deterministic rule)
ISD_ITERS = 300                 # Lee-Brickell iterations per step (fixed), p = 3, seeds 0, 1, ...


# ------------------------------------------------------------------------------------------- targets and keys
def reduced(gates, k):
    from classify import phase_polynomial
    P = phase_polynomial([frozenset(S) for S, _ in gates], [c for _, c in gates], set(range(k)), 3)
    return {T: CANON[len(T)] for T, v in P.items() if len(T) <= 3 and v % MODC[len(T)]}


def keyof(P):
    return tuple(sorted((tuple(sorted(T)), c) for T, c in P.items()))


def poly(key):
    return {frozenset(T): c for T, c in key}


def table(P, k):
    x = np.arange(1 << k); f = np.zeros(1 << k, dtype=np.int64)
    for T, c in P.items():
        m = sum(1 << q for q in T); f += c * ((x & m) == m)
    return f % 8


def mobius(tab, m):
    a = tab.copy()
    for i in range(m):
        b = 1 << i
        idx = np.arange(1 << m); sel = idx[(idx & b) != 0]
        a[sel] = (a[sel] - a[sel ^ b]) % 8
    return a


def restrict(P, k, basis):
    """Polynomial of x -> f(sum z_i basis_i) on len(basis) qubits, reduced mod diagonal Cliffords."""
    m = len(basis); tab = table(P, k); z = np.arange(1 << m); x = np.zeros(1 << m, dtype=np.int64)
    for i, v in enumerate(basis):
        x ^= ((z >> i) & 1) * v
    a = mobius(tab[x], m)
    out = {}
    for zz in range(1, 1 << m):
        T = frozenset(i for i in range(m) if zz >> i & 1)
        if len(T) <= 3 and a[zz] % MODC[len(T)]:
            out[T] = CANON[len(T)]
    return out


def essential_full(P, k):
    """True iff no CNOT frame leaves an output idle: for every direction a != 0, change basis so that a is the last
    basis vector and test whether some genuine monomial contains the last variable (exact, all 2^k - 1 directions)."""
    tab = table(P, k); z = np.arange(1 << k); last = 1 << (k - 1)
    small = [zz for zz in range(1 << k) if zz & last and bin(zz).count("1") <= 3]
    for a in range(1, 1 << k):
        piv = (a & -a).bit_length() - 1
        cols = [1 << i for i in range(k) if i != piv] + [a]
        x = np.zeros(1 << k, dtype=np.int64)
        for i, v in enumerate(cols):
            x ^= ((z >> i) & 1) * v
        g = mobius(tab[x], k)
        if not any(g[zz] % MODC[bin(zz).count("1")] for zz in small):
            return False
    return True


def canon_perm_map(key, k):
    """Canonical form under qubit permutations: WL colouring (3 rounds, permutation invariant), then the minimum
    relabelled key over every bijection that sends colour classes, in colour order, to consecutive positions (exact).
    Returns (canonical key, pos) with pos[q] the new label of qubit q."""
    P = poly(key)
    col = [tuple(sum(1 for T in P if q in T and len(T) == s) for s in (1, 2, 3)) for q in range(k)]
    for _ in range(3):
        new = []
        for q in range(k):
            ms = sorted((len(T), tuple(sorted(col[p] for p in T if p != q))) for T in P if q in T)
            new.append((col[q], tuple(ms)))
        col = new
    order = sorted(set(col)); cells = [[q for q in range(k) if col[q] == c] for c in order]
    best = None; bpos = None
    for parts in itertools.product(*[itertools.permutations(c) for c in cells]):
        seq = [q for p in parts for q in p]; pos = {q: i for i, q in enumerate(seq)}
        kk = tuple(sorted((tuple(sorted(pos[q] for q in T)), c) for T, c in key))
        if best is None or kk < best:
            best, bpos = kk, pos
    return best, bpos


def canon_perm(key, k):
    return canon_perm_map(key, k)[0]


def spec(key):
    T1 = [f"T{T[0]}" for T, c in key if len(T) == 1]
    T2 = [f"CS{''.join(map(str, T))}" for T, c in key if len(T) == 2]
    T3 = [f"CCZ{''.join(map(str, T))}" for T, c in key if len(T) == 3]
    return ",".join(T1 + T2 + T3)


# ------------------------------------------------------------------------------------------- coset and bounds
_SNF = {}


def coset(P, k):
    """tcount7.coset with the Smith form cached per k (same arithmetic)."""
    import tcount as TC
    if k not in _SNF:
        M = TC.parity_matrix(k); _SNF[k] = (M,) + tuple(TC.snf_z8(M))
    M, U, W, D = _SNF[k]
    Pr = {frozenset(q + 1 for q in T): c for T, c in P.items()}
    fv = TC._fvec(Pr, k); a0 = TC._solve_z8(U, W, D, fv)
    if TC._matvec_mod(M, a0) != [v % 8 for v in fv]:
        raise RuntimeError("not a level-3 diagonal gate")
    n = (1 << k) - 1
    basis = [sum(1 << r for r in range(n) if W[r][i] & 1) for i in range(n) if D[i] % 8 == 0]
    return sum(1 << r for r in range(n) if a0[r] & 1), basis, n


class Coset:
    def __init__(self, P, k):
        self.b, self.C0, self.nc = coset(P, k); self.k = k
        self.ech = []
        for v in self.C0:
            for e in self.ech:
                v = min(v, v ^ e)
            if v:
                self.ech.append(v)
        self.ech.sort(reverse=True)

    def contains(self, w):
        v = w ^ self.b
        for e in self.ech:
            v = min(v, v ^ e)
        return v == 0


def moments(word, k):
    ys = [y for y in range(1, 1 << k) if word >> (y - 1) & 1]
    def M(A):
        m = sum(1 << q for q in A)
        return sum(1 for y in ys if y & m == m) & 1
    return M


def rank2(rows):
    B = []
    for v in rows:
        for b in B:
            v = min(v, v ^ b)
        if v:
            B.append(v)
    return len(B)


def thm25(word, k):
    """Mazumder Thm 25: tau >= 2(k - d*) - r, from the moments of any residue (coset word)."""
    M = moments(word, k)
    Q = [sum((M({a}) if a == b else M({a, b})) << b for b in range(k)) for a in range(k)]
    r = rank2(Q)
    rows = [sum(M({a, p, q}) << a for a in range(k)) for p in range(k) for q in range(p, k)]
    dstar = k - rank2(rows)
    return 2 * (k - dstar) - r


def thm25P(P, k):
    """Thm 25 from the reduced polynomial: by Fact 15 the moments M_A (|A| <= 3) of every residue are the reduced
    coefficients (T: M_a, CS: M_ab, CCZ: M_abc), since a parity with even coefficient is a diagonal Clifford."""
    M = lambda A: 1 if frozenset(A) in P else 0
    Q = [sum((M({a}) if a == b else M({a, b})) << b for b in range(k)) for a in range(k)]
    rows = [sum(M({a, p, q}) << a for a in range(k)) for p in range(k) for q in range(p, k)]
    return 2 * rank2(rows) - rank2(Q)


def pure_cubic(P):
    return bool(P) and all(len(T) == 3 for T in P)


def lower_bounds(P, k):
    """[(value, source)] of proved lower bounds on tau(U)."""
    from isotropy_floor import floor
    import tcount as TC
    out = []
    C = Coset(P, k)
    t25 = thm25(C.b, k)
    assert t25 == thm25P(P, k)                             # moments from a coset word equal the reduced coefficients
    out.append((t25, "Mazumder Thm 25"))
    if pure_cubic(P):
        out.append((floor([sorted(T) for T in P], k), "isotropy floor (Mazumder Cor. 23)"))
    if k >= 7:
        hyps = []
        for h in range(1, 1 << k):
            piv = (h & -h).bit_length() - 1
            hyps.append([(1 << i) | ((h >> i & 1) << piv) for i in range(k) if i != piv])
        best25 = max(thm25P(R, k - 1) if R else 0 for R in (restrict(P, k, B) for B in hyps))
        out.append((best25, "Thm 25 on hyperplane restrictions (Lemmas 10, 11)"))
        subs = hyps if k == 7 else [[1 << i for i in S] for S in itertools.combinations(range(k), 6)]
        bestr = 0
        for B in subs:
            R = restrict(P, k, B)
            if R:
                Rr = {frozenset(q + 1 for q in T): c for T, c in R.items()}
                bestr = max(bestr, TC.tcount_rm(Rr, 6))
        out.append((bestr, "exact tau of restrictions to " + ("every hyperplane" if k == 7 else "28 coordinate 6-subsets")
                    + " (Lemmas 10, 11)"))
    return out


def readoffs(gates, k):
    """(plain, cleared) synthesis words of a factory: odd rotations' output parts; and after the CNOT map of the proof
    of Corollary 3, which clears the output part of r rotations with independent check parts."""
    odd = [(sum(1 << q for q in S if q < k), sum(1 << (q - k) for q in S if q >= k)) for S, c in gates if c % 2]
    plain = 0
    for a, b in odd:
        if a:
            plain ^= 1 << (a - 1)
    basis = []                                             # (reduced check vector, its K-image, pivot)
    for a, b in odd:
        acc = 0; v = b
        for bb, ka, p in basis:
            if v >> p & 1:
                v ^= bb; acc ^= ka
        if v:
            basis.append((v, a ^ acc, v.bit_length() - 1))
            basis.sort(key=lambda t: -t[2])
    cleared = 0
    for a, b in odd:
        acc = 0; v = b
        for bb, ka, p in basis:
            if v >> p & 1:
                v ^= bb; acc ^= ka
        assert v == 0
        a2 = a ^ acc
        if a2:
            cleared ^= 1 << (a2 - 1)
    return plain, cleared, len(basis)


def bz_exact(C, cap, lo):
    """tcount7.exists_below2 loop (one thread): returns (tau or None, hi seen, log)."""
    import tcount7 as T
    log = []; hi = cap
    while lo < hi and hi <= CAPMAX:
        t = time.process_time()
        found, seen, proved = T.exists_below2(C.b, C.C0, C.nc, hi, 1)
        log.append(dict(cap=hi, found=found, seen=seen, proved=proved, cpu_s=round(time.process_time() - t, 1)))
        if found:
            hi = seen
        elif proved is not None and proved >= hi:
            return hi, hi, log
        else:
            return None, hi, log                           # incomplete decision: no claim
    return (hi if lo >= hi else None), hi, log


def tbounds(P, k, words, recorded_hi=None, exact=True, extra_lows=()):
    """lo, hi with sources, and the lightest explicit coset word. extra_lows: proved lower bounds on the same tau
    (e.g. computed on a qubit-permuted frame of the same target, tau being permutation invariant)."""
    from k78_tau_upper import isd
    C = Coset(P, k)
    lows = lower_bounds(P, k) + list(extra_lows); lo = max(v for v, _ in lows); lo_src = sorted({s for v, s in lows if v == lo})
    W = [(bin(w).count("1"), s, w) for w, s in words if C.contains(w)]
    bad_words = sum(1 for w, s in words if not C.contains(w))
    W.append((bin(C.b).count("1"), "coset base", C.b))
    wbest = min(W)
    his = [(wbest[0], wbest[1])]
    if recorded_hi is not None:
        his.append((recorded_hi, "target list tcount_hi"))
    word = wbest[2]; wt = wbest[0]
    isd_steps = []
    if lo < wt and (recorded_hi is None or wt > recorded_hi):
        goal = lo if recorded_hi is None else max(lo, recorded_hi)
        while wt > goal:
            w2, seen, its = isd(C.b, C.C0, C.nc, wt - 1, ISD_ITERS, 3, seed=len(isd_steps))
            ok = bool(w2 is not None and seen < wt and C.contains(w2))
            isd_steps.append([wt - 1, its, seen if ok else None])
            if not ok:
                break
            word, wt = w2, seen
        if isd_steps and wt < wbest[0]:
            his.append((wt, "ISD (k78_tau_upper.isd, p=3)"))
    hi = min(v for v, _ in his)
    bz = None
    if exact and k <= 7 and lo < hi:
        tau, hseen, log = bz_exact(C, hi, lo)
        bz = log
        if hseen < hi:
            his.append((hseen, "Brouwer-Zimmermann word (exists_below2)")); hi = hseen
            while wt > hi:                                  # recover an explicit word for synthillation
                w2, seen, its = isd(C.b, C.C0, C.nc, hi, 4 * ISD_ITERS, 3, seed=100 + len(isd_steps))
                isd_steps.append([hi, its, seen if (w2 is not None and C.contains(w2)) else None])
                if w2 is None or not C.contains(w2):
                    break
                word, wt = w2, seen
        if tau is not None:
            lo = tau; lo_src = ["exact (tcount7.exists_below2, complete)"]
    hi_src = [s for v, s in his if v == hi]
    return dict(lo=lo, lo_source=lo_src, lo_all=lows, hi=hi, hi_source=hi_src, word_weight=wt, word=word,
                bad_words=bad_words, isd_steps=isd_steps, bz=bz)


# ------------------------------------------------------------------------------------------- synthillation
def mu_lempel(P, k):
    """Smallest Lempel factor size (Campbell-Howard Sec. III B, Lempel 1975): mu* = rank Q, + 1 if Q != 0 has a zero
    diagonal, with Q_aa = T coefficient of a mod 2 and Q_ab = CS coefficient of ab / 2 mod 2 (the constraints of
    opt_B). Equals the first feasible mu of the upward opt_B search on 766 of 766 classes (known-answer check)."""
    Q = [sum(((P.get(frozenset([a]), 0) % 2) if a == b else ((P.get(frozenset([a, b]), 0) // 2) % 2)) << b
             for b in range(k)) for a in range(k)]
    r = rank2(Q)
    return 0 if r == 0 else (r + 1 if not any(Q[a] >> a & 1 for a in range(k)) else r)


def syn_hi(k, P, A):
    import synthillation_strongest as SS, det_parts as DP
    mu = mu_lempel(P, k)
    B0 = DP.opt_B_fresh(k, P, mu); B1 = DP.opt_B_fresh(k, P, mu + 1)
    assert B0 is not None                                  # Lempel: a factor with mu* columns exists
    lin = any(len(T) == 1 for T in P); c = sum(1 << i for i in range(k) if frozenset([i]) in P)
    best = None
    for an, Av in (("A", A), ("A+0", A + [0])):
        for bn, Bv in (("B", B0), ("B+0", B0 + [0]), ("B(mu+1)", B1)):
            if Bv is None:
                continue
            res = SS.build(k, P, Av, Bv, lin, c)
            if res and (best is None or (res[0], res[1]) < (best[0], best[1])):
                best = (res[0], res[1], f"{an},{bn},case {res[2]}", len(Av))
    return mu, lin, best


def delta_min(lin, pa, pb):
    """Least Delta of any matrix that keeps the Table I head of case (lin, col A mod 2, col B mod 2) and appends m
    columns with zero data part: the head's check rows on blocks A, B, B (and c c c c) as printed in
    synthillation_compare.T5/C14, then the fewest extra check columns such that every monomial on check rows only
    (|S| <= 3 = all of them, R <= 3) has an even count, a necessary condition of the borrowed identity (coeffs_ok).
    Exhaustive over all multisets of extra check patterns of size < 2^R. Returns Delta = (4 if lin) + m."""
    import synthillation_compare as sc
    if lin:
        case = {(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}[(pa, pb)]
        tm = [r[:7] for r in sc.C14[case]]; sizes = [pa, pb, pb, 1, 1, 1, 1]; fixed = 4
    else:
        tm = [r[:3] for r in sc.T5[(pa, pb)]]; sizes = [pa, pb, pb]; fixed = 0
    R = len(tm)
    cnt = {S: sum(sz for bi, sz in enumerate(sizes) if all(tm[r][bi] for r in range(R) if S >> r & 1))
           for S in range(1, 1 << R)}
    for m in range(1 << R):
        for sub in itertools.combinations(range(1, 1 << R), m):
            if all((cnt[S] + sum(1 for q in sub if q & S == S)) % 2 == 0 for S in cnt):
                return fixed + m
    raise RuntimeError("no completion")


DELTA = {}


def delta_table():
    if not DELTA:
        for lin in (True, False):
            DELTA[lin] = {(pa, pb): delta_min(lin, pa, pb) for pa in (0, 1) for pb in (0, 1)}
    return DELTA


def syn_lo(lo, mu, lin, pure):
    """min over col(A) >= lo, col(B) >= mu* of col(A) + 2 col(B) + Delta_min(parities); Delta_min from delta_min()
    (= 8, 10, 9, 11 with linear terms, 0, 2, 1, 3 without; the printed Table I values except case 2, whose printed
    5 trailing columns admit no valid completion, and case 6, whose printed Delta = 1 disagrees with its own matrix
    'A B B 0 0' and with Example IV.3, 18 = 6 + 2*5 + 2). Pure CCZ: lo + (lo mod 2) (cases 9/11; case 10 and cases
    5-8 with B are never smaller)."""
    if lo is None:
        return None
    if pure:
        return lo + lo % 2
    D = delta_table()[lin]
    return min(a + 2 * b + D[(a % 2, b % 2)] for a in (lo, lo + 1) for b in (mu, mu + 1))


def cols(word, nc):
    return [y for y in range(1, nc + 1) if word >> (y - 1) & 1]


# ------------------------------------------------------------------------------------------- factory collection
def verify(k, N, gates):
    import unified_search as us
    g = [(frozenset(S), c % 8) for S, c in gates]
    f = dict(k=k, r=N - k, gates=g, stage="x", params="")
    ok = bool(us.is_factory(f, 3) and us.distance(f, 3) == 2)
    return ok, sum(us.inputs(c, 3) for _, c in g), g


def job_factory(item):
    t0 = time.process_time()
    kind, k, N, n_rec, src, payload, layer = item
    if kind == "csv":
        from merge_large_unified import regenerate_circuit
        r = payload
        g = regenerate_circuit(3, k, N, r["params"], r["search_stage"], r.get("raw"))
        if not g:
            return dict(ok=False, why="not rebuilt", src=src, k=k, cpu=time.process_time() - t0)
        gates = [(sorted(S), c) for S, c in g]
    else:
        gates = payload
    ok, n, g = verify(k, N, gates)
    P = reduced(g, k); key = keyof(P)
    why = None
    if n_rec is not None and n != n_rec:
        ok, why = False, f"n {n} != recorded {n_rec}"
    if ok and layer is not None and key != keyof({frozenset(x): 4 for x in layer}):
        ok, why = False, "output differs from the record's layer"
    # the same factory with its outputs relabelled into the permutation-canonical frame (check qubits unchanged):
    # a qubit relabelling keeps n and N, and the relabelled circuit is re-verified and its output key checked
    ckey, pos = canon_perm_map(key, k)
    gc = [(frozenset(pos[q] if q < k else q for q in S), c) for S, c in g]
    okc, nc_, gc = verify(k, N, [(sorted(S), c) for S, c in gc])
    if ok and not (okc and nc_ == n and keyof(reduced(gc, k)) == ckey):
        ok, why = False, "relabelled factory failed"
    return dict(ok=ok, why=why, src=src, k=k, N=N, n=n, key=key, ckey=ckey,
                gates=[[sorted(S), c] for S, c in gc], gates_raw_frame=[[sorted(S), c] for S, c in g],
                cpu=time.process_time() - t0)


def collect_items():
    items = []
    for r in csv.DictReader(open(os.path.join(O, "large_unified_factories.csv"))):
        if r["l"] == "3" and r["k"] in ("7", "8"):
            items.append(("csv", int(r["k"]), int(r["N"]), int(r["n"]), f"larger sweep: {r['search_stage']} {r['params']}", r, None))
    for k in (7, 8):
        for f in sorted(glob.glob(os.path.join(O, "staged_B", f"l3_k{k}", "c*.json"))):
            for lab, b in json.load(open(f))["best"].items():
                items.append(("g", k, b["N"], b["n"], f"sweep B stage (i): {b.get('G')} {b.get('H')} ({os.path.basename(f)})", b["gates"], None))
    for f in sorted(glob.glob(os.path.join(O, "k78_pipeline", "all_isd", "*.json"))):
        for t in json.load(open(f))["results"]:
            b = t["best"]
            if b:
                items.append(("g", t["k"], b["N"], b["n"], f"target list ii_eg:{b['kind']} ({t['id']})", b["gates"], t["ccz"]))
    for k in (7, 8):
        for p in ("ii_dj", "ii_eg", "iii_f"):
            d = json.load(open(os.path.join(O, "k78_pipeline", "full", f"l3_k{k}", f"{p}.json")))
            for o, b in d["best"].items():
                items.append(("g", k, b["N"], b["n"], f"target list pass {p}: {b.get('stage')}", b["gates"], b.get("ccz")))
    return items


def recorded_target_list():
    """per reduced key: min recorded tcount_hi, and the layer ids."""
    rec = {}
    for f in glob.glob(os.path.join(O, "k78_pipeline", "all_isd", "*.json")):
        for t in json.load(open(f))["results"]:
            key = keyof({frozenset(x): 4 for x in t["ccz"]})
            hi = t["tcount_hi"]
            if key not in rec or hi < rec[key]["hi"]:
                rec[key] = dict(hi=hi, lo=t["tcount_lo"], id=t["id"], fp=t.get("fp"))
    return rec


# ------------------------------------------------------------------------------------------- per target
def job_target(arg):
    t0 = time.process_time()
    k, key, facs, rec, members = arg                      # key = permutation-canonical key; members = raw frame keys
    P = poly(key)
    if not essential_full(P, k):
        return dict(k=k, key=key, genuine=False, members=len(members), cpu=time.process_time() - t0)
    extra = []
    for mk in members:                                     # frame-dependent bounds (6-subsets, hyperplanes) per frame
        if mk != key:
            extra += [(v, s + " [other frame]") for v, s in lower_bounds(poly(mk), k)]
    facs = sorted(facs, key=lambda f: (f["n"], f["N"], f["order"]))
    best = facs[0]
    words = []
    for f in facs:
        g = [(frozenset(S), c) for S, c in f["gates"]]
        pl, cl, r = readoffs(g, k)
        words.append((pl, "read-off synthesis of a factory"))
        words.append((cl, "Corollary 3 cleared read-off (weight <= n - r)"))
    tb = tbounds(P, k, words, rec["hi"] if rec else None, extra_lows=extra)
    nc = (1 << k) - 1
    A = cols(tb["word"], nc)
    mu, lin, sb = syn_hi(k, P, A)
    pure = pure_cubic(P)
    slo = syn_lo(tb["lo"], mu, lin, pure)
    n = best["n"]; shi = sb[0] if sb else None
    lo, hi = tb["lo"], tb["hi"]
    if slo is not None and n < slo:
        verdict = "fewer"
    elif shi is not None and n == shi and lo == hi and tb["word_weight"] == hi:
        verdict = "same"
    elif shi is not None and n < shi:
        verdict = "fewer_than_built"
    elif shi is not None and n > shi:
        verdict = "more"
    else:
        verdict = "open"
    canon = key
    raw_n = {}
    for f in facs:
        raw_n[f["key"]] = min(raw_n.get(f["key"], 10 ** 9), f["n"])
    out = dict(k=k, key=[[list(T), c] for T, c in key], canonical_perm=[[list(T), c] for T, c in canon], spec=spec(canon),
               raw_frame_keys=len(members), raw_frame_best_n=sorted(raw_n.values()),
               ccz_layer=pure, genuine=True, n=n, N=best["N"], source=best["src"],
               all_sources=sorted({f["src"].split(" (")[0] for f in facs}), factories_for_target=len(facs),
               gates=best["gates"], lo=lo, lo_source=tb["lo_source"], lo_all=tb["lo_all"], hi=hi, hi_source=tb["hi_source"],
               synthesis_A_weight=tb["word_weight"], bad_readoff_words=tb["bad_words"], isd_steps=tb["isd_steps"], bz=tb["bz"],
               mu=mu, linear_terms=lin, syn_hi=shi, syn_hi_N=sb and sb[1], syn_hi_variant=sb and sb[2], syn_lo=slo,
               verdict=verdict, same_as_every_table1=bool(verdict == "same" and slo == n),
               on_bound=bool(lo == hi == n - 1), target_list=rec and dict(id=rec["id"], fp=rec["fp"]))
    if pure:
        from k78_orbits import fingerprint
        out["gl_fingerprint"] = fingerprint([list(T) for T, _ in key], k)
    out["cpu"] = time.process_time() - t0
    return out


# ------------------------------------------------------------------------------------------- known answers
def ka_exact7(item):
    t0 = time.process_time()
    r, words = item
    P = {frozenset(x): 4 for x in r["ccz"]}
    tb = tbounds(P, 7, words, r["synthesis"])
    return dict(id=r["id"], want=r["tau"], lo=tb["lo"], hi=tb["hi"], ok=tb["lo"] == tb["hi"] == r["tau"],
                bz=tb["bz"], cpu=time.process_time() - t0)


def ka_small(r):
    t0 = time.process_time()
    import class_table_l3 as CT, tcount as TC, tcount7 as T7
    k = r["k"]; g = [(frozenset(S), c) for S, c in r["gates"]]
    P = reduced(g, k)              # the factory's own output: a CNOT-frame image of the class spec, same exact tau
    pl, cl, rk = readoffs(g, k)
    C = Coset(P, k)
    lows = lower_bounds(P, k)
    exact = TC.tcount_rm({frozenset(q + 1 for q in T): c for T, c in P.items()}, k)
    f1 = T7.exists_below2(C.b, C.C0, C.nc, r["tau"], 1); f2 = T7.exists_below2(C.b, C.C0, C.nc, r["tau"] + 1, 1)
    tb = tbounds(P, k, [(pl, "read-off"), (cl, "cleared")], None)
    ok = (exact == r["tau"] and all(v <= r["tau"] for v, _ in lows) and C.contains(pl) and C.contains(cl)
          and bin(cl).count("1") <= r["n"] - rk and not f1[0] and f1[2] >= r["tau"] and f2[0]
          and tb["lo"] == tb["hi"] == r["tau"])
    return dict(spec=r["spec"], output=spec(keyof(P)), k=k, tau=r["tau"], exact_decoder=exact, lows=lows, readoff=bin(pl).count("1"),
                cleared=bin(cl).count("1"), bz_below_tau=f1, bz_below_tau1=f2, pipeline=[tb["lo"], tb["hi"]], ok=ok,
                cpu=time.process_time() - t0)


def ka_syn(r):
    t0 = time.process_time()
    import class_table_l3 as CT, synthillation_compare as sc
    k = r["k"]; P = CT.parse(r["spec"])
    A = sc.opt_A(k, P)
    mu, lin, sb = syn_hi(k, P, A)
    Pc = {T: CANON[len(T)] for T in P}
    perm = list(range(k)); random.Random(k * 1000 + len(P)).shuffle(perm)   # canonical form must not depend on labels
    lows = lower_bounds(Pc, k)
    tau = r["tcount"]
    slo = syn_lo(tau, mu, lin, pure_cubic(Pc))
    return dict(k=k, spec=r["spec"], tau=tau, A=len(A), mu=mu, mu_ok=mu == r["mu"], syn_hi=sb and sb[0], variant=sb and sb[2], syn_lo=slo,
                canon_ok=canon_perm(keyof(Pc), k) == canon_perm(keyof({frozenset(perm[q] for q in T): c for T, c in Pc.items()}), k),
                lows=[v for v, _ in lows], lows_ok=all(v <= tau for v, _ in lows), genuine=essential_full(Pc, k),
                cpu=time.process_time() - t0)


def ka_restrict(arg):
    """k=6 class embedded at k = 7 or 8 in a random CNOT frame (fixed seed): every lower bound, including the
    hyperplane / 6-subset restriction bounds that the k<=6 test never runs, must be <= the exact tau."""
    t0 = time.process_time()
    import class_table_l3 as CT
    r, k, seed = arg
    rnd = random.Random(seed); P = CT.parse(r["spec"])
    while True:
        B = [rnd.randrange(1, 1 << k) for _ in range(k)]
        if rank2(B) == k:
            break
    Q = restrict(P, k, B)
    lows = lower_bounds(Q, k)
    return dict(spec=r["spec"], k=k, seed=seed, tau=r["tcount"], lows=lows, ok=all(v <= r["tcount"] for v, _ in lows),
                restriction_bound_eq_tau=any(v == r["tcount"] for v, s in lows if "restriction" in s),
                cpu=time.process_time() - t0)


def workers_default():
    load = os.getloadavg()[0]
    return max(1, min(8, (os.cpu_count() or 16) - 5 - math.ceil(load)))


def _init():
    try:
        if os.nice(0) < 19:
            os.nice(19 - os.nice(0))
    except OSError:
        pass


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=0)
    ap.add_argument("--out", default=os.path.join(O, "checks", "k78_frontier.json"))
    a = ap.parse_args()
    W = a.workers or workers_default(); T0 = time.time(); cpu = defaultdict(float)
    print("workers", W, flush=True)
    res = dict(settings=dict(CAPMAX=CAPMAX, ISD_ITERS=ISD_ITERS, workers=W))
    with Pool(W, initializer=_init, maxtasksperchild=200) as pool:
        # ---- known answers: synthillation and bounds on the 766 classes (always run, by this code)
        CT = json.load(open(os.path.join(O, "class_table_l3.json")))
        SS = {(x["k"], x["spec"]): x for x in json.load(open(os.path.join(O, "synthillation_strongest_l3.json")))}
        out = list(pool.imap(ka_syn, sorted(CT, key=lambda r: -r["k"]), chunksize=1))
        cpu["known answers: synthillation, bounds (766)"] = sum(o["cpu"] for o in out)
        bn = lambda o: SS[(o["k"], o["spec"])]["best_n"]
        match = sum(o["syn_hi"] == bn(o) for o in out)
        lo_ok = sum(o["syn_lo"] <= bn(o) for o in out)
        res["known_answer_synthillation"] = dict(classes=len(out), syn_hi_equals_best_n=match, syn_lo_le_best_n=lo_ok,
                                                 syn_lo_equals_best_n=sum(o["syn_lo"] == bn(o) for o in out),
                                                 delta_min=dict((("linear" if l else "no linear") + f" A{a}B{b}", v)
                                                                for l, D in delta_table().items() for (a, b), v in D.items()),
                                                 lower_bounds_le_tau=sum(o["lows_ok"] for o in out),
                                                 essential_k=sum(o["genuine"] for o in out),
                                                 mu_lempel_equals_class_table=sum(o["mu_ok"] for o in out),
                                                 canonical_form_permutation_invariant=sum(o["canon_ok"] for o in out),
                                                 mismatches=[o for o in out if o["syn_hi"] != bn(o)][:20],
                                                 thm25_eq_tau=sum(o["lows"][0] == o["tau"] for o in out))
        print("KA synthillation:", {x: v for x, v in res["known_answer_synthillation"].items() if x != "mismatches"}, flush=True)
        # ---- known answers: T-count
        SB = [r for r in json.load(open(os.path.join(O, "staged_best.json")))["rows"] if r["l"] == 3]
        pick = SB[::38][:20]
        o2 = list(pool.imap(ka_small, pick, chunksize=1))
        cpu["known answers: 20 k<=6 classes"] = sum(o["cpu"] for o in o2)
        E = json.load(open(os.path.join(O, "k78_pipeline", "exact7_all_isd.json")))
        fac = {}
        for f in glob.glob(os.path.join(O, "k78_pipeline", "all_isd", "*.json")):
            for t in json.load(open(f))["results"]:
                fac[t["id"]] = t
        items = []
        for r in E:
            g = [(frozenset(S), c) for S, c in fac[r["id"]]["best"]["gates"]]
            pl, cl, _ = readoffs(g, 7)
            items.append((r, [(pl, "read-off"), (cl, "cleared")]))
        o3 = list(pool.imap_unordered(ka_exact7, items, chunksize=1))
        cpu["known answers: 142 exact7 layers"] = sum(o["cpu"] for o in o3)
        # restriction bounds (k >= 7 only): k=6 classes with T or CS terms, embedded at k=7 and k=8
        K6 = [r for r in CT if r["k"] == 6 and not all(x.startswith("CCZ") for x in r["spec"].split(","))]
        rk = [(r, k, 1000 * k + i) for k in (7, 8) for i, r in enumerate(K6[::max(1, len(K6) // NRESTRICT)][:NRESTRICT])]
        o4 = list(pool.imap_unordered(ka_restrict, rk, chunksize=1))
        cpu[f"known answers: restriction bounds ({len(o4)})"] = sum(o["cpu"] for o in o4)
        res["known_answer_tcount"] = dict(k_le_6=dict(classes=len(o2), ok=sum(o["ok"] for o in o2), rows=o2),
                                          exact7=dict(layers=len(o3), ok=sum(o["ok"] for o in o3), rows=o3,
                                                      bz_runs=sum(1 for o in o3 if o["bz"]), fails=[o for o in o3 if not o["ok"]]),
                                          restriction=dict(embeddings=len(o4), ok=sum(o["ok"] for o in o4),
                                                           by_k={str(k): [sum(o["ok"] for o in o4 if o["k"] == k), sum(1 for o in o4 if o["k"] == k)] for k in (7, 8)},
                                                           restriction_bound_reaches_tau=sum(o["restriction_bound_eq_tau"] for o in o4),
                                                           rows=o4))
        print("KA tcount: k<=6", sum(o["ok"] for o in o2), "of", len(o2), "; exact7", sum(o["ok"] for o in o3), "of", len(o3),
              "; restriction", sum(o["ok"] for o in o4), "of", len(o4), flush=True)
        json.dump(res, open(a.out[:-5] + "_known_answers.json", "w"), default=str)
        if (match != len(out) or lo_ok != len(out) or sum(o["ok"] for o in o2) != len(o2) or sum(o["ok"] for o in o3) != len(o3)
                or sum(o["ok"] for o in o4) != len(o4)):
            print("KNOWN-ANSWER CHECK FAILED: stopping before the main run", flush=True)
            return
        # ---- factories
        items = collect_items()
        prefixes = ["larger sweep", "sweep B stage (i)", "target list ii_eg", "target list pass"]
        F = list(pool.imap(job_factory, items, chunksize=4))
        cpu["factory verification"] = sum(f["cpu"] for f in F)
        fstats = Counter((f["k"], f["src"].split(":")[0], f["ok"]) for f in F)
        res["factories"] = dict(total=len(F), verified=sum(f["ok"] for f in F),
                                by_source={f"k={k} {s} {'ok' if ok else 'FAILED'}": c for (k, s, ok), c in sorted(fstats.items())},
                                failures=[dict(src=f["src"], why=f.get("why")) for f in F if not f["ok"]][:50])
        print("factories:", res["factories"]["verified"], "of", len(F), "verified", flush=True)
        groups = defaultdict(list)
        for f in F:
            if f["ok"]:
                f["order"] = [i for i, s in enumerate(prefixes) if f["src"].startswith(s)][0]
                groups[(f["k"], f["ckey"])].append(f)
        rec = recorded_target_list()
        args = []
        for (k, key), fs in groups.items():
            members = sorted({f["key"] for f in fs})
            rs = [rec[m] for m in members if m in rec]          # target-list records of any frame of this target
            args.append((k, key, fs, min(rs, key=lambda x: x["hi"]) if rs else None, members))
        args.sort(key=lambda x: (-(x[0] == 7 and not all(len(T) == 3 for T, _ in x[1])), x[0]))   # slow (BZ) first
        R = []
        for i, r in enumerate(pool.imap_unordered(job_target, args, chunksize=1)):
            R.append(r)
            if i % 200 == 0:
                print(i, "of", len(args), "targets", round(time.time() - T0), "s", flush=True)
        cpu["targets"] = sum(r["cpu"] for r in R)
    G = [r for r in R if r["genuine"]]
    tot = {}
    for k in (7, 8):
        Rk = [r for r in G if r["k"] == k]
        tot[str(k)] = dict(
            targets_raw_frame_keys=sum(r["raw_frame_keys"] for r in Rk),
            targets_up_to_permutation_before_drop=sum(1 for r in R if r["k"] == k),
            dropped_not_k_essential=sum(1 for r in R if r["k"] == k and not r["genuine"]),
            targets=len(Rk), targets_up_to_permutation=len(Rk),
            targets_with_frames_of_different_n=sum(1 for r in Rk if len(set(r["raw_frame_best_n"])) > 1),
            same_as_every_table1=sum(r["same_as_every_table1"] for r in Rk),
            same_built_only=sum(1 for r in Rk if r["verdict"] == "same" and not r["same_as_every_table1"]),
            ccz_layers=sum(r["ccz_layer"] for r in Rk),
            ccz_gl_fingerprints=len({r["gl_fingerprint"] for r in Rk if r["ccz_layer"]}),
            by_verdict=dict(Counter(r["verdict"] for r in Rk)),
            by_verdict_ccz=dict(Counter(r["verdict"] for r in Rk if r["ccz_layer"])),
            by_verdict_other=dict(Counter(r["verdict"] for r in Rk if not r["ccz_layer"])),
            on_bound=sum(r["on_bound"] for r in Rk), lo_eq_hi=sum(r["lo"] == r["hi"] for r in Rk),
            lo_none=sum(r["lo"] is None for r in Rk), syn_hi_none=sum(r["syn_hi"] is None for r in Rk),
            n_minus_syn_hi=dict(Counter(r["n"] - r["syn_hi"] for r in Rk if r["syn_hi"] is not None)),
            n_minus_syn_lo=dict(Counter(r["n"] - r["syn_lo"] for r in Rk if r["syn_lo"] is not None)),
            source_of_best=dict(Counter(r["source"].split(":")[0] for r in Rk)))
    for k in (7, 8):
        cl = defaultdict(list)
        for r in G:
            if r["k"] == k and r["ccz_layer"]:
                cl[r["gl_fingerprint"]].append(r)
        rows = []
        for fp, v in cl.items():                           # CNOTs on the outputs are free: the class best n serves every member
            b = min(v, key=lambda r: (r["n"], r["N"]))
            lo = max(r["lo"] for r in v); hi = min(r["hi"] for r in v); shi = min(r["syn_hi"] for r in v)
            slo = max(r["syn_lo"] for r in v)
            vd = ("fewer" if b["n"] < slo else "same" if (b["n"] == shi and lo == hi) else "fewer_than_built" if b["n"] < shi
                  else "more" if b["n"] > shi else "open")
            rows.append(dict(fingerprint=fp, targets=len(v), n=b["n"], N=b["N"], spec=b["spec"], lo=lo, hi=hi, syn_hi=shi,
                             syn_lo=slo, verdict=vd, on_bound=bool(lo == hi == b["n"] - 1)))
        tot[str(k)]["ccz_gl_classes"] = sorted(rows, key=lambda r: (r["n"], r["spec"]))
        tot[str(k)]["ccz_gl_classes_by_verdict"] = dict(Counter(r["verdict"] for r in rows))
        inv = defaultdict(list)                            # unproven GL invariant for the other targets
        for r in G:
            if r["k"] == k and not r["ccz_layer"]:
                inv[(r["mu"], r["linear_terms"], r["lo_all"][0][0], r["hi"], r["lo"])].append(r)
        tot[str(k)]["other_more_with_same_invariants_smaller_n"] = sum(
            1 for v in inv.values() for r in v if r["verdict"] == "more" and any(x["n"] < r["n"] for x in v))
    res["totals"] = tot
    res["cpu_s"] = {x: round(v, 1) for x, v in cpu.items()}; res["cpu_s_total"] = round(sum(cpu.values()), 1)
    res["wall_s"] = round(time.time() - T0, 1)
    res["targets"] = sorted(G, key=lambda r: (r["k"], r["n"], r["spec"]))
    res["dropped"] = [dict(k=r["k"], key=r["key"], raw_frame_keys=r["members"]) for r in R if not r["genuine"]]
    json.dump(res, open(a.out, "w"), default=str)
    write_md(res, a.out[:-5] + ".md")
    print(json.dumps(tot, indent=1), flush=True)
    print("CPU s", res["cpu_s"], "wall s", res["wall_s"], flush=True)


def write_md(res, path):
    L = ["# k = 7, 8 factories of sweep B against synthillation (searches/k78_frontier.py)", ""]
    ks = res["known_answer_synthillation"]; kt = res["known_answer_tcount"]; kr = kt["restriction"]
    L += ["## Known-answer checks (all run by this code in this run)",
          f"- synthillation, 766 classes (k<=6): syn_hi with exact A equals synthillation_strongest_l3.json best_n on "
          f"{ks['syn_hi_equals_best_n']} of {ks['classes']}; syn_lo <= best_n on {ks['syn_lo_le_best_n']} of {ks['classes']} "
          f"and syn_lo = best_n on {ks['syn_lo_equals_best_n']} of {ks['classes']}; every lower bound <= exact tau on "
          f"{ks['lower_bounds_le_tau']} of {ks['classes']}; k essential outputs on {ks['essential_k']} of {ks['classes']}; "
          f"permutation-canonical form unchanged by a random relabelling on {ks['canonical_form_permutation_invariant']} of {ks['classes']}.",
          f"- least Delta per Table I head (exhaustive, delta_min): {ks['delta_min']}.",
          f"- T-count code: {kt['k_le_6']['ok']} of {kt['k_le_6']['classes']} k<=6 classes of staged_best.json; "
          f"{kt['exact7']['ok']} of {kt['exact7']['layers']} layers of exact7_all_isd.json ({kt['exact7']['bz_runs']} needed the exact search).",
          f"- restriction lower bounds (run only at k >= 7): k=6 classes with T or CS terms in a random CNOT frame at k=7, 8: "
          f"every bound <= tau on {kr['ok']} of {kr['embeddings']} (k=7 {kr['by_k']['7'][0]} of {kr['by_k']['7'][1]}, "
          f"k=8 {kr['by_k']['8'][0]} of {kr['by_k']['8'][1]}); a restriction bound reaches tau on {kr['restriction_bound_reaches_tau']} "
          f"of {kr['embeddings']}. A test, not a proof (the proof is Lemmas 10, 11 of Mazumder arXiv:2610.01024).", ""]
    L += ["## Factories", f"{res['factories']['verified']} of {res['factories']['total']} re-verified (is_factory, distance 2, n), "
          "each also after relabelling its outputs into the permutation-canonical frame.", ""]
    L += ["## Totals per k (targets = reduced polynomials up to qubit permutation)", "",
          "| k | raw-frame keys | targets | dropped (< k essential) | kept | CCZ layers | fewer | same (every Table I) | "
          "same (built only) | fewer_than_built | more | open | on bound | lo = hi |", "|" + "---|" * 14]
    for k, t in res["totals"].items():
        v = t["by_verdict"]
        L.append(f"| {k} | {t['targets_raw_frame_keys']} | {t['targets_up_to_permutation_before_drop']} | {t['dropped_not_k_essential']} | "
                 f"{t['targets']} | {t['ccz_layers']} | {v.get('fewer', 0)} | {t['same_as_every_table1']} | {t['same_built_only']} | "
                 f"{v.get('fewer_than_built', 0)} | {v.get('more', 0)} | {v.get('open', 0)} | {t['on_bound']} | {t['lo_eq_hi']} |")
    L += ["", "Verdicts split: CCZ layers " + "; ".join(f"k={k} {t['by_verdict_ccz']}" for k, t in res["totals"].items()),
          "", "Other targets " + "; ".join(f"k={k} {t['by_verdict_other']}" for k, t in res["totals"].items()), "",
          "CCZ layers per GL(k,2) class (fingerprint classes; at k=7, 8 their number equals the Burnside count of "
          "outputs/checks/k78_class_count.json, 6 and 20, so they are the classes): " + "; ".join(
              f"k={k} {len(t['ccz_gl_classes'])} classes {t['ccz_gl_classes_by_verdict']}" for k, t in res["totals"].items()), "",
          "Other targets are not reduced under GL(k,2): 'more' targets sharing (mu*, linear terms, Thm 25, hi, lo) with a "
          "target of smaller n (possible CNOT-frame artefacts, not proved): " + "; ".join(
              f"k={k} {t['other_more_with_same_invariants_smaller_n']} of {t['by_verdict'].get('more', 0)}" for k, t in res["totals"].items()), "",
          f"CPU s: {res['cpu_s']} (total {res['cpu_s_total']}); wall {res['wall_s']} s.", "",
          "Verdicts. fewer: n < syn_lo, i.e. fewer inputs than every Campbell-Howard Thm 1 / Table I construction "
          "(n = col(A) + 2 col(B) + Delta, A a phase-polynomial synthesis of U, B a Lempel factor); this does not bound "
          "their general quasitransversal G-matrix protocols (Thm 4), a family that includes our own factories. "
          "same (every Table I): n = syn_hi = syn_lo with exact tau. same (built only): n = syn_hi with exact tau, but "
          "syn_lo < n. fewer_than_built: below the construction built here only (an upper-bound comparison). "
          "more: n > syn_hi. open: otherwise. T-count bounds hold for Hadamard-free {CNOT, T} circuits with clean ancillas."]
    open(path, "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
