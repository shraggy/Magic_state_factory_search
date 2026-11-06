"""Leading-order output error and first-order acceptance of the three measured distil-then-synthesize routes of
searches/route_compare_l3.py, by exact fault enumeration, against the catalogue factory, for the 766 l=3 classes, k<=6.

Noise model (the same as searches/error_coefficients.py): every input T state carries an independent Z fault with
probability p; a fault on the input of a parity rotation on support S applies Z_S after it. Clifford gates,
measurements and feedforward are perfect. Every sub-factory used here has distance 2, so a single fault is detected
(the run is rejected) and the output error is kappa p^2 + O(p^3), where kappa sums, over the pairs of faults inside one
sub-factory run that pass its checks, acceptance times (1 - fidelity of the final target state). Pairs in different
sub-factory runs are each caught by their own run (order p^2 rejection, never an accepted error at order p^2).
For the catalogue factory this kappa is the pair count c of error_coefficients.coefficients (every accepted pair with a
nonzero output part gives an orthogonal state), recomputed here for all 766.

Sub-factories and conversions, each simulated exactly:
  [[8,3,2]]       the CCZ012 factory of outputs/staged_best.json (4 qubits), check measured in X, accepted on +;
  CCZ gate        the |CCZ> state teleported onto the target (CNOT and Z measurement, Clifford correction);
  CS gate         [[8,3,2]] followed by a Y measurement of one output and a Clifford correction on each outcome
                  (Beverland et al. arXiv:1904.01124 Prop. A.4), then teleported;
  C2T             the catalysed |CCZ> -> 2|T> conversion of Gidney-Fowler arXiv:1812.01238, Fig. 10 (its Quirk circuit:
                  X^{-1/2} on q2, X_0 X_1 controlled on |0> of q2, Z_0 Z_1 controlled on |-> of q2, Z^{-1/4} on q2,
                  the same controlled Z_0 Z_1); the Z^{-1/4} is powered by the catalyst |T>. Outputs three |T>; one is
                  kept as the next catalyst. Catalyst treatment, two models (the paper, Sec. 3, says a catalyst error
                  poisons all later outputs, so the factory belongs at the end of a distillation chain):
                    ideal   : the catalyst entering every block is perfect ("once the catalyst is paid for"), the model
                              of the 4 tau input count; kappa counts errors on the delivered T states only;
                    carried : the catalyst is passed from block to block within one target, starting perfect; an error
                              on the catalyst a block passes on corrupts every later block (Z^{-1/4} fault = Z on all
                              three outputs, checked), so for every block but the last all three outputs must be right.
                  A detected fault also discards the catalyst (Gidney-Fowler Sec. 5); re-bootstrapping is not costed.
  [[6m+2,3m,2]]   m CCZ states from one check qubit: blocks of [[8,3,2]] rotations glued along a tree (block j adds the
                  pivots of its ancestors to its supports, sign alternating with depth); verified for m=1..9 by
                  unified_search.is_factory, the output monomials (exactly m disjoint CCZ) and n = 6m+2. It is one
                  construction of this size; its fault patterns may differ from Campbell-Howard's.
Routes (input counts as in route_compare_l3.py):
  cat       4 tau: ceil(tau/2) C2T blocks, two T states each (one if tau is odd), synthesized by CNOT and T gates.
            The T states go to distinct nonzero parities, so a Z pattern on them gives an orthogonal target state
            unless it misses every used T state; kappa is independent of the synthesis.
  lab       per gate of the class representative: CCZ, CS each from its own [[8,3,2]] (8 inputs), T gates in C2T
            blocks (4 each); n_tel = min(n_lab, n_cat), ties broken by the lower kappa.
  batch     all m CS and CCZ gates from one [[6m+2,3m,2]] (block j to gate j), T gates as in lab; n_tel_batch =
            min(n_batch, n_tel), ties by the lower kappa. kappa minimised over 6 variants (which output of a block
            the CS Y measurement takes, gate-to-block order forward or reversed), not over all assignments.
            One pair here can corrupt several gates on overlapping qubits; the target fidelity is computed exactly
            on the 2^k basis, averaged over the Y outcomes and the teleportation outcomes.
Acceptance, first order: a factory of n inputs accepts with 1 - n p, and its expected inputs per accepted output are
n + n^2 p + O(p^2). For a route, rej = sum of the inputs of the sub-factory runs it uses (P[some run rejects] = rej p),
and the expected inputs are n_route + ov p with ov = sum_j w_j n_j^2 (w_j = share of run j the target uses, 1/2 for the
spare T of an odd block).
Known-answer checks (asserted, recorded): [[8,3,2]] kappa 28 = C(8,2) by statevector; [[10,2,2]] 9 per output, 13 joint
(Meier-Eastin-Knill arXiv:1204.4221; the circuit of searches/catalytic_check.py rebuilt on this engine); [[15,1,3]] no
accepted harmful single or pair, 35 harmful triples; statevector pair counts = error_coefficients.coefficients on 5
classes; noiseless C2T gives |T>^3 with fidelity 1; the Z-pattern computation equals the statevector one for m=1.
Deterministic, single process.
    python searches/route_kappa_l3.py   -> outputs/checks/route_kappa_l3.json
"""
import os, sys, json, csv, glob, itertools, statistics as st
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
csv.field_size_limit(sys.maxsize)
import unified_search as us
from error_coefficients import coefficients, gates_of

I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Z = np.diag([1., -1.]); H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
S = np.diag([1, 1j]); T = np.diag([1, np.exp(1j * np.pi / 4)]); Tp = T @ np.ones(2) / np.sqrt(2)
P0 = np.diag([1., 0.]); Pm = np.array([[1, -1], [-1, 1]]) / 2
TOL = 1e-9


def frac(x):
    f = Fraction(float(x)).limit_denominator(1 << 12); assert abs(float(f) - x) < 1e-7, x; return f


# ------------------------------------------------------------------ dense statevector engine (qubit 0 most significant)
class Circ:
    def __init__(s, n):
        s.n = n; s.psi = np.zeros(1 << n, complex); s.psi[0] = 1; s.inputs = 0; s.faults = set()
        s.bits = np.array([[(x >> (n - 1 - q)) & 1 for q in range(n)] for x in range(1 << n)])

    def op(s, U, q):
        m = [I2] * s.n; m[q] = U; R = m[0]
        for x in m[1:]: R = np.kron(R, x)
        return R

    def u1(s, U, q): s.psi = s.op(U, q) @ s.psi

    def par(s, sup): return s.bits[:, sorted(sup)].sum(1) % 2 if sup else np.zeros(1 << s.n, int)

    def rot(s, sup, c):                                    # input rotation exp(2 pi i c/8 [odd parity]), faultable
        s.psi = s.psi * np.exp(2j * np.pi * c / 8 * s.par(sup))
        if s.inputs in s.faults: s.psi = s.psi * (1 - 2 * s.par(sup))
        s.inputs += 1

    def u1_input(s, U, q):                                 # single-qubit T-type gate counted as an input
        s.u1(U, q)
        if s.inputs in s.faults: s.u1(Z, q)
        s.inputs += 1

    def ctrl(s, P, q, U):                                  # U controlled on the projector P of qubit q
        A = s.op(P, q); s.psi = (A @ U + (np.eye(1 << s.n) - A)) @ s.psi

    def proj(s, q, v):                                     # project qubit q on the 1-qubit state v (kept, not removed)
        s.psi = s.op(np.outer(v, v.conj()), q) @ s.psi

    def reduced(s, keep, fixed):                           # amplitudes on `keep` after projecting the other qubits
        out = s.psi.reshape([2] * s.n)
        for q in sorted(fixed, reverse=True):
            out = np.tensordot(out, fixed[q].conj(), axes=([q], [0]))
        return out.reshape(-1)


PLUS = np.ones(2) / np.sqrt(2); PI = np.array([1, 1j]) / np.sqrt(2); MI = np.array([1, -1j]) / np.sqrt(2)


_CORR = {}


def clifford_correction(gates, k, N):
    """The diagonal Clifford that returns the checks to |+>: minus every monomial of the noiseless phase polynomial
    that touches a check (Clifford by the factory condition, asserted). Basis order: qubit 0 most significant."""
    key = (tuple(sorted((tuple(sorted(s_)), c) for s_, c in gates)), k, N)
    if key not in _CORR:
        assert us.is_factory(dict(k=k, r=N - k, gates=[(frozenset(s_), c) for s_, c in gates]), 3)
        ph = us.phases([(frozenset(N - 1 - q for q in s_), c) for s_, c in gates], N, 3)   # us: qubit q = bit q
        P = us.poly_from_table(list(ph), N, 8)
        bits = np.array([[(x >> (N - 1 - q)) & 1 for q in range(N)] for x in range(1 << N)])
        f = np.zeros(1 << N)
        for Tm, a in P.items():
            Tq = [N - 1 - b for b in Tm]
            if all(q < k for q in Tq): continue
            f += a * bits[:, Tq].prod(1)
        _CORR[key] = np.exp(-2j * np.pi * f / 8)
    return _CORR[key]


def factory_run(gates, k, N, faults):
    """Distil with the factory rotations on N qubits (outputs 0..k-1), apply the diagonal Clifford correction, project
    the checks on |+>. Returns (acceptance probability, unnormalised state on the outputs, inputs)."""
    c = Circ(N); c.faults = set(faults)
    for q in range(N): c.u1(H, q)
    for sup, cf in gates: c.rot(sup, cf)
    c.psi = c.psi * clifford_correction(gates, k, N)
    v = c.reduced(list(range(k)), {q: PLUS for q in range(k, N)})
    return float(np.linalg.norm(v) ** 2), v, c.inputs


def unit(v): return v / np.linalg.norm(v)


def pair_count_statevector(gates, k, N, order=2):
    """Number of accepted harmful fault sets of each size <= order, by statevector, reference = noiseless output."""
    acc0, v0, n = factory_run(gates, k, N, ())
    assert abs(acc0 - 1) < TOL, acc0
    ref = unit(v0); res = {}
    for w in range(1, order + 1):
        cnt = 0
        for F in itertools.combinations(range(n), w):
            a, v, _ = factory_run(gates, k, N, F)
            if a > TOL and 1 - abs(np.vdot(ref, unit(v))) ** 2 > TOL: cnt += 1
        res[w] = cnt
    return n, res


# ------------------------------------------------------------------ CCZ factory, then conversions (exact statevector)
def ccz_ref(): return np.array([1, 1, 1, 1, 1, 1, 1, -1], complex) / np.sqrt(8)


def ccz_factory_state(G8, faults):
    """[[8,3,2]] accepted state on 3 qubits, with the diagonal Clifford correction fixed by the noiseless run."""
    acc, v, _ = factory_run(G8, 3, 4, faults)
    _, v0, _ = factory_run(G8, 3, 4, ())
    corr = ccz_ref() / unit(v0)                            # diagonal correction: noiseless output -> |CCZ>
    assert np.allclose(np.abs(corr), 1)
    w = corr * v / max(np.linalg.norm(v), 1e-300)
    return acc, w                                          # w normalised when acc > 0


def c2t(state3, cat_fault=False):
    """Gidney-Fowler Fig. 10 on |CCZ> (3 qubits); returns the 3-qubit output (ideally |T>^3)."""
    c = Circ(3); c.psi = state3.copy(); c.faults = {0} if cat_fault else set()
    c.u1(H @ S.conj().T @ H, 2)                            # X^{-1/2}
    c.ctrl(P0, 2, c.op(X, 0) @ c.op(X, 1))                # X_0 X_1 controlled on |0> of q2
    zz = c.op(Z, 0) @ c.op(Z, 1)
    c.ctrl(Pm, 2, zz); c.u1_input(T.conj().T, 2); c.ctrl(Pm, 2, zz)   # Z^{-1/4} powered by the catalyst
    return c.psi


TTT = np.kron(np.kron(Tp, Tp), Tp)


def fourier3(psi):
    """alpha_v with psi = sum_v alpha_v Z^v |T>^3 (requires psi = E |T>^3, E diagonal unimodular; asserted)."""
    E = psi / TTT; assert np.allclose(np.abs(E), 1, atol=1e-9), "C2T error not a diagonal unitary"
    y = np.array([[(x >> (2 - q)) & 1 for q in range(3)] for x in range(8)])
    return {v: complex(np.mean(E * (-1.) ** (y @ np.array(v)))) for v in itertools.product((0, 1), repeat=3)}


def cs_from_ccz(state3, qm):
    """Y-measure qubit qm of a 3-qubit state, correct each outcome by the diagonal Clifford fixed by the noiseless run;
    returns [(prob, 2-qubit state)]. Reference: CS|++> on the two other qubits in increasing order."""
    ref = np.array([1, 1, 1, 1j]) / 2; out = []
    for v in (PI, MI):
        def red(s3):
            t = np.tensordot(s3.reshape(2, 2, 2), v.conj(), axes=([qm], [0])); return t.reshape(-1)
        r0 = red(ccz_ref()); corr = ref / unit(r0); assert np.allclose(np.abs(corr), 1)
        r = red(state3); p = float(np.linalg.norm(r) ** 2)
        out.append((p, corr * r / np.sqrt(p) if p > TOL else r))
    return out


# ------------------------------------------------------------------ the batch [[6m+2,3m,2]]
def batch_gates(m):
    c = 3 * m; rot = {}; B = lambda j: [3 * j, 3 * j + 1, 3 * j + 2]
    path = {0: frozenset()}; sign = {0: 1}
    for j in range(1, m):
        pb, ps = (j - 1) // 3, (j - 1) % 3; path[j] = path[pb] | {B(pb)[ps]}; sign[j] = -sign[pb]
    for j in range(m):
        for r in range(4):
            for Xs in itertools.combinations(B(j), r):
                Sp = frozenset(Xs) | {c} | path[j]; rot[Sp] = (rot.get(Sp, 0) + sign[j]) % 8
    return [(Sp, v) for Sp, v in rot.items() if v]


def batch_output_ok(g, m):
    k = 3 * m
    for t in (1, 2, 3):
        mod = {1: 2, 2: 4, 3: 8}[t]
        for Tt in itertools.combinations(range(k), t):
            s = sum(v for Sp, v in g if set(Tt) <= Sp) * ((-2) ** (t - 1))
            want = 4 if (t == 3 and Tt[0] % 3 == 0 and Tt == (Tt[0], Tt[0] + 1, Tt[0] + 2)) else 0
            if (s - want) % mod: return False
    return True


def batch_pairs(m):
    """Accepted pairs of the batch and their Z patterns on the outputs (all pairs pass: one check qubit, distance 2)."""
    g = batch_gates(m); assert sum(us.inputs(v, 3) for _, v in g) == 6 * m + 2
    assert us.is_factory(dict(k=3 * m, r=1, gates=g), 3) and batch_output_ok(g, m) and us.distance(dict(k=3 * m, r=1, gates=g), 3) == 2
    f = []
    for Sp, v in g:
        for _ in range(us.inputs(v, 3)): f.append((frozenset(Sp - set(range(3 * m))), frozenset(Sp & set(range(3 * m)))))
    assert all(ch for ch, _ in f)
    return [o1 ^ o2 for (c1, o1), (c2, o2) in itertools.combinations(f, 2) if c1 == c2]


# ------------------------------------------------------------------ gate errors as diagonal functions, target fidelity
def ccz_err_table():
    return {a: [(1.0, lambda y, a=a: (-1) ** sum(ai * yi for ai, yi in zip(a, y)))] for a in itertools.product((0, 1), repeat=3)}


def cs_err_table(qm):
    """For each Z pattern a on the CCZ state: [(prob, d)] with d the diagonal error (2 bits) on CS|++> after correction."""
    tab = {}; ref = np.array([1, 1, 1, 1j]) / 2
    for a in itertools.product((0, 1), repeat=3):
        y3 = np.array([[(x >> (2 - q)) & 1 for q in range(3)] for x in range(8)])
        s3 = ccz_ref() * (-1.) ** (y3 @ np.array(a))
        lst = []
        for p, v2 in cs_from_ccz(s3, qm):
            d = v2 / ref; assert np.allclose(np.abs(d), 1)
            lst.append((p, {(y >> 1 & 1, y & 1): complex(d[y]) for y in range(4)}))
        tab[a] = lst
    return tab


_ALT = {}


def gate_alts(k, sup, br):
    """[(prob, array over x in {0,1}^k)] for one gate's error: Y/branch outcomes times uniform teleportation shifts
    (the shifts reduce to one when the error is a Z-type sign pattern, which shifts only by a global phase)."""
    xs = np.array([[(x >> q) & 1 for q in range(k)] for x in range(1 << k)])
    alts = []
    for p, d in br:
        vals = lambda mm: np.array([d[tuple(int(xs[x, q]) ^ mq for q, mq in zip(sup, mm))] for x in range(1 << k)])
        base = vals((0,) * len(sup))
        pauli = all(np.allclose(vals(mm) / base, (vals(mm) / base)[0]) for mm in itertools.product((0, 1), repeat=len(sup)))
        shifts = [(0,) * len(sup)] if pauli else list(itertools.product((0, 1), repeat=len(sup)))
        for mm in shifts: alts.append((p / len(shifts), vals(mm)))
    return alts


def target_fidelity(k, items):
    """items: [alts of each faulted gate]. Exact E_{branches, teleport outcomes} |mean_x prod_j d_j(x)|^2."""
    n_comb = int(np.prod([len(a) for a in items])); assert n_comb <= 1 << 16, n_comb
    F = 0.0
    for combo in itertools.product(*items):
        p = np.prod([c[0] for c in combo]); D = np.prod([c[1] for c in combo], axis=0)
        F += p * abs(np.mean(D)) ** 2
    return F


# ------------------------------------------------------------------ main
def parse(spec):
    out = []
    for g in spec.split(","):
        nm = g.rstrip("0123456789"); out.append((nm, tuple(int(c) for c in g[len(nm):])))
    return out


def main():
    checks = {}
    rows = [r for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"] if r["l"] == 3 and r["k"] <= 6]
    assert len(rows) == 766
    by = {r["spec"]: r for r in rows}
    G8 = [(frozenset(s), c) for s, c in by["CCZ012"]["gates"]]

    # (1) [[8,3,2]] kappa 28 by statevector
    n8, cnt = pair_count_statevector(G8, 3, 4); assert n8 == 8 and cnt == {1: 0, 2: 28}, cnt
    checks["[[8,3,2]]"] = dict(inputs=n8, accepted_harmful_singles=cnt[1], accepted_harmful_pairs=cnt[2], expected=28)

    # (2) [[10,2,2]] (catalytic_check.py circuit) on this engine: 9 per output, 13 joint
    L0 = [(0, 2, 3), (0, 1, 3), (0, 3), (0, 1, 2, 3), (0, 1, 2), (0, 1), (0, 2)]
    def mek(F):
        c = Circ(4); c.faults = set(F)
        for q in range(4): c.u1(H, q)
        for s_ in L0: c.rot(set(s_), 1)
        c.rot({0}, 1); c.u1(H @ S @ H, 1)
        c.ctrl(np.diag([0., 1.]), 1, c.op(X, 2) @ c.op(X, 3)); zz = c.op(Z, 2) @ c.op(Z, 3)
        c.ctrl(Pm, 1, zz); c.u1_input(T.conj().T, 1); c.ctrl(Pm, 1, zz); c.u1_input(T, 1); c.u1(S.conj().T, 1)
        v = c.reduced([2, 3], {0: PLUS, 1: PLUS}); return float(np.linalg.norm(v) ** 2), v, c.inputs
    a0, v0, n10 = mek(()); assert n10 == 10 and abs(a0 - 1) < TOL and abs(abs(np.vdot(np.kron(Tp, Tp), unit(v0))) - 1) < TOL
    e = [0.0, 0.0, 0.0]; singles_bad = 0
    for i in range(n10):
        a, v, _ = mek((i,)); singles_bad += a > TOL and 1 - abs(np.vdot(np.kron(Tp, Tp), unit(v))) ** 2 > TOL
    for F in itertools.combinations(range(n10), 2):
        a, v, _ = mek(F)
        if a < TOL: continue
        M = unit(v).reshape(2, 2)
        e[0] += a * (1 - np.real(Tp.conj() @ (M @ M.conj().T) @ Tp)); e[1] += a * (1 - np.real(Tp.conj() @ (M.T @ M.conj()) @ Tp))
        e[2] += a * (1 - abs(np.vdot(np.kron(Tp, Tp), unit(v))) ** 2)
    assert singles_bad == 0 and [round(x, 6) for x in e] == [9, 9, 13], e
    checks["[[10,2,2]]"] = dict(inputs=n10, per_output=[round(e[0], 9), round(e[1], 9)], joint=round(e[2], 9), expected="9, 9, 13")

    # (3) [[15,1,3]]: 35 at third order
    f15 = glob.glob(os.path.join(ROOT, "outputs", "large_unified_l3_with_quirk.csv.zip"))
    import zipfile
    with zipfile.ZipFile(f15[0]) as zf:
        name = [x for x in zf.namelist() if x.endswith(".csv")][0]
        rd = csv.DictReader(line.decode() for line in zf.open(name))
        r15 = next(r for r in rd if r["l"] == "3" and r["n"] == "15" and r["k"] == "1" and r["d"] == "3")
    G15 = [(frozenset(s_), c) for s_, c in gates_of(r15)]
    n15, cnt15 = pair_count_statevector(G15, 1, int(r15["N"]), order=3)
    assert n15 == 15 and cnt15 == {1: 0, 2: 0, 3: 35}, cnt15
    checks["[[15,1,3]]"] = dict(inputs=n15, N=int(r15["N"]), accepted_harmful=cnt15, expected="0, 0, 35")

    # (4) statevector pair count = error_coefficients.coefficients on 5 classes
    five = []
    for spec in ("T0", "CS01", "CCZ012", "CCZ014,CCZ023", "T0,CCZ012"):
        r = by[spec]; g = [(frozenset(s_), c) for s_, c in r["gates"]]
        n_, c_ = coefficients(g, r["k"], 3); nn, cc = pair_count_statevector(g, r["k"], r["N"])
        assert n_ == nn == r["n"] and c_ == cc[2] == r["c"] and cc[1] == 0, (spec, n_, c_, nn, cc, r["c"])
        five.append(dict(spec=spec, n=r["n"], N=r["N"], c_file=r["c"], c_coefficients=c_, c_statevector=cc[2]))
    checks["factory_kappa_5_classes"] = five

    # sub-factory tables from [[8,3,2]] pairs, exact statevector
    pairs8 = list(itertools.combinations(range(8), 2))
    states8 = {}
    for F in pairs8:
        a, w = ccz_factory_state(G8, F); assert abs(a - 1) < TOL; states8[F] = w
    _, w0 = ccz_factory_state(G8, ()); assert abs(abs(np.vdot(ccz_ref(), w0)) - 1) < TOL
    k_ccz = sum(1 - abs(np.vdot(ccz_ref(), w)) ** 2 for w in states8.values())
    refcs = np.array([1, 1, 1, 1j]) / 2
    k_cs = {qm: sum(sum(p * (1 - abs(np.vdot(refcs, v)) ** 2) for p, v in cs_from_ccz(w, qm)) for w in states8.values()) for qm in range(3)}
    for qm in range(3): assert all(abs(abs(np.vdot(refcs, v)) - 1) < TOL for p, v in cs_from_ccz(w0, qm))
    out0 = c2t(w0); assert abs(abs(np.vdot(TTT, out0)) - 1) < TOL
    al0 = fourier3(c2t(w0, cat_fault=True)); assert abs(abs(al0[(1, 1, 1)]) - 1) < TOL   # catalyst fault: Z on all three
    A = [fourier3(c2t(w)) for w in states8.values()]
    def K_c2t(cat, used, carried):
        """sum over the 28 pairs of 1 - fidelity; used = delivered outputs the target uses; carried: catalyst must be right too."""
        tot = 0.0
        for al in A:
            need = set(used) | ({cat} if carried else set())
            tot += 1 - sum(abs(x) ** 2 for v, x in al.items() if all(v[q] == 0 for q in need))
        return tot
    c2t_tab = {}
    for cat in range(3):
        dl = [q for q in range(3) if q != cat]
        c2t_tab[cat] = dict(K2=K_c2t(cat, dl, False), K1=min(K_c2t(cat, [q], False) for q in dl),
                            K2_carried=K_c2t(cat, dl, True), K1_carried=min(K_c2t(cat, [q], True) for q in dl))
    cat_q = min(range(3), key=lambda c_: (c2t_tab[c_]["K2"], c2t_tab[c_]["K1"], c2t_tab[c_]["K2_carried"]))
    K2, K1, K2c = (frac(c2t_tab[cat_q][x]) for x in ("K2", "K1", "K2_carried"))
    KCCZ = frac(k_ccz); qm_best = min(range(3), key=lambda q: k_cs[q]); KCS = frac(k_cs[qm_best])
    checks["subfactory_tables"] = dict(kappa_CCZ_from_832=str(KCCZ), kappa_CS_from_832_Y={q: str(frac(k_cs[q])) for q in range(3)},
                                       c2t_by_catalyst_qubit={c_: {x: str(frac(v)) for x, v in t.items()} for c_, t in c2t_tab.items()},
                                       c2t_catalyst_qubit=cat_q, c2t_noiseless_fidelity=round(abs(np.vdot(TTT, out0)) ** 2, 12),
                                       catalyst_fault_gives_Z_on_all_three=True)

    def kT(nt, carried):
        """kappa of nt T states from C2T blocks within one target (T states to distinct nonzero parities)."""
        if nt == 0: return Fraction(0)
        B = (nt + 1) // 2; last = K1 if nt % 2 else K2
        return (B - 1) * (K2c if carried else K2) + last

    def acc_T(nt):                                         # (rej, ov) for nt T states from C2T blocks
        B = (nt + 1) // 2; return 8 * B, Fraction(64 * (nt // 2)) + (32 if nt % 2 else 0)

    # batch tables: Z-pattern computation; consistency with statevector at m = 1
    ccz_tab = ccz_err_table(); cs_tab = {qm: cs_err_table(qm) for qm in range(3)}
    pat8 = batch_pairs(1)
    def br_ccz(a): return [(1.0, {y: ccz_tab[a][0][1](y) for y in itertools.product((0, 1), repeat=3)})]
    pa = lambda o: tuple(1 if q in o else 0 for q in range(3))
    k_ccz_pat = sum(1 - target_fidelity(3, [gate_alts(3, (0, 1, 2), br_ccz(pa(o)))]) for o in pat8 if o)
    assert len(pat8) == 28 and abs(k_ccz_pat - float(KCCZ)) < 1e-9, (k_ccz_pat, KCCZ)
    for qm in range(3):
        kp = sum(1 - target_fidelity(2, [gate_alts(2, (0, 1), cs_tab[qm][pa(o)])]) for o in pat8 if o)
        assert abs(kp - k_cs[qm]) < 1e-9, (qm, kp, k_cs[qm])
    checks["pattern_equals_statevector_m1"] = True
    PAIRS = {m: batch_pairs(m) for m in range(1, 10)}
    checks["batch_6m2"] = {m: dict(n=6 * m + 2, accepted_pairs=len(PAIRS[m]), harmful_pairs_on_outputs=sum(1 for o in PAIRS[m] if o))
                           for m in PAIRS}

    PAT = {m: {} for m in PAIRS}
    for m in PAIRS:
        for o in PAIRS[m]:
            if o: PAT[m][o] = PAT[m].get(o, 0) + 1

    def kappa_batch(k, gates_cc, qm, order):
        m = len(gates_cc); gl = gates_cc if order == 0 else gates_cc[::-1]; tot = 0.0
        for o, mult in PAT[m].items():
            items = []
            for j, (nm, sup) in enumerate(gl):
                a = tuple(1 if 3 * j + q in o else 0 for q in range(3))
                if not any(a): continue
                key = (k, nm, sup, qm if nm == "CS" else -1, a)
                if key not in _ALT:
                    _ALT[key] = gate_alts(k, sup, br_ccz(a) if nm == "CCZ" else cs_tab[qm][a])
                items.append(_ALT[key])
            tot += mult * (1 - target_fidelity(k, items))
        return frac(tot)

    out = []
    for r in rows:
        k, n, tau = r["k"], r["n"], r["tau"]
        g = [(frozenset(s_), c) for s_, c in r["gates"]]
        n_c, c_c = coefficients(g, k, 3); assert n_c == n and c_c == r["c"]
        lab = parse(r["spec"]); nt = sum(nm == "T" for nm, _ in lab); cc = [(nm, s_) for nm, s_ in lab if nm != "T"]
        ncs = sum(nm == "CS" for nm, _ in cc); nccz = len(cc) - ncs; m = len(cc)
        # cat
        n_cat = 4 * tau; k_cat = kT(tau, False); k_cat_c = kT(tau, True); rej_cat, ov_cat = acc_T(tau)
        # per gate
        n_lab = 8 * m + 4 * nt; k_lab = nccz * KCCZ + ncs * KCS + kT(nt, False); k_lab_c = nccz * KCCZ + ncs * KCS + kT(nt, True)
        rT, oT = acc_T(nt); rej_lab, ov_lab = 8 * m + rT, Fraction(64 * m) + oT
        # batch
        if m:
            kb = min(kappa_batch(k, cc, qm, od) for qm in range(3) for od in (0, 1))
            n_b = 6 * m + 2 + 4 * nt; k_b = kb + kT(nt, False); k_b_c = kb + kT(nt, True)
            rej_b, ov_b = 6 * m + 2 + rT, Fraction((6 * m + 2) ** 2) + oT
        else:
            n_b = None
        cands = {"cat": (n_cat, k_cat, k_cat_c, rej_cat, ov_cat), "lab": (n_lab, k_lab, k_lab_c, rej_lab, ov_lab)}
        if m: cands["batch"] = (n_b, k_b, k_b_c, rej_b, ov_b)
        pick = lambda names: min(((cands[x][0], cands[x][1], x) for x in names if x in cands))[2]
        routes = {"cat": "cat", "tel": pick(["lab", "cat"]), "tel_batch": pick(["batch", "lab", "cat"])}
        routes["best"] = routes["tel_batch"]                 # cheapest of the three by inputs (route_margin.py), ties by kappa
        row = dict(k=k, orbit=r["orbit"], spec=r["spec"], n=n, tau=tau, syn_n=r.get("syn_n"), kappa_factory=c_c,
                   rej_factory=n, ov_factory=n * n)
        for rn, cn in routes.items():
            nn, kk, kc, rj, ov = cands[cn]
            row[rn] = dict(construction=cn, n=nn, kappa=str(kk), kappa_carried=str(kc), rej=rj, ov=str(ov))
        out.append(row)

    rc = {r["spec"]: r for r in json.load(open(os.path.join(ROOT, "outputs", "checks", "route_compare_l3.json")))["rows"]}
    for row in out:                                        # input counts must equal route_compare_l3.json
        q = rc[row["spec"]]
        assert row["cat"]["n"] == q["n_cat"] and row["tel"]["n"] == q["n_tel"] and row["tel_batch"]["n"] == q["n_tel_batch"], row["spec"]

    def cmp(rows_, rn, key="kappa"):
        f = e_ = mo = 0; ratios = []
        for r in rows_:
            kr = Fraction(r[rn][key]); kf = r["kappa_factory"]
            f += kf < kr; e_ += kf == kr; mo += kf > kr; ratios.append(float(kf / kr) if kr else None)
        rr = [x for x in ratios if x is not None]
        return dict(classes=len(rows_), factory_lower=f, equal=e_, factory_higher=mo,
                    median_ratio_factory_over_route=round(st.median(rr), 3) if rr else None,
                    min_ratio=round(min(rr), 3) if rr else None, max_ratio=round(max(rr), 3) if rr else None)

    def cmp_acc(rows_, rn):
        a = dict(rej_factory_lower=0, rej_equal=0, rej_factory_higher=0, ov_factory_lower=0, ov_equal=0, ov_factory_higher=0)
        for r in rows_:
            for key, fv, rv in (("rej", r["rej_factory"], r[rn]["rej"]), ("ov", Fraction(r["ov_factory"]), Fraction(r[rn]["ov"]))):
                a[f"{key}_factory_lower" if fv < rv else (f"{key}_equal" if fv == rv else f"{key}_factory_higher")] += 1
        return a

    subsets = {"all": out, "below_synthillation": [r for r in out if r["syn_n"] is not None and r["n"] < r["syn_n"]],
               "factory_fewer_inputs_than_best_route": [r for r in out if r["n"] < r["best"]["n"]]}
    subsets["below_synthillation_and_fewer_than_best_route"] = [r for r in subsets["below_synthillation"] if r["n"] < r["best"]["n"]]
    summary = {}
    for nm_, rr in subsets.items():
        summary[nm_] = dict(classes=len(rr))
        for rn in ("cat", "tel", "tel_batch", "best"):
            summary[nm_][rn] = dict(kappa_ideal_catalyst=cmp(rr, rn), kappa_carried_catalyst=cmp(rr, rn, "kappa_carried"),
                                    acceptance=cmp_acc(rr, rn))
    assert summary["all"]["classes"] == 766 and summary["below_synthillation"]["classes"] == 330
    assert summary["factory_fewer_inputs_than_best_route"]["classes"] == 659
    assert summary["below_synthillation_and_fewer_than_best_route"]["classes"] == 314
    summary["cat_per_T"] = dict(K2_two_T=str(K2), K1_one_T=str(K1), K2_carried=str(K2c),
                                note="kappa_cat = (B-1) K2 + (K1 if tau odd else K2), B = ceil(tau/2); carried: K2_carried for all blocks but the last")
    res = dict(source="outputs/staged_best.json, outputs/checks/route_compare_l3.json, outputs/large_unified_l3_with_quirk.csv.zip ([[15,1,3]] only)",
               model="Z faults with probability p on every input; kappa = sum over accepted pairs within one sub-factory run of (1 - target fidelity)",
               known_answer_checks=checks, summary=summary, rows=out)
    os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
    json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "route_kappa_l3.json"), "w"), indent=1, default=str)
    print(json.dumps(dict(checks={k_: v for k_, v in checks.items() if k_ != "batch_6m2"}, summary=summary), indent=1, default=str))


if __name__ == "__main__":
    main()
