"""Full-Clifford classes of the 766 l=3 output classes with k <= 6, and the win count against synthillation under them.

The paper's classes are orbits of |U> = U|+>^k under CNOT frames, Pauli X and free diagonal gates (the group that
maps diagonal states to diagonal states). The full Clifford group can merge some of them (CS12 CS13 ~ CCZ by Hadamards,
Howard-Campbell). This script finds the merges exactly.

Method (exact, no search heuristics). For a pure state psi on k qubits, a Clifford C has C psi = D|+>^k with D diagonal
iff C psi has amplitudes of equal modulus, iff <C psi|Z^b|C psi> = 0 for every b != 0, iff the Lagrangian subspace
L = C^dag <Z_1..Z_k> C of the Pauli group (as F_2^{2k}) lies, apart from 0, inside the zero set
{P : <psi|P|psi> = 0}. Two Cliffords that send the same L to <Z_i> differ by a Clifford that fixes <Z_i>: CNOT frames,
Pauli X and Z, and diagonal Cliffords, i.e. exactly the group the existing labeller quotients by. So the diagonal states
Clifford-equivalent to psi are, class by class, the labels of C_L psi over every Lagrangian L in the zero set. We
enumerate all such L (csrc/lagrangians.c, exhaustive DFS over reduced bases, checked against the total count
prod_{i=1..k}(2^i+1) of Lagrangians), build C_L from H, S, CZ gates, check that C_L psi is flat, read its phase
polynomial and label it with the class labeller of the cascade (cascade_class.classes_of). Classes A, B are merged iff
B is a label of A. Independent checks: the relation must be symmetric; each class must label itself; merged classes
must have the same multiset of |<P>| (Pauli spectrum, an exact Clifford invariant), whose number of distinct values
is a lower bound on the number of full-Clifford classes.

Then per full-Clifford class: n = min over member classes of the smallest factory (staged_best n), syn = min over
members of the synthillation n (staged_best syn_n = synthillation_strongest_l3.json best_n), and wins/same/more.
Run for the factory table the paper quotes (outputs/archive/v4_published/outputs/staged_best.json, 295 wins) and for
the current outputs/staged_best.json.

    python searches/full_clifford_classes.py [--workers W]   -> outputs/checks/full_clifford_classes.{json,md}
Workers: at most min(3, 16 - 5 - load); niced to 19.
"""
import os, sys, json, time, argparse, subprocess, collections, pickle
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
OUT = os.path.join(ROOT, "outputs", "checks")
CSRC = os.path.join(HERE, "csrc", "lagrangians.c"); EXE = os.path.join(HERE, "csrc", "lagrangians")
CACHE = os.path.join(ROOT, "outputs", "staged_A", "_cache")          # read only
SOURCES = {"paper_v4": "outputs/archive/v4_published/outputs/staged_best.json", "current": "outputs/staged_best.json"}


# ------------------------------------------------------------------ states and Pauli spectrum
def parse(spec):
    P = {}
    for g in spec.split(","):
        name = g.rstrip("0123456789"); qs = frozenset(int(c) for c in g[len(name):])
        P[qs] = (P.get(qs, 0) + {"T": 1, "CS": 2, "CCZ": 4}[name]) % 8
    return P


def state(P, k):
    x = np.arange(1 << k); f = np.zeros(1 << k, dtype=np.int64)
    for T, c in P.items():
        m = sum(1 << q for q in T); f += c * ((x & m) == m)
    return np.exp(1j * np.pi / 4 * (f % 8)) / 2 ** (k / 2)


def wht(a):
    a = a.copy(); n = len(a); h = 1
    while h < n:
        a = a.reshape(-1, 2, h); a = np.stack([a[:, 0] + a[:, 1], a[:, 0] - a[:, 1]], 1).reshape(-1); h *= 2
    return a


def spectrum(psi, k):
    """|<psi| X^a Z^b |psi>| as an array S[a, b]; the Pauli X^a Z^b is encoded u = (a << k) | b."""
    N = 1 << k; x = np.arange(N); E = np.zeros((N, N), complex)
    for a in range(N):
        E[a] = wht(np.conj(psi[x ^ a]) * psi)
    return np.abs(E)


def spec_key(S):
    return tuple(sorted(collections.Counter(np.round(S.ravel(), 9).tolist()).items()))


# ------------------------------------------------------------------ Lagrangians and the Clifford that diagonalizes
def build_exe():
    if not os.path.exists(EXE) or os.path.getmtime(EXE) < os.path.getmtime(CSRC):
        subprocess.run(["cc", "-O2", "-o", EXE, CSRC], check=True)


def lagrangians(zero, k):
    """Every Lagrangian L with L \\ {0} inside zero (bool array over u), as reduced bases (fresh process per call)."""
    inp = f"{k}\n" + "".join("1" if v else "0" for v in zero) + "\n"
    out = subprocess.run([EXE], input=inp, capture_output=True, text=True, check=True).stdout
    return [tuple(int(t) for t in line.split()) for line in out.splitlines() if line.strip()]


def rank2(rows):
    rows = list(rows); r = 0
    for bit in range(32):
        p = next((i for i in range(r, len(rows)) if rows[i] >> bit & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i] >> bit & 1: rows[i] ^= rows[r]
        r += 1
    return r


def diagonalizer(basis, k):
    """Gate list (('H',q) | ('S',q) | ('CZ',i,j)) of a Clifford C with C P C^dag in <Z_i> for every P in span(basis)."""
    m = (1 << k) - 1
    A = [(u >> k) & m for u in basis]; B = [u & m for u in basis]
    for s in sorted(range(1 << k), key=lambda s: (bin(s).count("1"), s)):
        A2 = [(a & ~s) | (b & s) for a, b in zip(A, B)]; B2 = [(b & ~s) | (a & s) for a, b in zip(A, B)]
        if rank2(A2) == k: break
    else:
        raise AssertionError("no Hadamard subset makes the X part full rank")
    gates = [("H", q) for q in range(k) if s >> q & 1]
    A, B = A2, B2
    for q in range(k):                                         # row-reduce the X part to the identity
        p = next(i for i in range(q, k) if A[i] >> q & 1)
        A[q], A[p] = A[p], A[q]; B[q], B[p] = B[p], B[q]
        for i in range(k):
            if i != q and A[i] >> q & 1: A[i] ^= A[q]; B[i] ^= B[q]
    assert all(A[i] == 1 << i for i in range(k))
    M = [[B[i] >> j & 1 for j in range(k)] for i in range(k)]
    assert all(M[i][j] == M[j][i] for i in range(k) for j in range(k)), "not isotropic"
    gates += [("CZ", i, j) for i in range(k) for j in range(i + 1, k) if M[i][j]]
    gates += [("S", i) for i in range(k) if M[i][i]]
    gates += [("H", q) for q in range(k)]
    return gates


def apply(psi, gates, k):
    x = np.arange(1 << k); psi = psi.copy()
    for g in gates:
        if g[0] == "H":
            bit = 1 << g[1]; a0 = psi[x & ~bit]; a1 = psi[x | bit]
            psi = np.where(x & bit, a0 - a1, a0 + a1) / np.sqrt(2)
        elif g[0] == "S":
            psi = np.where(x >> g[1] & 1, 1j * psi, psi)
        else:
            psi = np.where((x >> g[1] & 1) & (x >> g[2] & 1), -psi, psi)
    return psi


def phase_poly(phi, k):
    """Phase polynomial {support: coeff mod 8} of a flat state with 8th-root phases, or a reason it is not one."""
    if np.max(np.abs(np.abs(phi) - 2 ** (-k / 2))) > 1e-9:
        return None, "not flat"
    r = phi / phi[0]; f = np.angle(r) / (np.pi / 4)
    fi = np.round(f).astype(np.int64)
    if np.max(np.abs(f - fi)) > 1e-7:
        return None, "phases not 8th roots"
    from unified_search import poly_from_table
    P = poly_from_table((fi % 8).tolist(), k, 8)
    for T, c in P.items():
        if c % (1 << (len(T) - 1)) or (len(T) > 3 and c % 8):
            return None, "not level 3"
    return P, None


# ------------------------------------------------------------------ labeller (the cascade's exact classes)
_LAB = {}


def labeller(k):
    if k not in _LAB:
        import cascade_class as CC
        path = os.path.join(CACHE, f"l3_k{k}_classes.pkl")
        classes, ctx = pickle.load(open(path, "rb")) if os.path.exists(path) else CC.classes_of(3, k)
        C, off, lab, F = ctx
        from enumerate_targets import index_of_poly
        def label(P):
            # the labeller expects the AND-monomial (Newton) form: index_of_poly / digits_of_poly re-derive it from a table
            Pn = {T: c for T, c in P.items()}
            if lab is not None:
                return str(int(lab[index_of_poly(Pn, k, 3, C, off)]))
            a, b = F.label(F.digits_of_poly(Pn)); return f"{a}:{b}"
        _LAB[k] = (classes, label)
    return _LAB[k]


def and_form(P_newton):
    """poly_from_table gives Newton (AND-monomial) coefficients; parse() gives the same form. Identity, kept explicit."""
    return P_newton


def one(item):
    k, orbit, spec = item
    t = time.time()
    psi = state(parse(spec), k); S = spectrum(psi, k)
    zero = np.zeros(1 << (2 * k), bool); A, B = np.nonzero(S < 1e-9); zero[(A << k) | B] = True; zero[0] = False
    Ls = lagrangians(zero, k)
    classes, label = labeller(k)
    labels = collections.Counter(); bad = collections.Counter(); seen = {}
    for basis in Ls:
        phi = apply(psi, diagonalizer(basis, k), k)
        P, why = phase_poly(phi, k)
        if P is None:
            bad[why] += 1; continue
        key = tuple(sorted((tuple(sorted(T)), c) for T, c in P.items()))
        if key not in seen:
            seen[key] = label(and_form(P))
        labels[seen[key]] += 1
    return dict(k=k, orbit=str(orbit), spec=spec, nullity=k - int(round(np.log2((np.abs(S - 1) < 1e-9).sum()))),
                spectrum=[[v, c] for v, c in spec_key(S)], n_zero=int(zero.sum()), n_lagrangians=len(Ls),
                labels=dict(labels), failures=dict(bad), seconds=round(time.time() - t, 2))


# ------------------------------------------------------------------ checks with known answers
def selftest():
    build_exe()
    for k in range(1, 6):                                       # all Lagrangians of F_2^{2k}
        z = np.ones(1 << (2 * k), bool); z[0] = False
        want = int(np.prod([2 ** i + 1 for i in range(1, k + 1)]))
        got = len(lagrangians(z, k)); assert got == want, (k, got, want)
    # every Lagrangian is diagonalized by its gate list (k = 3: all 135)
    k = 3; z = np.ones(1 << 6, bool); z[0] = False
    for basis in lagrangians(z, k):
        g = diagonalizer(basis, k)
        for u in basis:                                         # conjugate the Pauli through the gates
            a, b = u >> k, u & 7
            for G in g:
                if G[0] == "H": q = G[1]; ab, bb = a >> q & 1, b >> q & 1; a ^= (ab ^ bb) << q; b ^= (ab ^ bb) << q
                elif G[0] == "S": b ^= (a >> G[1] & 1) << G[1]
                else: i, j = G[1], G[2]; b ^= ((a >> i & 1) << j) ^ ((a >> j & 1) << i)
            assert a == 0, (basis, g)
    return dict(lagrangian_counts_k1_5="prod(2^i+1) reproduced", diagonalizer_k3="all 135 Lagrangians sent to <Z_i>")


def uf_classes(res):
    parent = {}
    def f(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for r in res:
        a = (r["k"], r["orbit"]); f(a)
        for o in r["labels"]:
            parent[f(a)] = f((r["k"], o))
    groups = collections.defaultdict(list)
    for r in res:
        groups[f((r["k"], r["orbit"]))].append(r)
    return list(groups.values())


def recount(groups, rows):
    R = {(r["k"], str(r["orbit"])): r for r in rows}
    out = collections.Counter(); member_wins_kept = 0; member_wins = 0; by = collections.Counter(); gain = 0
    detail = []
    for g in groups:
        ms = [R[(r["k"], r["orbit"])] for r in g]
        n = min(m["n"] for m in ms); syn = min(m["syn_n"] for m in ms)
        out["fewer" if n < syn else "same" if n == syn else "more"] += 1
        if n < syn: by[syn - n] += 1; gain += syn - n
        for m in ms:
            if m["n"] < m["syn_n"]:
                member_wins += 1; member_wins_kept += m["n"] < syn
        detail.append(dict(members=[f'{m["k"]}:{m["orbit"]}:{m["spec"]}' for m in ms], n=n, syn_n=syn,
                           member_n=[m["n"] for m in ms], member_syn_n=[m["syn_n"] for m in ms]))
    diag = collections.Counter("fewer" if r["n"] < r["syn_n"] else "same" if r["n"] == r["syn_n"] else "more" for r in R.values())
    return dict(diagonal_classes=len(R), diagonal_counts=dict(diag), full_clifford_classes=len(groups),
                full_counts=dict(out), full_fewer_by=dict(by), full_inputs_saved=gain,
                diagonal_wins_that_survive=member_wins_kept, diagonal_wins=member_wins), detail


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--kmax", type=int, default=6, help="test run on k <= kmax only (writes *_test files)"); a = ap.parse_args()
    os.nice(19)
    load = os.getloadavg()[0]; cap = min(3, max(0, int(16 - 5 - load)))
    W = min(a.workers, cap) if a.workers else cap
    while W == 0:
        print(f"load {load:.1f}: no core free, waiting", flush=True); time.sleep(60)
        load = os.getloadavg()[0]; W = min(3, max(0, int(16 - 5 - load)))
    st = selftest(); print("selftest ok", st, flush=True)
    rows = {name: [r for r in json.load(open(os.path.join(ROOT, p)))["rows"] if r["l"] == 3] for name, p in SOURCES.items()}
    keyset = {name: {(r["k"], str(r["orbit"])) for r in rr} for name, rr in rows.items()}
    assert keyset["paper_v4"] == keyset["current"] and len(keyset["current"]) == 766
    if a.kmax < 6:
        rows = {nm: [r for r in rr if r["k"] <= a.kmax] for nm, rr in rows.items()}
    items = sorted({(r["k"], str(r["orbit"]), r["spec"]) for r in rows["current"]})
    items.sort(key=lambda it: -it[0])                          # big cells first
    t0 = time.time(); res = []
    with Pool(W) as p:
        for i, r in enumerate(p.imap_unordered(one, items, chunksize=1)):
            res.append(r)
            if i % 25 == 0:
                print(f"{i + 1}/{len(items)} {time.time() - t0:.0f}s load {os.getloadavg()[0]:.1f}", flush=True)
    res.sort(key=lambda r: (r["k"], r["orbit"]))
    byk = {(r["k"], r["orbit"]): r for r in res}
    # checks
    chk = {}
    chk["self_label"] = sum(r["orbit"] in r["labels"] for r in res)
    chk["labels_outside_class_list"] = sum(1 for r in res for o in r["labels"] if (r["k"], o) not in byk)
    chk["asymmetric_pairs"] = sum(1 for r in res for o in r["labels"] if r["orbit"] not in byk[(r["k"], o)]["labels"])
    chk["diagonalization_failures_not_flat"] = sum(r["failures"].get("not flat", 0) for r in res)
    # a Lagrangian frame whose flat state is not a level-3 diagonal state (phases not 8th roots, or a pair coefficient
    # not even, e.g. phase pi/4 x1 x2): a diagonal state Clifford-equivalent to the target but outside the level-3 class
    # list; the whole frame family (CNOT, X, diagonal Cliffords) stays outside, so it cannot merge two listed classes
    chk["frames_outside_level3"] = sum(r["failures"].get("not level 3", 0) + r["failures"].get("phases not 8th roots", 0) for r in res)
    assert chk["diagonalization_failures_not_flat"] == 0
    groups = uf_classes(res)
    chk["groups_with_unequal_spectra"] = sum(len({json.dumps(m["spectrum"]) for m in g}) > 1 for g in groups)
    chk["groups_not_cliques"] = sum(any(o["orbit"] not in m["labels"] for m in g for o in g) for g in groups)
    spectra = {(r["k"], json.dumps(r["spectrum"])) for r in res}
    chk["distinct_spectra_lower_bound"] = len(spectra)
    chk["nullity_equals_k"] = sum(r["nullity"] == r["k"] for r in res)
    find = lambda spec: next(r for r in res if r["spec"] == spec)
    ka = {"CS01,CS02 ~ CCZ012": find("CCZ012")["orbit"] in find("CS01,CS02")["labels"],
          "T0,CS01,CS02 ~ T0,CCZ012": find("T0,CCZ012")["orbit"] in find("T0,CS01,CS02")["labels"]}
    assert all(ka.values()), ka
    rec = {}; det = {}
    for name in SOURCES:
        rec[name], det[name] = recount(groups, rows[name])
    sizes = collections.Counter(len(g) for g in groups)
    out = dict(method=__doc__.split("\n\n")[1].replace("\n", " "), sources=SOURCES, selftest=st, known_answers=ka,
               checks=chk, full_clifford_classes=len(groups), group_sizes=dict(sorted(sizes.items())),
               full_classes_by_k=dict(sorted(collections.Counter(g[0]["k"] for g in groups).items())),
               recount=rec, groups=[[f'{m["k"]}:{m["orbit"]}:{m["spec"]}' for m in g] for g in groups if len(g) > 1],
               detail_paper_v4=det["paper_v4"], per_class=res, seconds=round(time.time() - t0, 1), workers=W)
    os.makedirs(OUT, exist_ok=True)
    suf = "" if a.kmax == 6 else f"_test_k{a.kmax}"
    json.dump(out, open(os.path.join(OUT, f"full_clifford_classes{suf}.json"), "w"), indent=1)
    write_md(out, suf)
    print(json.dumps(dict(checks=chk, recount=rec, full=len(groups)), indent=1))


def write_md(o, suf=""):
    L = ["# Full-Clifford classes of the 766 l=3 output classes (k <= 6)", "",
         "Script: `searches/full_clifford_classes.py` (rerun: `python searches/full_clifford_classes.py`). "
         "Data: `outputs/checks/full_clifford_classes.json`.", "",
         "## Method", "", o["method"], "",
         "## Checks", "", f"- Self-test: {o['selftest']}", f"- Known answers: {o['known_answers']}"]
    L += [f"- {k}: {v}" for k, v in o["checks"].items()]
    L += ["", "## Result", "",
          f"- The 766 classes (counted up to CNOT, Pauli X and diagonal Cliffords) fall into **{o['full_clifford_classes']} "
          f"classes under the full Clifford group** (exact: every Clifford frame in which a state is diagonal is enumerated). "
          f"Distinct Pauli spectra give the independent lower bound {o['checks']['distinct_spectra_lower_bound']}.",
          f"- Class sizes (number of diagonal classes per full-Clifford class: count): {o['group_sizes']}",
          f"- Full-Clifford classes by k: {o['full_classes_by_k']}", ""]
    for name, r in o["recount"].items():
        L += [f"### Factories from `{o['sources'][name]}`", "",
              f"- Diagonal classes: {r['diagonal_counts']} of {r['diagonal_classes']}.",
              f"- Full-Clifford classes (n = min over members, synthillation = min over members): {r['full_counts']} "
              f"of {r['full_clifford_classes']}; fewer by {r['full_fewer_by']}; inputs saved in total {r['full_inputs_saved']}.",
              f"- Of the {r['diagonal_wins']} diagonal classes where the factory beats synthillation, "
              f"{r['diagonal_wins_that_survive']} still beat the cheapest synthillation over their whole full-Clifford class.", ""]
    L += ["## Merged groups (members k:orbit:spec)", ""] + ["- " + " ~ ".join(g) for g in o["groups"]]
    open(os.path.join(OUT, f"full_clifford_classes{suf}.md"), "w").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
