"""Check of the [[10,2,2]] catalytic conversion of App. S2.1d (and [[11,1,2]]) by exact statevector simulation.

Circuits, on the four parent qubits q0..q3 at l=3 (theta = pi/8, rotation Z^{otimes Lambda}(theta) = exp(i theta (I - prod Z))):
  common start: H on every qubit, then the seven rotations of layer L_0 that touch an output (the [[8,3,2]] stage);
  [[10,2,2]]: T_0 (the eighth input), measure q0 in X; sqrt(X)_1; C_1(X_2 X_3); C^-_1(Z_2 Z_3) T^dag_1 C^-_1(Z_2 Z_3);
              T_1 S^dag_1; measure q1 in X; accept only on the trivial record (+, +). Output on q2 q3.
  [[11,1,2]]: (no T_0) sqrt(X)_1; C_1(X_2 X_3); the same conversion; T_1 T_2 T_3 S^dag_1 S^dag_2 S^dag_3; measure q1 q2 q3
              in X; accept on all +. Output on q0.
C^-_1 is a control on the X-basis state |-> of q1. Checks, written to outputs/checks/catalytic_check.json:
  1. noiseless: acceptance probability 1 and fidelity 1 with the stated output;
  2. every single Z fault on each of the n input rotations (incl. the conversion T^dag_1 and T_1) is caught (changes
     the accepted record) or harmless (output fidelity 1); every fault pair that passes the trivial record with a harmed
     output is listed (inputs indexed in circuit order: the 7 rotations of L_0, then T_0 for [[10,2,2]], T^dag_1, T_1, ...);
  3. rejection probability to first order: n p, with n the number of inputs, from the count of caught single faults;
  4. for [[10,2,2]], the p^2 coefficients of the output error: each fault pair weighted by its acceptance probability
     times its infidelity, for each output (marginal) and for the pair (joint, at least one output wrong). Known answer
     (Meier-Eastin-Knill, arXiv:1204.4221, Sec. III): u(p) = 9p^2 per output, u2(p) = 13p^2 for the pair.
    python searches/catalytic_check.py
"""
import os, json, itertools
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TH = np.pi / 8
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]]); Z = np.diag([1, -1]); H = np.array([[1, 1], [1, -1]]) / np.sqrt(2)
S = np.diag([1, 1j]); T = np.diag([1, np.exp(1j * np.pi / 4)]); Tp = T @ np.array([1, 1]) / np.sqrt(2)
Pp = np.ones((2, 2)) / 2; Pm = np.array([[1, -1], [-1, 1]]) / 2
N = 4
L0 = [(0, 2, 3), (0, 1, 3), (0, 3), (0, 1, 2, 3), (0, 1, 2), (0, 1), (0, 2)]


def op(U, q):
    m = [I2] * N; m[q] = U; R = m[0]
    for x in m[1:]:
        R = np.kron(R, x)
    return R


def rot(Sup):
    return np.diag([np.exp(2j * TH) if sum((x >> (N - 1 - q)) & 1 for q in Sup) % 2 else 1 for x in range(2 ** N)])


def ctrl(P, q, U):
    A = op(P, q); return A @ U + (np.eye(2 ** N) - A)


def fault(Sup):
    R = np.eye(2 ** N)
    for q in Sup:
        R = op(Z, q) @ R
    return R


def run(kind, fi=None):
    """Returns {record: (prob, state)} after the measurements; fi = index of the faulted input (Z on its support)."""
    psi = np.zeros(2 ** N, complex); psi[0] = 1
    for q in range(N):
        psi = op(H, q) @ psi
    inputs = []
    def r(Sup, U=None):
        nonlocal psi
        psi = (U if U is not None else rot(Sup)) @ psi
        if fi is not None and len(inputs) in (fi if isinstance(fi, tuple) else (fi,)):
            psi = fault(Sup) @ psi
        inputs.append(Sup)
    for s in L0:
        r(s)
    if kind == "10":
        r((0,))                                          # T_0, the eighth input of [[10,2,2]] only
    sx = H @ S @ H
    psi = op(sx, 1) @ psi
    X23 = op(X, 2) @ op(X, 3); Z23 = op(Z, 2) @ op(Z, 3); C = ctrl(Pm, 1, Z23)
    psi = ctrl(np.diag([0, 1]), 1, X23) @ psi                # C_1(X_2 X_3): control on |1> of q1
    psi = C @ psi; r((1,), op(T.conj().T, 1)); psi = C @ psi   # conversion T^dag_1, an input rotation
    r((1,), op(T, 1))                                     # T_1
    psi = op(S.conj().T, 1) @ psi
    meas = [0, 1] if kind == "10" else [1, 2, 3]
    if kind == "11":
        r((2,), op(T, 2)); r((3,), op(T, 3)); psi = op(S.conj().T, 2) @ op(S.conj().T, 3) @ psi
    out = {}
    for rec in itertools.product("+-", repeat=len(meas)):
        v = psi.copy()
        for q, b in zip(meas, rec):
            v = op(Pp if b == "+" else Pm, q) @ v
        out["".join(rec)] = (float(np.linalg.norm(v) ** 2), v)
    return out, len(inputs)


def fid(kind, v):
    """Fidelity of the (normalised) accepted state with the stated output on the unmeasured qubits."""
    v = v / np.linalg.norm(v); psi = v.reshape([2] * N)
    if kind == "10":                                      # outputs q2 q3, q0 q1 projected on |+>
        p = np.array([1, 1]) / np.sqrt(2); u = np.einsum("a,b,abcd->cd", p, p, psi).reshape(4)
        want = np.kron(Tp, Tp)
    else:                                                 # output q0, q1 q2 q3 projected on |+>
        p = np.array([1, 1]) / np.sqrt(2); u = np.einsum("b,c,d,abcd->a", p, p, p, psi)
        u = u / np.linalg.norm(u)                           # |T> up to the one single-qubit Clifford C0 fixed by
        return float(abs(np.vdot(C0[0] @ Tp, u)) ** 2)      # the noiseless run (a faulted run is compared with it)
    u = u / np.linalg.norm(u)
    return float(abs(np.vdot(want, u)) ** 2)


C0 = [I2]
res = {}
for kind, label, ninp in (("10", "[[10,2,2]]", 10), ("11", "[[11,1,2]]", 11)):
    ok_rec = "+" * (2 if kind == "10" else 3)
    out, n = run(kind)
    p_acc, v = out[ok_rec]
    if kind == "11":                                      # fix C0 once, from the noiseless output
        w = v / np.linalg.norm(v); u = np.einsum("b,c,d,abcd->a", *[np.array([1, 1]) / np.sqrt(2)] * 3, w.reshape([2] * N))
        u = u / np.linalg.norm(u); cl = [I2, S, S.conj().T, Z, X, X @ S, X @ S.conj().T, X @ Z]
        C0[0] = max(cl, key=lambda C: abs(np.vdot(C @ Tp, u)))
    rec = dict(inputs=n, accept_prob_noiseless=round(p_acc, 12), fidelity_noiseless=round(fid(kind, v), 12))
    caught, harmless, bad = [], [], []
    for i in range(n):
        o, _ = run(kind, i)
        pa, vv = o[ok_rec]
        if pa < 1e-9:
            caught.append(i)
        elif fid(kind, vv) > 1 - 1e-9:
            harmless.append(i)
        else:
            bad.append(i)
    pairs_bad = []                                        # fault pairs that pass the trivial record and harm the output
    for i, j in itertools.combinations(range(n), 2):
        o, _ = run(kind, (i, j)); pa, vv = o[ok_rec]
        if pa > 1e-9 and fid(kind, vv) < 1 - 1e-9:
            pairs_bad.append([i, j])
    rec.update(single_faults=n, caught=len(caught), harmless=len(harmless), undetected_harmful=len(bad),
               rejection_first_order=f"{len(caught)}p", fault_pairs=n * (n - 1) // 2,
               pairs_accepted_and_harmful=pairs_bad)
    if kind == "10":                                      # p^2 coefficients of the output error, weighted by acceptance
        p = np.array([1, 1]) / np.sqrt(2); e2 = e3 = ej = 0.0; w1 = w12 = 0
        for i, j in pairs_bad:
            o, _ = run(kind, (i, j)); pa, vv = o[ok_rec]
            u = np.einsum("a,b,abcd->cd", p, p, (vv / np.linalg.norm(vv)).reshape([2] * N)).reshape(4)
            u = u / np.linalg.norm(u); M = u.reshape(2, 2)
            e2 += pa * (1 - np.real(Tp.conj() @ (M @ M.conj().T) @ Tp))
            e3 += pa * (1 - np.real(Tp.conj() @ (M.T @ M.conj()) @ Tp))
            ej += pa * (1 - abs(np.vdot(np.kron(Tp, Tp), u)) ** 2)
            w1 += abs(pa - 1) < 1e-9; w12 += abs(pa - 0.5) < 1e-9
        rec.update(output_error_p2=dict(output_q2=round(float(e2), 9), output_q3=round(float(e3), 9), joint=round(float(ej), 9),
                                        pairs_accepted_prob_1=int(w1), pairs_accepted_prob_half=int(w12)))
        assert round(e2, 6) == 9 and round(e3, 6) == 9 and round(ej, 6) == 13, (e2, e3, ej)
    res[label] = rec
    print(label, rec, flush=True)
os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "catalytic_check.json"), "w"), indent=1)
