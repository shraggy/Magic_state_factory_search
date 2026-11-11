"""Check that the [[10,2,2]] factory of App. S2.3 is a borrowed identity with its two output-only gates removed.

The noiseless [[10,2,2]] circuit F (searches/catalytic_check.py, measurements dropped) maps |+>^4 to
|+>_0 |+>_1 T|+>_2 T|+>_3, so B = T^dag_2 T^dag_3 F returns |+>^4 to itself: a borrowed identity with
12 non-Clifford rotations (the 7 of L_0, T_0, T^dag_1, T_1, T^dag_2, T^dag_3) and Clifford gates. B is neither
diagonal nor the identity operator, so it acts as the identity only on the state it returns. Removing T^dag_2 T^dag_3,
the only gates of B on the outputs q2 q3 alone, leaves F, whose output is their inverse, |T>^2, with 10 inputs.
Its distance 2 is the fault simulation of catalytic_check.py (10 of 10 single faults caught).
    python searches/bi_1022_check.py   -> outputs/checks/bi_1022_check.json
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import catalytic_check as cc

N = 4; op = cc.op; T = cc.T; S = cc.S; X = cc.X; Z = cc.Z; H = cc.H


def factory_unitary():
    """The [[10,2,2]] circuit of catalytic_check.run('10') as one 16x16 unitary, measurements dropped."""
    U = np.eye(2 ** N, dtype=complex)
    for s in cc.L0:
        U = cc.rot(s) @ U
    U = cc.rot((0,)) @ U                                         # T_0
    U = op(H @ S @ H, 1) @ U                                     # sqrt(X)_1
    X23 = op(X, 2) @ op(X, 3); Z23 = op(Z, 2) @ op(Z, 3); C = cc.ctrl(cc.Pm, 1, Z23)
    U = cc.ctrl(np.diag([0, 1]), 1, X23) @ U                     # C_1(X_2 X_3)
    U = C @ U; U = op(T.conj().T, 1) @ U; U = C @ U             # C^-_1(Z_2Z_3) T^dag_1 C^-_1(Z_2Z_3)
    U = op(T, 1) @ U; U = op(S.conj().T, 1) @ U                  # T_1 S^dag_1
    return U


plus = np.ones(2 ** N) / 4
F = factory_unitary()
B = op(T.conj().T, 2) @ op(T.conj().T, 3) @ F
Tp = cc.Tp
target = np.kron(np.kron(np.ones(2) / np.sqrt(2), np.ones(2) / np.sqrt(2)), np.kron(Tp, Tp))
res = dict(
    script="searches/bi_1022_check.py",
    bi_overlap_abs=float(abs(np.vdot(plus, B @ plus))),                   # 1 => B|+>^4 = e^{i phi}|+>^4
    factory_output_overlap_abs=float(abs(np.vdot(target, F @ plus))),     # 1 => F|+>^4 = |+>|+>|T>|T>
    B_is_identity_operator=bool(np.allclose(B, B[0, 0] * np.eye(2 ** N))),
    B_is_diagonal=bool(np.allclose(B, np.diag(np.diag(B)))),
    non_clifford_rotations_in_B=len(cc.L0) + 1 + 2 + 2,
    removed=["T^dag_2", "T^dag_3"],
    inputs_after_removal=len(cc.L0) + 1 + 2,
    distance_check="outputs/checks/catalytic_check.json: 10 of 10 single faults caught",
)
res["is_borrowed_identity"] = abs(res["bi_overlap_abs"] - 1) < 1e-12
res["output_is_T_T"] = abs(res["factory_output_overlap_abs"] - 1) < 1e-12
os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "bi_1022_check.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
