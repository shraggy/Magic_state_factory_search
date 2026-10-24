"""Build and check the synthillation construction for an arbitrary l=3 target given as {frozenset(qubits): coeff mod 8}."""
import sys, itertools
import os
R = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, R)
import synthillation_compare as sc
import numpy as np
def rank2(M):
    M = M.copy() % 2; r = 0
    for c in range(M.shape[1]):
        p = [i for i in range(r, M.shape[0]) if M[i, c]]
        if not p: continue
        M[[r, p[0]]] = M[[p[0], r]]
        for i in range(M.shape[0]):
            if i != r and M[i, c]: M[i] ^= M[r]
        r += 1
    return r
def mu_of(k, target):
    Q = np.zeros((k, k), int)
    for T, c in target.items():
        T = sorted(T)
        if len(T) == 1 and c % 2: Q[T[0], T[0]] = 1
        if len(T) == 2 and (c // 2) % 2: Q[T[0], T[1]] = Q[T[1], T[0]] = 1
    r = rank2(Q)
    return 0 if r == 0 else (r if Q.diagonal().any() else r + 1)
def build(k, target, A):
    mu = mu_of(k, target); tau = len(A)
    lin = any(len(T) == 1 and c % 2 for T, c in target.items())
    if mu == 0:
        # cases 9-11 (B = 0): A with all-ones row; case 11 pads one zero column when tau is odd
        tm = [[1] + ([0] if tau % 2 else [])] if True else None
        blocks = [A] + ([[0]] if tau % 2 else [])
        cols = [v | (1 << k) for v in A] + ([1 << k] if tau % 2 else [])
        ok = sc.coeffs_ok(cols, k + 1, k, target, 3) and all(y >> k for y in cols)
        return dict(N=len(cols), n=k + 1, ok=ok, tau=tau, mu=0)
    B = sc.opt_B(k, target, mu)
    if B is None: return None
    if lin:
        tm = sc.C14[{(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}[(tau % 2, mu % 2)]]
        c = sum(1 << i for i in range(k) if target.get(frozenset([i]), 0) % 2); blocks = [A, B, B] + [[c]] * 4
    else:
        tm = sc.T5[(tau % 2, mu % 2)]; blocks = [A, B, B]
    blocks += [[0]] * (len(tm[0]) - len(blocks)); Rr = len(tm)
    cols = [v | sum(1 << (k + r) for r in range(Rr) if tm[r][bi]) for bi, blk in enumerate(blocks) for v in blk]
    ok = sc.coeffs_ok(cols, k + Rr, k, target, 3) and all(y >> k for y in cols)
    return dict(N=len(cols), n=k + Rr, ok=ok, tau=tau, mu=mu)
