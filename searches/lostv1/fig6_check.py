"""Checks behind the appendix figure of the odd-weight parent (v1 Fig. 6 redrawn; figure parent_odd).

    nice -n 19 python3 searches/lostv1/fig6_check.py   -> outputs/checks/fig6_check.json   (seconds, one core)

Every statement drawn in the figure, by exact state-vector simulation (numpy, own code) and, for the cuts, also with
unified_search.is_factory on the full 2^N phase table (brute_force=True) and unified_search.classify_factory:
 1. the five-qubit odd-weight parent (16 rotations of odd weight, pi/8 each) is the identity on every basis state;
 2. its four cuts are [[8,4,2]], [[12,3,2]], [[14,2,2]], [[15,1,3]]: checks return to |+> with probability 1;
 3. the [[12,3,2]] output equals T^dag_2 T^dag_3 T^dag_4 T^dag_{234}|+++> exactly (four rotations, parities of rank 3),
    and its T-count is 4 (outputs/staged_best.json row CS01,CS02, tau=4; and unified_search.zcount);
 4. one catalytic T rotation on the parity q2q3q4 turns it into (T^dag|+>)^{x3}, and S on each wire into |T>^{x3},
    exactly (also checked: T on q2, CNOT(q3->q2), CNOT(q4->q2), S on each wire);
 5. the [[8,4,2]] output: sqrt(X)^dag on q1, then C_1(X_2X_3X_4), then an X measurement of q1 gives the + outcome with
    probability 1 and leaves exactly the [[12,3,2]] output of 3 (all 24 single-qubit Cliffords on q1 also listed);
 6. the six-qubit circuit of the first submission: not the identity (its phase takes several values), and its cuts after
    q0 and after q0q1q2 are an [[8,5,2]] and a [[12,3,2]] whose output equals the one in 3 up to a global phase.
 7. exact outputs drawn in the redrawn figure (Fig. S1): every cut of the odd-weight parent leaves exactly the inverse
    of the rotations to its right; the six-qubit [[8,5,2]] output is exactly T_12 T_34 T_45 T_35 T_1234 T_1235 T_1245
    |+>^5 (wires q1..q5), its [[12,3,2]] output exactly T^dag_3 T^dag_4 T^dag_5 T^dag_345 |+>^3; and the chain of
    panel (c), on the wires q1..q5 of the six-qubit circuit: CNOT(q1->q2) CNOT(q3->q2) maps the [[8,5,2]] output to
    |+>_1 times the [[8,4,2]] output of the odd-weight parent on q2..q5; sqrt(X)^dag_2, C_2(X_3X_4X_5) and an X
    measurement of q2 (outcome + with probability 1) leave the [[12,3,2]] output on q3q4q5; T_345 and S_3 S_4 S_5 leave
    |T>^3. Also: CNOT(q1->q2) CNOT(q3->q4) CNOT(q3->q5) leave |+>_1|+>_3 times CCZ CZ CZ CZ Z Z Z|+++> on q2q4q5.
"""
import itertools, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); SE = os.path.dirname(HERE); ROOT = os.path.dirname(SE)
sys.path.insert(0, SE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us

W = np.exp(1j * np.pi / 4)
def bits(N):
    x = np.arange(2 ** N); return [(x >> (N - 1 - q)) & 1 for q in range(N)]   # q0 is the most significant bit
def diag(N, rots, sgn=1):
    b = bits(N); ph = np.zeros(2 ** N)
    for R in rots: ph = ph + sum(b[q] for q in R) % 2
    return W ** (sgn * ph)
def plus(N): return np.ones(2 ** N, complex) / np.sqrt(2 ** N)
def project_plus(psi, N, checks):
    """Project checks (a prefix q0..q_{j-1}) on |+>; return acceptance probability and normalised output."""
    j = len(checks); assert checks == list(range(j))
    out = psi.reshape(2 ** j, -1).sum(axis=0) / np.sqrt(2 ** j)
    p = np.vdot(out, out).real; return p, out / np.sqrt(p)
def same(a, b):                                             # equal up to a global phase
    return bool(abs(abs(np.vdot(a, b)) - 1) < 1e-9)
def kron(*m):
    r = np.eye(1)
    for x in m: r = np.kron(r, x)
    return r
I2 = np.eye(2); H = np.array([[1, 1], [1, -1]]) / np.sqrt(2); S = np.diag([1, 1j]); T = np.diag([1, W])
def one(g, q, n): return kron(*[g if i == q else I2 for i in range(n)])
def cnot(c, t, n):
    M = np.zeros((2 ** n, 2 ** n)); b = bits(n)
    for x in range(2 ** n):
        y = x ^ ((b[c][x]) << (n - 1 - t)); M[y, x] = 1
    return M
tstate = T @ np.array([1, 1]) / np.sqrt(2)
out = {}

# 1-2. five-qubit odd-weight parent
five = [S_ for w in (5, 3, 1) for S_ in itertools.combinations(range(5), w)]
out["odd_parent_rotations"] = len(five)
out["odd_parent_distinct_phases"] = int(len(set(np.round(np.angle(diag(5, five)), 9) % (2 * np.pi))))
cuts = {}
for j in range(1, 5):
    ch = list(range(j)); keep = [R for R in five if set(R) & set(ch)]
    p, o = project_plus(diag(5, keep) * plus(5), 5, ch)
    outs = list(range(j, 5)); m = {q: i for i, q in enumerate(outs + ch)}
    f = dict(k=len(outs), r=j, gates=[(frozenset(m[q] for q in R), 1) for R in keep], stage="x", params="")
    info = us.classify_factory(f, 3)
    cuts[f"[[{len(keep)},{len(outs)}]]"] = dict(accept_prob=round(p, 12), is_factory_brute=us.is_factory(f, 3, brute_force=True),
                                               d=info["d"], reduced=info["reduced"], zcount=info["zcount"])
    if j == 2: out12 = o
out["odd_parent_cuts"] = cuts

# 3. the [[12,3,2]] output, exactly
target12 = diag(3, [(0,), (1,), (2,), (0, 1, 2)], sgn=-1) * plus(3)     # wires q2 q3 q4 -> 0 1 2
out["out12_equals_Tdag2_Tdag3_Tdag4_Tdag234"] = same(out12, target12)
sb = json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))
row = [r for r in sb["rows"] if r["l"] == 3 and r["spec"] == "CS01,CS02"]
out["staged_best_CS01,CS02"] = dict(n=row[0]["n"], tau=row[0]["tau"], k=row[0]["k"]) if row else None

# 4. catalytic split: T_2 (one |T> consumed), CNOT(q3->q2), CNOT(q4->q2), S on each wire
U = kron(S, S, S) @ cnot(2, 0, 3) @ cnot(1, 0, 3) @ one(T, 0, 3)
out["catalytic_split_gives_T^3"] = same(U @ out12, kron(tstate, tstate, tstate).ravel())
# the form drawn: one T rotation on the parity q2q3q4 (consumes one |T>), then S on each wire
U2 = kron(S, S, S) @ np.diag(diag(3, [(0, 1, 2)]))
out["parity_T234_then_S_gives_T^3"] = same(U2 @ out12, kron(tstate, tstate, tstate).ravel())
out["parity_T234_alone_gives_Tdag^3"] = same(np.diag(diag(3, [(0, 1, 2)])) @ out12, kron(tstate.conj(), tstate.conj(), tstate.conj()).ravel())
# without the catalytic T no CNOT circuit and diagonal Clifford layer gives |T>^3 (T-count 4 > 3): recorded from item 3

# 5. [[8,4,2]] output -> [[12,3,2]] output by a Clifford on q1, C_1(X2X3X4), X measurement of q1
p8, out8 = project_plus(diag(5, [R for R in five if 0 in R]) * plus(5), 5, [0])     # wires q1..q4
cl = {}
gens = {"H": H, "S": S}
words = [""]
for L in range(1, 9): words += ["".join(w) for w in itertools.product("HS", repeat=L)]
seen = []
for wd in words:
    g = np.eye(2)
    for ch in wd: g = gens[ch] @ g
    if any(abs(abs(np.trace(g.conj().T @ h)) - 2) < 1e-9 for h in seen): continue
    seen.append(g)
    fan = cnot(0, 3, 4) @ cnot(0, 2, 4) @ cnot(0, 1, 4)
    st = fan @ one(g, 0, 4) @ out8
    st = one(H, 0, 4) @ st                                     # X measurement of q1: H then Z, accept 0 (= |+>)
    r = st.reshape(2, 8)[0]; pr = np.vdot(r, r).real
    if pr > 1 - 1e-9:
        r = r / np.sqrt(pr)
        cl[wd or "I"] = dict(accept=round(pr, 12), equals_out12=same(r, target12), equals_conj_out12=same(r, target12.conj()))
out["n_single_qubit_cliffords"] = len(seen)
out["842_reduction_deterministic"] = cl
SX = np.array([[1 + 1j, 1 - 1j], [1 - 1j, 1 + 1j]]) / 2                      # sqrt(X)
named = {}
for nm, g in (("sqrtX", SX), ("sqrtX_dag", SX.conj().T)):
    st = one(H, 0, 4) @ cnot(0, 3, 4) @ cnot(0, 2, 4) @ cnot(0, 1, 4) @ one(g, 0, 4) @ out8
    r = st.reshape(2, 8)[0]; pr = np.vdot(r, r).real; r = r / np.sqrt(pr)
    named[nm] = dict(accept=round(pr, 12), equals_out12=same(r, target12),
                     then_split_gives_T3=same(U2 @ r, kron(tstate, tstate, tstate).ravel()))
out["842_to_out12_named"] = named

# 6. six-qubit circuit of the first submission (v1 Fig. 6(c))
six = [(0, 1, 3, 4, 5), (0, 2, 3, 4, 5), (0, 1, 5), (0, 2, 5), (0, 1, 4), (0, 2, 4), (0, 1, 3), (0, 2, 3),
       (1, 2, 3, 4, 5), (1, 2, 5), (1, 2, 4), (1, 2, 3), (3, 4, 5)]
out["six_distinct_phases"] = int(len(set(np.round(np.angle(diag(6, six)), 9) % (2 * np.pi))))
out["six_full_is_BI_brute"] = us.is_factory(dict(k=0, r=6, gates=[(frozenset(R), 1) for R in six]), 3, brute_force=True)
sc = {}
for ch in ([0], [0, 1, 2]):
    keep = [R for R in six if set(R) & set(ch)]
    p, o = project_plus(diag(6, keep) * plus(6), 6, ch)
    outs = [q for q in range(6) if q not in ch]; m = {q: i for i, q in enumerate(outs + ch)}
    f = dict(k=len(outs), r=len(ch), gates=[(frozenset(m[q] for q in R), 1) for R in keep], stage="x", params="")
    info = us.classify_factory(f, 3)
    e = dict(n=len(keep), k=len(outs), accept_prob=round(p, 12), is_factory_brute=us.is_factory(f, 3, brute_force=True),
             d=info["d"], reduced=info["reduced"])
    if len(ch) == 3:
        e["output_equals_odd_parent_12_3_2_output"] = same(o, out12)
    sc[str(ch)] = e
out["six_cuts"] = sc

# 7. exact outputs drawn in the figure and the conversion chain of panel (c)
ex = {}
for j in range(1, 5):                                     # odd-weight parent: output = inverse of the rotations right of the cut
    ch = list(range(j)); keep = [R for R in five if set(R) & set(ch)]; right = [R for R in five if not set(R) & set(ch)]
    _, o = project_plus(diag(5, keep) * plus(5), 5, ch)
    rr = [tuple(q - j for q in R) for R in right]
    ex[f"odd_cut{j}_equals_inverse_of_rotations_right"] = same(o, diag(5 - j, rr, sgn=-1) * plus(5 - j))
_, o852 = project_plus(diag(6, [R for R in six if 0 in R]) * plus(6), 6, [0])      # wires q1..q5 -> 0..4
_, o1232 = project_plus(diag(6, [R for R in six if set(R) & {0, 1, 2}]) * plus(6), 6, [0, 1, 2])   # q3q4q5
P852 = [(0, 1), (2, 3), (3, 4), (2, 4), (0, 1, 2, 3), (0, 1, 2, 4), (0, 1, 3, 4)]  # 12,34,45,35,1234,1235,1245
ex["six_852_equals_T12_T34_T45_T35_T1234_T1235_T1245"] = same(o852, diag(5, P852) * plus(5))
ex["six_1232_equals_Tdag3_Tdag4_Tdag5_Tdag345"] = same(o1232, target12)
st = cnot(2, 1, 5) @ cnot(0, 1, 5) @ o852                   # CNOT(q1->q2), CNOT(q3->q2)
ex["852_cnots_q1q2_q3q2_give_plus1_x_842_output_of_odd_parent"] = same(st, np.kron(plus(1), out8))
st2 = cnot(2, 4, 5) @ cnot(2, 3, 5) @ cnot(0, 1, 5) @ o852
ccz_like = diag(3, [(0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)]) * plus(3)   # (-1)^{OR}: CCZ CZ CZ CZ Z Z Z
# direct check on q2q4q5 (q1 = wire 0, q3 = wire 2 must be |+>): reorder wires to (q1,q3,q2,q4,q5)
t = st2.reshape(2, 2, 2, 2, 2).transpose(0, 2, 1, 3, 4).ravel()
ex["852_cnots_q1q2_q3q4_q3q5_give_plus1_plus3_x_CCZclass_on_q2q4q5"] = same(t, np.kron(plus(2), ccz_like))
ccz_exact = np.ones(8, complex); ccz_exact[7] = -1; ccz_exact /= np.sqrt(8)
b3 = bits(3); czzz = np.array([(-1) ** (b3[0][x] * b3[1][x] + b3[0][x] * b3[2][x] + b3[1][x] * b3[2][x] + b3[0][x] + b3[1][x] + b3[2][x]) for x in range(8)])
ex["CCZclass_equals_CCZ_CZ12_CZ13_CZ23_Z1_Z2_Z3"] = same(ccz_like, czzz * ccz_exact)
SXd = SX.conj().T
st3 = one(H, 1, 5) @ cnot(1, 4, 5) @ cnot(1, 3, 5) @ cnot(1, 2, 5) @ one(SXd, 1, 5) @ st   # X measurement of q2
r = st3.reshape(2, 2, 8)                                 # (q1, q2, q3q4q5)
pr = np.vdot(r[:, 0, :], r[:, 0, :]).real
ex["chain_Xmeas_q2_accept_prob"] = round(pr, 12)
r0 = r[:, 0, :].ravel() / np.sqrt(pr)
ex["chain_after_Xmeas_equals_plus1_x_1232_output"] = same(r0, np.kron(plus(1), target12))
fin = U2 @ (r[:, 0, :] / np.sqrt(pr))[0] * np.sqrt(2)       # q1 is |+>: take its 0 component (times sqrt 2)
ex["chain_T345_then_S_gives_T^3"] = same(fin, kron(tstate, tstate, tstate).ravel())
out["exact_outputs_and_chain"] = ex
for k_, v in out.items(): print(k_, v)
json.dump(out, open(os.path.join(ROOT, "outputs", "checks", "fig6_check.json"), "w"), indent=1, default=str)
