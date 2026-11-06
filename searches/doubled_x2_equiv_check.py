"""Are the four l=4 rows marked x2 (T_1, CS_12, T_{1,2}, CCZ_123) the doubled l=3 factories of the same class?
For each pair, search invertible GF(2) maps M of the N qubits that fix every output qubit (M e_q = e_q for q < k) and
send each check column to a vector with a nonzero check part, such that M maps the l=3 supports onto the l=4 supports;
all l=3 coefficients are 1 and all l=4 coefficients 2, so the doubled circuit then agrees with the l=4 one rotation for
rotation after a CNOT circuit on the checks. 'identity_works' records whether no relabelling is needed at all.
    nice -n 19 python searches/doubled_x2_equiv_check.py -> outputs/checks/doubled_x2_equiv.json
"""
import json, os, itertools
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = {(r["l"], r["spec"]): r for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]}
pairs = {"T1": ("T0", "sT0,sT0"), "CS12": ("CS01", "CT01,CT01"), "T12": ("T0,T1", "sT0,sT0,sT1,sT1"),
         "CCZ123": ("CCZ012", "CCS012,CCS012")}
v = lambda S: sum(1 << q for q in S)
def app(M, x):
    y = 0
    for q in range(len(M)):
        if x >> q & 1: y ^= M[q]
    return y
res = {}
for name, (s3, s4) in pairs.items():
    a, b = R[(3, s3)], R[(4, s4)]; k, N = a["k"], a["N"]
    A = [v(S) for S, c in a["gates"]]; B = {v(S) for S, c in b["gates"]}
    coef = all(c == 1 for _, c in a["gates"]) and all(c == 2 for _, c in b["gates"])
    same_N = b["N"] == N and len(A) == len(B)
    maps = 0
    if same_N:
        for cols in itertools.product(range(1 << k, 1 << N), repeat=N - k):
            M = [1 << q for q in range(k)] + list(cols)
            if {app(M, x) for x in A} == B and len({app(M, x) for x in A}) == len(A): maps += 1
    res[name] = dict(l3=(a["n"], N, a["stage"]), l4=(b["n"], b["N"], b["stage"]), coefficients_1_and_2=coef,
                     identity_works=set(A) == B, maps_found=maps, equivalent=coef and maps > 0)
json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "doubled_x2_equiv.json"), "w"), indent=1)
print(json.dumps(res, indent=1))
