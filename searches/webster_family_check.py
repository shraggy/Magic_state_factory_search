"""All-weights family [[2^(l+1)-2,1,2]] against the embedded repetition codes of Webster, Quintavalle and Bartlett,
New J. Phys. 25, 103018 (2023), Example 4.1 (repetition code on d bits, embedding V = M^t_d: all vectors of weight
1..t; their text: 'the embedded code has a transversal logical phase gate at level t').
For t = l = 2, 3, 4 and d = t+1, map each embedded qubit (vector v of weight 1..t) to its column in check/logical
coordinates (v1+v2, ..., v_{d-1}+v_d, v_d) and compare with our all-weights rotation set on N = l+1 qubits with the
output singleton removed. Also checks that [[14,1,2]] is Bravyi-Haah [[14,2,2]] with output 2 deleted.
    python searches/webster_family_check.py -> outputs/checks/webster_family_check.json
"""
import json, os, itertools, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
res = {}
for t in (2, 3, 4):
    d = t + 1; cols = set()
    for v in range(1, 2 ** d):
        if bin(v).count("1") > t: continue
        b = [(v >> i) & 1 for i in range(d)]
        u = [b[i] ^ b[i + 1] for i in range(d - 1)] + [b[d - 1]]
        cols.add(sum(x << i for i, x in enumerate(u)))
    ours = set(range(1, 2 ** d)) - {1 << (d - 1)}
    res[f"l={t}"] = dict(n=len(cols), equal=cols == ours)
sb = json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))
R = {x["spec"]: x for x in sb["rows"] if x["l"] == 3}
g = collections.Counter()
for S, c in R["T0,T1"]["gates"]:
    T = tuple(q if q < 1 else q - 1 for q in S if q != 1)
    if T: g[T] += c
g = {S for S, c in g.items() if c % 8}
tb = {tuple(S) for S, c in R["T0"]["gates"]}
res["bh_14_2_2_one_output_kept_equals_14_1_2"] = any(
    {tuple(sorted(p[q] for q in S)) for S in g} == tb for p in itertools.permutations(range(4)))
json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "webster_family_check.json"), "w"), indent=1)
print(res)
