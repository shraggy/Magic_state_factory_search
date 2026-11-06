"""Distance of two-group circuits from the block symmetry.

In a two-group circuit every w_S-subset C of the checks carries, for each admissible pair (w_O, w_S), all C(k, w_O)
output subsets. By the pair condition (Discussion; searches/distance3_pair_check.py) a factory has d >= 3 iff no two
input rotations share a check part while differing on the outputs. So a two-group factory has d >= 3 iff for every w_S
the admissible output weights give at most one output subset in total, sum_{w_O: (w_O,w_S) in W, w_S >= 1} C(k, w_O) <= 1:
each check weight is paired with the empty output set or, at k = 1, with the single output. In particular d = 2 for
every two-group factory with k >= 2 that has an input rotation touching the outputs; d >= 3 needs k = 1 (or no coupling).
This script checks the prediction against unified_search.distance on every two-group factory of the main sweep.
    python searches/two_group_distance.py     -> outputs/checks/two_group_distance.json
"""
import os, sys, json
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search as us
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def predict_d3(k, gates):
    out = set(range(k)); byS = {}
    for S, c in gates:
        if not us.inputs(c, 3):
            continue
        wS = len(S - out); wO = len(S & out)
        byS.setdefault(wS, set()).add(wO)
    if any(wS == 0 for wS in byS):
        return False                                        # an input rotation on the outputs alone: d = 1
    return all(sum(comb(k, a) for a in A) <= 1 for A in byS.values())


res = dict(checked=0, agree=0, d3=0, d3_k=[], mismatches=[])
for L in (2, 3, 4):
    us.S_MAX = 7
    for f in us.stage_two_group(L, 6, 9, 1, one_sign=True):
        if not us.is_factory(f, L):
            continue
        d = us.distance(f, L)
        if d == 1:
            continue
        p = predict_d3(f["k"], f["gates"]); real = not (d == 2)
        res["checked"] += 1; res["agree"] += p == real
        if real:
            res["d3"] += 1; res["d3_k"].append(f["k"])
        if p != real and len(res["mismatches"]) < 10:
            res["mismatches"].append(dict(l=L, k=f["k"], r=f["r"], d=str(d), params=f.get("params", "")))
res["d3_k"] = sorted(set(res["d3_k"]))
os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "two_group_distance.json"), "w"), indent=1)
print({k: v for k, v in res.items() if k != "mismatches"}, res["mismatches"][:3])
