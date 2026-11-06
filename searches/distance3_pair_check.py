"""Distance beyond 2: a factory has d >= 3 iff no single rotation and no pair of rotations is an undetected fault that
changes the output. A faulty rotation flips the checks in its check part S \\ O and acts on the outputs by its output
part S & O, so a pair is undetected iff the two check parts are equal, and harmful iff the output parts differ.
Checked here against unified_search.distance on [[15,1,3]], [[14,2,2]] and random distance-2 factories.
    python searches/distance3_pair_check.py
"""
import sys, os, itertools, json, csv, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search as us
csv.field_size_limit(sys.maxsize)

def pair_ok(gates, k):
    out = set(range(k)); seen = {}
    for S, c in gates:
        if us.inputs(c, 3) == 0: continue
        chk = frozenset(S - out); op = frozenset(S & out)
        if not chk: return False
        if chk in seen and seen[chk] != op: return False
        seen.setdefault(chk, op)
    return True

def rm(N, out, removed):
    g = [(frozenset(S), 1) for w in range(1, N + 1, 2) for S in itertools.combinations(range(N), w)]
    return [(S, c) for S, c in g if S not in removed]

# Reed-Muller l=3 on N=5 odd weights; outputs first
f15 = rm(5, {0}, {frozenset({0})})                      # [[15,1,3]]: output q0
f14 = rm(5, {0, 1}, {frozenset({0}), frozenset({1})})   # [[14,2,2]]: outputs q0 q1
for name, g, k in (("[[15,1,3]]", f15, 1), ("[[14,2,2]]", f14, 2)):
    f = dict(stage="x", k=k, r=5 - k, gates=g, params="")
    print(name, "n", len(g), "is_factory", us.is_factory(f, 3), "distance", us.distance(f, 3), "pair condition", pair_ok(g, k))
R = list(csv.DictReader(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs", "l3_catalogue_k6.csv"))))
agree = tot = 0
for r in random.Random(3).sample(R, 80):
    k = int(r["k"]); g = [(frozenset(S), c) for S, c in json.loads(r["gates"])]
    f = dict(stage="x", k=k, r=int(r["N"]) - k, gates=g, params="")
    d = us.distance(f, 3); agree += (d >= 3) == pair_ok(g, k); tot += 1
print(f"catalogue factories: pair condition <=> d >= 3 on {agree} of {tot}")
