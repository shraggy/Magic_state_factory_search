"""Third-level targets inside the l=4 class lists (computed check, rerunnable).

For every l=3 class with k<=4 in outputs/staged_best.json, double each rotation coefficient (c -> 2c mod 16), so that
each T rotation becomes two sqrt(T) rotations, verify the circuit as an l=4 factory (unified_search.is_factory,
distance 2), label it with the exact l=4 class labeller, and compare its sqrt(T) count with the l=4 best of that class.
Also counts the l=4 classes whose best factory uses only even coefficients.
    nice -n 19 python searches/doubled_l3_at_l4.py  -> outputs/checks/doubled_l3_at_l4.json
"""
import json, sys, os
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification")); os.chdir(ROOT)
import unified_search as us, cascade_class as CC
d = json.load(open("outputs/staged_best.json"))
best4 = {(x["k"], x["orbit"]): x for x in d["rows"] if x["l"] == 4}
even4 = [(x["k"], x["orbit"]) for x in d["rows"] if x["l"] == 4 and all(c % 2 == 0 for _, c in x["gates"])]
rows = []
for k in range(1, 5):
    src = {1: "staged_small", 2: "staged_small", 3: "staged_l4fix", 4: "staged_A"}[k]
    cache = next(c for c in (f"outputs/{src}/_cache", "outputs/staged_A/_cache") if os.path.exists(f"{c}/l4_k{k}_classes.pkl"))
    classes, ctx = CC.classes_of(4, k, cache); lab = CC.labeller(4, k, ctx, classes)
    for x in [x for x in d["rows"] if x["l"] == 3 and x["k"] == k]:
        g = [(frozenset(S), 2 * c % 16) for S, c in x["gates"]]; f = dict(k=k, r=x["N"] - k, gates=g)
        ok = us.is_factory(f, 4) and us.distance(f, 4) == 2; n = sum(us.inputs(c, 4) for _, c in g)
        o = lab(f); b = best4[(k, o)]
        rows.append(dict(k=k, l3_spec=x["spec"], l3_n=x["n"], l4_orbit=o, l4_spec=b["spec"], l4_best_n=b["n"],
                         doubled_n=n, doubled_valid=ok, l4_best_even=(k, o) in even4))
s = dict(l3_classes_k_le_4=len(rows), l4_classes=len(best4), l4_best_even_only=len(even4),
         distinct_l4_classes_hit=len({(r["k"], r["l4_orbit"]) for r in rows}),
         all_doubled_valid=all(r["doubled_valid"] for r in rows),
         hit_classes_are_the_even_ones=sorted((r["k"], r["l4_orbit"]) for r in rows) == sorted(even4),
         doubled_equal=sum(r["doubled_n"] == r["l4_best_n"] for r in rows),
         doubled_fewer=[(r["l3_spec"], r["doubled_n"], r["l4_best_n"]) for r in rows if r["doubled_n"] < r["l4_best_n"]],
         doubled_more=sum(r["doubled_n"] > r["l4_best_n"] for r in rows))
json.dump(dict(source="outputs/staged_best.json", summary=s, rows=rows), open("outputs/checks/doubled_l3_at_l4.json", "w"), indent=1)
print(json.dumps(s))
