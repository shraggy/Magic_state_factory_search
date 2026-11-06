"""Corollary 3 (n >= tau(U) + 1) on every factory in the catalogues, with the count that meets it with equality.

Reads the l=3 per-class catalogue (outputs/l3_catalogue_k6.csv: 766 classes, exact T-count) and the larger sweep
(outputs/large_unified_factories.csv, rows with an exact Z(pi/2^l)-count). Reports, per file, how many factories are
checked, how many violate the bound (must be 0) and how many meet it with equality.
    python searches/catalogue_bound_check.py     -> outputs/checks/catalogue_bound_check.json
"""
import os, csv, sys, json
csv.field_size_limit(sys.maxsize)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
res = {}
R = list(csv.DictReader(open(os.path.join(O, "l3_catalogue_k6.csv"))))
res["l3_catalogue_k6"] = dict(checked=len(R), violate=sum(int(r["n"]) < int(r["tcount"]) + 1 for r in R),
                             equality=sum(int(r["n"]) == int(r["tcount"]) + 1 for r in R))
R = [r for r in csv.DictReader(open(os.path.join(O, "large_unified_factories.csv"))) if r["zcount"].strip().isdigit()]
res["large_unified_factories"] = dict(checked=len(R), violate=sum(int(r["n"]) < int(r["zcount"]) + 1 for r in R),
                                     equality=sum(int(r["n"]) == int(r["zcount"]) + 1 for r in R),
                                     by_level={l: sum(1 for r in R if r["l"] == l) for l in ("2", "3", "4")},
                                     by_distance={"2": sum(r["d"] == "2" for r in R), ">=3": sum(r["d"] != "2" for r in R)},
                                     equality_distances=sorted({r["d"] for r in R if int(r["n"]) == int(r["zcount"]) + 1}))
json.dump(res, open(os.path.join(O, "checks", "catalogue_bound_check.json"), "w"), indent=1)
print(json.dumps(res))
