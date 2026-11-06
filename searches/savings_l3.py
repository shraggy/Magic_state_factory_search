"""Relative size of the saving over synthillation at l=3 (k<=6), from outputs/staged_best.json (n, syn_n), with
syn_n checked against best_n of outputs/synthillation_strongest_l3.json.
    python searches/savings_l3.py -> outputs/checks/savings_l3.json
"""
import json, os, statistics as st
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
r = [x for x in json.load(open(os.path.join(O, "staged_best.json")))["rows"] if x["l"] == 3]
s = {(x["k"], x["spec"]): x["best_n"] for x in json.load(open(os.path.join(O, "synthillation_strongest_l3.json")))}
assert all(s[(x["k"], x["spec"])] == x["syn_n"] for x in r)
f = [x for x in r if x["n"] < x["syn_n"]]; pc = [(x["syn_n"] - x["n"]) / x["syn_n"] for x in f]
out = dict(classes=len(r), fewer=len(f), equal=sum(x["n"] == x["syn_n"] for x in r), more=sum(x["n"] > x["syn_n"] for x in r),
           syn_inputs_total=sum(x["syn_n"] for x in r), search_inputs_total=sum(x["n"] for x in r),
           saving_total=sum(x["syn_n"] - x["n"] for x in f), saving_pct_of_total=100 * sum(x["syn_n"] - x["n"] for x in f) / sum(x["syn_n"] for x in r),
           per_class_pct_min=100 * min(pc), per_class_pct_max=100 * max(pc), per_class_pct_median=100 * st.median(pc),
           fewer_at_k6=sum(x["k"] == 6 for x in f), fewer_at_k5=sum(x["k"] == 5 for x in f))
json.dump(out, open(os.path.join(O, "checks", "savings_l3.json"), "w"), indent=1); print(out)
