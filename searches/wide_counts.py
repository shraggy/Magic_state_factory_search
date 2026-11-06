"""Counts quoted in the paper that no other script writes: the wide-range sweep (two-group and disjoint/gluing circuits
over N<=18, k<=12 at l=2,3 and N<=16, k<=9 at l=4) from outputs/large_unified_factories.csv, and the n of the smallest
factory found for the k=8 CCZ-layer classes not on the bound, from outputs/k78_pipeline/full/l3_k8/iii_f.json.
A row is a Clifford-padded variant when essential_dim < k (App. "Clifford-padded factories").
    nice -n 19 python3 -I searches/wide_counts.py   -> outputs/checks/wide_counts.json"""
import csv, json, os, collections
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rows = list(csv.DictReader(open(os.path.join(R, "outputs/large_unified_factories.csv"))))
core = [r for r in rows if r["essential_dim"] == r["k"]]
out = dict(rows=len(rows), padded=len(rows) - len(core), per_level={})
for l in ("2", "3", "4"):
    S = [r for r in core if r["l"] == l]
    out["per_level"][l] = dict(
        circuits=len(S), d2=sum(r["d"] == "2" for r in S), d_ge3=sum(r["d"] != "2" for r in S),
        d_ge3_single_output=sum(r["d"] != "2" and r["k"] == "1" for r in S),
        by_stage=dict(collections.Counter(r["search_stage"] for r in S)),
        distinct_n_k_output=len({(r["n"], r["k"], r["reduced"]) for r in S}))
big = max(rows, key=lambda r: int(r["n"]))
out["largest"] = dict(l=int(big["l"]), n=int(big["n"]), k=int(big["k"]), d=big["d"], c=int(big["c"]))
d = json.load(open(os.path.join(R, "outputs/k78_pipeline/full/l3_k8/iii_f.json")))
off = [v["n"] for key, v in d["best"].items() if key not in d["on_bound_classes"]]
out["k8"] = dict(classes=len(d["best"]), on_bound=len(d["on_bound_classes"]),
                 n_hist_off_bound={str(n): c for n, c in sorted(collections.Counter(off).items())})
os.makedirs(os.path.join(R, "outputs/checks"), exist_ok=True)
json.dump(out, open(os.path.join(R, "outputs/checks/wide_counts.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
