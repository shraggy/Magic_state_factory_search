"""How far the catalogue factory lies below the cheapest of the three measured distil-then-synthesize routes
(outputs/checks/route_compare_l3.json: catalysed 4 tau, per-gate teleport, batched teleport), per l=3 class with k<=6.
Reports, overall, by k and on the 330 classes below synthillation: classes where the factory needs fewer inputs (N of M),
the margin in inputs (median, histogram) and as a share of the route, and the totals over all classes.
    python searches/route_margin.py   -> outputs/checks/route_margin.json
"""
import os, json, statistics as st, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")


def main():
    R = json.load(open(os.path.join(O, "checks", "route_compare_l3.json")))["rows"]
    best = lambda r: min(r["n_cat"], r["n_tel"], r["n_tel_batch"])
    def stats(S):
        pos = [best(r) - r["n"] for r in S if best(r) > r["n"]]
        rel = [(best(r) - r["n"]) / best(r) for r in S if best(r) > r["n"]]
        return dict(classes=len(S), factory_fewer=len(pos), median_margin=st.median(pos) if pos else None,
                    median_share_of_route=round(st.median(rel), 3) if rel else None,
                    margin_histogram=dict(sorted(collections.Counter(pos).items())),
                    median_margin_all=st.median(best(r) - r["n"] for r in S))
    W = [r for r in R if r["n"] < r["syn_n"]]
    res = dict(all=stats(R), below_synthillation=stats(W),
               by_k={k: stats([r for r in R if r["k"] == k]) for k in range(1, 7)},
               totals=dict(route=sum(best(r) for r in R), factory=sum(r["n"] for r in R),
                           cheaper_of_both=sum(min(r["n"], best(r)) for r in R)))
    t = res["totals"]; t["saving_vs_route"] = round(1 - t["cheaper_of_both"] / t["route"], 3)
    json.dump(res, open(os.path.join(O, "checks", "route_margin.json"), "w"), indent=1)
    print(json.dumps({k: res[k] for k in ("all", "below_synthillation", "totals")}, indent=1))


if __name__ == "__main__":
    main()
