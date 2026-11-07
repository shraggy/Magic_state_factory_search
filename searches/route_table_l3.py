"""Table of the factory against the three measured distil-then-synthesize routes (Table SIV, App. S4.3 of the Letter;
Table II, Sec. V E of the PRX Quantum version), for the 766 l=3 classes with k<=6, overall, on the 330 below
synthillation and by k. The paper prints the 'full_clifford' section; the top-level sections (each route on the diagonal
representative only, 659 / 47 / 60) are the weaker variant, kept as known answers.
Reads outputs/checks/route_kappa_l3.json (inputs, kappa, rejection per route; produced by searches/route_kappa_l3.py)
and outputs/checks/class_invariants_l3.json (product / entangled, single gate type). No new physics: it only counts.
'fewest' = the route with the fewest inputs ('best' in route_kappa_l3.json); 'lowkappa' = the route with the lowest
kappa at any input count ('min_kappa'). Per route (cat, tel, tel_batch) the input and kappa counts are also given.
Known-answer checks (asserted): the overall and 330 counts reproduce the summaries stored in route_kappa_l3.json and
outputs/checks/route_margin.json (659 / 47 / 60 inputs; 764 / 2 / 0 kappa; 756 / 3 / 7 against the lowest-kappa route;
696 / 30 / 40 rejection; 314 of 330), and the input counts per route reproduce outputs/checks/route_compare_l3.json.
Full-Clifford minimum ('full_clifford' section, added 2026-11-06): a route may use any Clifford, so for each class the
route is also taken as the cheapest (inputs, then kappa) over the members of its full-Clifford group
(outputs/checks/full_clifford_classes.json, 766 -> 622); a Clifford conversion is perfect in the noise model and
unitary, so kappa and rejection carry over unchanged. Known answer (asserted): CS01,CS02 -> 8 inputs, kappa 28, via the
CCZ012 member ([[8,3,2]]), since CS12 CS13|+++> = CX23 H2 CX12 S2 H2|CCZ> (arXiv:2610.31115 Eq. 3). The four-route
minimum also includes the lambda route of outputs/checks/lambda_route.json (route_n, 32 classes; inputs only).
    python searches/route_table_l3.py   -> outputs/checks/route_table_l3.json
"""
import os, json, collections, statistics as st
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); C = os.path.join(ROOT, "outputs", "checks")


def rel(a, b):
    return "factory_lower" if a < b else ("equal" if a == b else "factory_higher")


def tally(S):
    out = dict(classes=len(S))
    for key in ("best", "min_kappa", "cat", "tel", "tel_batch"):
        nm = {"best": "fewest", "min_kappa": "lowkappa"}.get(key, key)
        out[nm] = dict(inputs=dict(collections.Counter(rel(r["n"], r[key]["n"]) for r in S)),
                       kappa=dict(collections.Counter(rel(F(r["kappa_factory"]), F(r[key]["kappa"])) for r in S)),
                       rejection=dict(collections.Counter(rel(r["rej_factory"], r[key]["rej"]) for r in S)))
    rat = [F(r["kappa_factory"]) / F(r["best"]["kappa"]) for r in S]
    out["kappa_ratio_factory_over_fewest"] = dict(median=round(float(st.median(rat)), 3), min=round(float(min(rat)), 3),
                                                   max=round(float(max(rat)), 3))
    return out


def full_clifford(rows, inv):
    """Counts against the cheapest route over each full-Clifford group, and with the lambda route added."""
    G = json.load(open(os.path.join(C, "full_clifford_classes.json")))["groups"]
    by = {r["spec"]: r for r in rows}; grp = {}
    for g in G:
        mem = [x.rsplit(":", 1)[1] for x in g]; assert all(m in by for m in mem), g
        for m in mem: grp[m] = mem
    LR = {r["spec"]: r for r in json.load(open(os.path.join(C, "lambda_route.json")))["rows"]}
    best = lambda s, key: min((by[m][key] for m in grp.get(s, [s])), key=lambda b: (b["n"], F(b["kappa"])))
    lowk = lambda s: min((by[m]["min_kappa"] for m in grp.get(s, [s])), key=lambda b: (F(b["kappa"]), b["n"]))
    def lam(s, mem):
        v = [LR[m]["route_n"] for m in mem if m in LR and LR[m].get("route_n") is not None]
        return min(v) if v else None
    def n4(s, full):
        mem = grp.get(s, [s]) if full else [s]
        a = min(by[m]["best"]["n"] for m in mem); b = lam(s, mem)
        return a if b is None else min(a, b)
    def tal(S):
        bb = {r["spec"]: best(r["spec"], "best") for r in S}
        mar = [bb[r["spec"]]["n"] - r["n"] for r in S if r["n"] < bb[r["spec"]]["n"]]
        sh = [F(bb[r["spec"]]["n"] - r["n"], bb[r["spec"]]["n"]) for r in S if r["n"] < bb[r["spec"]]["n"]]
        return dict(classes=len(S),
                    inputs=dict(collections.Counter(rel(r["n"], bb[r["spec"]]["n"]) for r in S)),
                    kappa=dict(collections.Counter(rel(F(r["kappa_factory"]), F(bb[r["spec"]]["kappa"])) for r in S)),
                    kappa_lowest=dict(collections.Counter(rel(F(r["kappa_factory"]), F(lowk(r["spec"])["kappa"])) for r in S)),
                    rejection=dict(collections.Counter(rel(r["rej_factory"], bb[r["spec"]]["rej"]) for r in S)),
                    median_margin=float(st.median(mar)) if mar else None, margin_range=[min(mar), max(mar)] if mar else None,
                    median_share_of_route=round(float(st.median(sh)), 3) if sh else None,
                    kappa_ratio_median=round(float(st.median(F(r["kappa_factory"]) / F(bb[r["spec"]]["kappa"]) for r in S)), 3),
                    total_factory=sum(r["n"] for r in S), total_route=sum(bb[r["spec"]]["n"] for r in S),
                    inputs_four_routes_diagonal=dict(collections.Counter(rel(r["n"], n4(r["spec"], False)) for r in S)),
                    inputs_four_routes_full=dict(collections.Counter(rel(r["n"], n4(r["spec"], True)) for r in S)))
    W = [r for r in rows if r["n"] < r["syn_n"]]
    out = dict(source="outputs/checks/route_kappa_l3.json, outputs/checks/full_clifford_classes.json, outputs/checks/lambda_route.json",
               all=tal(rows), below_synthillation=tal(W), by_k={k: tal([r for r in rows if r["k"] == k]) for k in range(1, 7)})
    # per full-Clifford class (622): factory and route each the cheapest over the group
    seen, P = set(), []
    for r in rows:
        mem = tuple(sorted(grp.get(r["spec"], [r["spec"]])))
        if mem in seen: continue
        seen.add(mem); P.append((min(by[m]["n"] for m in mem), best(r["spec"], "best")["n"]))
    assert len(P) == 622
    out["per_full_clifford_class"] = dict(classes=len(P), inputs=dict(collections.Counter(rel(a, b) for a, b in P)))
    bb = {r["spec"]: best(r["spec"], "best") for r in rows}
    RF = [r for r in rows if bb[r["spec"]]["n"] < r["n"]]
    out["route_fewer_inputs"] = dict(classes=len(RF), by_k=dict(sorted(collections.Counter(r["k"] for r in RF).items())),
        max_tau=max(r["tau"] for r in RF), entangled=sum(not inv[r["spec"]]["product"] for r in RF),
        product=sum(inv[r["spec"]]["product"] for r in RF),
        margin_range=[min(r["n"] - bb[r["spec"]]["n"] for r in RF), max(r["n"] - bb[r["spec"]]["n"] for r in RF)],
        median_margin=float(st.median(r["n"] - bb[r["spec"]]["n"] for r in RF)),
        factory_kappa_lower=sum(F(r["kappa_factory"]) < F(bb[r["spec"]]["kappa"]) for r in RF),
        kappa_ratio_median=round(float(st.median(F(r["kappa_factory"]) / F(bb[r["spec"]]["kappa"]) for r in RF)), 3),
        new_by_merging=sorted(r["spec"] for r in RF if not r["best"]["n"] < r["n"]))
    out["route_lowered_by_merging"] = dict(
        all=sum(bb[r["spec"]]["n"] < r["best"]["n"] for r in rows),
        below_synthillation=sum(bb[r["spec"]]["n"] < r["best"]["n"] for r in W))
    out["lowered_by_lambda_route_after_merging"] = sum(n4(r["spec"], True) < bb[r["spec"]]["n"] for r in rows)
    out["k56_entangled"] = dict(collections.Counter(rel(r["n"], bb[r["spec"]]["n"]) for r in rows
                                                    if r["k"] in (5, 6) and not inv[r["spec"]]["product"]))
    out["k_le_4"] = dict(collections.Counter(rel(r["n"], bb[r["spec"]]["n"]) for r in rows if r["k"] <= 4))
    out["k56"] = dict(collections.Counter(rel(r["n"], bb[r["spec"]]["n"]) for r in rows if r["k"] in (5, 6)))
    dom, hk = [], []
    for r in rows:                                          # every route and every member of the full-Clifford group
        xs = [(m, by[m][key]) for m in grp.get(r["spec"], [r["spec"]]) for key in ("cat", "tel", "tel_batch", "min_kappa")]
        for m, x in xs:
            if x["n"] <= r["n"] and F(x["kappa"]) <= F(r["kappa_factory"]) and (x["n"] < r["n"] or F(x["kappa"]) < F(r["kappa_factory"])):
                dom.append((r["spec"], r["n"], r["kappa_factory"], m, x["construction"], x["n"], x["kappa"])); break
        m, x = min(xs, key=lambda t: (F(t[1]["kappa"]), t[1]["n"]))
        if F(x["kappa"]) < F(r["kappa_factory"]): hk.append((r["spec"], r["n"], r["kappa_factory"], m, x["construction"], x["n"], x["kappa"]))
    out["route_at_least_as_good_in_both"] = dom
    out["factory_higher_kappa_than_lowkappa_route"] = hk
    out["lowkappa_construction"] = dict(collections.Counter(lowk(r["spec"])["construction"] for r in rows))
    bm = {r["spec"]: min(((m, by[m]["best"]) for m in grp.get(r["spec"], [r["spec"]])), key=lambda t: (t[1]["n"], F(t[1]["kappa"])))
          for r in rows}
    out["kappa_carried"] = {nm: dict(collections.Counter(rel(F(r["kappa_factory"]), F(bm[r["spec"]][1]["kappa_carried"])) for r in S))
                            for nm, S in (("all", rows), ("below_synthillation", W))}
    out["second_order_acceptance"] = dict(max_ov_factory=max(r["ov_factory"] for r in rows),
                                          max_ov_fewest=float(max(F(bm[r["spec"]][1]["ov"]) for r in rows)),
                                          ov_vs_fewest=dict(collections.Counter(rel(F(r["ov_factory"]), F(bm[r["spec"]][1]["ov"])) for r in rows)))
    nb = [(r, by[bm[r["spec"]][0]]) for r in rows if bm[r["spec"]][1]["construction"] == "batch"
          and not by[bm[r["spec"]][0]]["batch_search"]["exhaustive"]]
    out["batch_kappa_sensitivity"] = dict(cheapest_is_nonexhaustive_batch=len(nb),
        max_kappa_ratio_factory_over_route_there=round(float(max(F(r["kappa_factory"]) / F(b["best"]["kappa"]) for r, b in nb)), 3))
    # known answer: CS01,CS02 is made by [[8,3,2]] and Cliffords (arXiv:2610.31115 Eq. 3)
    ka = best("CS01,CS02", "best"); assert (ka["n"], F(ka["kappa"])) == (8, 28), ka
    assert by["CS01,CS02"]["best"]["n"] == 14 and by["CS01,CS02"]["n"] == 12
    out["known_answer"] = dict(CS01_CS02=dict(route_n=ka["n"], route_kappa=ka["kappa"], factory_n=12,
                                              diagonal_route_n=by["CS01,CS02"]["best"]["n"], via="[[8,3,2]] for CCZ012"))
    return out


def main():
    R = json.load(open(os.path.join(C, "route_kappa_l3.json")))
    rows, summ = R["rows"], R["summary"]
    inv = {r["spec"]: r for r in json.load(open(os.path.join(C, "class_invariants_l3.json")))["rows"]}
    assert len(rows) == 766 and all(r["spec"] in inv for r in rows)
    W = [r for r in rows if r["n"] < r["syn_n"]]
    res = dict(source="outputs/checks/route_kappa_l3.json, outputs/checks/class_invariants_l3.json",
               all=tally(rows), below_synthillation=tally(W),
               by_k={k: tally([r for r in rows if r["k"] == k]) for k in range(1, 7)})
    # where a route needs fewer inputs than the factory
    RF = [r for r in rows if r["best"]["n"] < r["n"]]
    nT = lambda s: sum(g.startswith("T") for g in s.split(","))
    res["route_fewer_inputs"] = dict(
        classes=len(RF), by_k=dict(sorted(collections.Counter(r["k"] for r in RF).items())),
        max_tau=max(r["tau"] for r in RF), product=sum(inv[r["spec"]]["product"] for r in RF),
        entangled=sum(not inv[r["spec"]]["product"] for r in RF),
        single_gate_type=sum(bool(inv[r["spec"]]["single"]) for r in RF),
        with_T_gate_in_label=sum(nT(r["spec"]) > 0 for r in RF),
        factory_kappa_lower=sum(F(r["kappa_factory"]) < F(r["best"]["kappa"]) for r in RF),
        margin_histogram=dict(sorted(collections.Counter(r["n"] - r["best"]["n"] for r in RF).items())),
        classes_k_le_4=[(r["spec"], r["n"], r["best"]["n"], r["kappa_factory"], r["best"]["kappa"]) for r in RF if r["k"] <= 4])
    # entangled / product at k = 5, 6
    res["k56_by_type"] = {f"{k}_{t}": dict(collections.Counter(rel(r["n"], r["best"]["n"]) for r in rows
                                                                 if r["k"] == k and inv[r["spec"]]["product"] == (t == "product")))
                          for k in (5, 6) for t in ("entangled", "product")}
    # every class where some route is at least as good in inputs and kappa (and better in one)
    dom = []
    for r in rows:
        for key in ("cat", "tel", "tel_batch", "min_kappa"):
            x = r[key]; a, b = x["n"] <= r["n"], F(x["kappa"]) <= F(r["kappa_factory"])
            if a and b and (x["n"] < r["n"] or F(x["kappa"]) < F(r["kappa_factory"])):
                dom.append((r["spec"], r["n"], r["kappa_factory"], x["construction"], x["n"], x["kappa"])); break
    res["route_at_least_as_good_in_both"] = dom
    res["factory_higher_kappa_than_lowkappa_route"] = [(r["spec"], r["n"], r["kappa_factory"], r["min_kappa"]["construction"],
                                                        r["min_kappa"]["n"], r["min_kappa"]["kappa"]) for r in rows
                                                       if F(r["kappa_factory"]) > F(r["min_kappa"]["kappa"])]
    B = [r for r in rows if "batch_search" in r]
    nb = [r for r in B if not r["batch_search"]["exhaustive"] and r["best"]["construction"] == "batch"]
    res["batch_kappa_sensitivity"] = dict(
        rows_with_batch=len(B),
        max_relative_lowering_vs_6_variants=round(float(max(1 - F(r["batch_search"]["kappa"]) / F(r["batch_search"]["kappa_6variants"]) for r in B)), 4),
        cheapest_is_nonexhaustive_batch=len(nb),
        max_kappa_ratio_factory_over_route_there=round(float(max(F(r["kappa_factory"]) / F(r["best"]["kappa"]) for r in nb)), 3))
    RFk = [F(r["kappa_factory"]) / F(r["best"]["kappa"]) for r in RF]
    res["route_fewer_inputs"]["kappa_ratio_factory_over_route_median"] = round(float(st.median(RFk)), 3)
    res["second_order_acceptance"] = dict(max_ov_factory=max(r["ov_factory"] for r in rows),
                                          max_ov_fewest=float(max(F(r["best"]["ov"]) for r in rows)),
                                          ov_vs_fewest=dict(collections.Counter(rel(F(r["ov_factory"]), F(r["best"]["ov"])) for r in rows)))
    # known-answer checks against the stored summaries
    a = res["all"]
    assert a["fewest"]["inputs"] == {"factory_lower": 659, "equal": 47, "factory_higher": 60}, a["fewest"]["inputs"]
    sk = summ["all"]["best"]["kappa_ideal_catalyst"]
    assert (a["fewest"]["kappa"].get("factory_lower", 0), a["fewest"]["kappa"].get("equal", 0),
            a["fewest"]["kappa"].get("factory_higher", 0)) == (sk["factory_lower"], sk["equal"], sk["factory_higher"])
    sm = summ["all"]["min_kappa"]["kappa_ideal_catalyst"]
    assert (a["lowkappa"]["kappa"].get("factory_lower", 0), a["lowkappa"]["kappa"].get("equal", 0),
            a["lowkappa"]["kappa"].get("factory_higher", 0)) == (sm["factory_lower"], sm["equal"], sm["factory_higher"])
    sa = summ["all"]["best"]["acceptance"]
    assert (a["fewest"]["rejection"]["factory_lower"], a["fewest"]["rejection"]["equal"],
            a["fewest"]["rejection"]["factory_higher"]) == (sa["rej_factory_lower"], sa["rej_equal"], sa["rej_factory_higher"])
    M = json.load(open(os.path.join(C, "route_margin.json")))
    assert res["below_synthillation"]["fewest"]["inputs"]["factory_lower"] == M["below_synthillation"]["factory_fewer"] == 314
    for k in range(1, 7):
        assert res["by_k"][k]["fewest"]["inputs"].get("factory_lower", 0) == M["by_k"][str(k)]["factory_fewer"]
    RC = {r["spec"]: r for r in json.load(open(os.path.join(C, "route_compare_l3.json")))["rows"]}
    for r in rows:
        c = RC[r["spec"]]
        assert (r["cat"]["n"], r["tel"]["n"], r["tel_batch"]["n"]) == (c["n_cat"], c["n_tel"], c["n_tel_batch"])
    res["full_clifford"] = full_clifford(rows, inv)
    json.dump(res, open(os.path.join(C, "route_table_l3.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in res["full_clifford"].items() if k != "by_k"}, indent=0))
    for k in range(1, 7):
        t = res["full_clifford"]["by_k"][k]; print("full", k, t["classes"], t["inputs"], t["kappa"], t["kappa_lowest"], t["rejection"])
    print(json.dumps({k: res[k] for k in ("all", "below_synthillation", "route_fewer_inputs", "k56_by_type",
                                          "route_at_least_as_good_in_both", "factory_higher_kappa_than_lowkappa_route",
                                          "second_order_acceptance", "batch_kappa_sensitivity")}, indent=0))
    for k in range(1, 7):
        t = res["by_k"][k]; print(k, t["classes"], t["fewest"]["inputs"], t["fewest"]["kappa"], t["lowkappa"]["kappa"])


if __name__ == "__main__":
    main()
