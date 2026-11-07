"""Synthillation alone against the distil-then-synthesize routes, for the 766 l=3 classes with k<=6 (inputs only).

Question (mock Referee A, 2026-11-06): does the synthillation construction alone already beat the same routes on most
classes, and how much of the factory's input saving over the routes does it already give?

Route per class: exactly the 'full_clifford' route of searches/route_table_l3.py, i.e. the cheapest (inputs, then kappa)
'best' route of outputs/checks/route_kappa_l3.json over the members of the class's full-Clifford group
(outputs/checks/full_clifford_classes.json), with the lambda route (outputs/checks/lambda_route.json) also allowed; the
lambda route lowers no class after merging (route_table_l3.json, lowered_by_lambda_route_after_merging = 0), asserted.
Synthillation per class: syn_n of route_kappa_l3.json (= outputs/staged_best.json syn_n, the strongest Campbell-Howard
construction built by searches/class_table_l3.py). Reported two ways: per diagonal class (as the referee did), and,
on the same footing as the route, the minimum over the full-Clifford group ('syn_full'); and once per full-Clifford
class (622), with factory, synthillation and route all minimised over the class ('per_full_clifford_class_all_minimised').
Known answers (asserted): factory vs route 608 / 62 / 96; route total 27,130; factory total 21,364; factory never above
synthillation (330 lower, 436 equal).
    python searches/syn_vs_routes.py   -> outputs/checks/syn_vs_routes.json
"""
import os, json, collections
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); C = os.path.join(ROOT, "outputs", "checks")


def rel(a, b):
    return "lower" if a < b else ("equal" if a == b else "higher")


def main():
    rows = json.load(open(os.path.join(C, "route_kappa_l3.json")))["rows"]
    assert len(rows) == 766
    by = {r["spec"]: r for r in rows}; grp = {}
    for g in json.load(open(os.path.join(C, "full_clifford_classes.json")))["groups"]:
        mem = [x.rsplit(":", 1)[1] for x in g]
        for m in mem: grp[m] = mem
    LR = {r["spec"]: r for r in json.load(open(os.path.join(C, "lambda_route.json")))["rows"]}
    route = {}
    for r in rows:
        mem = grp.get(r["spec"], [r["spec"]])
        b = min((by[m]["best"] for m in mem), key=lambda x: (x["n"], F(x["kappa"])))["n"]
        lam = [LR[m]["route_n"] for m in mem if m in LR and LR[m].get("route_n") is not None]
        assert not lam or min(lam) >= b
        route[r["spec"]] = b
    syn = {r["spec"]: r["syn_n"] for r in rows}
    synf = {r["spec"]: min(by[m]["syn_n"] for m in grp.get(r["spec"], [r["spec"]])) for r in rows}
    fac = {r["spec"]: r["n"] for r in rows}
    S = [r["spec"] for r in rows]
    # known answers
    assert collections.Counter(rel(fac[s], route[s]) for s in S) == {"lower": 608, "equal": 62, "higher": 96}
    assert sum(route.values()) == 27130 and sum(fac.values()) == 21364
    assert collections.Counter(rel(fac[s], syn[s]) for s in S) == {"lower": 330, "equal": 436}

    def block(sv, name):
        joint = collections.Counter((rel(fac[s], route[s]), rel(sv[s], route[s])) for s in S)
        tot = sum(sv.values()); Rt = sum(route.values()); Ft = sum(fac.values())
        pos_f = sum(max(0, route[s] - fac[s]) for s in S); pos_s = sum(max(0, route[s] - sv[s]) for s in S)
        W = [s for s in S if fac[s] < syn[s]]
        return dict(
            synthillation=name,
            syn_vs_route=dict(collections.Counter(rel(sv[s], route[s]) for s in S)),
            joint_factory_vs_route__syn_vs_route={f"{a}|{b}": v for (a, b), v in sorted(joint.items())},
            totals=dict(factory=Ft, synthillation=tot, route=Rt),
            net_saving=dict(factory=Rt - Ft, synthillation=Rt - tot, share_by_synthillation=round((Rt - tot) / (Rt - Ft), 4)),
            saving_on_classes_where_lower=dict(factory=pos_f, synthillation=pos_s, share_by_synthillation=round(pos_s / pos_f, 4)),
            factory_below_synthillation=dict(classes=len(W), inputs_saved=sum(sv[s] - fac[s] for s in W) if sv is syn else
                                             sum(syn[s] - fac[s] for s in W),
                                             syn_vs_route_there=dict(collections.Counter(rel(sv[s], route[s]) for s in W))),
            factory_beats_route_only=sorted(s for s in S if fac[s] < route[s] and not sv[s] < route[s]))
    out = dict(source="outputs/checks/route_kappa_l3.json, outputs/checks/full_clifford_classes.json, outputs/checks/lambda_route.json",
               note="inputs only; route = route_table_l3.json full_clifford route; 'lower' = left item needs fewer inputs",
               factory_vs_route=dict(collections.Counter(rel(fac[s], route[s]) for s in S)),
               per_diagonal_class=block(syn, "syn_n per diagonal class"),
               full_clifford_min=block(synf, "syn_n minimised over the full-Clifford group, as the route is"))
    # equal footing: one entry per full-Clifford class (622), factory, synthillation and route each minimised over it
    reps = {}
    for s in S:
        reps.setdefault(tuple(sorted(grp.get(s, [s]))), s)
    G = list(reps)
    assert len(G) == 622, len(G)
    fg = {g: min(fac[m] for m in g) for g in G}; sg = {g: min(syn[m] for m in g) for g in G}; rg = {g: route[reps[g]] for g in G}
    assert all(rg[g] == min(route[m] for m in g) for g in G)
    Rg, Fg, Sg = sum(rg.values()), sum(fg.values()), sum(sg.values())
    out["per_full_clifford_class_all_minimised"] = dict(
        classes=len(G), note="factory, synthillation and route each minimised over the full-Clifford class",
        factory_vs_route=dict(collections.Counter(rel(fg[g], rg[g]) for g in G)),
        syn_vs_route=dict(collections.Counter(rel(sg[g], rg[g]) for g in G)),
        factory_vs_syn=dict(collections.Counter(rel(fg[g], sg[g]) for g in G)),
        joint_factory_vs_route__syn_vs_route={f"{a}|{b}": v for (a, b), v in sorted(collections.Counter(
            (rel(fg[g], rg[g]), rel(sg[g], rg[g])) for g in G).items())},
        totals=dict(factory=Fg, synthillation=Sg, route=Rg),
        net_saving=dict(factory=Rg - Fg, synthillation=Rg - Sg, share_by_synthillation=round((Rg - Sg) / (Rg - Fg), 4)))
    for b in (out["per_diagonal_class"], out["full_clifford_min"]):
        b["factory_beats_route_only_count"] = len(b["factory_beats_route_only"])
    json.dump(out, open(os.path.join(C, "syn_vs_routes.json"), "w"), indent=1)
    print(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "factory_beats_route_only"})
                      for k, v in out.items()}, indent=1))


if __name__ == "__main__":
    main()
