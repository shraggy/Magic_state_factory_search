"""Distil-then-synthesize routes against the cascade's factory, for the 766 l=3 classes with k<=6.

Two achievable distance-2 routes, both using Clifford gates, measurements and feedforward (so they lie outside the
unitary form of Corollary 3 and are not bounded by it):
  cat  : distil T states with the catalysed 8T->2T factory (Gidney-Fowler, arXiv:1812.01238; [[8,3,2]] followed by a
         catalysed CCZ->2T conversion, 4 T inputs per T state once the catalyst is paid for) and synthesize the target
         from tau of them by CNOT and T gates: n_cat = 4 tau.
  tel  : distil each gate of the representative's label separately with the smallest factory for it in this
         catalogue (CCZ: [[8,3,2]], 8 inputs; CS: [[12,2,2]], 12; T: the catalysed route, 4) and teleport the gates
         onto the outputs with Clifford corrections: n_lab = 8 #CCZ + 12 #CS + 4 #T, and n_tel = min(n_lab, n_cat).
         The label is one representative of the class, so n_tel is an upper bound for the class, not minimised over it.
Leading-order output error of the catalysed route: about 14 p^2 per T state (28 p^2 per CCZ->2T block, catalyst
errors neglected), so c_cat ~ 14 tau; this estimate is computed, unchecked.
Counts, against n of outputs/staged_best.json: route < n, = n, > n, over all 766 and over the classes with
n < n_syn (fewer inputs than the synthillation construction); c < 14 tau.
Deterministic; reads outputs/staged_best.json only.
    python searches/route_compare_l3.py   -> outputs/checks/route_compare_l3.json
Known-answer checks (asserted): the factory for CCZ012 has n=8 (route cat 28, tel 8); T0 has n=14 (cat 4).
"""
import os, json, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COST = {"CCZ": 8, "CS": 12, "T": 4}


def label_cost(spec):
    return sum(COST[g.rstrip("0123456789")] for g in spec.split(","))


def main():
    rows = [r for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]
            if r["l"] == 3 and r["k"] <= 6]
    assert len(rows) == 766
    out = []; S = collections.Counter()
    for r in rows:
        n, tau = r["n"], r["tau"]; ncat = 4 * tau; nlab = label_cost(r["spec"]); ntel = min(nlab, ncat)
        few = r.get("syn_n") is not None and n < r["syn_n"]
        for nm, v in (("cat", ncat), ("tel", ntel)):
            rel = "fewer" if v < n else ("equal" if v == n else "more")
            S[f"{nm}_{rel}"] += 1
            if few:
                S[f"{nm}_{rel}_among_fewer_than_syn"] += 1
        S["fewer_than_syn"] += few
        S["c_below_14tau"] += r["c"] < 14 * tau
        out.append(dict(k=r["k"], orbit=r["orbit"], spec=r["spec"], n=n, tau=tau, c=r["c"], syn_n=r.get("syn_n"),
                        n_cat=ncat, n_label=nlab, n_tel=ntel, c_cat_est=14 * tau))
    by = {r["spec"]: r for r in out}
    assert by["CCZ012"]["n"] == 8 and by["CCZ012"]["n_cat"] == 28 and by["CCZ012"]["n_tel"] == 8
    assert by["T0"]["n"] == 14 and by["T0"]["n_cat"] == 4
    fewer_cat = sorted((r["k"], r["spec"], r["n"], r["n_cat"]) for r in out if r["n_cat"] < r["n"])
    fewer_tel = sorted((r["k"], r["spec"], r["n"], r["n_tel"]) for r in out if r["n_tel"] < r["n"])
    eq_tel_syn = [(r["spec"], r["n"], r["n_tel"]) for r in out
                  if r["syn_n"] is not None and r["n"] < r["syn_n"] and r["n_tel"] == r["n"]]
    s295 = [r for r in out if r["syn_n"] is not None and r["n"] < r["syn_n"]]
    summary = dict(classes=len(out), **S,
                   min_margin_cat_among_fewer_than_syn=min(r["n_cat"] - r["n"] for r in s295),
                   min_margin_tel_among_fewer_than_syn=min(r["n_tel"] - r["n"] for r in s295),
                   cat_fewer_classes=fewer_cat, tel_fewer_classes=fewer_tel, tel_equal_classes_among_fewer_than_syn=eq_tel_syn)
    os.makedirs(os.path.join(ROOT, "outputs", "checks"), exist_ok=True)
    json.dump(dict(source="outputs/staged_best.json", summary=summary, rows=out),
              open(os.path.join(ROOT, "outputs", "checks", "route_compare_l3.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in summary.items() if not isinstance(v, list)}, indent=1))
    print("cat fewer:", fewer_cat); print("tel fewer:", fewer_tel); print("tel equal among fewer-than-syn:", eq_tel_syn)


if __name__ == "__main__":
    main()
