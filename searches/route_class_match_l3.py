"""Do the catalogue factory and the distil-then-synthesize routes make the same target class? For each of the 766 l=3
classes with k<=6 in outputs/staged_best.json, the noiseless output polynomial of the factory (its gate list, outputs
0..k-1) and the class label (the representative 'spec' that searches/route_compare_l3.py and route_kappa_l3.py
synthesize) are both mapped to their exact class with the cascade's labeller (searches/cascade_class.py: orbit tables
of enumerate_targets.py for k<=5, enumerate_targets_fiber.py at k=6), and compared with the row's 'orbit'.
The two outputs differ as states on most classes, but lie in one orbit under CNOT circuits, Pauli X and diagonal
Cliffords, so each is made from the other by free (Clifford) operations and the input and kappa comparisons compare
the same target class.
Known-answer check (asserted): CCZ012 and T0 map to their own rows' orbits; a perturbed polynomial (one extra T) maps
to a different orbit.
    nice -n 19 python searches/route_class_match_l3.py   -> outputs/checks/route_class_match_l3.json
Uses the pickled class tables in outputs/staged_A/_cache/ when present (rebuilt otherwise, up to about an hour).
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import cascade_class as CC
from cascade_class import index_of_poly
from classify import phase_polynomial
from coset_certify import parse_target


def main():
    rows = [r for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"] if r["l"] == 3 and r["k"] <= 6]
    assert len(rows) == 766
    res = dict(classes=0, factory_matches_orbit=0, label_matches_orbit=0, both=0, mismatches=[], by_k={})
    known = {}
    for k in range(1, 7):
        classes, (C, off, lab, fiber) = CC.classes_of(3, k, os.path.join(ROOT, "outputs", "staged_A", "_cache"))
        def L(P):
            if lab is not None: return str(int(lab[index_of_poly(P, k, 3, C, off)]))
            lb = fiber.label(fiber.digits_of_poly(P)); return lb if isinstance(lb, str) else f"{lb[0]}:{lb[1]}"
        kk = [r for r in rows if r["k"] == k]; good = 0
        for r in kk:
            g = r["gates"]
            a = L(phase_polynomial([frozenset(s) for s, _ in g], [c for _, c in g], set(range(k)), 3))
            b = L(parse_target(r["spec"], 3)); o = str(r["orbit"])
            res["classes"] += 1; res["factory_matches_orbit"] += a == o; res["label_matches_orbit"] += b == o
            if a == b == o: res["both"] += 1; good += 1
            else: res["mismatches"].append(dict(spec=r["spec"], orbit=o, factory=a, label=b))
            if r["spec"] in ("CCZ012", "T0"): known[r["spec"]] = (a == o)
        if k == 3:                                          # perturbation: CCZ012 times T0 is another class
            P = dict(parse_target("CCZ012", 3)); P[frozenset([0])] = 1
            known["CCZ012*T0_differs"] = L(P) != next(str(r["orbit"]) for r in kk if r["spec"] == "CCZ012")
        res["by_k"][k] = f"{good} of {len(kk)}"
    assert all(known.values()) and len(known) == 3, known
    res["known_answer_checks"] = known
    json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "route_class_match_l3.json"), "w"), indent=1)
    print(json.dumps({x: res[x] for x in ("classes", "factory_matches_orbit", "label_matches_orbit", "both", "by_k", "known_answer_checks")}))


if __name__ == "__main__":
    main()
