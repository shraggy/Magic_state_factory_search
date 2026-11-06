"""Distance-2 protocols of Gong et al. (arXiv:2608.09727, Table 1), each converted to T inputs and set against the
smallest factory the search finds for the same output class (outputs/staged_best.json).

Input conversions, the cheapest named in that paper: a CS state from 3 T (3T->CS), a CCZ state from 8 T (8T->CCZ
synthillation) or from 2 CS (Beverland et al., Fig. 13), a CCS state from its smallest l=4 factory here. Output classes
are labelled with the exact class labeller (TOF# = CCZ_{123}CCZ_{345} is the class of CCZ014,CCZ023).
    python searches/gong_compare.py   -> outputs/checks/gong_compare.json
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
sys.path[:0] = [HERE, os.path.join(ROOT, "classification")]
import cascade_class as CC
from enumerate_targets import index_of_poly
from coset_certify import parse_target
import argparse as _ap
_p = _ap.ArgumentParser(); _p.add_argument("--sb", default="staged_best.json", help="staged-best file under outputs/")
_p.add_argument("--out-dir", default=None, help="write every output here instead of outputs/ (and the manuscript)")
_A = _p.parse_args(); SB_FILE = _A.sb; OUTD = _A.out_dir

SB = {(r["l"], r["k"], r["spec"]): r for r in json.load(open(os.path.join(O, SB_FILE)))["rows"]}


def label3(spec, k):
    classes, ctx = CC.classes_of(3, k, os.path.join(O, "staged_A", "_cache"))
    return CC.labeller(3, k, ctx, classes), classes


def main():
    rows = []
    def add(name, nin, kind, t_in, l, k, spec, note=""):
        r = SB.get((l, k, spec))
        rows.append(dict(protocol=name, inputs=nin, input_type=kind, t_inputs=t_in, output=spec, l=l, k=k,
                         ours_n=r and r["n"], ours_N=r and r["N"], tau=r and r["tau"], synthillation_n=r and r["syn_n"],
                         on_bound=bool(r and r["tau"] is not None and r["n"] == r["tau"] + 1), note=note))
    for m, spec in ((1, "CS01"), (2, "CS03,CS12"), (3, "CS05,CS14,CS23")):
        add(f"({2*m+2})CS -> {m}CS", 2 * m + 2, "CS", 3 * (2 * m + 2), 3, 2 * m, spec)
    add("8CS -> TOF#", 8, "CS", 24, 3, 5, "CCZ014,CCZ023", "TOF# = CCZ_{123}CCZ_{345}, the class of CCZ014,CCZ023")
    add("8CCZ -> 2CCZ", 8, "CCZ", 64, 3, 6, "CCZ045,CCZ123", "64 T via 8T->CCZ; 16 CS (48 T) via 2CS->CCZ")
    add("8CCS -> 1CCS", 8, "CCS", None, 4, 3, "CCS012",
        "l=4: our factory takes sqrt(T) inputs; the protocol takes CCS inputs, each needing at least the "
        f"{SB[(4, 3, 'CCS012')]['n']} sqrt(T) inputs of our smallest CCS factory, so 8 of them exceed that factory")
    # check the TOF# class claim exactly
    lab, classes = label3("CCZ014,CCZ023", 5)
    C, off, lb, fib = CC.classes_of(3, 5, os.path.join(O, "staged_A", "_cache"))[1]
    o = str(int(lb[index_of_poly(parse_target("CCZ012,CCZ234", 3), 5, 3, C, off)]))
    assert classes[o]["spec"] == "CCZ014,CCZ023"
    json.dump(dict(source="Gong et al. arXiv:2608.09727 Table 1 (d=2 rows); outputs/" + SB_FILE, rows=rows),
              open(os.path.join(OUTD or os.path.join(O, "checks"), "gong_compare.json"), "w"), indent=1)
    for r in rows:
        print(f"{r['protocol']:16s} {str(r['t_inputs']):>4s} T  vs ours n={r['ours_n']} (tau {r['tau']}, synthillation {r['synthillation_n']}, bound {r['on_bound']})")


if __name__ == "__main__":
    main()
