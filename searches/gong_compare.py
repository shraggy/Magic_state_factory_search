"""Distance-2 protocols of Gong et al. (arXiv:2608.09727, Table 1), each converted to T inputs and set against the
smallest factory the search finds for the same output class (outputs/staged_best.json).

Input conversions at the protocol's own distance 2, without distillation (the strongest we know, AGENTS rule 5): a CS
state from 3 T (3T->CS, named in that paper), a CCZ state from 4 T (the measurement-assisted Toffoli* of Jones,
arXiv:1212.5069; Beverland et al., arXiv:1904.01124, Table 1 and Fig. 11) or from 7 T with CNOT and T gates
(tau(CCZ) = 7), a CCS state from 7 sqrt(T) (Beverland et al., Table 1; tau(CCS) = 7). Fed through the 8T->CCZ factory
the 8CCZ->2CCZ protocol uses 64 T, but the result is a distance-4 protocol (Remark 4.4 of that paper), so 64 is
recorded only in the note. The 4T->CCZ conversion uses a measurement and so lies outside the form of Corollary 3.
Output classes are labelled with the exact class labeller (TOF# = CCZ_{123}CCZ_{345} is the class of CCZ014,CCZ023).
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
    add("8CCZ -> 2CCZ", 8, "CCZ", 32, 3, 6, "CCZ045,CCZ123",
        "32 T via 4T->CCZ (Jones Toffoli*, with a measurement); 56 via CNOT+T (tau=7); 48 via 2CS->CCZ and 3T->CS; "
        "64 via 8T->CCZ is the distance-4 protocol of Haah-Hastings (Remark 4.4 of arXiv:2608.09727)")
    add("8CCS -> 1CCS", 8, "CCS", 56, 4, 3, "CCS012",
        "l=4: 7 sqrt(T) per CCS (Beverland et al. Table 1; tau(CCS)=7), 56 in all; our factory distils CCS from "
        f"{SB[(4, 3, 'CCS012')]['n']} sqrt(T) inputs")
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
