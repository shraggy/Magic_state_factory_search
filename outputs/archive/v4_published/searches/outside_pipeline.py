"""Factories smaller than the final pipeline's, found only by earlier runs that the pipeline does not reproduce.
The paper reports the final pipeline only (outputs/staged_best.json); these are kept, re-verified and explained here.

Sources compared per l=3 class with the final best (outputs/staged_best.json):
  sweep A, stage (iii) as first run  outputs/staged_A (k=5,6), outputs/staged_pc (k<=4); passes iii_r, iii_f
  staged_A_v2                 outputs/solver_v2/staged_A_v2 (first XOR free-SAT rebuild, seeded from sweep A's iii_r)
  staged_A_v3                 outputs/solver_v2/staged_A_v3 (first stage-(iii) chain from the gluing seed)
How they were obtained: each comes from restricted SAT (or free SAT seeded by it) whose allowed output parts used a
Lempel factor returned by synthillation_compare.opt_B in z3's global context. Which of the many valid factors that
returns depends on every earlier z3 call in the same process, so these runs searched a different, history-dependent
set of output parts. The final pipeline computes the parts in a fresh z3 context (searches/det_parts.py) and runs every
query in a fresh process, which makes it independent of history and of the x86_64 host but no longer reaches these.
Every factory is re-verified here: us.is_factory (monomial form and brute force on all 2^N states), distance 2, n, and
the exact class label of the l=3 labeller.
    python searches/outside_pipeline.py   -> outputs/outside_pipeline/smaller_than_pipeline.csv, README.md
"""
import os, sys, csv, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
import cascade_class as CC
import argparse as _ap
_p = _ap.ArgumentParser(); _p.add_argument("--sb", default="staged_best.json", help="staged-best file under outputs/")
_p.add_argument("--out-dir", default=None, help="write every output here instead of outputs/ (and the manuscript)")
_A = _p.parse_args(); SB_FILE = _A.sb; OUTD = _A.out_dir

SOURCES = [("sweep A, stage (iii) as first run (2026-11-01/02)", lambda k: os.path.join(O, "staged_A" if k >= 5 else "staged_pc", f"l3_k{k}"),
            ("iii_r", "iii_f")),
           ("first XOR free-SAT rebuild (staged_A_v2, 2026-11-03)", lambda k: os.path.join(O, "solver_v2", "staged_A_v2", f"l3_k{k}"),
            ("iii_r", "iii_f")),
           ("first stage-(iii) chain from the gluing seed (staged_A_v3, 2026-11-03)",
            lambda k: os.path.join(O, "solver_v2", "staged_A_v3", f"l3_k{k}"), None)]


def main():
    final = {(r["k"], str(r["orbit"])): r for r in json.load(open(os.path.join(O, SB_FILE)))["rows"] if r["l"] == 3}
    best = {}
    for name, d_of, passes in SOURCES:
        for k in range(1, 7):
            d = d_of(k)
            files = [os.path.join(d, f"{p}.json") for p in passes] if passes else glob.glob(os.path.join(d, "iii*.json"))
            for f in files:
                if not os.path.exists(f): continue
                for o, b in json.load(open(f))["best"].items():
                    key = (k, o)
                    if key in final and "gates" in b and b["n"] < final[key]["n"] and (key not in best or b["n"] < best[key][1]["n"]):
                        best[key] = (name, b, os.path.relpath(f, ROOT))
    rows = []
    for (k, o), (name, b, src) in sorted(best.items()):
        g = [(frozenset(S), c) for S, c in b["gates"]]; f = dict(k=k, r=b["N"] - k, gates=g)
        classes, ctx = CC.classes_of(3, k, os.path.join(O, "staged_A", "_cache"))
        label = CC.labeller(3, k, ctx, classes)
        n = sum(us.inputs(c, 3) for _, c in g)
        assert us.is_factory(f, 3) and us.distance(f, 3) == 2 and n == b["n"] and label(f) == o, (k, o)
        fr = final[(k, o)]
        rows.append(dict(k=k, orbit=o, spec=fr["spec"], tau=fr["tau"], n=b["n"], N=b["N"], final_n=fr["n"],
                         synthillation_n=fr["syn_n"], source=name, file=src, stage=b.get("stage", ""),
                         gates=json.dumps([[sorted(S), c] for S, c in g]),
                         quirk=us.quirk_url(g, b["N"], 3, k) if b["N"] <= 14 else ""))
    od = OUTD or os.path.join(O, "outside_pipeline"); os.makedirs(od, exist_ok=True)
    with open(os.path.join(od, "smaller_than_pipeline.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]) if rows else ["k"]); w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} classes where an earlier run has a smaller valid factory than the final pipeline")


if __name__ == "__main__":
    main()
