"""Best factory per output class from the staged cascade, every cell with an exact class list, in one file that the
paper's tables and figures read.

Cells and their final pass (cascade order i_tg, i_gs, i_gh, ii_dj, ii_eg, iii_r, iii_f; the last file present):
  l=2, k<=6 and l=4, k<=2   outputs/staged_small/       (searches/staged.py, one process per cell)
  l=3, k<=6                 outputs/staged_final_b3/    (stages (i)-(ii) of sweep A; stage (iii) per class the best of
                                                        the three deterministic runs v4, v5, v6, searches/build_staged_best3.py),
                                                        or the directory given by --final-dir (below)
  l=4, k=3                  outputs/staged_l4fix/       (searches/staged.py)
  l=4, k=4                  outputs/staged_A/           (sweep A, searches/staged_class.py)
  --extra DIR ...           later staged runs of the same cells, merged per class (none used for the paper: staged_best.json has extra=[])
For each class: n, N, the stage that gave the final factory, the first pass whose cumulative best reached that n
(ties to the earliest pass, as in Table I), the gates, the exact Z(pi/2^l)-count tau (l=3 class_table_l3.json, else
unified_search.zcount), the strongest synthillation n at l=3 (synthillation_strongest_l3.json) and the output-error
coefficient c (error_coefficients.coefficients). Every factory is re-verified (us.is_factory, distance 2, and the
exact class label of the cell's labeller) before it is written; a failure stops the script. Each l=4 row carries
third_level_target: true for the 22 classes whose target is a level-3 gate (see third_level_target).
    python searches/staged_best.py [--extra outputs/staged_A_rerun] -> outputs/staged_best.json
    python searches/staged_best.py --final-dir staged_final --out staged_best_v4.json
        the earlier published version: l=3 stage (iii) from the v4 run alone (outputs/staged_final,
        build_staged_final.py); its outputs are kept in outputs/archive/v4_published/
"""
import os, sys, json, argparse
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
import cascade_class as CC
from coset_certify import parse_target
from error_coefficients import coefficients

PASSES = ["i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg", "iii_r", "iii_f"]   # l=4 has no stage (iii)
CELLS = [("staged_small", 2, k) for k in range(1, 7)] + [("staged_final", 3, k) for k in range(1, 7)] + \
        [("staged_small", 4, 1), ("staged_small", 4, 2), ("staged_l4fix", 4, 3), ("staged_A", 4, 4)]


def passes(d):
    return {p: json.load(open(os.path.join(d, f"{p}.json"))) for p in PASSES if os.path.exists(os.path.join(d, f"{p}.json"))}


def cell(src, L, k, extra):
    d = os.path.join(O, src, f"l{L}_k{k}"); P = passes(d); last = P[[p for p in PASSES if p in P][-1]]
    best = {o: dict(b) for o, b in last["best"].items()}
    for x in extra:                                           # later runs: smaller n per class wins
        e = os.path.join(ROOT, x, f"l{L}_k{k}")
        for p, rec in (passes(e).items() if os.path.isdir(e) else []):
            for o, b in rec["best"].items():
                if o in best and b["n"] < best[o]["n"]:
                    best[o] = dict(b, stage=b["stage"] + f" ({os.path.basename(x)})")
    first = {}
    for o, b in best.items():                                 # first pass whose cumulative best reaches the final n
        first[o] = next((p for p in PASSES if p in P and o in P[p]["best"] and P[p]["best"][o]["n"] <= b["n"]),
                        b["stage"].split(":")[0])
    return last["specs"], best, first, {p: P[p].get("seconds") for p in P}


def third_level_target(spec, L=4):
    """True if the level-L target is a level-(L-1) diagonal gate: the coefficient on every set S of qubits is a
    multiple of 2^|S| (mod 2^L), i.e. the target is a product of squares of level-L gates (sT^2 = T, CT^2 = CS,
    CCS^2 = CCZ). At l=4, k<=4 these are the 22 classes of l=3 with k<=4 (outputs/checks/doubled_l3_at_l4.json);
    they are not fourth-level results, and the tables and figures leave them out of the l=4 counts."""
    return all(c % (1 << len(S)) == 0 for S, c in parse_target(spec, L).items())


def main(extra, final_dir="staged_final_b3", out="staged_best.json"):
    TAU3 = {(r["k"], r["spec"]): r["tcount"] for r in json.load(open(os.path.join(O, "class_table_l3.json")))}
    SYN = {(r["k"], r["spec"]): r["best_n"] for r in json.load(open(os.path.join(O, "synthillation_strongest_l3.json")))}
    rows, meta = [], {}
    for src, L, k in [(final_dir if s == "staged_final" else s, L, k) for s, L, k in CELLS]:
        specs, best, first, secs = cell(src, L, k, extra)
        cache = next((c for c in (os.path.join(O, src, "_cache"), os.path.join(O, "staged_A", "_cache"))
                      if os.path.exists(os.path.join(c, f"l{L}_k{k}_classes.pkl"))), os.path.join(O, src, "_cache"))
        classes, ctx = CC.classes_of(L, k, cache)
        label = CC.labeller(L, k, ctx, classes)
        assert set(best) == set(classes), (L, k, "classes without a factory", sorted(set(classes) - set(best))[:5])
        for o in sorted(best):
            b = best[o]; spec = specs[o]
            g = [(frozenset(S), c) for S, c in b["gates"]]; f = dict(k=k, r=b["N"] - k, gates=g)
            n = sum(us.inputs(c, L) for _, c in g)
            assert us.is_factory(f, L) and us.distance(f, L) == 2 and n == b["n"] and label(f) == o, (L, k, o)
            tau = TAU3.get((k, spec)) if L == 3 else us.zcount(parse_target(spec, L), k, L)
            try:
                c = coefficients(g, k, L)[1]
            except AssertionError:
                c = None
            row = dict(l=L, k=k, orbit=o, spec=spec, n=n, N=b["N"], stage=b["stage"], first=first[o], tau=tau,
                       syn_n=SYN.get((k, spec)) if L == 3 else None, c=c, gates=b["gates"])
            if L == 4:
                row["third_level_target"] = third_level_target(spec)
            rows.append(row)
        meta[f"l{L}_k{k}"] = dict(source=f"outputs/{src}/l{L}_k{k}", classes=len(classes), pass_seconds=secs)
        print(f"l={L} k={k}: {len(best)} classes from {src}", flush=True)
    json.dump(dict(cells=meta, extra=extra, rows=rows), open(os.path.join(O, out), "w"))
    t4 = sum(r.get("third_level_target", False) for r in rows); n4 = sum(r["l"] == 4 for r in rows)
    print(f"{len(rows)} classes -> {os.path.relpath(os.path.join(O, out), ROOT)}; l=4: {n4} listed, "
          f"{t4} third-level targets, {n4 - t4} genuine fourth-level")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--extra", nargs="*", default=[])
    ap.add_argument("--final-dir", default="staged_final_b3"); ap.add_argument("--out", default="staged_best.json")
    a = ap.parse_args(); main(a.extra, a.final_dir, a.out)
