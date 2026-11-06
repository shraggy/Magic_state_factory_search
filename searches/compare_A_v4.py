"""Compare the stage-(iii) cascade (merged by merge_A_v3.py into OUT) with outputs/staged_best.json, staged_A_v2 and
the history-dependent run (staged_A_v3), per l=3 class with k <= 6.
  per class: tau, staged_best n, A_v2 n, v3 n (buggy run, kept separate), v4 n (this run), which stage set v4's best
  counts: v4 vs staged_best lower/equal/higher, inputs saved, below tau+1 (must be 0); v4 vs A_v2; classes where v3
  (or A_v2) holds a smaller valid factory than v4, with that factory re-verified (is_factory monomial and brute force
  for N<=9, distance, exact label) so it can be reported separately
  the 9 classes above the earlier catalogue (outputs/staged_A/compare_l3.json worse_list): staged_best, v4, catalogue
    python searches/compare_A_v4.py OUT      -> OUT/compare.json, OUT/best_per_class.json
"""
import os, sys, json, csv
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
NINE = ["1892", "2320:771", "2320:3224", "2320:9752", "2320:16992", "2320:82457", "2320:82505", "2320:82512", "41984:4228"]
CELLS = [(1, "staged_pc"), (2, "staged_pc"), (3, "staged_pc"), (4, "staged_pc"), (5, "staged_A"), (6, "staged_A")]


def best_of(path):
    return json.load(open(path))["best"] if os.path.exists(path) else {}


def main(out):
    ref = {(r["k"], r["orbit"]): r for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"] if r["l"] == 3}
    cat = {}
    for k in (5, 6):
        p = os.path.join(ROOT, "outputs", "earlier_catalogues", f"l3_catalogue_k{k}.csv")
        if os.path.exists(p):
            for r in csv.DictReader(open(p)):
                cat[(k, frozenset(r["spec"].split(",")))] = int(r["n"])
    rows, smaller_elsewhere = [], []
    for k, _ in CELLS:
        v4 = best_of(os.path.join(out, f"l3_k{k}", "iii.json"))
        v2 = best_of(os.path.join(ROOT, "outputs", "solver_v2", "staged_A_v2", f"l3_k{k}", "iii_f.json"))
        v3 = best_of(os.path.join(ROOT, "outputs", "solver_v2", "staged_A_v3", f"l3_k{k}", "iii.json"))
        for o, b in v4.items():
            r = ref[(k, o)]
            row = dict(k=k, orbit=o, spec=r["spec"], tau=r["tau"], staged_best=r["n"], A_v2=v2.get(o, {}).get("n"),
                       v3=v3.get(o, {}).get("n"), v4=b["n"], N=b["N"], stage=b["stage"],
                       catalogue=cat.get((k, frozenset(r["spec"].split(",")))))
            rows.append(row)
            for name, src in (("A_v2", v2), ("v3", v3)):
                if o in src and src[o]["n"] < b["n"]:
                    smaller_elsewhere.append(dict(k=k, orbit=o, source=name, n=src[o]["n"], v4=b["n"], N=src[o]["N"],
                                                  gates=src[o]["gates"]))
    # re-verify the factories that beat v4 elsewhere
    if smaller_elsewhere:
        import staged as STG, staged_class as SC, unified_search as us
        for k, src in CELLS:
            mine = [x for x in smaller_elsewhere if x["k"] == k]
            if not mine:
                continue
            cache, classes = SC.setup(3, k, os.path.join(ROOT, "outputs", src))
            for x in mine:
                g = [(frozenset(S), c) for S, c in x["gates"]]; f = dict(k=k, r=x["N"] - k, gates=g)
                x["check"] = dict(valid=us.is_factory(f, 3), brute=us.is_factory(f, 3, brute_force=True) if x["N"] <= 9 else None,
                                  d=str(us.distance(f, 3)), n=sum(us.inputs(c, 3) for _, c in g),
                                  label_ok=STG.label_fix(f, 3, k, cache)[0] == x["orbit"])
                del x["gates"]
    cmp = lambda a, b: "lower" if a < b else "equal" if a == b else "higher"
    cnt = lambda key: {c: sum(1 for r in rows if r[key] is not None and cmp(r["v4"], r[key]) == c) for c in ("lower", "equal", "higher")}
    summ = dict(classes=len(rows), v4_vs_staged_best=cnt("staged_best"),
                inputs_saved_vs_staged_best=sum(max(0, r["staged_best"] - r["v4"]) for r in rows),
                below_tau_plus_1=sum(1 for r in rows if r["tau"] is not None and r["v4"] < r["tau"] + 1),
                at_tau_plus_1=sum(1 for r in rows if r["tau"] is not None and r["v4"] == r["tau"] + 1),
                v4_vs_A_v2=cnt("A_v2"), v4_vs_v3=cnt("v3"),
                smaller_elsewhere=len(smaller_elsewhere),
                nine=[{a: r[a] for a in ("k", "orbit", "staged_best", "A_v2", "v3", "v4", "catalogue")}
                      for r in rows if r["orbit"] in NINE])
    json.dump(dict(summary=summ, smaller_elsewhere=smaller_elsewhere, rows=rows), open(os.path.join(out, "compare.json"), "w"), indent=0)
    json.dump({f"l3_k{r['k']}|{r['orbit']}": r for r in rows}, open(os.path.join(out, "best_per_class.json"), "w"), indent=0)
    print(json.dumps(summ, indent=1)); print("smaller elsewhere:", [(x["k"], x["orbit"], x["source"], x["n"], x["v4"], x.get("check")) for x in smaller_elsewhere])


if __name__ == "__main__":
    main(sys.argv[1])
