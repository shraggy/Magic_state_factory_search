"""Compare the multi-factor stage-(iii) trial (OUT = outputs/solver_v2/staged_A_v5, merged by merge_A_v3.py) with
outputs/staged_best.json (the final pipeline), v4, A_v2, v3 and synthillation (outputs/synthillation_strongest_l3.json,
best_n), per l=3 class with k <= 6, and with the outside-pipeline list (outputs/outside_pipeline/smaller_than_pipeline.csv).
Read-only on everything but OUT/compare_v5.json and OUT/still_smaller_than_v5.csv.
  counts lower/equal/higher of v5 vs each reference; inputs saved vs staged_best; below tau+1 (must be 0)
  the 7 classes where v4 > published sweep A and the 35 smaller_elsewhere classes: how many v5 recovers
  vs synthillation: fewer / same / more inputs (classes with a synthillation value)
  still_smaller_than_v5: every earlier factory (staged_best, A_v2, v3, the outside-pipeline CSV) with n < v5,
  re-verified (is_factory monomial and brute force N<=9, distance != 1, exact class label, n), with its source.
    python searches/compare_A_v5.py OUT
"""
import os, sys, json, csv
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
CELLS = [(1, "staged_pc"), (2, "staged_pc"), (3, "staged_pc"), (4, "staged_pc"), (5, "staged_A"), (6, "staged_A")]
csv.field_size_limit(10 ** 9)


def best(p):
    return json.load(open(p))["best"] if os.path.exists(p) else {}


def main(out):
    S = os.path.join(ROOT, "outputs", "solver_v2")
    sb = {(r["k"], r["orbit"]): r for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"] if r["l"] == 3}
    syn = {(r["k"], r["spec"]): r.get("best_n") for r in json.load(open(os.path.join(ROOT, "outputs", "synthillation_strongest_l3.json")))}
    cmp4 = json.load(open(os.path.join(S, "staged_A_v4", "compare.json")))
    ts = json.load(open(os.path.join(S, "testset_v5.json")))
    cands = []                                            # (k, orbit, n, N, gates, source, how)
    rows = []
    for k, _ in CELLS:
        v5 = best(os.path.join(out, f"l3_k{k}", "iii.json"))
        refs = dict(v4=best(os.path.join(S, "staged_A_v4", f"l3_k{k}", "iii.json")),
                    v5=best(os.path.join(S, "staged_A_v5", f"l3_k{k}", "iii.json")) if os.path.basename(out.rstrip("/")) != "staged_A_v5" else {},
                    A_v2=best(os.path.join(S, "staged_A_v2", f"l3_k{k}", "iii_f.json")),
                    v3=best(os.path.join(S, "staged_A_v3", f"l3_k{k}", "iii.json")))
        for o, b in v5.items():
            r = sb[(k, o)]
            row = dict(k=k, orbit=o, spec=r["spec"], tau=r["tau"], staged_best=r["n"], run=b["n"], N=b["N"], stage=b["stage"],
                       synthillation=syn.get((k, r["spec"])), **{a: refs[a].get(o, {}).get("n") for a in refs})
            if row.get("v5") is not None and row.get("v4") is not None:
                row["min_v4_v5"] = min(row["v4"], row["v5"])
            rows.append(row)
            if r["n"] < b["n"] and r.get("gates"):
                cands.append((k, o, r["n"], r["N"], r["gates"], "staged_best.json (final pipeline)", r["stage"]))
            for a in refs:
                if o in refs[a] and refs[a][o]["n"] < b["n"]:
                    cands.append((k, o, refs[a][o]["n"], refs[a][o]["N"], refs[a][o]["gates"], a, refs[a][o]["stage"]))
    v5n = {(r["k"], r["orbit"]): r["run"] for r in rows}
    op = os.path.join(ROOT, "outputs", "outside_pipeline", "smaller_than_pipeline.csv")
    if os.path.exists(op):
        for r in csv.DictReader(open(op)):
            key = (int(r["k"]), r["orbit"])
            if key in v5n and int(r["n"]) < v5n[key]:
                cands.append((key[0], key[1], int(r["n"]), int(r["N"]), json.loads(r["gates"]), "outside_pipeline: " + r["source"], r["stage"]))
    # one entry per (class, smallest n), keep all sources that reach it
    import staged as STG, staged_class as SC, unified_search as us
    still = {}
    for c in cands:
        key = (c[0], c[1])
        if key not in still or c[2] < still[key]["n"]:
            still[key] = dict(k=c[0], orbit=c[1], n=c[2], N=c[3], gates=c[4], sources=[c[5]], how=[c[6]])
        elif c[2] == still[key]["n"] and c[5] not in still[key]["sources"]:
            still[key]["sources"].append(c[5]); still[key]["how"].append(c[6])
    for k, src in CELLS:
        mine = [x for x in still.values() if x["k"] == k]
        if not mine:
            continue
        cache, classes = SC.setup(3, k, os.path.join(ROOT, "outputs", "staged_A"))
        for x in mine:
            g = [(frozenset(S_), c) for S_, c in x["gates"]]; f = dict(k=k, r=x["N"] - k, gates=g)
            ok = dict(valid=us.is_factory(f, 3), brute=us.is_factory(f, 3, brute_force=True) if x["N"] <= 9 else None,
                      d=str(us.distance(f, 3)), n=sum(us.inputs(c, 3) for _, c in g), label_ok=STG.label_fix(f, 3, k, cache)[0] == x["orbit"])
            x["check"] = ok; x["ok"] = ok["valid"] and ok["brute"] is not False and ok["d"] != "1" and ok["label_ok"] and ok["n"] == x["n"]
    cmpc = lambda a, b: "lower" if a < b else "equal" if a == b else "higher"
    cnt = lambda key: {c: sum(1 for r in rows if r.get(key) is not None and cmpc(r["run"], r[key]) == c) for c in ("lower", "equal", "higher")}
    earlier = {}
    for x in cmp4["smaller_elsewhere"]:
        earlier[(x["k"], x["orbit"])] = min(earlier.get((x["k"], x["orbit"]), 99), x["n"])
    rec7 = sum(1 for k, o in ts["higher"] if v5n[(k, o)] <= sb[(k, o)]["n"])
    rec35 = sum(1 for k, o in ts["smaller_elsewhere"] if v5n[(k, o)] <= earlier[(k, o)])
    ws = [r for r in rows if r["synthillation"] is not None]
    summ = dict(classes=len(rows), vs_staged_best=cnt("staged_best"),
                inputs_saved_vs_staged_best=sum(max(0, r["staged_best"] - r["run"]) for r in rows),
                inputs_lost_vs_staged_best=sum(max(0, r["run"] - r["staged_best"]) for r in rows),
                below_tau_plus_1=sum(1 for r in rows if r["tau"] is not None and r["run"] < r["tau"] + 1),
                vs_v4=cnt("v4"), vs_A_v2=cnt("A_v2"), vs_v3=cnt("v3"),
                vs_v5=cnt("v5") if any(r.get("v5") is not None for r in rows) else None,
                vs_min_v4_v5=cnt("min_v4_v5") if any(r.get("min_v4_v5") is not None for r in rows) else None,
                recovered_7=f"{rec7} of {len(ts['higher'])}", recovered_35=f"{rec35} of {len(ts['smaller_elsewhere'])}",
                vs_synthillation=dict(classes=len(ws), fewer=sum(r["run"] < r["synthillation"] for r in ws),
                                      same=sum(r["run"] == r["synthillation"] for r in ws), more=sum(r["run"] > r["synthillation"] for r in ws)),
                still_smaller_than_v5=len(still), still_smaller_verified=sum(x["ok"] for x in still.values()))
    json.dump(dict(summary=summ, rows=rows, still_smaller=sorted(still.values(), key=lambda x: (x["k"], x["orbit"]))),
              open(os.path.join(out, "compare_" + os.path.basename(out.rstrip("/")).replace("staged_A_", "") + ".json"), "w"), indent=0)
    with open(os.path.join(out, "still_smaller_than_" + os.path.basename(out.rstrip("/")).replace("staged_A_", "") + ".csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["k", "orbit", "spec", "v5_n", "earlier_n", "N", "sources", "how", "verified", "gates"])
        for x in sorted(still.values(), key=lambda x: (x["k"], x["orbit"])):
            w.writerow([x["k"], x["orbit"], sb[(x["k"], x["orbit"])]["spec"], v5n[(x["k"], x["orbit"])], x["n"], x["N"],
                        " | ".join(x["sources"]), " | ".join(x["how"]), x["ok"], json.dumps(x["gates"])])
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
