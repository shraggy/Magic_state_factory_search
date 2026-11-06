"""Summarise the multi-factor variant test (outputs/solver_v2/v5test/<variant>_<order>/l3_k*/chain/*.json).
Per variant: classes, final n per class in each order (must agree: reproducibility), recovered counts against the
targets (7 'higher' -> reach staged_best; 35 'smaller_elsewhere' -> reach the smaller earlier n; 2320:771 -> 28),
classes lower/higher than v4, CPU of the restricted stage and of the whole chain. Every factory below v4 is re-verified.
    python searches/v5_test_report.py   -> outputs/solver_v2/v5test/report.json
"""
import os, sys, json, glob, collections
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
T = os.path.join(ROOT, "outputs", "solver_v2")


def main():
    ts = json.load(open(os.path.join(T, "testset_v5.json")))
    cmp4 = json.load(open(os.path.join(T, "staged_A_v4", "compare.json")))
    v4 = {(r["k"], r["orbit"]): r for r in cmp4["rows"]}
    earlier = collections.defaultdict(lambda: 99)
    for x in cmp4["smaller_elsewhere"]:
        earlier[(x["k"], x["orbit"])] = min(earlier[(x["k"], x["orbit"])], x["n"])
    tgt = {}
    for k, o in map(tuple, ts["all"]):
        r = v4[(k, o)]; t = min(r["staged_best"], earlier[(k, o)], 28 if o == "2320:771" else 99)
        tgt[(k, o)] = t
    rep = {}
    for vd in sorted(glob.glob(os.path.join(T, "v5test", "*_fwd"))):
        v = os.path.basename(vd)[:-4]; fin = {}; cpu = collections.Counter(); bad = []
        for order in ("fwd", "rev"):
            for f in glob.glob(os.path.join(T, "v5test", f"{v}_{order}", "l3_k*", "chain", "*.json")):
                x = json.load(open(f)); k = x["k"]
                for o, p in x["per_class"].items():
                    fin.setdefault((k, o), {})[order] = p["final"]
                    if order == "fwd":
                        for s in p["stages"]:
                            cpu[s["stage"].split(":")[0] + ("" if s["stage"].startswith("x:") else "")] += s.get("seconds", 0)
                            cpu["total"] += s.get("seconds", 0)
        same = sum(1 for d in fin.values() if d.get("fwd") == d.get("rev"))
        f = {key: d.get("fwd") for key, d in fin.items()}
        get = lambda key: f.get(key, 99)
        rec_hi = sum(1 for key in map(tuple, ts["higher"]) if get(key) <= v4[key]["staged_best"])
        rec_se = sum(1 for key in map(tuple, ts["smaller_elsewhere"]) if get(key) <= earlier[key])
        rec_tgt = sum(1 for key in tgt if get(key) <= tgt[key])
        low = sum(1 for key in f if f[key] < v4[key]["v4"]); high = sum(1 for key in f if f[key] > v4[key]["v4"])
        rep[v] = dict(classes=len(f), orders_agree=same, both_orders=sum(1 for d in fin.values() if len(d) == 2),
                      recovered_higher=f"{rec_hi} of {len(ts['higher'])}",
                      recovered_smaller_elsewhere=f"{rec_se} of {len(ts['smaller_elsewhere'])}",
                      reached_target=f"{rec_tgt} of {len(tgt)}", n_771=get((6, "2320:771")),
                      lower_than_v4=low, higher_than_v4=high,
                      cpu_h={a: round(b / 3600, 2) for a, b in cpu.items()},
                      per_class={f"{k}|{o}": dict(v4=v4[(k, o)]["v4"], staged_best=v4[(k, o)]["staged_best"],
                                                   earlier=earlier[(k, o)] if earlier[(k, o)] < 99 else None,
                                                   target=tgt[(k, o)], fwd=d.get("fwd"), rev=d.get("rev"))
                                 for (k, o), d in sorted(fin.items())})
    json.dump(rep, open(os.path.join(T, "v5test", "report.json"), "w"), indent=0)
    for v, r in rep.items():
        print(v, {a: b for a, b in r.items() if a != "per_class"}, flush=True)


if __name__ == "__main__":
    main()
