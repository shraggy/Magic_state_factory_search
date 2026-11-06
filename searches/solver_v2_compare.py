"""Compare the solver_v2 runs (xor_free.py variants, syn_eg.py) with sweep A, per class at l=3 k=5,6.

For every variant directory under OUT/l3_k{k}/ (xor_*, syn_eg): classes done, total and per-class seconds, the best
n it credits to each class (any factory it found, labelled with the class it deposits), and per class:
  ref    = n in outputs/staged_best.json (published sweep A best)
  seed   = the seed cap the variant started from
  own    = best n from the variant (seeded runs report only factories below the seed)
  final  = min(seed, own): the pass's result if it replaced its stage
Counts: classes where own < ref (lowered), own == ref, own > ref or none (did not reach ref). CP-SAT reference: the
per-class seconds of outputs/staged_A/l3_k{k}/free_bool_w5_d600/ and free_coef_1.5e9_unseeded/.
Every factory with own < ref is re-verified here: us.is_factory (monomial and brute force), us.distance != 1, and the
exact class label (staged.label_fix) equal to the class.
    python searches/solver_v2_compare.py OUT          -> OUT/compare.json, OUT/lowered.json
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))


def per_class_secs(pat):
    s = {}
    for f in glob.glob(pat):
        for o, p in json.load(open(f))["per_class"].items():
            s[o] = p.get("seconds", 0.0)
    return s


def main(out, verify=True):
    ref = {}
    for r in json.load(open(os.path.join(ROOT, "outputs", "staged_best.json")))["rows"]:
        if r["l"] == 3 and r["k"] in (5, 6):
            ref[(r["k"], r["orbit"])] = r["n"]
    rep, lowered = {}, []
    for k in (5, 6):
        A = os.path.join(ROOT, "outputs", "staged_A", f"l3_k{k}")
        cp = per_class_secs(os.path.join(A, "free_bool_w5_d600", "*.json"))
        z3 = per_class_secs(os.path.join(A, "free_coef_1.5e9_unseeded", "*.json"))
        for vd in sorted(glob.glob(os.path.join(out, f"l3_k{k}", "*"))):
            if not os.path.isdir(vd):
                continue
            v = os.path.basename(vd); done, secs, own, per = set(), {}, {}, {}
            for f in glob.glob(os.path.join(vd, "*.json")):
                x = json.load(open(f)); done |= set(x["classes"])
                for o, p in x["per_class"].items():
                    secs[o] = p.get("seconds", 0.0); per[o] = p
                for r in x["results"]:
                    o = str(r["orbit"])
                    if o not in own or r["n"] < own[o]["n"]:
                        own[o] = r
            ran = [o for o in done if not per.get(o, {}).get("skipped")]
            low = sorted(o for o in own if (k, o) in ref and own[o]["n"] < ref[(k, o)])
            eq = sorted(o for o in own if (k, o) in ref and own[o]["n"] == ref[(k, o)])
            # final n of the pass = min(seed, own); classes (searched) where that stays above the published n
            fin = {o: min([x for x in (per[o].get("seed"), own.get(o, {}).get("n")) if x is not None], default=None)
                   for o in ran}
            miss = sorted(o for o in ran if fin[o] is None or fin[o] > ref[(k, o)])
            own_miss = sorted(o for o in ran if (k, o) in ref and (o not in own or own[o]["n"] > ref[(k, o)])
                              and per[o].get("seed") is not None and per[o]["seed"] > ref[(k, o)])
            rep[f"k{k}/{v}"] = dict(classes_done=len(done), searched=len(ran), skipped=len(done) - len(ran),
                                    cpu_seconds=round(sum(secs.values()), 1),
                                    cpsat_bool_seconds_same_classes=round(sum(cp.get(o, 0) for o in ran), 1),
                                    z3_coef_seconds_same_classes=round(sum(z3.get(o, 0) for o in ran), 1),
                                    max_class_seconds=max(secs.values(), default=0),
                                    lowered=len(low), equal=len(eq), final_above_ref=len(miss),
                                    final_above_ref_list={o: [ref[(k, o)], fin[o]] for o in miss},
                                    seed_above_ref_not_recovered=len(own_miss),
                                    lowered_list={o: [ref[(k, o)], own[o]["n"], own[o]["N"]] for o in low})
            for o in low:
                lowered.append(dict(own[o], k=k, orbit=o, variant=v, ref=ref[(k, o)]))
    if verify and lowered:
        import staged as STG, staged_class as SC, unified_search as us
        for k in (5, 6):
            mine = [r for r in lowered if r["k"] == k]
            if not mine:
                continue
            cache, classes = SC.setup(3, k, os.path.join(ROOT, "outputs", "staged_A"))
            for r in mine:
                g = [(frozenset(S), c) for S, c in r["gates"]]; f = dict(k=k, r=r["N"] - k, gates=g)
                r["check"] = dict(is_factory=us.is_factory(f, 3), brute=us.is_factory(f, 3, brute_force=True),
                                  distance=str(us.distance(f, 3)), n=sum(us.inputs(c, 3) for _, c in g),
                                  label=STG.label_fix(f, 3, k, cache)[0], spec=classes[r["orbit"]]["spec"])
                r["check"]["ok"] = (r["check"]["is_factory"] and r["check"]["brute"] and r["check"]["distance"] != "1"
                                    and r["check"]["label"] == r["orbit"] and r["check"]["n"] == r["n"])
    json.dump(rep, open(os.path.join(out, "compare.json"), "w"), indent=1)
    json.dump(lowered, open(os.path.join(out, "lowered.json"), "w"), indent=0)
    for v, x in rep.items():
        print(v, {a: b for a, b in x.items() if a not in ("lowered_list", "final_above_ref_list")}, x["lowered_list"],
              "above:", x["final_above_ref_list"])
    best = {}
    for r in lowered:
        key = (r["k"], r["orbit"])
        if key not in best or r["n"] < best[key]["n"]:
            best[key] = r
    print(f"{len(best)} classes lowered below staged_best by some variant; verified "
          f"{sum(1 for r in best.values() if r.get('check', {}).get('ok'))} of {len(best)}")


if __name__ == "__main__":
    main(sys.argv[1], verify="--noverify" not in sys.argv)
