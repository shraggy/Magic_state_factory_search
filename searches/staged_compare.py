"""Known-answer check of the staged cascade: per class, the final best n of outputs/staged/l{L}_k{k}/ against every
earlier source (cascade v4 and v3 runs, and at l=3 the merged catalogue outputs/l3_catalogue_k6.csv and the
group-symmetric run outputs/checks/orbit_group_l3.json). Reports, per cell, N of M classes equal / better / worse,
and lists every worse class with the source that beats it.
    python searches/staged_compare.py --L 3 --k 1 2 3   -> outputs/staged/compare_l{L}.json
"""
import os, sys, json, csv, argparse
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")


def refs(L, k):
    R = {}                                                      # spec -> (n, source)
    def put(spec, n, src):
        if n and (spec not in R or n < R[spec][0]):
            R[spec] = (n, src)
    for d in ("cascade_v4", "cascade_main_v3", "cascade"):
        p = os.path.join(O, d, f"l{L}_k{k}.json")
        if os.path.exists(p):
            s = json.load(open(p))["summary"]
            if s.get("classes"):
                for o, b in s["best"].items():
                    put(b["spec"], b["n"], d)
    if L == 3:
        for r in csv.DictReader(open(os.path.join(O, "l3_catalogue_k6.csv"))):
            if int(r["k"]) == k:
                put(r["spec"], int(r["n"]), "catalogue:" + r["source"])
        for x in json.load(open(os.path.join(O, "checks", "orbit_group_l3.json")))["classes"]:
            if x["k"] == k and x.get("og_n"):
                put(x["spec"], x["og_n"], "orbit_group")
    return R


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--L", type=int); ap.add_argument("--k", type=int, nargs="+")
    ap.add_argument("--out", default=os.path.join(O, "staged")); a = ap.parse_args()
    res = {}
    for k in a.k:
        last = None
        for p in ("i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg", "iii_r", "iii_f"):
            f = os.path.join(a.out, f"l{a.L}_k{k}", f"{p}.json")
            if os.path.exists(f):
                last = json.load(open(f)); lastp = p
        if last is None:
            continue
        R = refs(a.L, k); eq = bt = wr = miss = 0; worse = []
        for o, spec in last["specs"].items():
            b = last["best"].get(o); r = R.get(spec)
            if b is None:
                miss += 1; worse.append(dict(spec=spec, staged=None, ref=r)); continue
            if r is None or b["n"] < r[0]:
                bt += 1
            elif b["n"] == r[0]:
                eq += 1
            else:
                wr += 1; worse.append(dict(spec=spec, staged=b["n"], stage=b["stage"], ref=r))
        res[k] = dict(classes=len(last["specs"]), equal=eq, better=bt, worse=wr, missing=miss, last_pass=lastp, worse_list=worse)
        print(f"l={a.L} k={k}: {len(last['specs'])} classes | equal {eq}, better {bt}, worse {wr}, missing {miss} (last pass {lastp})")
        for w in worse[:10]:
            print("   ", w)
    json.dump(res, open(os.path.join(a.out, f"compare_l{a.L}.json"), "w"), indent=1)
