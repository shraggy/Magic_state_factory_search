"""Stage (iii) at l=3, k<=6 with different SAT layers: the same cascade (stages (i)-(ii) of sweep A, the gluing best as
seed, history-independent parts, a fresh process per query, CryptoMiniSat at 2e6 conflicts, x86_64), run three times with
a different list of SAT layers. Each layer starts from the running best n and passes only that n on (chain_iii2.run), so
a layer that lowers no class can be dropped without changing any later layer.
  v4  outputs/solver_v2/staged_A_v4   xs:None, xs:6, x:5, x:6, x:None            (its last form, z3 syn, lowered none)
  v5  outputs/solver_v2/staged_A_v5   xsb:8:None, xsb:8:6, x:5, x:6, x:None
  v6  outputs/solver_v2/staged_A_v6   xs:None, xs:6, x:5, x:6, x:None, xsb:8:None, x2:5, x2:6, x2:None
xs = restricted SAT, one synthesis; xsb:8 = restricted SAT, up to 8 syntheses; x / x2 = free SAT; the number is the
support-width limit. Per variant: layers, CPU-h per form, classes with fewer inputs than the strongest synthillation
construction (synthillation_strongest_l3.json) and inputs fewer in total; also the per-class best over the three runs.
Also the same counts for stages (i)-(ii) alone (no SAT), and checks that the layers of v6 that lowered no class never changed any class's n.
    python searches/sat_variants.py   -> outputs/checks/sat_variants.json
"""
import os, json, glob
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
RUNS = {v: os.path.join(O, "solver_v2", f"staged_A_{v}") for v in ("v4", "v5", "v6")}
FORM = lambda st: "restricted" if st.startswith("xs") else ("free" if st.startswith("x") else "other")


def final_n(d):
    out = {}
    for k in range(1, 7):
        for o, b in json.load(open(os.path.join(d, f"l3_k{k}", "iii.json")))["best"].items():
            out[(k, o)] = int(b["n"])
    return out


def main():
    syn = {}
    for r in json.load(open(os.path.join(O, "staged_best.json")))["rows"]:
        if r["l"] == 3: syn[(r["k"], str(r["orbit"]))] = r["syn_n"]
    res, finals = {}, {}
    for v, d in RUNS.items():
        s = json.load(open(os.path.join(d, "summary.json")))
        chain = next(x["chain"] for x in s.values() if isinstance(x, dict) and "chain" in x)
        per = {}
        for x in s.values():
            if isinstance(x, dict) and "per_stage" in x:
                for st, y in x["per_stage"].items():
                    p = per.setdefault(st, dict(lowered=0, inputs_saved=0, cpu_h=0.0))
                    p["lowered"] += y["lowered"]; p["inputs_saved"] += y["inputs_saved"]; p["cpu_h"] += y["cpu_s"] / 3600
        sat = [st for st in chain if FORM(st) != "other"]
        f = finals[v] = final_n(d); assert set(f) == set(syn) and len(f) == 766
        cpu = {fm: round(sum(per[st]["cpu_h"] for st in sat if FORM(st) == fm), 1) for fm in ("restricted", "free")}
        res[v] = dict(layers=sat, n_layers=len(sat), n_restricted=sum(FORM(st) == "restricted" for st in sat),
                      n_free=sum(FORM(st) == "free" for st in sat), cpu_h=cpu, cpu_h_total=round(sum(cpu.values()), 1),
                      per_layer={st: dict(per[st], cpu_h=round(per[st]["cpu_h"], 1)) for st in chain},
                      fewer=sum(f[c] < syn[c] for c in f), inputs_fewer=sum(syn[c] - f[c] for c in f if f[c] < syn[c]))
    best = {c: min(finals[v][c] for v in finals) for c in syn}
    res["best_of_three"] = dict(fewer=sum(best[c] < syn[c] for c in best),
                                inputs_fewer=sum(syn[c] - best[c] for c in best if best[c] < syn[c]),
                                lower_than_v6=sorted(f"{k}:{o}" for (k, o) in best if best[(k, o)] < finals["v6"][(k, o)]))
    pre = {}                                                  # stages (i)-(ii) only: the seed of stage (iii), shared by all runs
    for k in range(1, 7):
        for ps in ("i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg"):
            for o, b in json.load(open(os.path.join(RUNS["v6"], f"l3_k{k}", f"{ps}.json")))["best"].items():
                pre[(k, o)] = min(pre.get((k, o), 10 ** 9), int(b["n"]))
    assert set(pre) == set(syn)
    res["without_sat"] = dict(fewer=sum(pre[c] < syn[c] for c in syn), same=sum(pre[c] == syn[c] for c in syn),
                              more=sum(pre[c] > syn[c] for c in syn),
                              inputs_fewer=sum(syn[c] - pre[c] for c in syn if pre[c] < syn[c]),
                              wins_created_by_sat=sum(best[c] < syn[c] <= pre[c] for c in syn),
                              wins_deepened_by_sat=sum(best[c] < pre[c] < syn[c] for c in syn),
                              more_fixed_by_sat=sum(pre[c] > syn[c] >= best[c] for c in syn))
    idle = [st for st, p in res["v6"]["per_layer"].items() if p["lowered"] == 0 and FORM(st) != "other"]
    tot = changed = 0
    for fn in glob.glob(os.path.join(RUNS["v6"], "l3_k*", "chain", "*.json")):
        for v in json.load(open(fn))["per_class"].values():
            for st in v["stages"]:
                if st["stage"] in idle and "skipped" not in st:
                    tot += 1; changed += st["n"] != st["seed"]
    res["v6_idle_layers"] = dict(layers=idle, runs=tot, changed_n=changed)
    json.dump(res, open(os.path.join(O, "checks", "sat_variants.json"), "w"), indent=1)
    for v in ("v4", "v5", "v6"):
        r = res[v]; print(v, r["n_layers"], r["cpu_h"], r["cpu_h_total"], r["fewer"], r["inputs_fewer"])
    print("best of three", res["best_of_three"]); print("without SAT", res["without_sat"]); print("v6 idle", res["v6_idle_layers"])


if __name__ == "__main__":
    main()
