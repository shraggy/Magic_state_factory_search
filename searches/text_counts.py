"""Counts quoted in the paper's prose that no other script writes (entangled / product classes among those with fewer
inputs than synthillation, N = 8 or 9, gap to the bound >= 17, equal n with the same N as synthillation, per-form
lowerings of stage (iii), examples), from outputs/staged_best.json (the per-class best of the stage-(iii) runs v4, v5, v6,
key "b3") and, for comparison, the archived v4-only build outputs/archive/v4_published/outputs/staged_best.json ("old").
Per-form lowerings in the winning run's own order: outputs/staged_final_b3/l3_k*/attribution.json.
    python searches/text_counts.py   -> outputs/checks/text_counts.json
"""
import json, os, statistics, collections, glob
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = R + "/outputs"
OLD = "archive/v4_published/outputs/staged_best.json"
def supports(spec):
    out = []
    for t in spec.split(","):
        nm = t.rstrip("0123456789"); out.append(frozenset(int(c) for c in t[len(nm):]))
    return out
def components(spec):
    S = supports(spec); par = {}
    def f(x):
        while par.setdefault(x, x) != x: x = par[x]
        return x
    for s in S:
        s = sorted(s)
        for q in s[1:]: par[f(q)] = f(s[0])
        f(s[0])
    return len({f(q) for s in S for q in s})
def types(spec): return {t.rstrip("0123456789") for t in spec.split(",")}
def proofs():
    """(k, orbit) -> max N on which free_coef* z3 query is unsat (stage-(iii)-independent proofs)."""
    mx = {}
    for src in ("staged_A", "staged_pc"):
        for f in glob.glob(os.path.join(O, src, "l3_k*", "free_coef*", "*.json")):
            d = json.load(open(f)); k = int(d["k"])
            for o, v in d["per_class"].items():
                m = max([s["N"] for s in v.get("steps") or [] if s.get("last") == "unsat"], default=0)
                mx[(k, o)] = max(mx.get((k, o), 0), m)
    return mx
PR = proofs()
def stats(name):
    SB = json.load(open(os.path.join(O, name))); rows = [r for r in SB["rows"] if r["l"] == 3]
    fin = {(r["k"], str(r["orbit"])): r for r in rows}
    fewer = [r for r in rows if r["n"] < r["syn_n"]]
    ent = [r for r in fewer if components(r["spec"]) == 1]
    ex = {}
    for sp in ("CCZ156,CCZ234,CS12", "CS01,CS23,CS45", "CS01,CS23,T4", "CCZ015,CCZ234,CS01"):
        ex[sp] = [(r["k"], r["n"], r["syn_n"], r["N"]) for r in rows if r["spec"] == sp]
    # stage (iii) "what the stage adds": classes lowered by iii_r below ii_eg, of those fewer than syn, by iii_f below iii_r
    def cum(p):
        out = {}
        for key, m in SB["cells"].items():
            if not key.startswith("l3_"): continue
            k = int(key.split("_k")[1]); f = os.path.join(R, m["source"], f"{p}.json")
            for o, b in json.load(open(f))["best"].items(): out[(k, o)] = b["n"]
        return out
    eg, ir, if_ = cum("ii_eg"), cum("iii_r"), cum("iii_f")
    low_r = [key for key in fin if ir[key] < eg[key]]; low_f = [key for key in fin if if_[key] < ir[key]]
    syn = {key: fin[key]["syn_n"] for key in fin}
    onN = [key for key in fin if PR.get(key, 0) >= fin[key]["N"] and fin[key]["n"] != fin[key]["tau"] + 1]
    return dict(file="outputs/" + name, classes=len(rows), fewer=len(fewer), same=sum(r["n"] == r["syn_n"] for r in rows),
                more=sum(r["n"] > r["syn_n"] for r in rows), inputs_fewer_total=sum(r["syn_n"] - r["n"] for r in fewer),
                fewer_entangled=len(ent), fewer_product=len(fewer) - len(ent),
                fewer_product_specs=sorted(r["spec"] for r in fewer if components(r["spec"]) > 1),
                mixed_types=sum(len(types(r["spec"])) >= 2 for r in rows),
                N_8_or_9=sum(r["N"] in (8, 9) for r in rows), N_hist=dict(sorted(collections.Counter(r["N"] for r in rows).items())),
                gap_ge17_of_761=sum(r["n"] - r["tau"] - 1 >= 17 for r in rows if r["n"] != r["tau"] + 1),
                on_bound=sum(r["n"] == r["tau"] + 1 for r in rows),
                iii_r_lowered_below_ii_eg=len(low_r), iii_r_lowered_and_fewer_than_syn=sum(ir[key] < syn[key] for key in low_r),
                iii_f_lowered_below_iii_r=len(low_f), iii_f_lowered_all_fewer_than_syn=all(if_[key] < syn[key] for key in low_f),
                iii_f_lowered_fewer_than_syn=sum(if_[key] < syn[key] for key in low_f),
                iii_distinct_lowered=len(set(low_r) | set(low_f)),
                smallest_on_own_N_by_free_coef=len(onN),
                proofs_N7=sum(PR.get(key, 0) >= 7 for key in fin), proofs_N6=sum(PR.get(key, 0) == 6 for key in fin),
                examples=ex)
SYN = {(x["k"], x["spec"]): x for x in json.load(open(O + "/synthillation_strongest_l3.json"))}
def extra(name, att_dir=None):
    rows = [r for r in json.load(open(os.path.join(O, name)))["rows"] if r["l"] == 3]
    eq = [r for r in rows if r["n"] == r["syn_n"]]
    d = dict(equal_n_same_N=sum(SYN[(r["k"], r["spec"])]["best_N"] == r["N"] for r in eq),
             equal_n_one_more_N=sum(r["N"] - SYN[(r["k"], r["spec"])]["best_N"] == 1 for r in eq))
    if att_dir:
        fin = {(r["k"], str(r["orbit"])): r for r in rows}; A = {}
        for k in range(1, 7): A.update({(k, o): v for o, v in json.load(open(f"{O}/{att_dir}/l3_k{k}/attribution.json")).items()})
        rr = {x for x, v in A.items() if v["saved_restricted"] > 0}; ff = {x for x, v in A.items() if v["saved_free"] > 0}
        d.update(form_restricted_lowered=len(rr), form_restricted_lowered_end_fewer=sum(fin[x]["n"] < fin[x]["syn_n"] for x in rr),
                 form_free_lowered=len(ff), form_free_lowered_end_fewer=sum(fin[x]["n"] < fin[x]["syn_n"] for x in ff),
                 form_lowered_union=len(rr | ff), form_lowered_both=len(rr & ff))
    return d
out = dict(old=stats(OLD), b3=stats("staged_best.json"),
           old_extra=extra(OLD), b3_extra=extra("staged_best.json", "staged_final_b3"),
           note="entangled = support hypergraph of the class representative's gates is connected on all k outputs; "
                "pass counts are differences of consecutive cumulative pass files (for b3 see the limits in build_staged_best3.py)")
json.dump(out, open(os.path.join(O, "checks", "text_counts.json"), "w"), indent=1)
for k in out["old"]:
    if k not in ("examples", "fewer_product_specs", "N_hist"): print(f"{k:42s} {out['old'][k]!s:>30s} {out['b3'][k]!s:>30s}")
print(out["old_extra"], out["b3_extra"]); print(out["old"]["examples"], out["b3"]["examples"]); print(out["old"]["fewer_product_specs"]); print(out["b3"]["fewer_product_specs"])
