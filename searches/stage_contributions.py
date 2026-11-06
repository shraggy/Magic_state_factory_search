"""Table I: contribution of each stage of the staged cascade at l=3 over the 766 output classes with k <= 6, and the
per-cell stage table of the Supplement (Table SII), both from outputs/staged_best.json and the per-pass files it names.

Table I, one row per stage form, in cascade order: reached = classes the form alone gives a factory for (its own
finds, before the best of earlier stages is applied; for i_tg, i_gs, i_gh, ii_dj the cell sweep files, for later
passes every factory recorded for the class); best = classes whose final n it reaches first (ties to the earliest);
fewer = of the classes with fewer inputs than synthillation, those it reaches first.
Table SII, one row per cell: classes, classes reached after (i), (ii), (iii), and the CPU time of each stage.
At l=3 the cells come from outputs/staged_final_b3 (build_staged_best3.py: stage (iii) per class the best of the three
deterministic runs v4, v5, v6) and the own finds of stage (iii) from the per-class records of all three runs (xs, xsb =
restricted, x, x2 = free). The stage-(iii) seconds of Table SII are those of v4, as build_staged_best3.py writes them.
If the staged-best file given by --sb has its l=3 cells in outputs/staged_final (the earlier v4-only build, archived in
outputs/archive/v4_published/), only v4's records are read and the outputs carry the suffix _v4 (no manuscript copy).
    python searches/stage_contributions.py   -> outputs/stage_contributions_l3.json, outputs/checks/stage_table.{json,tex},
                                                Magic_State_Factories_tex/MSD_manuscript/stage_table.tex
    python searches/stage_contributions.py --out-dir DIR
        -> DIR/stage_contributions_l3.json, DIR/stage_table.{json,tex} (nothing under outputs/ or the manuscript)
"""
import os, sys, json, glob, collections, argparse
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
_a = argparse.ArgumentParser(); _a.add_argument("--sb", default="staged_best.json")
_a.add_argument("--out-dir", default=None, help="write stage_contributions_l3.json, stage_table.{json,tex} here (no suffix, no manuscript copy)")
_A = _a.parse_args()
SB_NAME = _A.sb; OUTD = _A.out_dir
SB = json.load(open(os.path.join(O, SB_NAME)))
B3 = os.path.basename(os.path.dirname(SB["cells"]["l3_k1"]["source"])) == "staged_final_b3"   # else the v4-only build
SUF = "" if B3 else "_v4"
P_CONTRIB = os.path.join(OUTD, "stage_contributions_l3.json") if OUTD else os.path.join(O, f"stage_contributions_l3{SUF}.json")
P_TJSON = os.path.join(OUTD, "stage_table.json") if OUTD else os.path.join(O, "checks", f"stage_table{SUF}.json")
P_TTEX = os.path.join(OUTD, "stage_table.tex") if OUTD else os.path.join(O, "checks", f"stage_table{SUF}.tex")
III_FORM = lambda st: "iii_r" if st.startswith("xs") else ("iii_f" if st.startswith("x") else None)
III_RUNS = ["v4", "v5", "v6"] if B3 else ["v4"]
FORMS = [("i_tg", "(i) two-group"), ("i_gs", "(i) group-symmetric"), ("i_gh", "(i) check groups"),
         ("ii_dj", "(ii) disjoint"), ("ii_eg", "(ii) entangled"), ("iii_r", "(iii) restricted"), ("iii_f", "(iii) free")]
STAGE = {"i_tg": "i", "i_gs": "i", "i_gh": "i", "ii_dj": "ii", "ii_eg": "ii", "iii_r": "iii", "iii_f": "iii"}


def own_finds(L, k):
    """{pass: {orbit: smallest n that pass itself found}} for one cell."""
    src = os.path.join(ROOT, SB["cells"][f"l{L}_k{k}"]["source"]); out = collections.defaultdict(dict)
    def put(p, o, n):
        if n is not None and (o not in out[p] or n < out[p][o]): out[p][o] = n
    for p in ("i_tg", "i_gs", "i_gh", "ii_dj"):
        f = os.path.join(src, "cell", f"{p}.json")
        if os.path.exists(f):
            for o, b in json.load(open(f))["best"].items(): put(p, o, b["n"])
    for f in glob.glob(os.path.join(src, "cls", "*.json")) + glob.glob(os.path.join(src, "free_*", "*.json")):
        res = json.load(open(f))["results"]
        for r in res if isinstance(res, list) else eval(res):
            put(r["stage"].split(":")[0], str(r["orbit"]), r["n"])
    if L == 3 and os.path.basename(src) == f"l3_k{k}" and "staged_final" in src:   # stage (iii) of the final pipeline
        assert os.path.basename(os.path.dirname(src)) == ("staged_final_b3" if B3 else "staged_final"), src
        for run in III_RUNS:                                                            # per-class records of each form
            chain = os.path.join(O, "solver_v2", f"staged_A_{run}", f"l3_k{k}", "chain")
            for f in glob.glob(os.path.join(chain, "*.json")):
                for o, v in json.load(open(f))["per_class"].items():
                    for st in v["stages"]:
                        if III_FORM(st["stage"]) and st.get("n") is not None and "skipped" not in st:
                            put(III_FORM(st["stage"]), o, st["n"])
        a = os.path.join(ROOT, "outputs", "staged_A", f"l3_k{k}") if k >= 5 else os.path.join(ROOT, "outputs", "staged_pc", f"l3_k{k}")
        for p in ("i_tg", "i_gs", "i_gh", "ii_dj"):                                      # stages (i)-(ii): sweep A's own
            f = os.path.join(a, "cell", f"{p}.json")
            if os.path.exists(f):
                for o, b in json.load(open(f))["best"].items(): put(p, o, b["n"])
        for f in glob.glob(os.path.join(a, "cls", "*.json")):
            res = json.load(open(f))["results"]
            for r in res if isinstance(res, list) else eval(res):
                if not r["stage"].startswith("iii"): put(r["stage"].split(":")[0], str(r["orbit"]), r["n"])
        if not glob.glob(os.path.join(a, "cls", "*.json")):                               # k <= 4: one-process passes
            for p in ("i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg"):
                f = os.path.join(src, f"{p}.json")
                if os.path.exists(f):
                    d = json.load(open(f))
                    for o in d["new"]: put(p, o, d["best"][o]["n"])
                    for o, (_, b) in d["lowered"].items(): put(p, o, b)
        return out
    if not glob.glob(os.path.join(src, "cls", "*.json")):     # one-process cells: the 'new' and 'lowered' of each pass
        for p, _ in FORMS:
            f = os.path.join(src, f"{p}.json")
            if not os.path.exists(f): continue
            d = json.load(open(f))
            for o in d["new"]: put(p, o, d["best"][o]["n"])
            for o, (_, b) in d["lowered"].items(): put(p, o, b)
    return out


def main():
    rows = [r for r in SB["rows"] if r["l"] == 3]
    first = collections.Counter(r["first"] for r in rows)
    fewer = collections.Counter(r["first"] for r in rows if r["n"] < r["syn_n"])
    reach = collections.Counter(); match = collections.Counter()      # match: the form's own n equals the final n
    fin = {(r["k"], str(r["orbit"])): r["n"] for r in rows}
    for k in range(1, 7):
        for p, d in own_finds(3, k).items():
            reach[p] += len(d); match[p] += sum(n == fin[(k, o)] for o, n in d.items())
    T1 = [dict(form=f, label=lab, reached=reach[f], matches_final=match[f], best=first[f], fewer=fewer[f]) for f, lab in FORMS]
    tot = dict(classes=len(rows), fewer=sum(r["n"] < r["syn_n"] for r in rows), same=sum(r["n"] == r["syn_n"] for r in rows),
               more=sum(r["n"] > r["syn_n"] for r in rows), at_bound=sum(r["n"] == r["tau"] + 1 for r in rows),
               fewer_by=dict(collections.Counter(r["syn_n"] - r["n"] for r in rows if r["n"] < r["syn_n"])),
               fewer_by_k=dict(collections.Counter(r["k"] for r in rows if r["n"] < r["syn_n"])),
               gap_ge17=sum(r["n"] - r["tau"] - 1 >= 17 for r in rows))
    json.dump(dict(table=T1, totals=tot, source=f"outputs/{SB_NAME}", iii_runs=III_RUNS),
              open(P_CONTRIB, "w"), indent=1)
    # Table SII
    T2 = []
    for key, m in SB["cells"].items():
        L, k = int(key[1]), int(key.split("_k")[1])
        P = {}
        for p, _ in FORMS:
            if p == "iii_r" and L != 3:                       # restricted SAT is part of the pipeline at l=3 only
                continue
            f = os.path.join(ROOT, m["source"], f"{p}.json")
            if os.path.exists(f): P[p] = json.load(open(f))
        after = {st: max([P[p]["reached"] for p in P if STAGE[p] == st], default=None) for st in ("i", "ii", "iii")}
        cpu = {st: round(sum(P[p]["seconds"] or 0 for p in P if STAGE[p] == st) / 3600, 2) for st in ("i", "ii", "iii")}
        T2.append(dict(l=L, k=k, classes=m["classes"], after=after, cpu_hours=cpu,
                       cpu_sum=any(P[p].get("seconds_are_cpu_sum") for p in P), source=m["source"]))
    json.dump(dict(rows=T2, source=f"outputs/{SB_NAME}"), open(P_TJSON, "w"), indent=1)
    h = lambda x: "--" if x is None else (f"{x:.1f}" if x >= 0.05 else "$<0.1$")
    fmt = lambda x: "--" if x is None else f"{x:,}".replace(",", "{,}")
    tex = [r"\begin{tabular}{c c r r r r r r r}", r"\toprule",
           r" & & & \multicolumn{3}{c}{classes reached after} & \multicolumn{3}{c}{CPU hours} \\",
           r"\cmidrule(lr){4-6}\cmidrule(lr){7-9}", r"$l$ & $k$ & classes & (i) & (ii) & (iii) & (i) & (ii) & (iii) \\", r"\midrule"]
    for r in sorted(T2, key=lambda r: (r["l"], r["k"])):
        a, c = r["after"], r["cpu_hours"]
        tex.append(f"{r['l']} & {r['k']} & {fmt(r['classes'])} & {fmt(a['i'])} & {fmt(a['ii'])} & {fmt(a['iii'])} & "
                   f"{h(c['i'])} & {h(c['ii'])} & {h(c['iii'])} \\\\")
    tex += [r"\bottomrule", r"\end{tabular}"]
    body = "% generated by searches/stage_contributions.py\n" + "\n".join(tex) + "\n"
    open(P_TTEX, "w").write(body)
    if B3 and not OUTD:
        open(os.path.join(ROOT, "Magic_State_Factories_tex", "MSD_manuscript", "stage_table.tex"), "w").write(body)   # the copy main.tex inputs
    print(json.dumps(T1, indent=0)); print(tot)


if __name__ == "__main__":
    main()
