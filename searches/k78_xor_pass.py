"""Merge k78_xor.py chunk files (one or more variants) into the stage (iii) pass file iii_f.json (staged.py schema),
seeded with the gluing pass ii_eg.json: per class keep the smaller (n, N) (credit stays with the earlier pass on a tie).
    python searches/k78_xor_pass.py DEST RUNDIR [RUNDIR ...]  -> DEST/l3_k{7,8}/iii_f.json
"""
import os, sys, json, glob
dest, runs = sys.argv[1], sys.argv[2:]
R = []
for d in runs:
    for f in sorted(glob.glob(os.path.join(d, "c*of*.json"))):
        x = json.load(open(f))
        for r in x["results"]:
            r["variant"] = os.path.basename(d.rstrip("/")); R.append(r)
for k in (7, 8):
    prev = json.load(open(os.path.join(dest, f"l3_k{k}", "ii_eg.json")))
    best = {o: dict(b) for o, b in prev["best"].items()}; pn = {o: b["n"] for o, b in best.items()}
    rows, cpu = [], 0.0
    for r in R:
        if r["k"] != k:
            continue
        c = sum(u["cpu_s"] for u in r["runs"]); cpu += c
        rows.append(dict(orbit=r["fp"], id=r["id"], variant=r["variant"], cap=r["cap"], skipped=r.get("skipped"),
                         runs=[dict(N=u["N"], status=u["status"], n=u["n"], verified=u.get("verified"), cpu_s=u["cpu_s"]) for u in r["runs"]],
                         n=r["best"] and r["best"]["n"], cpu_s=round(c, 1)))
        b = r["best"]; o = r["fp"]
        if b and (o not in best or (b["n"], b["N"]) < (best[o]["n"], best[o]["N"])):
            old = best.get(o, {})
            best[o] = dict(old, n=b["n"], N=b["N"], stage="iii_f:xor:" + r["variant"], gates=b["gates"], ccz=r["ccz"], id=r["id"])
            if old.get("tau") and old["tau"][1] is not None:
                best[o]["class_on_bound"] = bool(old.get("tau_exact") and b["n"] == old["tau"][1] + 1)
    new = sorted(o for o in best if o not in pn); low = {o: [pn[o], best[o]["n"]] for o in best if o in pn and best[o]["n"] < pn[o]}
    rec = dict(L=3, k=k, pass_="iii_f", settings=dict(variants=runs), seconds=round(cpu, 1), seconds_are_cpu_sum=True,
               classes=prev["classes"], reached=len(best), new=new, lowered=low, inputs_saved=sum(a - b for a, b in low.values()),
               prev="ii_eg", best=best, rows=rows, specs=prev["specs"],
               on_bound_classes=sorted(o for o, b in best.items() if b.get("class_on_bound")))
    json.dump(rec, open(os.path.join(dest, f"l3_k{k}", "iii_f.json"), "w"))
    print(f"k={k}: {len(rows)} target runs, reached {len(best)}/{prev['classes']}, +{len(new)} new, lowered {len(low)} by "
          f"{rec['inputs_saved']}; {rec['seconds']} CPU-s", flush=True)
