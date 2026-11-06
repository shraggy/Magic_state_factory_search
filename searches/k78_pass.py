"""Merge the k78_eg.py chunk files into a pipeline pass file with the schema of staged.py (best per class, new, lowered,
inputs_saved, seconds, settings, rows, specs), so it can seed later passes at k = 7, 8 the way ii_eg seeds iii_r/iii_f
at k <= 6. "Class" here is the GL(k,2)-invariant fingerprint id of outputs/k78_pipeline/targets_all.json
(searches/k78_orbits.py): layers with different fingerprints are in different classes; two fingerprint ids are never
merged, but one id could in principle hold more than one class (not proved injective at k = 7, 8).
The previous pass is sweep B stage (i): its CCZ-layer outputs (outputs/k78_pipeline/staged_B_ccz_best.json, best per
fingerprint, re-verified). The pass keeps per class the smallest (n, N) over B and every target of the class, and
records per class the T-count interval [isotropy floor, best synthesis] of the target that gave it.
    python searches/k78_pass.py RUNDIR [DEST]    -> DEST (default outputs/k78_pipeline/full)/l3_k{7,8}/ii_eg.json
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); P = os.path.join(ROOT, "outputs", "k78_pipeline")
run = sys.argv[1]; DEST = sys.argv[2] if len(sys.argv) > 2 else os.path.join(P, "full")
TA = json.load(open(os.path.join(P, "targets_all.json"))); fps = TA["fingerprints"]
fid_of = {}
for t in TA["targets"]:
    fid_of[(t["k"], json.dumps(sorted(map(sorted, t["ccz"]))))] = t["fp"]
R, settings, files = [], None, sorted(glob.glob(os.path.join(run, "c*of*.json")))
M = int(files[0].split("of")[1].split(".")[0]) if files else 0
for f in files:
    d = json.load(open(f)); R += d["results"]; settings = d["settings"]
B = json.load(open(os.path.join(P, "staged_B_ccz_best.json")))
from k78_orbits import fingerprint
from isotropy_floor import floor as iso_floor
for k in (7, 8):
    prevb = {}
    for r in B[str(k)]:
        fid = next((v for v, x in fps.items() if x["k"] == k and x["fp"] == r["fp"]), None)
        if fid and r["verified"] and (fid not in prevb or (r["n"], r["N"]) < (prevb[fid]["n"], prevb[fid]["N"])):
            prevb[fid] = dict(n=r["n"], N=r["N"], stage=f"B:i_gs:{r['G']}x{r['H'][2:]}", gates=r["gates"], ccz=r["ccz"], tau=[r["floor"], None])
    best = {o: dict(v) for o, v in prevb.items()}; rows = []; cpu = 0.0
    for r in R:
        if r["k"] != k or not r["best"]:
            continue
        o = r.get("fp") or fid_of[(k, json.dumps(sorted(map(sorted, r["ccz"]))))]
        b = r["best"]; cpu += sum(r["cpu_s"].values())
        rows.append(dict(orbit=o, id=r.get("id"), n=b["n"], N=b["N"], stage="ii_eg:" + b["kind"], tau=[r["tcount_lo"], r["tcount_hi"]],
                         on_bound=r["on_bound"]))
        cur = best.get(o)
        if cur is None or (b["n"], b["N"]) < (cur["n"], cur["N"]):        # Book semantics: strictly fewer (n, N)
            best[o] = dict(n=b["n"], N=b["N"], stage="ii_eg:" + b["kind"], gates=b["gates"], ccz=r["ccz"], id=r.get("id"),
                           tau=[r["tcount_lo"], r["tcount_hi"]], on_bound=r["on_bound"])
    # T-count interval: per target (fingerprints are not proved injective, so bounds are never combined across
    # targets). For a class, the bound of the target that gave its best factory, or for a B factory the best interval
    # among this pass's targets with the same layer (else the isotropy floor alone)
    tau_t = {json.dumps(sorted(map(sorted, r["ccz"]))): (r["tcount_lo"], r["tcount_hi"], r.get("tcount_exact", r["tcount_lo"] == r["tcount_hi"]))
             for r in R if r["k"] == k}
    for o, b in best.items():
        lo, hi, ex = tau_t.get(json.dumps(sorted(map(sorted, b["ccz"]))), (b["tau"][0], None, False))
        b["tau"] = [lo, hi]; b["tau_exact"] = bool(hi is not None and lo == hi)
        b["isotropy_floor"] = iso_floor(b["ccz"], k)
        b["tau_lower_from"] = "tcount7.exists_below2 (exact, k=7)" if lo > b["isotropy_floor"] else "isotropy floor"
        b["class_on_bound"] = bool(b["tau_exact"] and b["n"] == hi + 1)
    new = sorted(o for o in best if o not in prevb)
    lowered = {o: [prevb[o]["n"], best[o]["n"]] for o in best if o in prevb and best[o]["n"] < prevb[o]["n"]}
    classes = sorted(v for v, x in fps.items() if x["k"] == k)
    rec = dict(L=3, k=k, pass_="ii_eg", settings=dict(settings or {}, chunks=M, synthesis="frontier_k7.bound + isd"),
               seconds=round(cpu, 1), seconds_are_cpu_sum=True, classes=len(classes), reached=len(best), new=new, lowered=lowered,
               inputs_saved=sum(a - b for a, b in lowered.values()), prev="sweep B stage (i), CCZ-layer outputs",
               best=best, rows=rows, specs={o: "fingerprint " + o for o in classes},
               targets=sum(1 for r in R if r["k"] == k), targets_with_factory=len(rows),
               on_bound_classes=sorted(o for o, b in best.items() if b["class_on_bound"]))
    os.makedirs(os.path.join(DEST, f"l3_k{k}"), exist_ok=True)
    json.dump(rec, open(os.path.join(DEST, f"l3_k{k}", "ii_eg.json"), "w"))
    print(f"k={k}: {rec['targets_with_factory']}/{rec['targets']} targets with a factory; classes reached {len(best)}/{len(classes)} "
          f"(+{len(new)} not in B), lowered {len(lowered)} by {rec['inputs_saved']}; on the bound: {len(rec['on_bound_classes'])} classes; "
          f"{rec['seconds']} CPU-s", flush=True)
    for o in classes:
        b = best.get(o)
        print("  ", o, b and (b["n"], b["N"], b["stage"], "tau", b["tau"], b["tau_lower_from"][:8], "bound" if b["class_on_bound"] else "", "B", prevb.get(o, {}).get("n")), flush=True)
