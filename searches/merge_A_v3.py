"""Merge the stage-(iii) cascade runs (chain_iii.py) into a complete sweep-A directory in the staged.py schema.

For each l=3 cell k=1..6: copy i_tg .. ii_eg unchanged (outputs/staged_pc for k<=4, outputs/staged_A for k=5,6) and
write one pass file per chain stage, cumulative in chain order, named iii_<stage>.json (e.g. iii_xs_None_2e6.json),
plus iii.json = the final best after stage (iii) (the role of iii_f.json; also written as iii_f.json so check_A2 reads
it). Each pass file: best per class after that stage, lowered/new against the previous pass, inputs_saved, CPU seconds
of that stage (sum over classes), and `first`: classes whose final best was first reached at this stage.
Factories credited to another class than searched for are included, as in staged_class.do_merge.
    python searches/merge_A_v3.py OUT        -> OUT/l3_k{k}/*.json, OUT/summary.json
"""
import os, sys, json, glob, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
PRE = ["i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg"]
CELLS = [(1, "staged_pc"), (2, "staged_pc"), (3, "staged_pc"), (4, "staged_pc"), (5, "staged_A"), (6, "staged_A")]


def save(p, o):
    json.dump(o, open(p + ".tmp", "w")); os.replace(p + ".tmp", p)


def main(out):
    summ = {}
    for k, src in CELLS:
        s = os.path.join(ROOT, "outputs", src, f"l3_k{k}"); d = os.path.join(out, f"l3_k{k}")
        for p in PRE:
            shutil.copyfile(os.path.join(s, p + ".json"), os.path.join(d, p + ".json"))
        eg = json.load(open(os.path.join(s, "ii_eg.json"))); specs = eg["specs"]
        files = sorted(glob.glob(os.path.join(d, "chain", "c*of*.json")))
        chain, res, per, done = None, [], {}, set()
        for f in files:
            x = json.load(open(f)); chain = chain or x["chain"]; assert x["chain"] == chain
            res += x["results"]; per.update(x["per_class"]); done |= set(x["classes"])
        missing = sorted(set(specs) - done)
        best = {o: dict(b) for o, b in eg["best"].items()}; prev_file = eg
        stage_first = {}
        for st in chain:
            prevn = {o: b["n"] for o, b in best.items()}
            for r in res:
                if r["stage"] != "iii:" + st:
                    continue
                o = str(r["orbit"]); cur = best.get(o)
                if o in specs and (cur is None or (r["n"], r["N"]) < (cur["n"], cur["N"])):
                    best[o] = dict(n=r["n"], N=r["N"], stage=r["stage"], gates=r["gates"])
            low = {o: [v, best[o]["n"]] for o, v in prevn.items() if best[o]["n"] < v}
            new = sorted(o for o in best if o not in prevn)
            cpu = sum(sum(z.get("seconds", 0.0) for z in p["stages"] if z["stage"] == st) for p in per.values())
            name = "iii_" + st.replace(":", "_")
            save(os.path.join(d, name + ".json"), dict(L=3, k=k, pass_=name, settings=dict(stage=st, chain=chain),
                 seconds=round(cpu, 1), seconds_are_cpu_sum=True, classes=len(specs), reached=len(best), new=new,
                 lowered=low, inputs_saved=sum(a - b for a, b in low.values()), best=best, rows=[], steps=None,
                 specs=specs, missing_class_files=missing))
            stage_first[st] = dict(lowered=len(low), inputs_saved=sum(a - b for a, b in low.values()), cpu_s=round(cpu, 1))
        fin = dict(L=3, k=k, pass_="iii", settings=dict(chain=chain), classes=len(specs), reached=len(best),
                   lowered={o: [eg["best"][o]["n"], best[o]["n"]] for o in eg["best"] if best[o]["n"] < eg["best"][o]["n"]},
                   best=best, rows=[], steps=None, specs=specs, missing_class_files=missing,
                   seconds=round(sum(v["cpu_s"] for v in stage_first.values()), 1), seconds_are_cpu_sum=True)
        fin["inputs_saved"] = sum(a - b for a, b in fin["lowered"].values())
        save(os.path.join(d, "iii.json"), fin); save(os.path.join(d, "iii_f.json"), dict(fin, pass_="iii_f"))
        firsts = {}
        for o, b in best.items():
            firsts[b["stage"].split(":")[0] if not b["stage"].startswith("iii:") else b["stage"]] = \
                firsts.get(b["stage"].split(":")[0] if not b["stage"].startswith("iii:") else b["stage"], 0) + 1
        summ[f"l3_k{k}"] = dict(classes=len(specs), missing=len(missing), chain=chain, per_stage=stage_first,
                                final_source=firsts, lowered_below_ii_eg=len(fin["lowered"]), inputs_saved=fin["inputs_saved"])
        print(f"k={k}: {summ[f'l3_k{k}']}", flush=True)
    save(os.path.join(out, "summary.json"), summ)


if __name__ == "__main__":
    main(sys.argv[1])
