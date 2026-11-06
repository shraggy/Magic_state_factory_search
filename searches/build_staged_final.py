"""The l=3 cells of the cascade with the final stage (iii), in the pass layout staged_best.py and stage_contributions.py
read. Source: outputs/solver_v2/staged_A_v4 (searches/chain_iii2.py: history-independent synthesis parts, every XOR
query in a fresh process; x86_64), whose passes i_tg..ii_eg are copied from sweep A unchanged and whose stage (iii) ran
in this order from the gluing best: restricted XOR SAT with supports unlimited then <= 6 (xs:None, xs:6), free XOR SAT
with supports <= 5, <= 6, unlimited (x:5, x:6, x:None), then z3 synthesis SAT (syn), which lowered no class and is not
part of the final pipeline. Each form's pass file holds the cumulative best after it, so dropping the last form
changes nothing. Written here:
  l3_k{K}/{i_tg,i_gs,i_gh,ii_dj,ii_eg}.json  copied
  l3_k{K}/iii_r.json   restricted SAT: the cumulative best after xs:6 (lowered/new/seconds summed over xs:None, xs:6)
  l3_k{K}/iii_f.json   free SAT: the cumulative best after x:None (summed over x:5, x:6, x:None)
  l3_k{K}/forms.json   per form: lowered, inputs saved, CPU s (for Table SII and the README)
    python searches/build_staged_final.py   -> outputs/staged_final/l3_k{1..6}/
"""
import os, json, shutil
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
SRC = os.path.join(O, "solver_v2", "staged_A_v4"); DST = os.path.join(O, "staged_final")
COPY = ["i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg"]
STAGES = {"iii_r": ["iii_xs_None_2e6", "iii_xs_6_2e6"], "iii_f": ["iii_x_5_2e6", "iii_x_6_2e6", "iii_x_None_2e6"]}


def main():
    for k in range(1, 7):
        s, d = os.path.join(SRC, f"l3_k{k}"), os.path.join(DST, f"l3_k{k}")
        os.makedirs(d, exist_ok=True)
        for p in COPY:
            shutil.copy(os.path.join(s, f"{p}.json"), os.path.join(d, f"{p}.json"))
        forms, prev = {}, json.load(open(os.path.join(s, "ii_eg.json")))["best"]
        for st, parts in STAGES.items():
            P = [json.load(open(os.path.join(s, f"{p}.json"))) for p in parts]
            last = P[-1]
            lowered = {o: [prev[o]["n"], b["n"]] for o, b in last["best"].items() if o in prev and b["n"] < prev[o]["n"]}
            new = [o for o in last["best"] if o not in prev]
            rec = dict(last, pass_=st, settings=dict(forms=parts, solver="CryptoMiniSat 5.17.0 via searches/xor_sat2.py",
                                                     driver="searches/chain_iii2.py", source=os.path.relpath(s, ROOT)),
                       seconds=round(sum(x["seconds"] for x in P), 1), seconds_are_cpu_sum=True,
                       lowered=lowered, new=new, inputs_saved=sum(a - b for a, b in lowered.values()))
            json.dump(rec, open(os.path.join(d, f"{st}.json"), "w"))
            for x in P:
                forms[x["pass_"]] = dict(lowered=len(x["lowered"]), inputs_saved=x["inputs_saved"], cpu_s=x["seconds"])
            prev = last["best"]
        syn = json.load(open(os.path.join(s, "iii_syn.json")))
        assert all(syn["best"][o]["n"] == prev[o]["n"] for o in prev), k          # z3 synthesis SAT adds nothing
        forms["iii_syn (not in the pipeline)"] = dict(lowered=len(syn["lowered"]), cpu_s=syn["seconds"])
        json.dump(forms, open(os.path.join(d, "forms.json"), "w"), indent=1)
        print(f"l=3 k={k}: {len(prev)} classes; " + ", ".join(f"{f} {v['lowered']}" for f, v in forms.items()), flush=True)


if __name__ == "__main__":
    main()
