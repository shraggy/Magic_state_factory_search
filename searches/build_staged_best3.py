"""The l=3 cells of the cascade with stage (iii) taken per class as the best of the three deterministic stage-(iii) runs,
in the pass layout of outputs/staged_final, so that staged_best.py and stage_contributions.py read it unchanged.

Sources (searches/chain_iii2.py, same stages (i)-(ii) of sweep A and the same gluing seed ii_eg; each layer starts from
the running best n of the class and passes only n on):
  v4  outputs/solver_v2/staged_A_v4   xs:None, xs:6, x:5, x:6, x:None            (z3 syn after it lowered none, not used)
  v5  outputs/solver_v2/staged_A_v5   xsb:8:None, xsb:8:6, x:5, x:6, x:None
  v6  outputs/solver_v2/staged_A_v6   xs:None, xs:6, x:5, x:6, x:None, xsb:8:None, x2:5, x2:6, x2:None
xs, xsb = restricted SAT (fixed output parts); x, x2 = free SAT. The passes i_tg..ii_eg are identical in the three runs
(checked here, n and gates) and are copied from v4, as build_staged_final.py does.

Per class: the final factory is the one with the smallest n over the three runs' iii.json; ties to the smallest N, then
v6, then v4, then v5. Credit (rule of this file): the final n is credited to the form of the layer that first reached it
in the winning run's own order (restricted or free); if no layer lowered the seed, to the pass of stages (i)-(ii) that
staged_best.py finds. Since a layer only records an n strictly below its seed, this credit is expressed in the pass
layout by
  iii_r  per class the smallest n among the seed and the restricted layers' own finds of the winning run, with the
         gates of the pass file after that layer (the seed's gates if no restricted layer lowered it)
  iii_f  per class the final factory
so that staged_best.py's rule 'first pass whose cumulative best reaches the final n' gives iii_r exactly for the
restricted-credited classes and iii_f for the free-credited ones (asserted below). Example: in v6 a class lowered by x:5
and then further by xsb:8 is credited to restricted; one lowered by xs:None and then further by x:5 to free.
Two limits of this credit, both recorded so that nothing has to be read off the pass layout:
  (a) The 'lowered' and 'inputs_saved' of iii_r / iii_f are differences between consecutive cumulative passes, not what
      each form did: in v6 the forms interleave (xs, x, xsb, x2), so a class that free SAT lowered and restricted SAT
      lowered further sits wholly in iii_r (11 of 766 classes, 22 inputs). What each form saved, following the winning
      run's own layer order, is in form_lowered / form_inputs_saved of iii_r and iii_f and in attribution.json
      (saved_restricted, saved_free); per run and per layer in forms.json. Quote those, not 'lowered' / 'inputs_saved'.
  (b) The credit depends on the tie-break between runs: most classes tie on n and N across the runs, and the run that
      wins the tie (v6 first, the version the pipeline description keeps) decides the form. attribution.json gives, per
      class, the credited form in every run that ties with the winner (forms_tied), and summary.json the totals for
      'credited by the winner', 'restricted in every tied run' and 'restricted in some tied run'.
Seconds of iii_r / iii_f are those of v4 (the minimum version, which the time table quotes), as in staged_final; the
CPU seconds of every layer of every run are in forms.json and in seconds_by_run. Written here:
  outputs/staged_final_b3/l3_k{K}/{i_tg,i_gs,i_gh,ii_dj,ii_eg}.json  copied from v4
  outputs/staged_final_b3/l3_k{K}/iii_r.json, iii_f.json  as above; every best entry also carries run, layer, form
  outputs/staged_final_b3/l3_k{K}/attribution.json  per class: run, layer, form, n, N, seed n, n of v4, v5, v6
  outputs/staged_final_b3/l3_k{K}/forms.json  per run, per layer: lowered, inputs saved, CPU s
  outputs/staged_final_b3/summary.json  over k=1..6: credit counts (winner / every tied run / some tied run) and the
                                        per-form savings in the winning run's order
    python searches/build_staged_best3.py   -> outputs/staged_final_b3/l3_k{1..6}/
"""
import os, json, glob, shutil
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
RUNS = {v: os.path.join(O, "solver_v2", f"staged_A_{v}") for v in ("v4", "v5", "v6")}
PREF = {"v6": 0, "v4": 1, "v5": 2}                            # tie-break after n and N
DST = os.path.join(O, "staged_final_b3")
COPY = ["i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg"]
FORM = lambda st: "restricted" if st.startswith("xs") else ("free" if st.startswith("x") else None)
PASSFILE = lambda st: "iii_" + st.replace(":", "_")


def load_run(v, k):
    d = os.path.join(RUNS[v], f"l3_k{k}"); per, chain = {}, None
    for f in glob.glob(os.path.join(d, "chain", "*.json")):
        c = json.load(open(f)); per.update(c["per_class"]); assert chain in (None, c["chain"]); chain = c["chain"]
    layers = [st for st in chain if FORM(st)]
    P = {st: json.load(open(os.path.join(d, PASSFILE(st) + ".json"))) for st in layers}
    final = json.load(open(os.path.join(d, "iii.json")))
    return dict(dir=d, per=per, layers=layers, P=P, final=final)


def credit(pc, seed_n):
    """(form that first reached the final n or None, {form: inputs saved}, {form: 1 if it lowered the class}) in run order."""
    first, cur, saved, low = None, seed_n, {"restricted": 0, "free": 0}, {"restricted": 0, "free": 0}
    for st in pc["stages"]:
        if not FORM(st["stage"]) or "skipped" in st:
            continue
        if st["n"] < cur:
            fm = FORM(st["stage"]); saved[fm] += cur - st["n"]; low[fm] = 1; cur = st["n"]
            if cur == pc["final"] and first is None: first = fm
    return first, saved, low


def main():
    tot = dict(classes=0, credit_winner={}, restricted_every_tied=0, restricted_some_tied=0, free_every_tied=0,
               free_some_tied=0, tie_dependent=0, form_lowered={"restricted": 0, "free": 0},
               form_inputs_saved={"restricted": 0, "free": 0}, pass_lowered={}, pass_inputs_saved={})
    for k in range(1, 7):
        R = {v: load_run(v, k) for v in RUNS}
        s4, d = R["v4"]["dir"], os.path.join(DST, f"l3_k{k}")
        os.makedirs(d, exist_ok=True)
        for p in COPY:
            a = json.load(open(os.path.join(s4, f"{p}.json")))["best"]
            for v in ("v5", "v6"):
                b = json.load(open(os.path.join(R[v]["dir"], f"{p}.json")))["best"]
                assert set(a) == set(b) and all(a[o]["n"] == b[o]["n"] and a[o]["gates"] == b[o]["gates"] for o in a), (k, p, v)
            shutil.copy(os.path.join(s4, f"{p}.json"), os.path.join(d, f"{p}.json"))
        seed = json.load(open(os.path.join(s4, "ii_eg.json")))
        best_r, best_f, attr, fsave, flow = {}, {}, {}, {}, {}
        for o in sorted(seed["best"]):
            sb = seed["best"][o]
            win = min(RUNS, key=lambda v: (R[v]["final"]["best"][o]["n"], R[v]["final"]["best"][o]["N"], PREF[v]))
            W = R[win]; fb = W["final"]["best"][o]; n_fin = fb["n"]; pc = W["per"][o]
            assert pc["seed"] == sb["n"] and pc["final"] == n_fin, (k, o, win)
            first, rmin, cur = None, (sb["n"], None), sb["n"]           # rmin: (n, layer) over seed + restricted finds
            for st in pc["stages"]:
                if not FORM(st["stage"]) or "skipped" in st:
                    continue
                n = st["n"]
                if n < cur:                                            # this layer's own find
                    if FORM(st["stage"]) == "restricted" and n < rmin[0]: rmin = (n, st["stage"])
                    if n == n_fin and first is None: first = st["stage"]
                    cur = n
            assert cur == n_fin, (k, o, win)
            form = FORM(first) if first else None
            rb = dict(sb) if rmin[1] is None else dict(W["P"][rmin[1]]["best"][o])
            assert rb["n"] == rmin[0] and rb["n"] <= sb["n"] and n_fin <= rb["n"], (k, o)
            assert (rb["n"] == n_fin) == (form in (None, "restricted")), (k, o, form)   # first-pass credit = form credit
            tag = dict(run=win, layer=first, form=form)
            f_w, saved, low = credit(pc, sb["n"]); assert f_w == form, (k, o)
            tied = {v: credit(R[v]["per"][o], sb["n"])[0] for v in RUNS
                    if (R[v]["final"]["best"][o]["n"], R[v]["final"]["best"][o]["N"]) == (n_fin, fb["N"])}
            ft = set(tied.values()); key = form or "no SAT lowering"
            tot["classes"] += 1; tot["credit_winner"][key] = tot["credit_winner"].get(key, 0) + 1
            tot["restricted_every_tied"] += ft == {"restricted"}; tot["restricted_some_tied"] += "restricted" in ft
            tot["free_every_tied"] += ft == {"free"}; tot["free_some_tied"] += "free" in ft; tot["tie_dependent"] += len(ft) > 1
            for fm in saved:
                tot["form_inputs_saved"][fm] += saved[fm]; tot["form_lowered"][fm] += low[fm]
            fsave[o], flow[o] = saved, low
            best_r[o] = dict(rb, run=win if rmin[1] else None, layer=rmin[1], form="restricted" if rmin[1] else None)
            best_f[o] = dict(fb, **tag)
            attr[o] = dict(tag, n=n_fin, N=fb["N"], seed=sb["n"], spec=seed["specs"][o],
                           saved_restricted=saved["restricted"], saved_free=saved["free"], forms_tied=tied,
                           **{f"n_{v}": R[v]["final"]["best"][o]["n"] for v in RUNS})
        forms, secs = {}, {}
        for v, W in R.items():
            forms[v] = {st: dict(lowered=len(W["P"][st]["lowered"]), inputs_saved=W["P"][st]["inputs_saved"],
                                 cpu_s=W["P"][st]["seconds"]) for st in W["layers"]}
            secs[v] = {fm: round(sum(W["P"][st]["seconds"] for st in W["layers"] if FORM(st) == fm), 1)
                       for fm in ("restricted", "free")}
        prev = seed["best"]
        for name, B, fm in (("iii_r", best_r, "restricted"), ("iii_f", best_f, "free")):
            lowered = {o: [prev[o]["n"], B[o]["n"]] for o in B if B[o]["n"] < prev[o]["n"]}
            rec = dict(L=3, k=k, pass_=name, classes=len(B), reached=len(B), new=[], lowered=lowered,
                       inputs_saved=sum(a - b for a, b in lowered.values()), best=B, rows=[], steps=None,
                       form_lowered=sum(flow[o][fm] for o in flow), form_inputs_saved=sum(fsave[o][fm] for o in fsave),
                       note="lowered / inputs_saved: difference to the previous cumulative pass; form_lowered / "
                            "form_inputs_saved: what this form did in the winning run's layer order (quote these)",
                       specs=seed["specs"], missing_class_files=[], seconds=secs["v4"][fm], seconds_are_cpu_sum=True,
                       seconds_by_run={v: secs[v][fm] for v in RUNS},
                       settings=dict(rule="per-class best of v4, v5, v6 (n, then N, then v6, v4, v5); credit to the form "
                                          "of the layer that first reached the final n in the winning run",
                                     form=fm, runs={v: os.path.relpath(R[v]["dir"], ROOT) for v in RUNS},
                                     layers={v: [st for st in R[v]["layers"] if FORM(st) == fm] for v in RUNS},
                                     solver="CryptoMiniSat 5.17.0 via searches/xor_sat2.py", driver="searches/chain_iii2.py",
                                     builder="searches/build_staged_best3.py", seconds="v4 (minimum version)"))
            json.dump(rec, open(os.path.join(d, f"{name}.json"), "w"))
            tot["pass_lowered"][name] = tot["pass_lowered"].get(name, 0) + len(lowered)
            tot["pass_inputs_saved"][name] = tot["pass_inputs_saved"].get(name, 0) + rec["inputs_saved"]
            prev = B
        json.dump(attr, open(os.path.join(d, "attribution.json"), "w"), indent=1)
        json.dump(forms, open(os.path.join(d, "forms.json"), "w"), indent=1)
        cnt = {}
        for a in attr.values():
            key = f"{a['run']}:{a['layer']}" if a["layer"] else "no SAT lowering"; cnt[key] = cnt.get(key, 0) + 1
        print(f"l=3 k={k}: {len(attr)} classes; credit {dict(sorted(cnt.items()))}", flush=True)
    tot["tie_break"] = "n, then N, then v6, v4, v5"
    json.dump(tot, open(os.path.join(DST, "summary.json"), "w"), indent=1); print(json.dumps(tot))


if __name__ == "__main__":
    main()
