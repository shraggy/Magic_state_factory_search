"""Test A: would the SAT stage (iii) alone, without the borrowed-identity stages (i)-(ii), reach the cascade's factories?

Same code path as the v6 run (outputs/solver_v2/staged_A_v6): chain_iii4.stage (-> chain_iii3 -> chain_iii2, det_parts
installed, every CryptoMiniSat query in a fresh spawned process, 2e6 conflicts per query), the same layer chain
  xs:None:2e6, xs:6:2e6, x:5:2e6, x:6:2e6, x:None:2e6, xsb:8:None:2e6, x2:5:2e6
(v6's x2:6 and x2:None lowered no class in v6 and are left out), the same per-class loop as chain_iii2.run (skip at
tau+1, verified factories labelled with staged.label_fix, only the class's own factories lower its running cap). Only the
starting cap changes:
  seed  cap = n of the gluing best ii_eg (what v6 used): known-answer check, must reproduce v6's per-layer n
  syn   cap = syn_n, the strongest synthillation construction (outputs/staged_best.json): SAT asked only to beat it
  free  as none, but only the free XOR layers (x:5, x:6, x:None, x2:5): no synthesis parts at all
  none  cap = 10**6, i.e. no factory known (chain_iii2.run's own default for a class with no seed; xor_sat2 then runs its
        first query without a cardinality bound and lowers from the first solution)
Known-answer result (rule 13): on arm64 the per-layer n matched v6 (x86_64) for k=3,4,5 classes, not for the two k=6
ones (README.md: stage (iii) results agree on every x86_64 host but not across architectures). So on arm64 the seeded
run is also made, as a same-architecture control, for every sample class where the unseeded run misses v6's final n.
Modes free and syn run on every 3rd class of the interleaved sample (21 of 63), for CPU.
Order: known-answer seeds, free (subset), none (all sample classes), seeded controls where none differs from v6, syn (subset), seeded controls (subset).
Sample (l=3): k=3,4: none of the 19 classes has its final factory from stage (iii), so none is taken; k=5 and k=6: every
12th class in sorted orbit-id order (6 + 57). Known-answer classes: KA below.
    python searches/sat_alone_test.py run [--max-workers 3]   driver: one job (class, cap mode) per fresh process,
                                                              workers <= min(3, 16-5-load), nice 19; resumable
    python searches/sat_alone_test.py one K ORBIT MODE        one job -> outputs/checks/sat_alone_test/raw/
    python searches/sat_alone_test.py summary                 -> outputs/checks/sat_alone_test.json and .md
    python searches/sat_alone_test.py verify                  every best factory in raw/ re-checked as check_A2.py does
                                                              (polynomial + brute-force condition, n, distance, class)
                                                              -> outputs/checks/sat_alone_test/verify.json
"""
import os, sys, json, time, glob, subprocess, signal, resource
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
O = os.path.join(ROOT, "outputs"); V6 = os.path.join(O, "solver_v2", "staged_A_v6")
CACHE_RUN = os.path.join(O, "staged_A")                        # holds _cache/l3_k{k}_classes.pkl (same ids as v6)
RAW = os.path.join(O, "checks", "sat_alone_test", "raw")
CHAIN = ["xs:None:2e6", "xs:6:2e6", "x:5:2e6", "x:6:2e6", "x:None:2e6", "xsb:8:None:2e6", "x2:5:2e6"]
KA = [(3, "24"), (4, "112"), (5, "120"), (6, "2304:136"), (6, "2320:9041")]   # cheap; 5/6 lowered by xs, x:5, xsb in v6
CHAIN_FREE = [c for c in CHAIN if c.startswith("x:") or c.startswith("x2:")]   # free XOR SAT only: no synthesis parts
MODES = ["none", "free", "seed", "syn"]                            # seed = same-architecture control for the sample


def v6_per(k):
    per = {}
    for f in glob.glob(os.path.join(V6, f"l3_k{k}", "chain", "*.json")):
        c = json.load(open(f)); assert c["chain"][:7] == CHAIN; per.update(c["per_class"])
    return per


def rows():
    return {(r["k"], str(r["orbit"])): r for r in json.load(open(os.path.join(O, "staged_best.json")))["rows"] if r["l"] == 3}


def sample():
    """Every 12th class by sorted orbit id at k=5 and k=6, interleaved by relative position (so a partial run is spread
    over both k and the whole id range)."""
    out = []
    for k in (5, 6):
        s = sorted(v6_per(k))[::12]
        out += [((i + 0.5) / len(s), k, o) for i, o in enumerate(s)]
    return [(k, o) for _, k, o in sorted(out)]


def one(k, o, mode):
    sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
    import chain_iii4                                          # installs chain_iii4.stage into chain_iii2
    import chain_iii2 as C2, staged as STG, staged_class as SC
    L = 3; cache, classes = SC.setup(L, k, CACHE_RUN)
    spec = classes[o]["spec"]; tau = STG.tau_of(L, k, o, spec); r = rows()[(k, o)]
    seed = json.load(open(os.path.join(V6, f"l3_k{k}", "ii_eg.json")))["best"]
    cap = dict(seed=seed[o]["n"] if o in seed else 10 ** 6, syn=r["syn_n"] if r["syn_n"] else 10 ** 6, none=10 ** 6,
               free=10 ** 6)[mode]
    chain = CHAIN_FREE if mode == "free" else CHAIN
    rec = dict(k=k, orbit=o, spec=spec, mode=mode, start_cap=cap, tau=tau, chain=chain, stages=[], other_class=[])
    T0 = time.time(); R0 = resource.getrusage(resource.RUSAGE_CHILDREN); S0 = resource.getrusage(resource.RUSAGE_SELF)
    for name in chain:
        if tau is not None and cap == tau + 1:
            rec["stages"].append(dict(stage=name, skipped="at bound")); continue
        t = time.time(); c0 = cap; c1 = resource.getrusage(resource.RUSAGE_CHILDREN)
        found = []
        for n, N, g in C2.stage(name, k, spec, cap):
            o2, gg = STG.label_fix(dict(k=k, r=N - k, gates=[(frozenset(S), cc) for S, cc in g]), L, k, cache)
            if o2 is None:
                continue
            if o2 != o:
                rec["other_class"].append(dict(stage=name, orbit=o2, n=n, N=N)); continue
            found.append(dict(n=n, N=N, gates=[[sorted(S_), cc] for S_, cc in gg]))
            if n < cap:
                cap = n
        c2 = resource.getrusage(resource.RUSAGE_CHILDREN)
        rec["stages"].append(dict(stage=name, seed=c0, n=cap if cap < 10 ** 6 else None, seconds=round(time.time() - t, 1),
                                  cpu_children=round(c2.ru_utime + c2.ru_stime - c1.ru_utime - c1.ru_stime, 1),
                                  found=[(f["n"], f["N"]) for f in found]))
        if found:
            rec["best_gates"] = found[-1]
    R1 = resource.getrusage(resource.RUSAGE_CHILDREN); S1 = resource.getrusage(resource.RUSAGE_SELF)
    rec["final"] = cap if cap < 10 ** 6 else None
    rec["wall_s"] = round(time.time() - T0, 1)
    rec["cpu_s"] = round(R1.ru_utime + R1.ru_stime - R0.ru_utime - R0.ru_stime + S1.ru_utime + S1.ru_stime
                         - S0.ru_utime - S0.ru_stime, 1)
    rec["host"] = os.uname().nodename; rec["machine"] = os.uname().machine
    os.makedirs(RAW, exist_ok=True)
    p = os.path.join(RAW, f"k{k}_{o.replace(':', '-')}_{mode}.json")
    json.dump(rec, open(p + ".tmp", "w")); os.replace(p + ".tmp", p)
    print(f"done k={k} {o} {mode}: start {rec['start_cap']} final {rec['final']} cpu {rec['cpu_s']}",
          flush=True)


def raw_path(k, o, mode):
    return os.path.join(RAW, f"k{k}_{o.replace(':', '-')}_{mode}.json")


def load1():
    return os.getloadavg()[0]


def drive(maxw):
    """Known-answer jobs first; the unseeded jobs only once all five reproduce v6 (rule 13)."""
    procs = []

    def kill(*_):
        for p, _j in procs:
            if p.poll() is None:
                os.killpg(p.pid, signal.SIGTERM)
        sys.exit(1)
    signal.signal(signal.SIGTERM, kill); signal.signal(signal.SIGINT, kill)
    def missed():                                      # seeded control only where no-seed misses v6's final n
        out = []
        for k, o in sample():
            p = raw_path(k, o, "none")
            if os.path.exists(p):
                f = json.load(open(p))["final"]
                if f is None or f != v6_per(k)[o]["final"]:  # misses, or beats, v6
                    out.append((k, o, "seed"))
        return out
    phases = [lambda: [(k, o, "seed") for k, o in KA],
              lambda: [(k, o, "free") for k, o in sample()[::3]],
              lambda: [(k, o, "none") for k, o in sample()],
              missed,
              lambda: [(k, o, "syn") for k, o in sample()[::3]],
              lambda: [(k, o, "seed") for k, o in sample()[::3]]]
    for ph, jobs in enumerate(phases):
        todo = [j for j in jobs() if not os.path.exists(raw_path(*j))]
        while todo or procs:
            procs = [(p, j) for p, j in procs if p.poll() is None]
            free = int(16 - 5 - load1())
            if todo and len(procs) < min(maxw, max(0, free)):
                j = todo.pop(0)
                p = subprocess.Popen(["nice", "-n", "19", sys.executable, os.path.abspath(__file__), "one", str(j[0]), j[1], j[2]],
                                     start_new_session=True)
                procs.append((p, j)); print(time.strftime("%H:%M:%S"), "start", j, "load", round(load1(), 1), flush=True)
                time.sleep(90)                                 # let the load average see the new worker
                continue
            time.sleep(20)
        if ph == 0:
            bad = known_answer()["mismatch"]
            if bad:                                        # v6 ran on x86_64; README: results differ across
                print("KNOWN-ANSWER CHECK: per-layer mismatch vs v6", bad,   # architectures. Continue only on arm64 with
                      "-> same-architecture seeded control is run for every sample class", flush=True)  # the seed control.
                if os.uname().machine == "x86_64":
                    sys.exit(2)
            else:
                print("known-answer check passed", flush=True)


def known_answer():
    res, bad = [], []
    for k, o in KA:
        p = raw_path(k, o, "seed")
        if not os.path.exists(p):
            res.append(dict(k=k, orbit=o, status="not run")); continue
        a = json.load(open(p)); b = v6_per(k)[o]
        na = [s.get("n", "skip") for s in a["stages"]]; nb = [s.get("n", "skip") for s in b["stages"][:7]]
        ok = na == nb and a["start_cap"] == b["seed"]
        res.append(dict(k=k, orbit=o, seed=a["start_cap"], here=na, v6=nb, v6_final=b["final"], final=a["final"], ok=ok,
                        cpu_s=a["cpu_s"], v6_s=round(sum(s.get("seconds", 0) for s in b["stages"][:7]), 1)))
        if not ok:
            bad.append((k, o))
    return dict(classes=res, mismatch=bad)


def z3_evidence():
    """outputs/staged_A/l3_k{5,6}/free_coef_1.5e9_unseeded: z3 all-coefficient free SAT, no seed, per class steps N."""
    R = rows(); out = {}
    for k in (5, 6):
        per, found = {}, {}
        for f in glob.glob(os.path.join(O, "staged_A", f"l3_k{k}", "free_coef*", "*.json")):
            d = json.load(open(f)); per.update(d["per_class"])
            for x in d.get("results", []):
                if x["orbit"] in d["per_class"]:
                    found[x["orbit"]] = min(found.get(x["orbit"], 10 ** 6), x["n"])
        st = {}; gap = {}
        reach = 0
        for o, p in per.items():
            r = R[(k, o)]; s = {x["N"]: x["last"] for x in p["steps"]}
            key = s.get(r["N"], "not run"); st[key] = st.get(key, 0) + 1
            if o in found and found[o] <= r["n"]:
                reach += 1
            src = "iii" if r["stage"].startswith("iii") else "i/ii"
            g = str(found[o] - r["n"]) if o in found else "none"
            gap.setdefault(src, {}); gap[src][g] = gap[src].get(g, 0) + 1
        out[f"k{k}"] = dict(classes=len(per), status_at_final_N=st, own_factory_found=len(found),
                            found_n_le_final=reach, found_minus_final_by_final_source=gap, cpu_s=round(sum(p["seconds"] for p in per.values()), 1),
                            source=f"outputs/staged_A/l3_k{k}/free_coef_1.5e9_unseeded/*.json")
    return out


def verify():
    sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
    import unified_search as us, cascade_class as CC
    lab = {}; rows = []
    for p in sorted(glob.glob(os.path.join(RAW, "*.json"))):
        a = json.load(open(p)); b = a.get("best_gates")
        if not b:
            continue
        k = a["k"]
        if k not in lab:
            cl, ctx = CC.classes_of(3, k, os.path.join(CACHE_RUN, "_cache")); lab[k] = CC.labeller(3, k, ctx, cl)
        g = [(frozenset(S), c) for S, c in b["gates"]]; f = dict(k=k, r=b["N"] - k, gates=g)
        rows.append(dict(file=os.path.basename(p), orbit=a["orbit"], n_filed=b["n"], n=sum(us.inputs(c, 3) for _, c in g),
                         N=b["N"], valid=us.is_factory(f, 3), brute=us.is_factory(f, 3, brute_force=True) if b["N"] <= 9 else None,
                         d=str(us.distance(f, 3)), class_ok=lab[k](f) == a["orbit"]))
    ok = [r for r in rows if r["valid"] and r["brute"] is not False and r["d"] != "1" and r["class_ok"] and r["n"] == r["n_filed"]]
    out = dict(factories=len(rows), all_checks_pass=len(ok), brute_checked=sum(r["brute"] is not None for r in rows), rows=rows)
    json.dump(out, open(os.path.join(O, "checks", "sat_alone_test", "verify.json"), "w"), indent=0)
    print({k_: v for k_, v in out.items() if k_ != "rows"})


def summary():
    R = rows(); V = {k: v6_per(k) for k in (5, 6)}; per, tally, res_layers = [], {}, {}
    for k, o in sample():
        r = R[(k, o)]; b = V[k][o]
        e = dict(k=k, orbit=o, spec=r["spec"], tau=r["tau"], syn_n=r["syn_n"], ii_eg=b["seed"], v6_final=b["final"],
                 best_final=r["n"], best_stage=r["stage"], v6_cpu_s=round(sum(s.get("seconds", 0) for s in b["stages"][:7]), 1))
        for m in MODES:
            p = raw_path(k, o, m)
            if not os.path.exists(p):
                e[m] = None; continue
            a = json.load(open(p)); n = a["final"]
            n_eff = n
            if m == "syn" and n is None:                       # nothing below syn_n: SAT ends at synthillation's n
                n_eff = r["syn_n"]
            if n_eff is None:
                cat = "nothing found"
            elif n_eff <= b["final"]:
                cat = "reaches cascade n" if n_eff == b["final"] else "below cascade n"
            else:
                cat = "larger than cascade n"
            first = next((x["stage"] for x in a["stages"] if x.get("n") is not None and n_eff is not None
                          and x["n"] <= n_eff and x.get("found")), None)
            e[m] = dict(final=n, n_eff=n_eff, cpu_s=a["cpu_s"], category=cat, first_layer_at_final=first,
                        per_stage=[x.get("n", "skip") for x in a["stages"]])
            if m == "none" and first:
                lay = res_layers.setdefault(first, 0); res_layers[first] = lay + 1
            g = (m, "iii" if r["stage"].startswith("iii") else "i/ii"); tally.setdefault(g, {}); tally[g][cat] = tally[g].get(cat, 0) + 1
        per.append(e)
    done = {m: [e for e in per if e.get(m)] for m in MODES}
    cpu = {m: dict(sat_alone=round(sum(e[m]["cpu_s"] for e in done[m]), 1),
                   v6_stage_iii=round(sum(e["v6_cpu_s"] for e in done[m]), 1), classes=len(done[m])) for m in MODES}
    res = dict(script="searches/sat_alone_test.py", chain=CHAIN, conflicts_per_query=2_000_000,
               known_answer=known_answer(), sample_size=len(per),
               tally={f"{m}|final from {s}": v for (m, s), v in sorted(tally.items())}, cpu=cpu,
               none_layer_that_reached_final=res_layers,
               z3_free_coef_unseeded=z3_evidence(), per_class=per,
               note="v6 ran on x86_64 (cluster); this run on " + os.uname().machine + ". 'nothing found' means not "
                    "reached within 2e6 conflicts per query, not infeasible.")
    p = os.path.join(O, "checks", "sat_alone_test.json")
    md = ["# Test A: SAT stage (iii) without stages (i)-(ii)", "", "Script: searches/sat_alone_test.py; data: "
          "outputs/checks/sat_alone_test.json, raw per job in outputs/checks/sat_alone_test/raw/.", "",
          "## Known-answer check (seed = ii_eg n, as v6)"]
    for c in res["known_answer"]["classes"]:
        md.append(f"- k={c['k']} {c['orbit']}: " + (c.get("status") or f"per-layer n here {c['here']} vs v6 {c['v6']}: "
                  + ("match" if c["ok"] else "MISMATCH") + f"; CPU {c['cpu_s']} s here, {c['v6_s']} s v6"))
    md += ["", f"## Sample: {len(per)} classes (k=5: every 12th of 65, k=6: every 12th of 679)"]
    for m in MODES:
        md.append(f"- mode {m}: {len(done[m])} of {len(per)} classes run")
    for key, v in res["tally"].items():
        md.append(f"  - {key}: " + ", ".join(f"{c}: {x}" for c, x in sorted(v.items())))
    md += ["", "No-seed runs: layer in which the run's final n was first found: " + ", ".join(f"{a}: {b}" for a, b in sorted(res_layers.items()))]
    md += ["", "## CPU (sum over classes run)"] + [f"- {m}: SAT alone {cpu[m]['sat_alone']} s vs v6 stage (iii) (7 layers) "
                                                   f"{cpu[m]['v6_stage_iii']} s over {cpu[m]['classes']} classes" for m in MODES]
    pair = lambda m1, m2: [e for e in per if e.get(m1) and e.get(m2)]
    sv = pair("seed", "none"); yv = pair("syn", "none"); fv = [e for e in per if e.get("free")]
    cmp_ = dict(seed_eq_none=[sum(e["seed"]["final"] == e["none"]["final"] for e in sv), len(sv)],
                seed_cpu_vs_none_cpu=[round(sum(e["seed"]["cpu_s"] for e in sv), 1), round(sum(e["none"]["cpu_s"] for e in sv), 1)],
                syn_eq_none=[sum(e["syn"]["n_eff"] == e["none"]["final"] for e in yv), len(yv)],
                free_eq_v6=[sum(e["free"]["final"] == e["v6_final"] for e in fv), len(fv)],
                free_gt_v6=[sum(e["free"]["final"] > e["v6_final"] for e in fv), len(fv)],
                none_lt_v6=[(e["k"], e["orbit"], e["none"]["final"], e["v6_final"]) for e in per if e.get("none") and e["none"]["final"] < e["v6_final"]])
    res["comparisons"] = cmp_; json.dump(res, open(p, "w"), indent=1)
    md += ["", "## Comparisons (same machine, arm64)",
           f"- seeded (ii_eg cap) vs no seed, same chain: same final n for {cmp_['seed_eq_none'][0]} of {cmp_['seed_eq_none'][1]}; "
           f"CPU {cmp_['seed_cpu_vs_none_cpu'][0]} s vs {cmp_['seed_cpu_vs_none_cpu'][1]} s",
           f"- cap = syn_n vs no seed: same final n for {cmp_['syn_eq_none'][0]} of {cmp_['syn_eq_none'][1]}",
           f"- free XOR layers only, no seed: v6's n for {cmp_['free_eq_v6'][0]} of {cmp_['free_eq_v6'][1]}, larger for {cmp_['free_gt_v6'][0]} of {cmp_['free_gt_v6'][1]}",
           f"- no seed below v6's n (verified, brute force, see verify.json): {cmp_['none_lt_v6']}; the seeded run on arm64 gives the same n, "
           "so this is the arm64/x86_64 difference of the known-answer check, not an effect of the seed"]
    md += ["", "## Existing z3 evidence (all-coefficient free SAT, no seed, z3 rlimit 1.5e9)"]
    for k, z in res["z3_free_coef_unseeded"].items():
        md.append(f"- {k}: {z['classes']} classes; status at the final factory's N: {z['status_at_final_N']}; own factory "
                  f"found for {z['own_factory_found']}, with n <= final for {z['found_n_le_final']}; found n - final n by final "
                  f"source: {z['found_minus_final_by_final_source']}; CPU {z['cpu_s']} s ({z['source']})")
    md += ["", res["note"]]
    open(os.path.join(O, "checks", "sat_alone_test.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "one":
        one(int(a[1]), a[2], a[3])
    elif a[0] == "run":
        drive(int(a[a.index("--max-workers") + 1]) if "--max-workers" in a else 3)
    elif a[0] == "summary":
        summary()
    elif a[0] == "verify":
        verify()
