"""Cause check for the known-answer mismatch of sat_alone_test.py (rule 8): v6 (x86_64) lowered k=6 class 2320:9041 from
30 to 28 in layer xsb:8:None:2e6; on arm64 the same layer at 2e6 conflicts did not. Reruns that one layer from cap 30 at
higher conflict limits (fresh process per query, as in v6) and reports whether n=28 is reached.
    python searches/sat_alone_probe.py   -> outputs/checks/sat_alone_test/probe_2320-9041.json
"""
import os, sys, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))

if __name__ == "__main__":
    import chain_iii4, chain_iii2 as C2, staged_class as SC, staged as STG
    cache, cl = SC.setup(3, 6, os.path.join(ROOT, "outputs", "staged_A")); o = "2320:9041"; spec = cl[o]["spec"]; out = []
    for c in ("2e6", "4e6", "8e6"):
        t = time.time(); r = C2.stage("xsb:8:None:" + c, 6, spec, 30)
        lab = [(STG.label_fix(dict(k=6, r=N - 6, gates=[(frozenset(S), cc) for S, cc in g]), 3, 6, cache)[0], n, N) for n, N, g in r]
        out.append(dict(conflicts=c, found=lab, seconds=round(time.time() - t, 1), machine=os.uname().machine)); print(out[-1], flush=True)
        json.dump(out, open(os.path.join(ROOT, "outputs", "checks", "sat_alone_test", "probe_2320-9041.json"), "w"), indent=1)
        if any(x[0] == o and x[1] <= 28 for x in lab):
            break
