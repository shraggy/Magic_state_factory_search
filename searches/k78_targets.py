"""Target list for the k=7, 8 gluing pass (searches/k78_eg.py): every CCZ layer we hold, deduplicated, each tagged with
its GL(k,2)-invariant fingerprint (searches/k78_orbits.py; equal layers up to CNOTs on the outputs have equal
fingerprints, so distinct fingerprints are distinct output classes; the converse is not proved). Sources:
  paper   outputs/k78_pipeline/targets_paper62.json (the 57 on-bound layers and the 5 Table SIII(b) rows)
  cand    outputs/k{7,8}_candidates/cand_*.json and best.json (random-layer candidate search)
  B       outputs/k78_pipeline/staged_B_ccz.json (sweep B stage (i) outputs whose phase polynomial is a CCZ layer;
          written by searches/k78_staged_B_ccz.py on the cluster)
  orbit   outputs/k78_pipeline/orbits/k{7,8}_s*.json (one random representative per fingerprint)
Layers with fewer than k essential outputs are dropped.
    python searches/k78_targets.py    -> outputs/k78_pipeline/targets_all.json
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); sys.path.insert(0, HERE)
from k78_orbits import fingerprint, genuine
O = os.path.join(ROOT, "outputs")
src = []
for t in json.load(open(os.path.join(O, "k78_pipeline", "targets_paper62.json"))):
    src.append((t["k"], t["ccz"], "paper:" + t["id"]))
for k in (7, 8):
    d = os.path.join(O, f"k{k}_candidates")
    for f in sorted(glob.glob(d + "/cand_*.json")) + [d + "/best.json"]:
        for n, v in json.load(open(f)).items():
            for c in (v if isinstance(v, list) else [v]):
                src.append((k, c["ccz"], "cand:" + os.path.basename(f)))
    for r in json.load(open(os.path.join(O, "k78_pipeline", "staged_B_ccz.json")))[str(k)]:
        src.append((k, r["ccz"], f"B:{r['file']}:{r['G']}:{r['H']}:n{r['n']}"))
    for f in sorted(glob.glob(os.path.join(O, "k78_pipeline", "orbits", f"k{k}_s*.json"))):
        for L in json.load(open(f))["result"][str(k)]["reps"]:
            src.append((k, L, "orbit:" + os.path.basename(f)))
seen, out, fpid = {}, [], {}
for k, L, s in src:
    key = (k, tuple(sorted(tuple(sorted(x)) for x in L)))
    if key in seen:
        out[seen[key]]["src"].append(s); continue
    if not genuine([list(x) for x in key[1]], k):
        continue
    fp = fingerprint([list(x) for x in key[1]], k)
    fid = fpid.setdefault((k, fp), f"k{k}f{sum(1 for kk, _ in fpid if kk == k)}")
    seen[key] = len(out)
    out.append(dict(id=f"t{len(out)}", k=k, ccz=[list(x) for x in key[1]], fp=fid, src=[s]))
for t in out:                                         # paper targets keep their reference n
    for s in t["src"]:
        if s.startswith("paper:"):
            t["paper_id"] = s[6:]
json.dump(dict(fingerprints={f"{v}": dict(k=k, fp=fp) for (k, fp), v in fpid.items()}, targets=out),
          open(os.path.join(O, "k78_pipeline", "targets_all.json"), "w"))
from collections import Counter
print(len(out), "targets;", Counter(t["k"] for t in out), "; fingerprints:", Counter(k for k, _ in fpid))
