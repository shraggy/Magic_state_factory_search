"""Assemble sweep A with the XOR free-SAT pass as iii_f, in a new directory (default outputs/solver_v2/staged_A_v2).

For each l=3 cell: copy the unchanged pass files i_tg .. iii_r (from outputs/staged_pc for k<=4, outputs/staged_A for
k=5,6) and write iii_f.json with searches/xor_merge.py from the XOR runs in XORDIR (the variants given, all seeded
from the copied iii_r). l=2 and l=4 cells have no free-SAT change and are not copied (sweep A's own files stand).
The published directories are only read; run where their class caches exist (_cache/l3_k{k}_classes.pkl),
otherwise the labeller would build and write one there.
    python searches/build_A_v2.py XORDIR [OUT] [--vars x_w5_c2e6,x_w6_c2e6,...]
"""
import os, sys, json, shutil
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import xor_merge as XM

PRE = ["i_tg", "i_gs", "i_gh", "ii_dj", "ii_eg", "iii_r"]
CELLS = [(1, "staged_pc"), (2, "staged_pc"), (3, "staged_pc"), (4, "staged_pc"), (5, "staged_A"), (6, "staged_A")]


def main(xordir, out, vars_):
    for k, src in CELLS:
        pk = os.path.join(ROOT, "outputs", src, "_cache", f"l3_k{k}_classes.pkl")
        assert os.path.exists(pk), f"missing {pk}: would be built inside a published directory"
    for k, src in CELLS:
        s = os.path.join(ROOT, "outputs", src, f"l3_k{k}"); d = os.path.join(out, f"l3_k{k}"); os.makedirs(d, exist_ok=True)
        for p in PRE:
            shutil.copyfile(os.path.join(s, p + ".json"), os.path.join(d, p + ".json"))
        # xor_merge reads seeds from SEEDDIR/l3_k{k}/iii_r.json and runs from OUT/l3_k{k}/xor_*; link the runs in
        xd = os.path.join(xordir, f"l3_k{k}")
        for v in sorted(os.listdir(xd) if os.path.isdir(xd) else []):
            if v.startswith("xor_") and (vars_ is None or v[4:] in vars_):
                dst = os.path.join(d, v)
                if not os.path.exists(dst):
                    shutil.copytree(os.path.join(xordir, f"l3_k{k}", v), dst)
        XM.merge(k, os.path.join(ROOT, "outputs", src), out, vars_)
    json.dump(dict(xordir=xordir, variants=sorted(vars_) if vars_ else "all", cells=[f"l3_k{k}" for k, _ in CELLS],
                   copied=PRE, note="iii_f = xor_merge over the variants; earlier passes copied unchanged"),
              open(os.path.join(out, "BUILD.json"), "w"), indent=1)


if __name__ == "__main__":
    a = sys.argv[1:]; v = None
    if "--vars" in a:
        i = a.index("--vars"); v = set(a[i + 1].split(",")); a = a[:i] + a[i + 2:]
    main(a[0], a[1] if len(a) > 1 else os.path.join(ROOT, "outputs", "solver_v2", "staged_A_v2"), v)
