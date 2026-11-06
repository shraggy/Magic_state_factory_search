"""Check the Lagrangian enumerator of searches/full_clifford_classes.py (csrc/lagrangians.c) against the total number
prod_{i=1}^k (2^i+1) of Lagrangian subspaces of F_2^{2k}, for k = 1..6 (the self-test of full_clifford_classes.py
covers k = 1..5; k = 6 takes about 10 CPU minutes).
    python searches/check_lagrangian_counts.py   -> outputs/checks/lagrangian_counts.json
"""
import os, sys, json, subprocess, time
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
EXE = os.path.join(HERE, "csrc", "lagrangians"); CSRC = os.path.join(HERE, "csrc", "lagrangians.c")
if not os.path.exists(EXE) or os.path.getmtime(EXE) < os.path.getmtime(CSRC):
    subprocess.run(["cc", "-O2", "-o", EXE, CSRC], check=True)
rows = []
for k in range(1, 7):
    want = 1
    for i in range(1, k + 1): want *= 2 ** i + 1
    t = time.time(); p = subprocess.Popen([EXE], stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True)
    p.stdin.write(f"{k}\n0" + "1" * ((1 << (2 * k)) - 1) + "\n"); p.stdin.close()
    got = sum(1 for line in p.stdout if line.strip()); p.wait()
    rows.append(dict(k=k, enumerated=got, prod_2i_plus_1=want, ok=got == want, seconds=round(time.time() - t, 1)))
    print(rows[-1], flush=True)
out = dict(script="searches/check_lagrangian_counts.py", rows=rows, all_ok=all(r["ok"] for r in rows))
json.dump(out, open(os.path.join(ROOT, "outputs", "checks", "lagrangian_counts.json"), "w"), indent=1)
print("all_ok", out["all_ok"])
