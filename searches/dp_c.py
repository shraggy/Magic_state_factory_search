"""Stage (i) DP in C with a one-byte table: a drop-in for orbit_group.dp_sweep reaching syndrome spaces of 2^25..2^34.

orbit_group.dp_sweep keeps an int64 table of 2^bits entries plus one packed take-bit array per cell (bits <= 24 by
default). csrc/xor_dp.c keeps one uint8 per state (saturating, 255 = not reachable below 255 inputs) and recovers an
optimum from the final table alone, so memory is 2^bits bytes (2^30 = 1 GiB, 2^32 = 4 GiB). Cells of cost >= 255
cannot be in any optimum below 255 and are dropped; so the outputs reported are exactly those of dp_sweep with
min n <= 254 (all sweep-A/B factories have n < 100). Same return value as dp_sweep: ((out, cells, cost, Dorb), info).
    python searches/dp_c.py --check K       # agreement with dp_sweep (min n per output) on every G x H combo, bits <= 24
"""
import os, sys, subprocess, time, json
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import orbit_group as OG

SRC = os.path.join(HERE, "csrc", "xor_dp.c")
BIN = os.environ.get("XOR_DP_BIN", os.path.join(os.path.expanduser("~"), "msf_scratch", "xor_dp"))


def binary():
    if not os.path.exists(BIN) or os.path.getmtime(BIN) < os.path.getmtime(SRC):
        os.makedirs(os.path.dirname(BIN), exist_ok=True)
        subprocess.check_call(["cc", "-O3", "-march=native", "-o", BIN, SRC])
    return BIN


def dp_sweep_c(k, r, orbs, bits_max=32):
    cells, cost, H, D, Dorb = OG.build(k, r, orbs)          # OG.build is swapped by orbit_group2.with_cells for G x H
    basis, h, combos = OG.rowbasis(H, D)
    bits = len(basis)
    info = dict(cells=len(cells), cond_rows=len(H), out_rows=len(D), h=h, d=bits - h, bits=bits, impl="c-uint8")
    if bits > bits_max:
        return None, info
    cols = [sum(((basis[b] >> ci) & 1) << b for b in range(bits)) for ci in range(len(cells))]
    cheapest = {}
    for ci in range(len(cells)):
        if cols[ci] and cost[ci] < 255 and (cols[ci] not in cheapest or cost[ci] < cost[cheapest[cols[ci]]]):
            cheapest[cols[ci]] = ci
    use = sorted(cheapest.values()); info["cells_used"] = len(use)
    inp = f"{bits} {h} {len(use)}\n" + "".join(f"{cols[ci]} {cost[ci]}\n" for ci in use)
    t = time.time()
    p = subprocess.run([binary()], input=inp, capture_output=True, text=True, check=True)
    info["c_seconds"] = round(time.time() - t, 2)
    out = {}
    for line in p.stdout.split("\n"):
        if not line:
            continue
        a = list(map(int, line.split())); tt, n, x = a[0], a[1], sorted(use[j] for j in a[2:])
        s = tt << h
        acc = 0
        for ci in x:
            acc ^= cols[ci]
        assert acc == s and sum(cost[ci] for ci in x) == n, "C traceback inconsistent"
        dig = tuple(int(bin(s & cm).count("1") & 1) for cm in combos)
        out[dig] = (n, x)
    info["reachable_outputs"] = len(out)
    return (out, cells, cost, Dorb), info


def check(k, bits_max=24):
    import orbit_group2 as OG2
    tot = agree = 0; rows = []
    for name, G, orbs, of, _ in OG.group_menu(k):
        for r in range(1, 5):
            for hn, Hg, corbs, cof in OG2.check_groups(r):
                t = time.time(); a, ia = OG2.with_cells(k, r, orbs, corbs, OG.dp_sweep, bits_max); ta = time.time() - t
                if a is None:
                    continue
                t = time.time(); b, ib = OG2.with_cells(k, r, orbs, corbs, dp_sweep_c, bits_max); tb = time.time() - t
                A = {d: v[0] for d, v in a[0].items() if v[0] <= 254}; Bm = {d: v[0] for d, v in b[0].items()}
                ok = A == Bm; tot += 1; agree += ok
                rows.append(dict(G=name, H=hn, r=r, bits=ia["bits"], outputs=len(A), same=ok,
                                 numpy_s=round(ta, 2), c_s=round(tb, 2)))
                if not ok:
                    print("MISMATCH", rows[-1], flush=True)
    print(f"k={k}: dp_sweep_c equals dp_sweep (min n per output) on {agree} of {tot} G x H combinations", flush=True)
    return rows


if __name__ == "__main__":
    if sys.argv[1] == "--check":
        k = int(sys.argv[2]); rows = check(k)
        d = os.path.join(ROOT, "outputs", "solver_v2"); os.makedirs(d, exist_ok=True)
        json.dump(rows, open(os.path.join(d, f"dp_c_check_k{k}.json"), "w"), indent=0)
