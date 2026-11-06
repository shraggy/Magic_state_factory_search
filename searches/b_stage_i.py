"""B, stage (i) at l=3 beyond the exact class lists: every runnable group-symmetric sweep (orbit_group.dp_sweep over
G x H cells, G a group of the menu on the k outputs, H a check group on r <= rmax checks, syndrome space <= 2^bits_max),
each optimum expanded to rotations and verified (us.is_factory, us.distance != 1). Outputs that a CNOT circuit reduces
to fewer than k outputs are dropped (essential_dim < k); the rest are labelled by their Clifford-reduced form
(us.classify_factory: identity frame plus up to 4,000 random CNOT frames, an upper bound on the class beyond k = 4,
as in the larger sweep). Split over tasks by (G, H, r) combination.
    python searches/b_stage_i.py --list K OUT                 # -> OUT/l3_k{K}_combos.json (runnable combinations)
    python searches/b_stage_i.py --one K I M OUT              # combinations with index % M == I -> OUT/l3_k{K}/c{I}of{M}.json
    python searches/b_stage_i.py --merge K OUT                # best n per reduced output -> OUT/l3_k{K}_best.json
"""
import os, sys, json, time, glob
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import orbit_group as OG, orbit_group2 as OG2, unified_search as us
L, RMAX, BITS = 3, 4, 24


def combos(k):
    out = []
    for name, G, orbs, of, _ in OG.group_menu(k):
        for r in range(1, RMAX + 1):
            for hn, H, corbs, cof in OG2.check_groups(r):
                cells, cost, Hr, D, Dorb = OG2.build2(k, r, orbs, corbs)
                basis, h, cm = OG.rowbasis(Hr, D)
                if len(basis) <= BITS:
                    out.append(dict(G=name, H=hn, r=r, bits=len(basis)))
    return out


def find(k, name, hn, r):
    G = [m for m in OG.group_menu(k) if m[0] == name][0]
    Hh = [h for h in OG2.check_groups(r) if h[0] == hn][0]
    return G[2], Hh[2]


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "--list":
        k, out = int(sys.argv[2]), sys.argv[3]; os.makedirs(out, exist_ok=True)
        c = combos(k); json.dump(c, open(os.path.join(out, f"l3_k{k}_combos.json"), "w")); print(len(c), "combinations")
    elif m == "--one":
        k, I, M, out = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
        C = json.load(open(os.path.join(out, f"l3_k{k}_combos.json")))
        dst = os.path.join(out, f"l3_k{k}", f"c{I}of{M}.json"); os.makedirs(os.path.dirname(dst), exist_ok=True)
        best, stats = {}, []
        for i, c in enumerate(C):
            if i % M != I:
                continue
            t = time.time(); orbs, corbs = find(k, c["G"], c["H"], c["r"])
            res, info = OG2.with_cells(k, c["r"], orbs, corbs, OG.dp_sweep, BITS)
            nout = nkept = 0
            if res:
                out_, cells, cost, Dorb = res
                for dig, (n, x) in out_.items():
                    if not x:
                        continue
                    nout += 1
                    g = OG2.gates2(k, c["r"], orbs, corbs, cells, x); f = dict(k=k, r=c["r"], gates=g)
                    if not (us.is_factory(f, L) and us.distance(f, L) != 1):
                        continue
                    info2 = us.classify_factory(f, L)
                    if info2 is None or info2["essential_dim"] != k:
                        continue
                    nkept += 1; red = info2["reduced"]
                    if red not in best or (n, k + c["r"]) < (best[red]["n"], best[red]["N"]):
                        best[red] = dict(n=n, N=k + c["r"], G=c["G"], H=c["H"], frame_exact=info2["frame_exact"],
                                         gates=[[sorted(S), cc] for S, cc in g])
            stats.append(dict(c, outputs=nout, genuine=nkept, seconds=round(time.time() - t, 1)))
            print(json.dumps(stats[-1]), flush=True)
        json.dump(dict(k=k, chunk=[I, M], host=os.uname().nodename, stats=stats, best=best), open(dst + ".tmp", "w"))
        os.replace(dst + ".tmp", dst)
    elif m == "--merge":
        k, out = int(sys.argv[2]), sys.argv[3]; best = {}; secs = 0.0; files = glob.glob(os.path.join(out, f"l3_k{k}", "*.json"))
        for f in files:
            x = json.load(open(f)); secs += sum(s["seconds"] for s in x["stats"])
            for red, b in x["best"].items():
                if red not in best or (b["n"], b["N"]) < (best[red]["n"], best[red]["N"]):
                    best[red] = b
        json.dump(dict(k=k, files=len(files), cpu_seconds=round(secs, 1), outputs=len(best), best=best),
                  open(os.path.join(out, f"l3_k{k}_best.json"), "w"))
        print(f"k={k}: {len(best)} distinct genuine reduced outputs from {len(files)} files, {round(secs/3600,2)} CPU-h")
