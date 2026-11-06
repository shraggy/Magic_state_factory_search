"""Catalogue of the search beyond the class lists, in the column layout of export_csv.py (sweep A), with a Quirk link
for every factory. Three parts of the one cascade, by regime:
  l=3, k=7..9, stage (i) group-symmetric   outputs/staged_B/l3_k{K}/*.json          (searches/b_stage_i.py; k=8,9 partial)
  l=3, k=7..9, stage (ii) disjoint gluing  outputs/k78_pipeline/full/l3_k{K}/ii_dj.json  (searches/k789_ii_dj.py)
  l=3, k=7,8, CCZ layers, stages (ii) entangled gluing and (iii) free SAT (XOR) on a target list
                                           outputs/k78_pipeline/all_isd/*.json, full/l3_k{7,8}/iii_f.json
  l=2,3,4, k<=12, N<=18, stages (i) two-group and (ii) disjoint gluing (the wide range)
                                           outputs/large_unified_with_quirk.csv.zip (searches/large_unified.py)
There is no class list here, so `clifford_class` holds the Clifford-reduced output label (an upper bound on the class:
labels depend on the CNOT frame) and, for CCZ layers, the GL(k,2) fingerprint and the target id. T-count is exact where
known (`zcount`), else `tcount_lo`-`tcount_hi`. Every factory with gates is re-verified (us.is_factory, distance != 1,
n recomputed) before it is written; the wide-range rows are copied (verified when that run wrote them).
    python searches/export_beyond_csv.py   -> outputs/beyond_class_lists_with_quirk.csv.zip
"""
import os, sys, csv, json, glob, io, zipfile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(ROOT, "classification"))
import unified_search as us
csv.field_size_limit(1 << 30)                 # wide-range Quirk links exceed csv's default field limit

QUIRK_NMAX = 14
COLS = ["regime", "clifford_class", "l", "k", "n", "c", "N", "d", "search_stage", "params", "reduced", "raw", "zcount",
        "tcount_lo", "tcount_hi", "on_bound", "essential_dim", "frame_exact", "gates", "quirk"]


def checked(k, N, gates, n):
    g = [(frozenset(S), c) for S, c in gates]
    f = dict(k=k, r=N - k, gates=g)
    assert us.is_factory(f, 3) and us.distance(f, 3) != 1, (k, N, n)
    assert sum(us.inputs(c, 3) for _, c in g) == n, (k, N, n)
    return g


def row(regime, label, k, n, N, stage, params, gates, lo=None, hi=None, reduced="", frame_exact=""):
    g = checked(k, N, gates, n)
    exact = lo is not None and lo == hi
    return dict(regime=regime, clifford_class=label, l=3, k=k, n=n, c="", N=N, d=2, search_stage=stage, params=params,
                reduced=reduced, raw=reduced, zcount=lo if exact else "", tcount_lo="" if lo is None else lo,
                tcount_hi="" if hi is None else hi, on_bound="" if lo is None else bool(exact and n == lo + 1),
                essential_dim=k, frame_exact=frame_exact, gates=json.dumps([[sorted(S), c] for S, c in g]),
                quirk=us.quirk_url(g, N, 3, k) if N <= QUIRK_NMAX else "")


def main():
    rows = []
    for k in (7, 8, 9):                                   # stage (i), group-symmetric, best per reduced output
        best = {}
        for f in glob.glob(os.path.join(O, "staged_B", f"l3_k{k}", "*.json")):
            for red, v in json.load(open(f))["best"].items():
                if red not in best or (v["n"], v["N"]) < (best[red]["n"], best[red]["N"]):
                    best[red] = v
        for red, v in sorted(best.items(), key=lambda x: (x[1]["n"], x[0])):
            rows.append(row("l=3 k=7-9", red, k, v["n"], v["N"], "i_gs", f"G={v['G']},{v['H']}", v["gates"],
                            reduced=red, frame_exact=v["frame_exact"]))
        dj = json.load(open(os.path.join(O, "k78_pipeline", "full", f"l3_k{k}", "ii_dj.json")))["best"]
        for red, v in sorted(dj.items()):                  # stage (ii), disjoint gluing
            rows.append(row("l=3 k=7-9", red, k, v["n"], v["N"], "ii_dj", v["stage"].split(":", 1)[1], v["gates"],
                            reduced=red))
    for f in sorted(glob.glob(os.path.join(O, "k78_pipeline", "all_isd", "*.json"))):    # CCZ layers, every target
        for t in json.load(open(f))["results"]:
            b = t["best"]
            spec = ",".join("CCZ" + "".join(map(str, x)) for x in t["ccz"])
            rows.append(row("l=3 k=7,8 CCZ layers", f"{t['fp']}:{t['id']}", t["k"], b["n"], b["N"], "ii_eg",
                            b["kind"], b["gates"], t["tcount_lo"], t["tcount_hi"], reduced=spec))
    zf = zipfile.ZipFile(os.path.join(O, "large_unified_with_quirk.csv.zip"))      # the wide range, copied
    name = [x for x in zf.namelist() if x.endswith(".csv")][0]
    wide = list(csv.DictReader(io.TextIOWrapper(zf.open(name), encoding="utf-8")))
    for r in wide:
        rows.append(dict(regime="wide range", clifford_class=r["reduced"], l=r["l"], k=r["k"], n=r["n"], c=r.get("c", ""),
                         N=r["N"], d=r["d"], search_stage=r["search_stage"], params=r["params"], reduced=r["reduced"],
                         raw=r["raw"], zcount=r["zcount"], tcount_lo=r["zcount"], tcount_hi=r["zcount"], on_bound="",
                         essential_dim=r["essential_dim"], frame_exact=r["frame_exact"], gates="", quirk=r.get("quirk", "")))
    buf = io.StringIO(); w = csv.DictWriter(buf, COLS, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
    out = os.path.join(O, "beyond_class_lists_with_quirk.csv.zip")
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("beyond_class_lists_with_quirk.csv", buf.getvalue())
    by = {}
    for r in rows: by[r["regime"]] = by.get(r["regime"], 0) + 1
    print(f"{len(rows)} rows ({by}); {sum(bool(r['quirk']) for r in rows)} with a Quirk link -> {out}")


if __name__ == "__main__":
    main()
