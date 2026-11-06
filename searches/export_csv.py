"""The main-sweep catalogues from the staged cascade (outputs/staged_best.json), in the column layout of the earlier
files so that readers of them keep working:
  outputs/smallest_per_class.csv.zip      one row per output class (l=2, k<=6; l=3, k<=6; l=4, k<=4), with a Quirk link
                                          for circuits of up to 14 qubits
  outputs/all_factories.csv               the same factories, without the Quirk column (the earlier file also held
                                          non-smallest factories; the staged run keeps only the best per class)
  outputs/all_factories_with_quirk.csv.zip  all_factories.csv with the Quirk column
  outputs/l3_catalogue_k6.csv             the 766 l=3 rows (k, spec, tcount, n, N, c, source, syn_strongest_n, gates, quirk, verified)
Columns that described Clifford padding (padded_from, clifford_padded, masks) are left empty: the staged cascade does
not pad. third_level_target = yes on the 22 l=4 rows whose target is a level-3 gate (staged_best.py): third-level
targets, kept for completeness and not counted as fourth-level classes. Every row's factory was re-verified by searches/staged_best.py; here n, N and the class are copied.
    python searches/export_csv.py
"""
import os, sys, csv, json, zipfile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); O = os.path.join(ROOT, "outputs")
sys.path.insert(0, HERE)
import unified_search as us
import argparse as _ap
_p = _ap.ArgumentParser(); _p.add_argument("--sb", default="staged_best.json", help="staged-best file under outputs/")
_p.add_argument("--out-dir", default=None, help="write every output here instead of outputs/")
_A = _p.parse_args(); SB_FILE = _A.sb; OUTD = _A.out_dir

QUIRK_NMAX = 14
STAGE = {"i_tg": "two-group", "i_gs": "group-symmetric", "i_gh": "group-symmetric (check groups)",
         "ii_dj": "gluing", "ii_eg": "entangled-gluing", "iii_r": "restricted-SAT", "iii_f": "free-SAT",
         "iii_xs": "restricted-SAT", "iii_x": "free-SAT"}
COLS = ["clifford_class", "l", "k", "n", "c", "N", "d", "search_stage", "padded_from", "params", "reduced", "raw",
        "zcount", "essential_dim", "clifford_padded", "masks", "reference", "third_level_target"]


def main():
    SB = json.load(open(os.path.join(O, SB_FILE))); D = OUTD or O
    rows = []
    for r in sorted(SB["rows"], key=lambda r: (r["l"], r["k"], r["n"], r["N"], r["spec"])):
        g = [(frozenset(S), c) for S, c in r["gates"]]
        q = us.quirk_url(g, r["N"], r["l"], r["k"]) if r["N"] <= QUIRK_NMAX else ""
        rows.append(dict(clifford_class=f"L{r['l']}-k{r['k']}-{r['orbit']}", l=r["l"], k=r["k"], n=r["n"],
                         c="" if r["c"] is None else r["c"], N=r["N"], d=2, search_stage=STAGE[r["first"]],
                         padded_from="", params=r["stage"], reduced=r["spec"], raw=r["spec"],
                         zcount="" if r["tau"] is None else r["tau"], essential_dim=r["k"], clifford_padded="",
                         masks="", reference="", third_level_target="yes" if r.get("third_level_target") else "", quirk=q, gates=json.dumps(r["gates"]), syn_n=r["syn_n"], first=r["first"]))
    def write(path, cols, rs):
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, cols, extrasaction="ignore"); w.writeheader(); w.writerows(rs)
    tmp = os.path.join(D, "smallest_per_class.csv")      # 63 MB: shipped zipped (GitHub file limit)
    write(tmp, COLS + ["quirk"], rows)
    with zipfile.ZipFile(tmp + ".zip", "w", zipfile.ZIP_DEFLATED) as z:
        z.write(tmp, "smallest_per_class.csv")
    os.remove(tmp)
    write(os.path.join(D, "all_factories.csv"), COLS, rows)
    tmp = os.path.join(D, "all_factories_with_quirk.csv")
    write(tmp, COLS + ["quirk"], rows)
    with zipfile.ZipFile(tmp + ".zip", "w", zipfile.ZIP_DEFLATED) as z:
        z.write(tmp, "all_factories_with_quirk.csv")
    os.remove(tmp)
    l3 = [dict(k=r["k"], spec=r["reduced"], tcount=r["zcount"], n=r["n"], N=r["N"], c=r["c"], source=r["params"],
               syn_strongest_n=r["syn_n"], gates=r["gates"], quirk=r["quirk"], verified="yes") for r in rows if r["l"] == 3]
    write(os.path.join(D, "l3_catalogue_k6.csv"), ["k", "spec", "tcount", "n", "N", "c", "source", "syn_strongest_n",
                                                    "gates", "quirk", "verified"], l3)
    print(f"{len(rows)} rows ({len(l3)} at l=3); {sum(bool(r['quirk']) for r in rows)} with a Quirk link")


if __name__ == "__main__":
    main()
