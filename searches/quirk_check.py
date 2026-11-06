"""Check that every Quirk link in outputs/smallest_per_class.csv.zip contains exactly its factory: the first m parity-rotation
columns after the H layer (m = one column per power-of-two unit of each rotation, unified_search._units), decoded to
(support, coefficient mod 2^l), equal the gate list of the same class in outputs/staged_best.json. The columns after
them (Clifford corrections on the checks, output decomposition) are not compared here.
    python searches/quirk_check.py
"""
import os, sys, csv, json, math, cmath, urllib.parse, collections
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); O = os.path.join(ROOT, "outputs")
csv.field_size_limit(sys.maxsize)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import unified_search as us
STD = {"i": 1j, "-i": -1j, "√i": cmath.exp(1j * math.pi / 4), "√-i": cmath.exp(-1j * math.pi / 4), "Z": -1}


def factory_part(url, L, N, m):
    d = json.loads(urllib.parse.unquote(url.split("circuit=")[1])); custom = {g["id"]: g for g in d.get("gates", [])}
    acc, mod = collections.Counter(), 1 << L
    for c in d["cols"][1:1 + m]:
        S = tuple(q for q, g in enumerate(c) if g == "zpar")
        z = STD[c[N]] if c[N] in STD else complex(custom[c[N]]["matrix"].split(",")[0].strip("{}").replace("i", "j"))
        acc[S] = (acc[S] + round(cmath.phase(z) / (2 * math.pi) * mod)) % mod
    return sorted((S, c) for S, c in acc.items() if c)


def main():
    SB = {f"L{r['l']}-k{r['k']}-{r['orbit']}": r for r in json.load(open(os.path.join(O, "staged_best.json")))["rows"]}
    import io, zipfile
    with zipfile.ZipFile(os.path.join(O, "smallest_per_class.csv.zip")) as z:
        rows = list(csv.DictReader(io.TextIOWrapper(z.open("smallest_per_class.csv"), encoding="utf-8")))
    bad = []
    for r in rows:
        b = SB[r["clifford_class"]]; L = b["l"]
        g = sorted((tuple(sorted(S)), c % (1 << L)) for S, c in b["gates"] if c % (1 << L))
        m = sum(len(us._units(c, L)) for _, c in b["gates"])
        if g != factory_part(r["quirk"], L, b["N"], m):
            bad.append(r["clifford_class"])
    res = dict(links=len(rows), match=len(rows) - len(bad), mismatch=bad)
    json.dump(res, open(os.path.join(O, "checks", "quirk_check.json"), "w"), indent=1)
    print(f"{res['match']} of {res['links']} Quirk links contain exactly their factory")


if __name__ == "__main__":
    main()
