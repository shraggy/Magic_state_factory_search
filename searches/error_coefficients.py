"""Leading-order output error and acceptance of the Fig. 1 factories and of the prior construction for each target.

Noise model (as in the paper): every input magic state carries an independent Z fault with probability p; a
fault on an input of the parity rotation on support S applies Z_S. The factory accepts when every check qubit
returns the trivial X outcome, i.e. when the check parts of the faults cancel. All factories here have
distance 2, so acceptance is 1 - n p + O(p^2) and the probability that the accepted output state is wrong is
c p^2 + O(p^3), where c counts the pairs of faults whose check parts cancel and whose output parts do not.
For "distil the parts separately" the targets' parts come from independent factories, so n and c add.
Library: coefficients(gates, k, L) -> (n, c), used by staged_best.py and check_A.py.
"""
import csv, json, math, os, sys, urllib.parse
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
csv.field_size_limit(sys.maxsize)
import unified_search as us
ROOT = os.path.dirname(HERE)


def decode(url):
    """Factory rotations (support, coefficient units of pi/2^L) from a Quirk link of quirk_url()."""
    d = json.loads(urllib.parse.unquote(url.split("circuit=")[1])); custom = {g["id"]: g for g in d.get("gates", [])}
    cols = d["cols"]; b = next((i for i, c in enumerate(cols) if "Bloch" in c), len(cols))
    def ang(g):
        if g in ("i", "-i"): return (1 if g == "i" else -1) / 4
        if g in ("√i", "√-i"): return (1 if g == "√i" else -1) / 8
        z = complex(custom[g]["matrix"].split(",")[0].strip("{}").replace("i", "j")); return math.atan2(z.imag, z.real) / (2 * math.pi)
    return [(frozenset(q for q, g in enumerate(c) if g == "zpar"), ang([x for x in c if x not in (1, "zpar")][0]))
            for c in cols[:b] if "zpar" in c]


def gates_of(row):
    L = int(row["l"]); mod = 1 << L; g = {}
    for S, a in decode(row["quirk"]): g[S] = (g.get(S, 0) + round(a * mod)) % mod      # e^{2 pi i a} = a 2^L units of pi/2^L
    return [(S, c) for S, c in g.items() if c]


def coefficients(gates, k, L):
    """(n, c): input count and number of undetected weight-2 fault sets that corrupt the output."""
    out = set(range(k)); faults = []
    for S, c in gates:
        for _ in range(us.inputs(c, L)):                  # one fault location per input magic state
            faults.append((frozenset(S - out), frozenset(S & out)))
    assert all(ch for ch, _ in faults), "a single fault goes undetected: distance 1"
    bad = sum(1 for (c1, o1), (c2, o2) in combinations(faults, 2) if c1 == c2 and o1 != o2)
    return len(faults), bad


def syn_gates(names, mu):
    """The synthillation construction for a target, as parity rotations (one T each), outputs first."""
    import synthillation_compare as sc
    k, target = sc.tgt(names); A = sc.opt_A(k, target); tau = len(A)
    lin = any(len(T) == 1 for T in target); B = sc.opt_B(k, target, mu)
    if lin:
        tm = sc.C14[{(0, 0): 1, (0, 1): 2, (1, 0): 3, (1, 1): 4}[(tau % 2, mu % 2)]]
        c = sum(1 << i for i in range(k) if frozenset([i]) in target); blocks = [A, B, B] + [[c]] * 4
    else:
        tm = sc.T5[(tau % 2, mu % 2)]; blocks = [A, B, B]
    blocks += [[0]] * (len(tm[0]) - len(blocks)); R = len(tm)
    cols = [v | sum(1 << (k + r) for r in range(R) if tm[r][bi]) for bi, blk in enumerate(blocks) for v in blk]
    assert sc.coeffs_ok(cols, k + R, k, target, 3)
    return k, [(frozenset(q for q in range(k + R) if y >> q & 1), 1) for y in cols]


def set_partitions(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for p in set_partitions(rest):
        yield [[first]] + p
        for i in range(len(p)):
            yield p[:i] + [[first] + p[i]] + p[i + 1:]
