"""Compare the four-qubit parent of Fig. 1(a) with Fig. 16(a) of Litinski, Quantum 3, 205 (2019).
Supports and signs of his 15 pi/8 rotations transcribed by hand from his Fig. 16(a) (qubits top to bottom 0..3):
four single-qubit -pi/8, then ZIIZ+, ZZIZ-, ZIZZ-, ZZZZ+, IZZZ-, IZIZ+, IIZZ+, ZZZI-, ZZII+, ZIZI+, IZZI+.
Checks (a) that his supports are all 15 nonempty subsets of 4 qubits, as in our parent, (b) that both his circuit and
ours (all signs +) are the identity up to a global phase on all 16 basis states, and (c) that the two differ only by
the signs of the odd-weight rotations, a Clifford change at l=3.
    python searches/litinski_fig16a_check.py -> outputs/checks/litinski_fig16a_check.json
"""
import json, os, itertools
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lit = [("ZIII", -1), ("IZII", -1), ("IIZI", -1), ("IIIZ", -1), ("ZIIZ", 1), ("ZZIZ", -1), ("ZIZZ", -1), ("ZZZZ", 1),
       ("IZZZ", -1), ("IZIZ", 1), ("IIZZ", 1), ("ZZZI", -1), ("ZZII", 1), ("ZIZI", 1), ("IZZI", 1)]
sup = lambda s: frozenset(i for i, c in enumerate(s) if c == "Z")
L = {sup(s): g for s, g in lit}
ours = {frozenset(S): 1 for w in range(1, 5) for S in itertools.combinations(range(4), w)}
def phases(circ):      # phase of Z^{S}(sign*pi/8) on |x>, in units of pi/8: -sign*(-1)^{x.S} up to global
    return {sum(g * (1 if sum(x[q] for q in S) % 2 else -1) for S, g in circ.items()) % 16
            for x in itertools.product((0, 1), repeat=4)}
res = dict(supports_equal=set(L) == set(ours), n_rotations=len(lit),
           litinski_identity=len(phases(L)) == 1, ours_identity=len(phases(ours)) == 1,
           differ_only_on_odd_weight_signs=all((L[S] == 1) == (len(S) % 2 == 0) for S in L))
json.dump(res, open(os.path.join(ROOT, "outputs", "checks", "litinski_fig16a_check.json"), "w"), indent=1)
print(res)
