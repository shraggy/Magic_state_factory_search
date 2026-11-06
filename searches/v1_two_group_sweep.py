"""Reproduce the first-submission (arXiv:2606.28518v1) two-group sweep statistics quoted in the response letter.
Sweep: l in {2,3,4}, k in 1..7, n = k+1..11, s_total, s_O in 1..7, s_S = 1; every 2^m sign pattern on check-only gates.
Uses searches/Two_group.py unchanged. Expected: tuples 7203, sign patterns 107940, valid 21920
(l=2: 15600, l=3: 4678, l=4: 1642), groups (l,k,N,d,output) 689, of which non-trivial 535."""
import json, sys, os
from collections import Counter
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import Two_group as T
tuples = signs = 0; rows = []
for l in (2, 3, 4):
    for k in range(1, 8):
        for n in range(k + 1, 12):
            for st in range(1, 8):
                for so in range(1, 8):
                    p, Wt, Wo = T.build_gate_set(k, n - k, st, so)
                    if not p:
                        continue
                    tuples += 1; signs += 2 ** sum(1 for x in p if x[0] == 0)
                    for sd in T.find_valid_signs(l, n, k, p):
                        N = T.extract_factory(k, n - k, sd, Wo, Wt)
                        d = T.factory_distance(k, n - k, sd)
                        o, _ = T.classify_output(l, k, n - k, sd, set(Wo) & set(Wt))
                        rows.append((l, k, N, d, o))
out = {"tuples": tuples, "sign_patterns": signs, "valid_configurations": len(rows),
       "valid_per_l": {str(a): b for a, b in sorted(Counter(r[0] for r in rows).items())},
       "groups_l_k_N_d_output": len(set(rows)),
       "groups_nontrivial": len({r for r in rows if r[4] != 'trivial'})}
print(json.dumps(out, indent=1))
