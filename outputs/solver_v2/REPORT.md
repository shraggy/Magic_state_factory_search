# Solver efficiency (draft appendix): sweep A, l=3, k<=6

Saved by the main agent from the solver agent's final hand-off (2026-11-03); the agent could not write report files.
Status: COMPLETE (final hand-off 2026-11-03). Every number below is computed; the headline counts (110/656/0, 232 inputs)
were re-derived independently by the main agent from staged_A_v2/ against outputs/staged_best.json, with every lowered
factory re-verified (is_factory, distance 2, n).

## Method
At l=3, free SAT reduces to choosing a set of coefficient-1 rotations subject to GF(2) parity rows. Sweep A encoded
each row as an integer equation, `sum = 2z + b`, and solved it with CP-SAT (`bool_sat.py`, supports of at most 5
qubits, 600 deterministic s per query). It also ran z3 over coefficients (step limit 1.5e9). The new form,
`searches/xor_sat2.py` (md5 fad4496cfce184c216c810eb89504f54), passes each row to CryptoMiniSat 5.17.0 as a native
XOR clause, which the solver handles by Gauss-Jordan elimination. It encodes n < cap as a sequential counter and keeps
one solver per (class, N), lowering the cap through assumptions. The limit is a number of conflicts per query, on one
thread. Seeds (the published iii_r best: `outputs/staged_A` for k=5,6, `outputs/staged_pc` for k<=4; no `:unseeded`
or `:seed=` flag), N = k+1..k+3, the bound skip and the exact labelling are those of sweep A.

## Checks
- Known answers: n = 8, 12, 14, 12 for CCZ012, CS01, T0,T1 and CCZ014,CCZ023, each proved smallest.
- Every factory found is re-verified: is_factory monomial and brute force, distance 2, exact class label, n
  recomputed (38 of 38 rows of `lowered_partial.json`). None is below n = tau+1.
- Reproducible per binary build: the instance is byte-identical on laptop and cluster; answers at a fixed limit were
  identical on 6 cluster nodes, Intel 6258R and AMD EPYC 7763 (job 1513410); the laptop's arm64 build can give a
  different n at the same limit.

## Speed
The free-SAT pass at k=6 cost 5,587,568 CPU-s in sweep A (`outputs/staged_A/l3_k6/iii_f.json`). On the 676 searched
k=6 classes, the CP-SAT form alone took 656,649 s. The XOR form at 2e5 conflicts and supports <= 5 took 11,393 s,
with a longest class of 89 s: 58x less than CP-SAT alone. On all 766 classes, all six XOR variants together took
494,800 s, against 6,232,284 s for sweep A's free SAT: 12.6x less.

## Result (`outputs/solver_v2/acceptance_iii_f.json`, `searches/xor_acceptance.py`)
| | better | equal | worse | CPU-s |
|---|---|---|---|---|
| union of all 6 variants | 110 | 656 | 0 | 516,010 |
| sweep A iii_f, same classes | | | | 6,232,284 |

12.1x less CPU. Per variant (better / equal / worse, CPU-s):

| variant | 2e5 | 2e6 |
|---|---|---|
| w5 | 28 / 735 / 3, 12,287 | 97 / 669 / 0, 137,994 |
| w6 | 6 / 753 / 7, 15,555 | 39 / 726 / 1, 149,247 |
| wNone | 3 / 755 / 8, 18,005 | 15 / 749 / 2, 182,924 |

- w5 at 2e6 alone passes (0 worse) at 45x less CPU than sweep A's free SAT.
- Most gains are at N = k+3 (N=9 at k=6), where CP-SAT returned UNKNOWN on 676 of 679 k=6 classes.

## Sweep A rebuilt (`outputs/solver_v2/staged_A_v2/`, `searches/build_A_v2.py`)
Passes i_tg..iii_r copied unchanged; iii_f is the merge of the 6 variants (`searches/xor_merge.py`). Gain over iii_r:
6 classes at k=5 (14 inputs), 117 at k=6 (248 inputs), none at k<=4. Against staged_best.json: 110 lower (232 inputs),
656 equal, 0 higher, none below n = tau+1. Checks (`searches/check_A2.py`, `staged_A_v2/checks/`): 766 of 766 valid,
exact class, n, distance 2, c; brute force 759 of 759 checked, 7 k=6 factories with N > 9 skipped by design.

## Step 2: seeded from ii_eg (separate comparison, not merged; `outputs/solver_v2/step2_compare.json`)
53 rows in 45 distinct k=6 classes go below staged_A_v2 (mostly by 2; 2352:3258 34 -> 30), all at N=9, all verified.
3 of the 45 are in the 9-class investigation (2320:16992, 82505, 82512).

| run (k=6, 679 classes) | CPU-s | below staged_best | below staged_A_v2 |
|---|---|---|---|
| z3 synthesis SAT, seeded at ii_eg n | 699,856 | 35 | 12 |
| XOR synthesis-restricted, wNone, 2e6 | 116,952 | 87 | 40 |
| XOR free, wNone, 2e6 | 142,457 | 13 | 1 |

At k=5 each run goes below staged_best on 2 classes, none below staged_A_v2. The seed changes what is found; the XOR
form of synthesis SAT is 6x cheaper than z3 synthesis and finds the most.

## Stage (i), optional, parked
`searches/csrc/xor_dp.c` with `searches/dp_c.py` agrees with dp_sweep on 40 of 40 combinations at k=4 and 42 of 42 at
k=5 (`dp_c_check_k{4,5}.json`); 8x faster at 24 bits; a 28-bit k=7 combination in 6.8 s. In sweep B the cost is
labelling (about 18 s per output in `us.classify_factory`), not the DP (2-4 s).

## The limit
The exact class lists stop near k=6 at l=3 and k=4 at l=4: at l=3 the output group has dimension
C(k,1)+C(k,2)+C(k,3) (2^41 elements at k=6, 2^63 at k=7), and the CNOT (affine) group has order |GL(k,2)| 2^k (about
2^40 at k=6, 2^54 at k=7). A faster solver shortens the search per class but cannot produce a class list beyond that.

## Files
Scripts (new, in searches/): xor_sat2.py, xor_free.py, xor_merge.py, solver_v2_compare.py, xor_acceptance.py,
syn_eg.py, build_A_v2.py, check_A2.py, dp_c.py, csrc/xor_dp.c. Laptop outputs: outputs/solver_v2/{calib, probe,
lowered_partial.json, dp_c_check_k4.json, dp_c_check_k5.json, acceptance_iii_f.json, step2_compare.json, full/,
staged_A_v2/}; cluster copies in ~/ccc_projects/msf_sat/outputs/solver_v2/.
