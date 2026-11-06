# Test A: SAT stage (iii) without stages (i)-(ii)

Script: searches/sat_alone_test.py; data: outputs/checks/sat_alone_test.json, raw per job in outputs/checks/sat_alone_test/raw/.

## Known-answer check (seed = ii_eg n, as v6)
- k=3 24: per-layer n here [12, 12, 12, 12, 12, 12, 12] vs v6 [12, 12, 12, 12, 12, 12, 12]: match; CPU 3.0 s here, 11.2 s v6
- k=4 112: per-layer n here [16, 16, 16, 16, 16, 16, 16] vs v6 [16, 16, 16, 16, 16, 16, 16]: match; CPU 99.5 s here, 184.2 s v6
- k=5 120: per-layer n here [22, 22, 22, 22, 22, 22, 22] vs v6 [22, 22, 22, 22, 22, 22, 22]: match; CPU 306.6 s here, 681.5 s v6
- k=6 2304:136: per-layer n here [22, 22, 22, 22, 22, 22, 22] vs v6 [24, 24, 22, 22, 22, 22, 22]: MISMATCH; CPU 640.5 s here, 1356.5 s v6
- k=6 2320:9041: per-layer n here [30, 30, 30, 30, 30, 30, 30] vs v6 [30, 30, 30, 30, 30, 28, 28]: MISMATCH; CPU 1077.2 s here, 1591.0 s v6

## Sample: 63 classes (k=5: every 12th of 65, k=6: every 12th of 679)
- mode none: 63 of 63 classes run
- mode free: 21 of 63 classes run
- mode seed: 8 of 63 classes run
- mode syn: 21 of 63 classes run
  - free|final from i/ii: larger than cascade n: 1, reaches cascade n: 8
  - free|final from iii: larger than cascade n: 5, reaches cascade n: 7
  - none|final from i/ii: reaches cascade n: 32
  - none|final from iii: below cascade n: 3, reaches cascade n: 28
  - seed|final from i/ii: reaches cascade n: 3
  - seed|final from iii: below cascade n: 3, reaches cascade n: 2
  - syn|final from i/ii: reaches cascade n: 9
  - syn|final from iii: below cascade n: 1, reaches cascade n: 11

No-seed runs: layer in which the run's final n was first found: x2:5:2e6: 1, x:5:2e6: 9, xs:None:2e6: 48, xsb:8:None:2e6: 5

## CPU (sum over classes run)
- none: SAT alone 57057.0 s vs v6 stage (iii) (7 layers) 121268.1 s over 63 classes
- free: SAT alone 7517.2 s vs v6 stage (iii) (7 layers) 38030.8 s over 21 classes
- seed: SAT alone 6797.0 s vs v6 stage (iii) (7 layers) 12983.8 s over 8 classes
- syn: SAT alone 18608.0 s vs v6 stage (iii) (7 layers) 38030.8 s over 21 classes

## Comparisons (same machine, arm64)
- seeded (ii_eg cap) vs no seed, same chain: same final n for 8 of 8; CPU 6797.0 s vs 6960.3 s
- cap = syn_n vs no seed: same final n for 21 of 21
- free XOR layers only, no seed: v6's n for 15 of 21, larger for 6 of 21
- no seed below v6's n (verified, brute force, see verify.json): [(6, '2352:1674', 30, 32), (6, '2352:3240', 30, 32), (6, '2352:8944', 30, 32)]; the seeded run on arm64 gives the same n, so this is the arm64/x86_64 difference of the known-answer check, not an effect of the seed

## Existing z3 evidence (all-coefficient free SAT, no seed, z3 rlimit 1.5e9)
- k5: 65 classes; status at the final factory's N: {'unknown': 63, 'not run': 1, 'unsat': 1}; own factory found for 65, with n <= final for 35; found n - final n by final source: {'i/ii': {'0': 28, '4': 5, '2': 8, '8': 3}, 'iii': {'4': 5, '2': 9, '0': 7}}; CPU 554808.1 s (outputs/staged_A/l3_k5/free_coef_1.5e9_unseeded/*.json)
- k6: 679 classes; status at the final factory's N: {'unsat': 3, 'unknown': 672, 'not run': 4}; own factory found for 679, with n <= final for 56; found n - final n by final source: {'i/ii': {'0': 44, '16': 53, '20': 25, '14': 40, '10': 23, '18': 50, '24': 4, '12': 33, '22': 7, '4': 3, '8': 6, '26': 2, '6': 4, '2': 3, '34': 1}, 'iii': {'14': 63, '18': 88, '16': 81, '12': 26, '20': 46, '10': 20, '0': 12, '22': 20, '8': 6, '36': 1, '4': 1, '24': 8, '2': 7, '30': 1, '6': 1}}; CPU 4930919.2 s (outputs/staged_A/l3_k6/free_coef_1.5e9_unseeded/*.json)

v6 ran on x86_64 (cluster); this run on arm64. 'nothing found' means not reached within 2e6 conflicts per query, not infeasible.
