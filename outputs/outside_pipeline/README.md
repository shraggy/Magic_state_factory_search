# Factories found outside the final pipeline

The paper reports only factories that the final pipeline finds (`outputs/staged_best.json` for the cells with a class
list; `outputs/staged_B`, `outputs/k78_pipeline` and `outputs/large_unified_*` beyond them). This folder keeps the
valid factories that earlier runs found and the final pipeline does not reproduce.

`smaller_than_pipeline.csv` (written by `searches/outside_pipeline.py`): 4 `l=3`, `k=6` classes where a run outside
the final pipeline has a smaller verified factory than the final pipeline, each `n=30` on `N=9` against `n=32`:

| orbit | target | T-count | synthillation `n` | source |
|---|---|---|---|---|
| `2320:19024` | `T4,CS03,CS12,CCZ023,CCZ025,CCZ045,CCZ124,CCZ135` | 13 | 34 | `staged_A_v3` (`iii_x_6_2e6.json`, stage `iii:xs:6`) |
| `2352:1674` | `T4,T5,CS03,CS12,CCZ012,CCZ035,CCZ045,CCZ134` | 11 | 32 | Test A, `outputs/checks/sat_alone_test/raw/k6_2352-1674_none.json` |
| `2352:3240` | `T4,T5,CS03,CS12,CCZ012,CCZ015,CCZ024,CCZ134` | 13 | 34 | Test A, `.../k6_2352-3240_none.json` |
| `2352:8944` | `T4,T5,CS03,CS12,CCZ023,CCZ034,CCZ045,CCZ134,CCZ135` | 13 | 34 | Test A, `.../k6_2352-8944_none.json` |

For `2352:1674` the final pipeline's `n=32` equals synthillation, so this class is not among the 330 below
synthillation; the factory here (30 against 32) would make a 331st, but it stays outside the pipeline and outside
the headline count. Columns: class, target, T-count, `n` and `N` of the factory, `n` of the final pipeline,
synthillation `n`, the run and file it came from, its gates and a Quirk link (for `N <= 14`).

**Test A rows.** `searches/sat_alone_test.py` (2026-11-05, arm64 laptop, CryptoMiniSat, 2e6 conflicts per query,
seven of the nine SAT layers of the third run, `N <= k+3`) ran the SAT chain without stages (i) and (ii). On these 3
of its 63 sampled classes it ends at 30, below the x86_64 pipeline's 32; the arm64 run started from the stage-(ii)
factory also gives 30, and on 2 of the 5 known-answer classes the per-layer values differ between arm64 and x86_64
(`outputs/checks/sat_alone_test.json`, `known_answer.mismatch`). The lower values are consistent with that
arm64/x86_64 difference; these classes have not been rerun on x86_64 without a seed, so the cause is not
established. All Test A factories pass the brute-force check (`outputs/checks/sat_alone_test/verify.json`) and are
re-verified by `searches/outside_pipeline.py`.

**History.** With the final pipeline taken from the v4 stage-(iii) run alone (archived in
`outputs/archive/v4_published/outputs/outside_pipeline/`), this file had 35 classes. The final pipeline now takes, per
class, the best of the three deterministic stage-(iii) runs v4, v5, v6 (`searches/build_staged_best3.py`,
`outputs/staged_final_b3`): it reaches the earlier `n` on 34 of the 35; on the remaining one it stays 2 inputs above (32 against 30).
All 35 now need fewer inputs than synthillation in the final pipeline (35 of 35; checked against
`outputs/staged_best.json`).

**How they were obtained.** Each of the 35 came from restricted SAT (or free SAT seeded by it) in one of these runs:

| run | classes (now) | classes (v4-only pipeline) | where |
|---|---|---|---|
| sweep A, stage (iii) as first run (2026-11-01/02) | 0 | 6 | `outputs/staged_A`, `outputs/staged_pc` |
| first XOR free-SAT rebuild (2026-11-03) | 0 | 2 | `outputs/solver_v2/staged_A_v2` |
| first stage-(iii) chain from the gluing seed (2026-11-03) | 1 | 27 | `outputs/solver_v2/staged_A_v3` |

In those runs the allowed output parts of restricted SAT came from `synthillation_compare.opt_B`, which solves the
Lempel condition in z3's global context and returns one of many valid factors; which one depended on every earlier z3
call in the same process (four successive calls on one class gave four different sets). The runs therefore searched a
history-dependent set of output parts, and in sweep A a shared z3 context also linked check-pattern SAT to the
synthesis SAT run after it. The final pipeline computes the parts in a fresh z3 context (`searches/det_parts.py`) and
runs every query in a fresh process (`searches/chain_iii2.py`): its answers are the same in any order and on every
x86_64 host tested, but it does not reach this one.

**Checks.** Every factory here was re-verified by `searches/outside_pipeline.py`: the factory condition (monomial form and
brute force on all `2^N` states), distance 2, `n`, and the exact class label.

A deterministic stage that tries several Lempel factors per class (`searches/multi_parts.py`, `searches/chain_iii3.py`)
is under test; if it reaches this factory it moves into the pipeline and leaves this folder.
