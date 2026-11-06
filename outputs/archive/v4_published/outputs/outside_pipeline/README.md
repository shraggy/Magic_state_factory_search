# Factories found outside the final pipeline

The paper reports only factories that the final pipeline finds (`outputs/staged_best.json` for the cells with a class
list; `outputs/staged_B`, `outputs/k78_pipeline` and `outputs/large_unified_*` beyond them). This folder keeps the
valid factories that earlier runs found and the final pipeline does not reproduce.

`smaller_than_pipeline.csv` (written by `searches/outside_pipeline.py`): 35 `l=3`, `k=6` classes where an earlier run
has a smaller factory than the final pipeline, by 2 inputs for 34 and 4 for 1. All 35 have fewer inputs than
synthillation; for 17 of them the final pipeline also has fewer, for the other 18 it has the same. Columns: class,
target, T-count, `n` and `N` of the earlier factory, `n` of the final pipeline, synthillation `n`, the run and file it
came from, its gates and a Quirk link (for `N <= 14`).

**How they were obtained.** Each comes from restricted SAT (or free SAT seeded by it) in one of these runs:

| run | classes | where |
|---|---|---|
| sweep A, stage (iii) as first run (2026-11-01/02) | 6 | `outputs/staged_A`, `outputs/staged_pc` |
| first XOR free-SAT rebuild (2026-11-03) | 2 | `outputs/solver_v2/staged_A_v2` |
| first stage-(iii) chain from the gluing seed (2026-11-03) | 27 | `outputs/solver_v2/staged_A_v3` |

In those runs the allowed output parts of restricted SAT came from `synthillation_compare.opt_B`, which solves the
Lempel condition in z3's global context and returns one of many valid factors; which one depended on every earlier z3
call in the same process (four successive calls on one class gave four different sets). The runs therefore searched a
history-dependent set of output parts, and in sweep A a shared z3 context also linked check-pattern SAT to the
synthesis SAT run after it. The final pipeline computes the parts in a fresh z3 context (`searches/det_parts.py`) and
runs every query in a fresh process (`searches/chain_iii2.py`): its answers are the same in any order and on every
x86_64 host tested, but it does not reach these 35.

**Checks.** Every factory was re-verified by `searches/outside_pipeline.py`: the factory condition (monomial form and
brute force on all `2^N` states), distance 2, `n`, and the exact class label.

A deterministic stage that tries several Lempel factors per class (`searches/multi_parts.py`, `searches/chain_iii3.py`)
is under test; if it reaches these factories they move into the pipeline and leave this folder.
