# Searches

A unified search, SAT-free except for an optional last-resort SAT stage, built from three complementary numerical searches for borrowed-identity
distillation factories, plus a closed-form family with no search needed. See the [paper](https://arxiv.org/abs/2606.28518) for
the theory (borrowed-identity condition, Theorems 2/3/5); this file covers how to run each script
and what its parameters mean.

## The Clifford hierarchy, briefly

Each script is parametrized by a level `l`: it searches for factories that consume noisy copies of
the level-`l` magic state `|θ⟩ := |0⟩ + e^(2iθ)|1⟩`, `θ = π/2^l`, and Clifford operations, to
produce `k` higher-fidelity copies at distance `d`.

| `l` | weight-1 gate | weight-2 gate | weight-3 gate | weight-4 gate |
|---|---|---|---|---|
| 2 | `S`  | `CZ`  | —     | —      |
| 3 | `T`  | `CS`  | `CCZ` | —      |
| 4 | `√T` | `CT`  | `CCS` | `CCCZ` |

`l = 3` is where `T` lives; `l = 2` is Clifford; `l ≥ 4` are higher levels. All three scripts below
are written generically in `l`; the Clifford hierarchy level isn't hardcoded anywhere.

## The searches

| Search | Files | Purpose |
|---|---|---|
| **Unified** | `unified_search.py` | The recommended entry point. Three SAT-free stages plus an optional fourth (`--sat`): (i) two-group circuits, (ii) output blocks sharing one check register (2-adic linear solve, the synthillation family), (iii) CNOT relabeling of the output block, which adds output wires at the same `N`. Every circuit is re-verified by brute-force phase simulation, classified up to Clifford equivalence, and given an exact `Z(pi/2^l)` count and a Quirk link. Writes `../outputs/all_factories.csv` (every distinct factory) and `../outputs/smallest_per_class.csv` (smallest per output class), each row with its search stage (`two-group`, `gluing`, `padding`, `SAT`) and Quirk link, plus `search_timing.json`, for whatever range is passed on the command line (`--out-dir` to write elsewhere). The source of the paper's Tables SII(a-c). Output qubits are labeled from 1 (`T_{1}T_{2}`), as in the paper. |
| **Two-group** (App. S5) | `Two_group.py`, `two_group.ipynb` | Stage (i). Partitions the `N` circuit qubits into `k` output qubits and `N-k` check qubits, with gate weights compatible with the output/check symmetry (three skip parameters `s_total`, `s_O`, `s_S=1`). Recovers all 21,920 valid `[[n,k,2]]` (parameter-tuple, sign) configurations for `l∈{2,3,4}`, `k≤7`, `N≤11`, in well under a second. |
| **Block gluing** (App. S5) | `symfree_search.py` | Stage (ii). Solves the borrowed-identity condition directly for a chosen output block-size partition. Reaches multi-block and mixed-output states, one `CS` block and one `CCZ` block from the same factory, that the symmetric ansatz cannot produce. At `l=4` it charges `pi/8` rotations as inputs (they are not Clifford at that level). |
| **Sequential / malleable** (App. S6) | `sequential.ipynb` | Builds a circuit one step at a time, as in the paper's Definition: step `j` anchors qubit `q_j` and applies every gate on `q_j` plus a `(w-1)`-subset of the later qubits, for `w` in `W_j` (the anchor counts toward the weight). Validity is checked by brute-force phase simulation. This is the construction behind the paper's malleable circuits. |
| Exact minimality (App. S7) | `coset_certify.py` | Not a search: the exact smallest `n` for a target on `N` circuit qubits at `l=3`, as the minimum weight of a coset of punctured `RM(N-4,N)`, with no symmetry assumed. Exact while the coset dimension is small (`N<=6`); `capped_exists` decides "is there a factory with `n<n*`" at `N=7`. |
| Capped exact prover | `capped_certify.py`, `csrc/capped_minweight.c` | Proves "no factory with `n<n*` on `N` qubits" by a Brouwer–Zimmermann decision over disjoint information sets, with a compiled inner loop (vectors up to 256 coordinates). Source of the `N=7` minimality results; build with `cc -O3 -o searches/csrc/capped_minweight searches/csrc/capped_minweight.c -lpthread`. |
| SAT certifier | `sat_certify.py` | z3 version of the same question over all parity rotations, with lex-leader symmetry breaking on check qubits. Used by the optional `--sat` stage of `unified_search.py`. |
| Synthillation comparison | `synthillation_compare.py` | Builds the synthillation template for every `l=3` target with `k<=6` and checks it with the borrowed-identity condition. Source of the synthillation comparison in Fig. 1 and `outputs/landscape_l3.json`. |
| Large-`k` sweep | `large_k_sweep.py` | Block gluing for every `l=3` target on `k=7,8` outputs (about 20 min). Writes `outputs/l3_k78_targets.{csv,json}`. |
| Synthillation for any target | `synthillation_any_target.py` | The same synthillation build-and-check as `synthillation_compare.py`, for an arbitrary target (used for `k=7,8`). |
| Large sweep | `large_sweep.py` | Two-group stage to large `N` (default `N<=14`): one verified factory per `(l, k, n)`, resumable; writes `outputs/large_sweep.csv` (Quirk links for `N<=14`). Shows how far the closed form reaches. |
| Figures | `plot_new_factories.py`, `plot_literature.py`, `plot_stages.py`, `plot_landscape.py`, `plot_supp_figures.py` | Fig. 1, Fig. E1, Fig. S1, `figures/landscape_l3.png`, and the two-group background points, from the CSV/JSON files in `outputs/`. |
| Landscape figure | `plot_landscape.py` | Fig. E1, from `outputs/landscape_l3.json` and `outputs/landscape_l4.json`. |
| Figures | `plot_supp_figures.py` | Regenerates Fig. 1, Fig. S1 and Fig. S2 from the notebook searches. Two-group distances are computed by enumerating fault sets (`Two_group.factory_distance`). |
| Closed-form family (Theorem 3) | `distance_check.py` | Not a search: verifies the true GF(2) distance of the analytic `s=1`/`s=2` symmetric family (`[[2^(l+1)-2,1,2]]`, `[[2^(l+1)-1,1,3]]`) by brute force against the claimed values. |

## Running them

```
python3 searches/unified_search.py   # writes ../outputs/{all_factories,smallest_per_class}.csv + search_timing.json; ~40 min at the full range
python3 searches/Two_group.py        # library module; run via two_group.ipynb, or import directly
jupyter notebook searches/two_group.ipynb
jupyter notebook searches/sequential.ipynb   # ~30 s; writes ../figures/sequential.png
python3 searches/symfree_search.py   # writes ../figures/block_search_levels.png, reads ../outputs/factory_catalogue_l{2,3,4}.csv
python3 searches/distance_check.py
```

To widen any search past its current swept range, pass `--kmax`, `--nmax`, `--rmax` to
`unified_search.py`, or edit the loop bounds (`L_RANGE`, `n_max`, `k_max`) in the others. The
condition being checked doesn't change, only how far it's swept. Most of the unified search's run
time at `l<=3` is the output classifier's search over CNOT frames, not the search itself (seconds for
two-group, minutes for gluing and padding at `l=4`);
`search_timing.json` records the two separately.

## What the unified search does not reach

At distance 2 every gate touching a check qubit is detected, so a factory is fixed by what it
deposits on the output block. The three stages above cover every output expressible as a block
partition after a CNOT relabeling. Outputs that are genuinely entangled beyond that need a SAT
search; at `l=3` the known ones are `[[16,4,2]]`, `[[16,5,2]]` (CS-type) and `[[12,6,2]]`
(CCZ-type). The exact `Z(pi/2^l)` count is not computed at `k=6` for `l=2,4` (or `k=7` at `l=3`).

## Sequential search: no new factories, but the malleability engine

Running `sequential.ipynb` over `l∈{2,3,4}`, `N≤8`, skip parameters `s_j∈{1,2,3,4}` per step and
both signs, it recovers the paper's malleable chain `[[8,3,2]]→CCZ`, `[[12,2,2]]→CS`,
`[[14,1,2]]→T` at `l=3`, `N=4`, and its analogues at `l=2` and `l=4`. Every factory it finds is
also found, at equal or smaller `n`, by the unified search, and sign flips never produce new
parameters. So it is not a source of records; its role is the malleability it makes explicit, one
parent circuit yielding several output types. The catalytic endpoints `[[11,1,2]]→T` and
`[[10,2,2]]→T,T` need an extra Hadamard step the notebook does not implement.

An earlier version of this notebook omitted the anchor from each gate's weight, so it never
generated the gate on the anchor alone and could not reproduce the chain above; that is fixed.

**Caveat if you extend this search:** a term with coefficient nonzero mod `2^l` can still be
Clifford if that coefficient is even mod `2^size`. Check genuineness with
`classification/classify.py`, not a bare "coefficient ≠ 0" test (see that folder's README for why).
