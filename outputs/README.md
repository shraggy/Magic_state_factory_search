# Outputs

## Catalogues

`factory_catalogue_l2.csv`, `factory_catalogue_l3.csv`, `factory_catalogue_l4.csv` — one file per
Clifford-hierarchy level, each row a distinct `(N, k, decomposition)` recovered by the two-group
search, the symmetry-free search, or both, classified by `classification/classify.py`.

| column | meaning |
|---|---|
| `N`, `k`, `d` | the factory: `N` input magic states, `k` output qubits, distance `d` (always 2 here) |
| `degree` | Clifford-reduced genuine degree of the deposited state, minimized over output CNOT frames |
| `essential_dim` | how many of the `k` output qubits actually carry genuine content in the best frame; `essential_dim < k` means the factory is *padded* |
| `degree_note` | set only when `essential_dim`/`degree` is an upper bound (random-frame search, `k≥5`) rather than exact |
| `t_count` *(l=3 only)* | exact minimal T-count (Reed–Muller coset decoder); not a column at all in the `l=2`/`l=4` files, since it's never defined there |
| `t_count_note` *(l=3 only)* | why `t_count` is blank for this row, when it is |
| `decomposition` | qubit-indexed gate content, e.g. `CS01+CS02` (one entangled 3-qubit block) vs `CS01+CS23` (two independent 2-qubit blocks) — see `classification/README.md` |
| `n` | total circuit qubits (`k` outputs + `n-k` checks) |
| `s_total`, `s_O` | **two-group rows only** — the parameters that regenerate this row: `Two_group.py`'s `build_gate_set(k, n-k, s_total, s_O)` (`s_S` is fixed at 1) |
| `parts` | **symmetry-free rows only** — the output block-size partition that regenerates this row: `symfree_search.py`'s `solve_binding(parts, n-k, l)` |
| `search` | which search(es) found this exact `(N, k, decomposition)` |

**Every row carries its own recipe.** For a two-group row, `(l, n, k, s_total, s_O)` regenerates it;
for a symmetry-free row, `(l, parts, n-k)` does. Concretely:

```
# a two-group row: l=3, n=4, k=2, s_total=1, s_O=1  ([[12,2,2]] CS)
python3 classification/export_circuit.py two-group --l 3 --n 4 --k 2 --s_total 1 --s_O 1

# a symmetry-free row: l=3, parts=3+2, n-k=2  ([[18,5,2]] CS+CCZ)
python3 classification/export_circuit.py symfree --l 3 --parts 3,2 --checks 2
```

`parts` is specifically the list of output-block sizes the symmetry-free search was asked to
realize *before* it solves for the checks — e.g. `parts=3,2` means "find a factory whose 5 outputs
split into one weight-≤3 block and one weight-≤2 block" (`symfree_search.py`'s `solve_binding`
argument, from `_partitions(k, L)`). It has no two-group equivalent because two-group's ansatz is
always fully symmetric across all `k` outputs at once — there's no block structure to choose.

There is also a separate, distance-**3** closed-form family (Theorem 3) at every level, not in
these CSVs: the quantum Reed–Muller code `[[2^(l+1)-1,1,3]]` (the Steane code `[[7,1,3]]` at `l=2`,
the classic `15→1` factory `[[15,1,3]]` at `l=3`) — verified by `searches/distance_check.py`.

## Raw sweep logs

`output_sequential.txt`, `output_two_group.txt` — captured stdout from running the sequential and
two-group searches directly (not filtered/classified the way the CSVs are); useful for seeing the
full parameter sweep, including the two-group log's own "recovery check" against named literature
factories.

## Regenerating

There's no single canned rebuild script yet. The build loops are the same ones in
`Two_group.py`'s `find_valid_signs`/`build_gate_set` and `symfree_search.py`'s `search`/
`solve_binding` — reconstruct each valid factory's gate list, call
`classification.classify.classify(Gf, cf, O, l)`, and write out the distinct `(N, k, decomposition)`
rows.

## Distance-2 factories a SAT/exhaustive search finds that these two don't

A prior audit of this repository, cross-checked against a SAT/exhaustive classification, flagged
five specific `l=3` distance-2 factories as reachable only by SAT search. Re-checked here with
`classify.py` (same `essential_dim`/`degree`/`decomposition` signature that audit used):

`[[12,3,2]]a` (essential_dim 3, degree 2, `CS01+CS02`) turns out to be a **false negative**, not a
gap — two-group finds it exactly (`l=3 n=5 k=3 s_total=2 s_O=1`, `t_count=4`). It isn't a very good
factory either way: `[[12,2,2]]`, which we also find, gets the same degree-2 content at the same
`N=12` with one fewer output qubit and one fewer T-state (`t_count=3`) — `[[12,3,2]]a`'s third
qubit is entangled but doesn't add anything `[[12,2,2]]` doesn't already provide.

`[[16,6,2]]` (essential_dim 6, degree 3, a single connected 6-qubit `CCZ`-type block) isn't worth
chasing: the symmetry-free search already gets the same degree-3 content on 6 outputs more cheaply
(`N=14`) by depositing it as two *independent* `CCZ` blocks instead of one entangled one. `16 > 14`,
so even reaching the entangled version wouldn't be an improvement.

The other three are genuine gaps — real structure neither search here reaches:
- `[[16,4,2]]` (essential_dim 4, degree 2, `CS`) beats this repo's best single-block `CS` factory
  at `k=4` (`N=18`) by 2.
- `[[16,5,2]]` (essential_dim 5, degree 2, `CS`) beats this repo's best single-block `CS` factory
  at `k=5` (`N=105`) by nearly 7×.
- `[[12,6,2]]` (essential_dim 5, degree 3, `CCZ`, padded to `k=6`) beats this repo's best
  single-block `CCZ` factory at `k=6` (`N=124`) by over 10×, and even beats the `2×CCZ` direct-sum
  factory above (`N=14`) by 2.

Closing that gap needs a genuinely different search, not more sweeping of the ranges already
covered here.
