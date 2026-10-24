# Outputs

## The catalogues

Two files hold every factory, each row with a [Quirk](https://algassert.com/quirk) link. Every link is drawn
the same way: all qubits start in `|0>` and get `H`; each factory gate is one multi-qubit `Z(pi/2^l)` parity
rotation; the Clifford corrections return the checks to `|+>` (first row of Bloch spheres); then the
output decomposition is undone and a final `H` returns every qubit to `|0>` (second row of Bloch spheres).

| file | rows | what it lists |
|---|---|---|
| `smallest_per_class.csv` | 186 | the **smallest factory for each output class** (Clifford-reduced decomposition, per `k`), at `l=2,3,4`. These are the factories quoted in the paper; the paper's Tables SII(a-c) give representatives. |
| `all_factories.csv` | 1,564 | **every distinct factory** the unified search finds (all `n`), at `l=2,3,4`. |

`factory_catalogue_l{2,3,4}.csv` are `smallest_per_class.csv` split by level, with plain-text decompositions
and literature references.

### Clifford-equivalent factories are listed separately

In `all_factories.csv`, factories whose outputs are equivalent up to a Clifford circuit share the same
`clifford_class` label (e.g. `L3-007`), and the file is sorted so that each class is contiguous. They are
still **listed as separate rows**, because they are different factories: different `n`, `k` or circuit
size `N`, and different output gates. Two columns tell them apart:

- `raw` is the gate the factory actually deposits on its outputs; `reduced` is its Clifford-reduced class.
- `clifford_padded = yes` means the output is Clifford-equivalent to a state on fewer than `k` qubits
  (it delivers the same magic as a smaller factory, on more wires).

A padded factory adds no distillation power: the smaller factory followed by a CNOT circuit does the same.
It is worth having where that CNOT circuit is expensive, for instance in architectures where logical
Clifford gates are synthesized from multi-step Pauli measurements (such as bivariate-bicycle codes). There
the padded factory delivers the target state directly, spending only magic-state inputs it needs anyway. For
that reason no Clifford-equivalent factory is dropped from `all_factories.csv`.

Every row is distance `d=2` (the distance-3 Reed-Muller factories are not swept here), built only from
`Z(pi/2^l)` parity rotations, and re-verified by brute-force simulation of its phase table.

### How the files were generated (and how to regenerate them for any range)

One command runs the search and writes **both** files, each row with its search stage and Quirk link, plus
`search_timing.json` with the wall-clock time of every stage:

```
python3 searches/unified_search.py --levels 2 3 4 --kmax 6 --nmax 11 --rmax 4 --relabel 4 --sat --sat-time 60
```

Change any parameter and the same two CSVs are produced for that range (use `--out-dir` to keep them apart
from the published ones). `--levels` picks `l`; `--kmax` the number of outputs; `--nmax` the circuit size in
the two-group stage; `--rmax` the check qubits in the gluing stage; `--relabel` the extra output wires
added by Clifford padding; `--sat` adds the optional SAT stage.

The stages: (i) two-group circuits (`Two_group.py`, `s_total, s_O <= 7`; all signs `+1` at `l<=3`, every
check-only sign pattern at `l=4`); (ii) block gluing, output blocks sharing one check register
(`symfree_search.py`'s 2-adic linear solve); (iii) Clifford padding, CNOT relabeling of the output block at
the same `N`; (iv) optional SAT. Classifying each output (Clifford class, distance, z-count) is bookkeeping,
timed separately in `search_timing.json`, and is not needed to find a factory.

## Columns

| column | meaning |
|---|---|
| `clifford_class` | output class label, e.g. `L3-007`; rows sharing it have Clifford-equivalent outputs |
| `l`, `k`, `n`, `N`, `d` | level, output qubits, input magic states, circuit qubits, distance |
| `search_stage` | the stage that found the factory: `two-group`, `gluing`, `padding` or `SAT` |
| `padded_from` | for `padding` rows, the stage of the circuit that was relabeled (`two-group` or `gluing`) |
| `params` | stage parameters (`s_total,s_O` for two-group; block partition for gluing; target for SAT) |
| `reduced` | the output up to Clifford equivalence, as a product of generators |
| `raw` | the gate the circuit actually deposits on its outputs |
| `zcount` | exact minimal number of `Z(pi/2^l)` rotations to prepare the output unitarily (blank where not computed) |
| `essential_dim` | output qubits carrying genuine content after Clifford reduction |
| `clifford_padded` | `yes` when `essential_dim < k` |
| `masks` | for `padding` rows, the CNOT relabeling (one bit-mask per new output wire) |
| `quirk` | Quirk link of the circuit, drawn as described above |
