# Magic State Factory Search

Code and data for:

> S. Singh, C. Gidney, C. Jones, **"Borrowed Identities: Malleable Distillation Factories and a
> Unified Numerical Search,"** [arXiv:2606.28518](https://arxiv.org/abs/2606.28518) (2026).

The paper introduces a *borrowed-identity* condition for magic-state distillation factories: a
circuit only has to act as the identity on **one** input state, not on an entire codespace. This
is strictly weaker than the usual transversal-gate requirement, holds at every level of the
Clifford hierarchy, and turns factory-finding into a search problem that this repository
implements and sweeps.

## The Clifford hierarchy, and which state lives where

For `θ = π/2^l` (`l ≥ 1`), the **level-`l` magic state** is the (unnormalized) state
`|θ⟩ := |0⟩ + e^(2iθ)|1⟩`, i.e. the +1 eigenstate of the diagonal phase gate
`Z(2θ) = diag(1, e^(2iθ))`. More generally, the weight-`w` rotation

```
Z^⊗w(θ) = exp( iθ · (I − ⊗_{j=1}^w Z_j) )
```

applies phase `e^(2iθ)` to every basis state with an odd number of 1s among `w` designated
qubits. These gates are exactly the diagonal elements of the **`l`-th level of the Clifford
hierarchy**, `C_l` — the hierarchy defined recursively by `C_1` = Pauli group, `C_2` = Clifford
group (the normalizer of `C_1`), and `C_l` the normalizer of `C_{l-1}`. Concretely:

| `l` | `θ = π/2^l` | weight-1 gate | weight-2 gate | weight-3 gate | weight-4 gate |
|---|---|---|---|---|---|
| 2 | π/4  | `S`  | `CZ`  | —     | —      |
| 3 | π/8  | `T`  | `CS`  | `CCZ` | —      |
| 4 | π/16 | `√T` | `CT`  | `CCS` | `CCCZ` |

**`l = 3` is where the `T` state lives** (`e^(iπ/4)` is the non-Clifford phase of `T`); `l = 2` is
still inside the Clifford group (`S`, `CZ`); `l ≥ 4` are higher, more expensive levels. A
distillation factory consumes `N` noisy copies of a level-`l` magic state and Clifford operations
only, and outputs `k` higher-fidelity copies at distance `d` — this repo searches for `[[N,k,d]]`
factories at `l = 2, 3, 4`, and the framework is written to extend to any `l`.

## The three searches

This repository implements three complementary numerical searches for borrowed-identity
factories, corresponding to Sections III/IV and Appendix H of the paper. (There is also a
**closed-form family**, Theorem 3 of the paper, obtained analytically for fully symmetric
circuits — no search needed; [`distance_check.py`](distance_check.py) verifies its true GF(2)
distance against the claimed `[[2^(l+1)-2,1,2]]` / `[[2^(l+1)-1,1,3]]` values.)

| Search | Files | Purpose |
|---|---|---|
| **Two-group** (Sec. IV) | [`Two_group.py`](Two_group.py), [`two_group.ipynb`](two_group.ipynb) | The main systematic search. Partitions the `n` circuit qubits into `k` output qubits and `n-k` check qubits and only allows gate weights compatible with the output/check symmetry (parametrized by three "skip" values `s_total`, `s_O`, `s_S`). This reduces the borrowed-identity condition to a linear check mod `2^(l+1)`, so the sweep is fast: it recovers all **21,920** valid `[[N,k,2]]` factories for `l ∈ {2,3,4}`, `k ≤ 7`, `n ≤ 11`, including entangled- and multi-output families no single earlier framework reached at once. |
| **Symmetry-free / targeted-output** (App. H) | [`symfree_search.py`](symfree_search.py) | An existence-proof search that drops the two-group symmetry ansatz and solves the borrowed-identity condition directly (via a linear solve mod `2^l`, minimized over the kernel) for a *chosen* output partition. It reaches factories the symmetric ansatz structurally cannot, e.g. multi-block and **mixed-output** states (one part of the output is a `CS` state, another a `CCZ` state, from the *same* factory). |
| **Sequential / malleable** (Sec. III) | [`sequential.ipynb`](sequential.ipynb) | Builds a circuit up **one step at a time**: each step adds an anchor qubit plus a family of higher-weight gates, and checks the disentangling condition (Theorem 2) is independent of the earlier check-qubit's value at every step. This is the constructive procedure behind the paper's *malleable circuits*: because the check condition holds at every intermediate step, a single parent circuit can be stopped early (or have gates removed) to yield a different `(N,k)` endpoint with a different output type — the output magic-state type becomes a **compile-time choice**, not something fixed by the circuit's design. |

### Sequential search: no new factories, but it *is* the malleability engine

Running the (bug-fixed — see below) sequential search exhaustively over `l ∈ {2,3,4}`, `n ≤ 8`,
skip parameters `s_j ∈ {1,2,3,4}` per step, and both sign choices, and checking every result's
*genuine* output degree (the highest-weight output term whose coefficient survives **odd** mod
`2^l` — an even coefficient means that term is secretly Clifford, not real magic; see the caveat
below) gives a blunt answer: **every single non-Clifford result it finds is already found by the
two-group search, at the same `N`.** Concretely, its only genuine hits are the `k=2`, degree-1
family `[[6,2,2]]` (`l=2`, the smallest H-code), `[[14,2,2]]` (`l=3`, the smallest Bravyi–Haah
`k=2` triorthogonal code), and `[[30,2,2]]` (`l=4`, the same family's next level) — all three
already independently recovered by the two-group search. Zero new `(N,k)` pairs, zero new output
types.

That is expected, not a bug: the sequential ansatz is a strict special case of what the two-group
/ symmetry-free searches already cover (a single "anchor + higher-weight block" per step), so it
cannot beat them on record-setting `(N,k,d)`. Its value is different — it is the actual
*construction* the paper uses to demonstrate malleability (Sec. III.c, "Malleable circuits"): take
one parent circuit and show that different early-termination points, or the addition of a
catalytic conversion step, give genuinely distinct distillation factories from the *same* circuit.
The paper works two hand-built examples this way (App. G): a 4-qubit parent at `l=3` whose
successive steps/conversions yield `[[8,3,2]] → CCZ`, `[[12,2,2]] → CS`, `[[14,1,2]] → T`, and (via
an added catalytic step, outside what this notebook's plain sequential condition implements)
`[[11,1,2]] → T` and `[[10,2,2]] → T,T` — the last two use non-diagonal (Hadamard) gates and are
the only genuinely *new* factories the sequential/malleable idea produces in the paper;
`[[11,1,2]]` is introduced there for the first time. Running `sequential.ipynb` as coded does not
reach these (it only implements the diagonal, no-catalysis condition), but it is the tool that
verifies the underlying step-by-step chain is valid at every intermediate point — which is exactly
the property malleability requires.

**Caveat (important if you extend the sequential search):** a term with a nonzero-mod-`2^l`
coefficient can still be a *Clifford* term if that coefficient is even — it just means the "magic"
at that weight is secretly a lower-level (cheaper) rotation. A naive "coefficient != 0" degree
check will misreport such states as having real output degree/entanglement when there is none.
The genuine-degree check above (require an *odd* surviving coefficient) is what the numbers in
this README use; a looser check (nonzero-mod-`2^l` instead of odd) inflates the sequential
search's apparent non-Clifford hit count from 36 to 128 (out of 180 raw results) — all 36 genuine
hits reduce to just the 3 ties above.

## Catalogue: every distance-2 factory found, by level

The factories themselves live in three CSV files, not in this README:
[`factory_catalogue_l2.csv`](factory_catalogue_l2.csv),
[`factory_catalogue_l3.csv`](factory_catalogue_l3.csv),
[`factory_catalogue_l4.csv`](factory_catalogue_l4.csv) -- one file per Clifford-hierarchy level,
each row a distinct `(N, k, decomposition)` recovered by the two-group search, the symmetry-free
search, or both. Columns:

| column | meaning |
|---|---|
| `N`, `k`, `d` | the factory: `N` input magic states, `k` output qubits, distance `d` (always 2 here) |
| `degree` | the Clifford-reduced genuine degree of the deposited state (minimized over output CNOT frames -- see [Classification](#classification-exact-t-count-and-clifford-reduced-degree) below) |
| `essential_dim` | how many of the `k` output qubits actually carry genuine (non-Clifford) content in the best frame found; `essential_dim < k` means the factory is *padded* -- some output qubits are free |
| `degree_note` | set only when `essential_dim`/`degree` is a documented upper bound from a random-frame search rather than an exact `GL(k,2)` enumeration (`k >= 5` here) |
| `t_count` | the **exact** minimal T-count of the deposited state (Reed-Muller minimum-weight-coset decoder); only defined at `l = 3` |
| `decomposition` | the deposited state's named gate content per separable component, e.g. `CS+CCZ` for a mixed-output factory |
| `n` | total circuit qubits (`k` outputs + `n-k` checks) |
| `s_total`, `s_O` | **two-group only**: the skip parameters passed to `Two_group.py`'s `build_gate_set(k, n-k, s_total, s_O)` that produced this row (`s_S` is fixed at 1) |
| `parts` | **symmetry-free only**: the output block-size partition passed to `symfree_search.py`'s `solve_binding(parts, ncheck, l)` that produced this row |
| `search` | which search(es) found this exact `(N, k, decomposition)` |

So every row carries back its own recipe: for a two-group row, `(l, n, k, s_total, s_O)` regenerates
it via `Two_group.py`; for a symmetry-free row, `(l, parts, n-k)` regenerates it via
`symfree_search.py`. See [Pipeline](#pipeline-regenerating-the-catalogues-and-explicit-circuits)
below for the exact commands.

There is also a separate, distance-**3** closed-form family (Theorem 3, `s=1` or `s=2` symmetric
solutions) at every level, not in the CSVs: the quantum Reed-Muller code `[[2^(l+1)-1,1,3]]`, which
specializes to the **Steane code** `[[7,1,3]]` at `l=2` and the classic **15->1 T-distillation**
factory `[[15,1,3]]` (Bravyi-Kitaev) at `l=3` -- verified by [`distance_check.py`](distance_check.py).

## Classification: exact T-count and Clifford-reduced degree

`classify.py` is the shared classifier both catalogue builders call. It is adapted from Shraddha
Singh's own classification code in her `sj-magic-state-factory-searches` repository (used here with
her permission), combining two independent pieces:

- **Exact minimal T-count** ([`tcount.py`](tcount.py), copied verbatim): the Amy-Mosca /
  Reed-Muller minimum-weight-coset decoder. Given the deposited phase polynomial, it finds the true
  minimum number of T-gates needed, correctly identifying that different-looking gate strings can
  be the same physical resource (e.g. `T1 T2 CS12` is Clifford-equivalent to a single `T`, not 3
  separate gates). Exact and self-tested; only defined at `l = 3`.
- **Clifford-reduced degree and essential dimension** (adapted from her `audit_flag2.py`): minimizes
  `(essential_dim, degree)` over every output `CNOT` frame (`GL(k,2)`) -- a monomial of size `r` is
  *genuine* only if its coefficient is nonzero mod `2^r`, not merely nonzero, and not merely odd
  (both weaker tests give wrong answers -- the odd test alone misclassifies the real `[[4,2,2]]`
  Iceberg code as Clifford, since its genuine `CZ` coefficient is 2). Exact for `k <= 4`
  (`|GL(4,2)| = 20,160`); for `k >= 5` it is a documented upper bound from up to 60,000 random
  frames (flagged in `degree_note`, and it can only ever *overstate* degree, never understate it).
  This generalizes her `l=3`-only script to any level `l` used here.

`decompose()` then reads the best frame's genuine monomials as connected components (qubits tied
together by a shared genuine monomial are one component) and names each component by its weight --
this is what produces mixed labels like `CS+CCZ` for a factory whose output splits into independent
pieces.

## Pipeline: regenerating the catalogues and explicit circuits

1. **Sweep + classify.** The catalogue-build scripts sweep `Two_group.py`'s
   `(l, n, k, s_total, s_O)` and `symfree_search.py`'s `(l, parts, ncheck)` parameter spaces,
   reconstruct each valid factory's actual gate list, call `classify.classify(Gf, cf, O, l)` from
   `classify.py`, and keep the distinct `(N, k, decomposition)` rows -- this is exactly what produced
   `factory_catalogue_l{2,3,4}.csv` above. There is no single canned "rebuild" script yet; the sweep
   loops are the same ones in `Two_group.py`'s `find_valid_signs`/`build_gate_set` and
   `symfree_search.py`'s `search`/`solve_binding` -- wrap them with a call to `classify.classify` in
   place of (or alongside) each script's own looser degree check to reproduce a catalogue row.
2. **Get one factory's explicit circuit.** [`export_circuit.py`](export_circuit.py) takes a
   catalogue row's own parameters and reconstructs the circuit as an explicit binary matrix (rows =
   the `n` wires, columns = the `N` gates -- the same column convention as Shraddha Singh's
   `master_catalog` catalogue: one qubit-support per `pi/2^l` parity rotation):
   ```
   python3 export_circuit.py two-group --l 3 --n 4 --k 2 --s_total 1 --s_O 1   # [[12,2,2]] CS
   python3 export_circuit.py symfree   --l 3 --parts 3,2 --checks 2           # [[18,5,2]] CS+CCZ
   ```
   Each prints the gate list, the 0/1 matrix, and the row's classification (degree, essential
   dimension, decomposition, exact T-count) as an independent check that the reconstructed circuit
   matches the catalogue entry it came from.

## An active, extensible tool

This repository is meant to keep being used, not just to reproduce the paper's figures. The
borrowed-identity condition is written generically in `l` (see `is_sequential_valid`,
`is_borrowed_identity`, and the `f_coefficient`/`sigma` helpers behind Theorem 5), so finding a
**`Z(π/2^l)`-to-any-`D_l`** magic-state factory at a level or parameter range not yet swept is a
matter of widening `L_RANGE`/`n_max`/`k_max` in the relevant script, not writing new theory. The
`k ≤ 7`, `n ≤ 11` ceiling in the two-group search and the `n ≤ 8` ceiling in the sequential search
are properties of the swept range, not the framework — both scale polynomially in `n` and `k`, so
a distance-2 factory at larger `k` is reachable by widening the sweep.

## What isn't found here

A prior audit of this repository (`SEQUENTIAL_SEARCH_AUDIT.md`, checked against a SAT/exhaustive
classification) flagged five specific `l=3` distance-2 factories as reachable only by SAT search,
not by any construction here. Re-checked against the rebuilt catalogues above (`degree`,
`essential_dim`, `decomposition` from `classify.py`, matching that audit's own classifier):

| factory | signature (essential_dim, degree, decomposition) | found here? |
|---|---|---|
| `[[12,3,2]]a` | 3, 2, `CS` | **found** — two-group, `l=3 n=5 k=3 s_total=2 s_O=1` (`factory_catalogue_l3.csv`) |
| `[[16,4,2]]` | 4, 2, `CS` | not found |
| `[[16,5,2]]` | 5, 2, `CS` | not found |
| `[[12,6,2]]` | 5, 3, `CCZ` | not found |
| `[[16,6,2]]` | 6, 3, `CCZ` | not found |

So one of the five turns out to be a false negative from the earlier, looser classification, not a
genuine gap: `[[12,3,2]]a` is in the two-group catalogue once classified correctly. The other four
are genuinely absent from both searches here and, as far as I've checked, require either a SAT
solver or an exhaustive enumeration to reach.
