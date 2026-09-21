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

The tables below give, for each `(l, k, output degree)` found by *any* of the three searches, the
smallest `N` (i.e. the best known factory of that shape) among the three, which search found it,
and — where I could verify it — its status in the prior literature. All entries below are
distance `d = 2`.

Legend for the "output" column: `S`/`T`/`√T` = weight-1 (single-qubit) gate at that level;
`CZ`/`CS`/`CT` = weight-2; `CCZ`/`CCS` = weight-3; `CCCZ` = weight-4 (see hierarchy table above).

### l = 2 (θ = π/4, S/CZ level)

| `[[N,k,d]]` | output | search | known in literature |
|---|---|---|---|
| `[[6,1,2]]`  | `S`  | two-group | this work (all-weights construction, Sec. II) |
| `[[6,2,2]]`  | `S`  | two-group *(tied by sequential)* | H-code, Jones 2013 |
| `[[4,2,2]]`  | `CZ` | two-group | Iceberg code |
| `[[8,3,2]]`  | `S`  | symmetry-free (3×`S`) | — |
| `[[4,3,2]]`  | `CZ` | two-group | — |
| `[[8,4,2]]`  | `S`  | two-group | — |
| `[[6,4,2]]`  | `CZ` | two-group | — |
| `[[10,5,2]]` | `S`  | symmetry-free (5×`S`) | — |
| `[[10,5,2]]` | `CZ` | symmetry-free (2×`CZ`+`S`) | — |
| `[[10,6,2]]` | `S`  | two-group | — |
| `[[8,6,2]]`  | `CZ` | two-group | — |
| `[[24,7,2]]` | `S`  | two-group | — |
| `[[8,7,2]]`  | `CZ` | two-group | — |

### l = 3 (θ = π/8, T/CS/CCZ level — where T lives)

| `[[N,k,d]]` | output | search | known in literature |
|---|---|---|---|
| `[[14,1,2]]`  | `T`   | two-group | this work (all-weights construction, Sec. II) |
| `[[14,2,2]]`  | `T`   | two-group *(tied by sequential)* | Bravyi–Haah 2012, `k=2` triorthogonal instance |
| `[[12,2,2]]`  | `CS`  | two-group | this work (Sec. II derivative); independently matches the `CS(4)→CS(1)` AG-code factory (A. Gong, QEC 2026 oral presentation) |
| `[[8,3,2]]`   | `CCZ` | two-group | "the cube": Bombín–Martín-Delgado 2007 3D color code |
| `[[20,3,2]]`  | `T`   | symmetry-free (3×`T`) | — |
| `[[18,4,2]]`  | `CS`  | symmetry-free (2×`CS`) | **newly found** — this paper's `T`-to-`CS` synthillation family, `[[6m+6,2m,2]]` |
| `[[8,4,2]]`   | `CCZ` | two-group | — |
| `[[26,5,2]]`  | `T`   | symmetry-free (5×`T`) | — |
| `[[26,5,2]]`  | `CS`  | symmetry-free (2×`CS`+`T`) | — |
| `[[18,5,2]]`  | `CCZ` | symmetry-free (mixed `CCZ`+`CS`) | **newly found** — one `CS` and one `CCZ` output from one factory; not in the synthillation literature |
| `[[26,6,2]]`  | `T`   | two-group | Bravyi–Haah 2012, `k=6` instance |
| `[[24,6,2]]`  | `CS`  | symmetry-free (3×`CS`) | — |
| `[[14,6,2]]`  | `CCZ` | symmetry-free (2×`CCZ`) | Campbell–Howard `T`-to-`CCZ` synthillation family, `[[6m+2,3m,2]]`, `m=2` |
| `[[56,7,2]]`  | `T`   | two-group | — |
| `[[108,7,2]]` | `CS`  | two-group | — |
| `[[132,7,2]]` | `CCZ` | two-group | — |

### l = 4 (θ = π/16, √T/CT/CCS/CCCZ level)

| `[[N,k,d]]` | output | search | known in literature |
|---|---|---|---|
| `[[30,1,2]]`  | `√T`   | two-group | this work (all-weights construction, Sec. II) |
| `[[30,2,2]]`  | `√T`   | two-group *(tied by sequential)* | the `l=4` instance of the same closed-form `[[2^(l+1)-2,2,2]]` family that gives the H-code (`l=2`) and Bravyi–Haah's `k=2` code (`l=3`) — not independently named at `l=4` |
| `[[28,2,2]]`  | `CT`   | two-group | — |
| `[[44,3,2]]`  | `√T`   | symmetry-free (3×`√T`) | — |
| `[[44,3,2]]`  | `CT`   | symmetry-free (`CT`+`√T`) | — |
| `[[24,3,2]]`  | `CCS`  | two-group | — |
| `[[44,4,2]]`  | `√T`   | two-group | — |
| `[[42,4,2]]`  | `CT`   | symmetry-free (2×`CT`) | — |
| `[[24,4,2]]`  | `CCS`  | two-group | — |
| `[[16,4,2]]`  | `CCCZ` | two-group | — |
| `[[66,5,2]]`  | `√T`   | symmetry-free (5×`√T`) | — |
| `[[60,5,2]]`  | `CT`   | symmetry-free (2×`CT`+`√T`) | — |
| `[[46,5,2]]`  | `CCS`  | symmetry-free (`CCS`+`CT`) | — |
| `[[16,5,2]]`  | `CCCZ` | two-group | — |
| `[[58,6,2]]`  | `√T`   | two-group | — |
| `[[58,6,2]]`  | `CT`   | symmetry-free (3×`CT`) | — |
| `[[42,6,2]]`  | `CCS`  | symmetry-free (2×`CCS`) | — |
| `[[42,6,2]]`  | `CCCZ` | symmetry-free (`CCCZ`+`CT`) | — |
| `[[120,7,2]]` | `√T`   | two-group | — |
| `[[308,7,2]]` | `CCS`  | two-group | — |
| `[[344,7,2]]` | `CCCZ` | two-group | — |

Rows marked "this work" are recovered by the searches here but aren't (as far as I've checked)
independently named in prior literature; rows with no entry in the last column simply haven't
been cross-referenced yet — absence of a citation is not a claim of novelty. The full, unfiltered
dump (every valid parameter tuple, not just the smallest-`N` representative per output type) is in
[`factory_catalogue.csv`](factory_catalogue.csv) (two-group) and
[`output_two_group.txt`](output_two_group.txt) / [`output_sequential.txt`](output_sequential.txt)
(raw sweep logs).

There is also a separate, distance-**3** closed-form family (Theorem 3, `s=1` or `s=2` symmetric
solutions) at every level: the quantum Reed–Muller code `[[2^(l+1)-1,1,3]]`, which specializes to
the **Steane code** `[[7,1,3]]` at `l=2` and the classic **15→1 T-distillation** factory
`[[15,1,3]]` (Bravyi–Kitaev) at `l=3`.

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

The three searches above are all built around distance-2, and (with the exception of the
diagonal-only sequential construction) they only reach circuits built from parity-phase (diagonal
Clifford-hierarchy) gates. Concretely, out of scope for this repository as it stands:

- **Non-CSS / catalytic factories that need Hadamard gates.** `[[10,2,2]]` (Meier–Eastin–Knill,
  using a `[[4,2,2]]` code as an inner block) and `[[11,1,2]]` (introduced in the paper itself)
  both require a catalytic conversion step with a non-diagonal (Hadamard) gate. None of the three
  scripts here implement that extension — it's flagged in the paper as important future work
  ("Adding Hadamard gates to the extended borrowed-identity formalism ... is equally important").
- **Distance ≥ 3, multi-output factories.** This repo's searches target `d=2` (plus the
  single-output, closed-form `d=3` Reed–Muller chain above). A concurrent, independent line of
  work — H. Jacinto, X. Valcarce, V. Barizien, É. Gouzien, N. Sangouard,
  ["Exploring the landscape of compact magic-state distillation factories"](https://arxiv.org/abs/2606.07734)
  — uses a **SAT solver** to derive no-go theorems and new smallest-qubit-count protocols
  specifically for `d ≥ 3` (e.g. new `d=4` and `d=5` T-to-T / T-to-CCZ protocols on 8–11 qubits).
  That regime is outside what any of the three searches in this repo attempt.
- **The full Nezami–Haah classification.** S. Nezami and J. Haah,
  ["Classification of small triorthogonal codes,"](https://arxiv.org/abs/2107.09684) Phys. Rev. A
  106, 012437 (2022), enumerates *every* affine-equivalence class of the relevant Reed–Muller
  polynomials at `l=3` — a space of ~5×10^12 candidates for `N ≤ 38`, well beyond what a targeted
  search like the ones here would ever enumerate directly. The two-group search recovers that
  catalogue's distance-2 entries *at the largest `k` for each `N`*; entries at smaller `k` for the
  same `N` (if they exist and are useful) are not reproduced by anything in this repo, and would
  need either a similarly exhaustive enumeration or a SAT-based search to find.
