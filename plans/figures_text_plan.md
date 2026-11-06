# Plan: edits to Fig. 1, Fig. 2 and the main text (your comments since 10:32, grouped)

Status key: **[draft]** shown in a render above, not yet in the paper; **[todo]** not started.
Nothing below is in main.tex yet; the PDF still has the 09:49 figures and captions.

## A. Your comments, grouped

| # | Comment | Group |
|---|---|---|
| 1 | Does the paper rely heavily on "d=2 circuits are optimal compilations with checks"? | Text: framing |
| 2 | d=2 is the cultivation regime; only our follow-up work does d>2, AI-assisted, because a systematic search is hard | Text: Discussion |
| 3 | What does "circuit picture after Ref. 28" mean? Litinski does not draw four-qubit circuits | Fig. 1 caption |
| 4 | Show the small [[10,2,2]] snippet (two added T, one check) if there is space | Fig. 1 |
| 5 | q0–q3 start in \|+> and return to \|+> after the grey bar; that is the point of a BI | Fig. 1 |
| 6 | What is Z_S, why S? Drop the header, the caption is enough | Fig. 1 |
| 7 | Captions 4–5 sentences, never longer than the figure | Both captions |
| 8 | Fig. 2(b): "9 of 9, 766 of 766" says nothing; one line "l=2: 9, l=3: 766, …"; add l=2,3,4 to (b), grey only for l=3 synthillation | Fig. 2(b) |
| 9 | Write n ≥ T-count+1 instead of "bound", move it off the y axis | Fig. 2(a) |
| 10 | Are the crosses showing that they do not reach the bound? | Fig. 2(a) |
| 11 | k=8: "upper bound only" or lower bound? Something in contrast to "exact" | Fig. 2(a) |

## B. Fig. 1 (conceptual; single column, ≤ 3.4 × 1.85 in)

1. **[draft] \|+> in and out.** Every wire labelled \|+> on the left and \|+> on the right of the circuit (comment 5). Done in the current render.
2. **[draft] No header** (comment 6). The "Z_S(π/8)" line is removed. The figure uses no symbol the caption does not define; dots mark the qubits each π/8 rotation acts on.
3. **[draft] [[10,2,2]] branch** (comment 4; 11:0x: actual gates instead of "C", |+> instead of X, π/8 rotations
   drawn widest first in every layer). Drawn below the [[8,3,2]] cut on q1–q3: C, T†, C, T on q1, X measurement on q1, output \|T>⊗2 on q2 q3, with the label "[[10,2,2]], 8+2 inputs (C: Clifford)". It matches the gate list in App. S3. Open point: the arrow from the [[8,3,2]] label to the branch still crosses the "\|CCZ>" label; I'll route it from the cut line instead.
4. **[todo] Cut labels.** Two rows ([[n,k,2]] above, output state below, in the cut's colour); the current render has them. Check that [[14,1,2]] fits inside the right margin at print size.
5. **[todo] Litinski credit (comment 3).** You're right: Litinski draws [[15,1,3]], [[20,4,2]] and [[14,2,2]] as identity circuits, not this four-qubit parent, so "after Ref. 28" is misleading. Replace it with a sentence in the text, not the caption: "identity circuits of this kind were drawn for [[15,1,3]], [[14,2,2]] and [[20,4,2]] in Ref. [28]". Remove it from the caption.
6. **[todo] Caption (comment 7), 4 sentences, draft:**
   > "One borrowed identity, four factories. Each column is a π/8 rotation on the qubits marked; the 15 columns place every qubit in 8 rotations, every pair in 4 and every triple in 2, so all four qubits return to \|+> (Theorem 2). Cutting after a block and measuring the qubits to its left leaves the columns to its right, a synthesis of the output: \|CCZ>, \|CS> or \|T>. Two more T rotations and one more check on the \|CCZ> output give [[10,2,2]]."

## C. Fig. 2 (results; single column)

1. **[draft] Bound label (comment 9).** "n ≥ T-count+1" next to the dashed line, inside the plot. In the current render it still overlaps the 12 tick label; move it to x ≈ 10.5, y ≈ 15, to the right of the y axis and above the line.
2. **[draft] Crosses (comment 10).** Yes: the crosses are the two k=7 candidates, n = 20 and 22, whose *exact* T-count is 17, below n−1 = 19 and 21, so they do not reach the bound. Legend text: "× k=7: exact T-count, below the bound".
3. **[draft] k=8 bars (comment 11).** Each dotted bar runs from a proved *lower* bound (isotropy floor, 17) to a heuristic *upper* bound (18, 20, 22 at n = 20, 22, 24). The upper bound reaches n−1 at n = 20, 22 and 24, but the lower bound does not, so we can't tell whether they meet it. Legend text, in contrast to "exact": "⋮ k=8: lower to upper bound only".
4. **[draft] Panel (b), all levels (comment 8).** One bar per row:
   - l=2 (9 classes) and l=4, k≤4 (4,167) as hatched "no synthillation" bars;
   - l=3 split into k≤4 (22), k=5 (65) and k=6 (679), grey = same inputs as synthillation, orange = fewer;
   - labels on the right: "22: 0 fewer", "65: 11 fewer", "679: 106 fewer", "9", "4,167".
   - The separate coverage sentence is deleted; the bars carry it.
   - Open choice: you suggested one line "l=2: 9, l=3: 766, l=4: 4,167"; the bars already show this, so I'd drop the line. Tell me if you want both.
5. **[todo] Caption (comment 7), 5 sentences, draft:**
   > "Smallest factory found for each of the 766 output classes at l=3 with k ≤ 6, at the exact T-count of its target. (a) Every factory meets n ≥ T-count+1 (Corollary 3, dashed); the five on it are the smallest possible, open circles in the synthillation paper, diamonds found by the search, labelled with their output-error coefficient c. Orange: fewer inputs than every synthillation construction; ring: also beats every route that converts distilled states by stabilizer operations. Crosses: k=7 candidates with exact T-count below the bound; dotted bars: k=8, lower to upper bound only. (b) Classes with a factory at each level, and at l=3 the share needing fewer inputs than synthillation (none needs more)."

## D. Main text

1. **[todo] Framing (comment 1).** Keep "a synthesis of the target", never "optimal compilation", except for the five factories at the bound. Add one sentence after Corollary 3: "Most factories are far from this: the median smallest factory uses 17 more inputs than the bound and three check qubits (`outputs/l3_best_per_class.json`), so the checks cost synthesis, and the search measures how much." This answers the referee concern that d=2 is "just compilation".
2. **[todo] Why d=2 is tractable, and d>2 is not (comment 2).** In the Discussion, replace "where the checks can no longer be symmetric" with:
   > "At distance 2 a fault is caught whenever its rotation touches a check, so the factory condition stays linear in the coefficients; at d ≥ 3 no two rotations may share a check pattern while differing on the outputs, a condition on supports rather than coefficients, which the linear stages cannot impose. Our follow-up work reaches d ≥ 3 at l = 3 with group-symmetry restrictions and AI-assisted search~[jain2026symmetry]; a systematic search there is still missing."

   **Check before adding:** the pair condition on [[15,1,3]] (none) and [[14,2,2]] (some), with a short script.
   **Your point (11:0x):** say explicitly what going beyond d=2 requires, and that this paper supplies it: the follow-up
   AI-assisted d>2 search *used* the borrowed-identity condition of this paper and added the pair condition above;
   it did not replace it. Draft: "Beyond distance 2 the condition of Theorem 2 still applies, and one more is needed:
   no two rotations may share a check pattern while differing on the outputs. Our follow-up work imposes both, with
   group-symmetry restrictions and AI-assisted search, to reach d ≥ 3 at l = 3 [jain2026symmetry]; a systematic search
   for the second condition is the open problem." 
3. **[todo] Litinski sentence (B5)** in the introduction's Litinski paragraph, so the credit is accurate.
4. **[todo] Recount** after all edits: words (APS rule), abstract characters, End Matter pages, and caption lengths, each no taller than its figure.

## E. Order

1. Fig. 1 fixes (B3–B6).
2. Fig. 2 fixes (C1, C5).
3. Captions in main.tex.
4. Pair-condition check, then text edits D1–D3.
5. Compile, counts and page check.
6. Render pages 1–3 and look at them.

Roughly 45 minutes.

## F. Decisions for you

1. Fig. 2(b): drop the one-line coverage text (my suggestion), or keep it?
2. The D1 sentence (median 17 above the bound, three checks): include it?
3. The D2 d>2 sentence: include it, after the pair-condition check?
