## Raw Diff

No feedback

## Fixture Contract

No feedback

## 6-Cycle Clique Expectations

No feedback

## Support-Extension Logic

Why do I need `_graph_support` specifically? is it used somewhere else? if not, why not just encode this information in json file and use it here?

continue but make a note that we need to rewrite this unittest

Resolution: rewrote the support-extension regression to remove `_graph_support` and dropped activated support from the fixture entirely. `wang_magron_example_3_support_extension` now keeps only observable edge expectations, while the test hard-codes `activated_supp` locally as the two semantic support monomials generated from `yz` and `yz * x`. `seed_edges` remain in the fixture only to compute the added-edge delta.

## Why NoElimination Is Dense

This test needs to be re-written too. reference uses maximal chordal extension, i.e maximalelimination() method here.

Resolution: rewrote the cycle regression to target `MaximalElimination()` semantics via `wang_magron_example_3_maximal_chordal`. The test now asserts the exact maximal-elimination clique and keeps the `MF()` checks at the summary-stat level.

## Why Seed Edges Can Induce New Edges

No feedback

Verification: `julia --project -e 'using Pkg; Pkg.test()'` passed with `Pass 2309  Broken 2  Total 2311`.
