# Expected Value Extraction: `wang-magron-2020-table-9`

## Source

Section 6.2, Table 9 of Wang and Magron (2020): trace minimization for the nc Broyden banded function, selected row n = 10.

## Problem Formulation

Section 6.2 studies unconstrained trace minimization for the noncommutative Broyden banded function. On page 24, the paper defines

```math
f_{Bb}(x) = \sum_{i=1}^n \left(2X_i + 5X_i^3 + 1 - \sum_{j \in J_i}(X_j + X_j^2)\right)^\star
\left(2X_i + 5X_i^3 + 1 - \sum_{j \in J_i}(X_j + X_j^2)\right),
```

with `J_i = {j | j != i, max(1, i - 5) <= j <= min(n, i + 1)}`. Table 9 selects the row `n = 10`; the repository mirrors that row in `docs/src/examples/literate/wang_magron_2021_example_6_2.jl` with `TABLE_9_REF` and reproduces it via `run_trace_case(broyden_banded_objective, 10, 3; ts_algo = MMD())`.

## How To Obtain The Expected Values

Run this snippet in the current project's native framework/tooling to reproduce or inspect the assertion target.

```julia
using Test
include("docs/src/examples/literate/wang_magron_2021_example_6_2.jl")

observed = run_trace_case(broyden_banded_objective, TABLE_9_REF.n, 3; ts_algo=MMD())
@test observed.objective ≈ TABLE_9_REF.objective atol = 1e-4
@test TABLE_9_REF.paper_sparse_mb == 29
```

## Why These Values

The selected row identifier n = 10 and the paper reference values mb = 29 and objective = 0 are deterministic and are the stable literature targets echoed by the in-repo example constants. Reasoning: Table 9 reports those values directly, and docs/src/examples/literate/wang_magron_2021_example_6_2.jl binds the same row as TABLE_9_REF.

## Why Not Others

Table 9 also reports runtime, but runtime varies by machine and solver configuration, so it does not define correctness. Dense-column dashes are memory-limit outcomes, not deterministic assertion targets for this repo.
