# Expected Value Extraction: `wang-magron-2020-table-10`

## Source

Section 6.2, Table 10 of Wang and Magron (2020): trace minimization for the nc Broyden tridiagonal function, selected row n = 20.

## Problem Formulation

Section 6.2 studies unconstrained trace minimization for the noncommutative Broyden tridiagonal function. On page 24, the paper defines

```math
f_{Bt}(x) = (3X_1 - 2X_1^2 - 2X_2 + 1)^\star(3X_1 - 2X_1^2 - 2X_2 + 1)
+ \sum_{i=2}^{n-1}(3X_i - 2X_i^2 - X_{i-1} - 2X_{i+1} + 1)^\star(3X_i - 2X_i^2 - X_{i-1} - 2X_{i+1} + 1)
+ (3X_n - 2X_n^2 - X_{n-1} + 1)^\star(3X_n - 2X_n^2 - X_{n-1} + 1),
```

with no constraints. Table 10 selects the row `n = 20`; the repository mirrors that row in `docs/src/examples/literate/wang_magron_2021_example_6_2.jl` with `TABLE_10_REF` and reproduces it via `run_trace_case(broyden_tridiagonal_objective, 20, 2; ts_algo = MMD())`. Reasoning: page 24 gives the objective family, page 29 gives the selected trace-benchmark row, and the Literate example encodes the same case directly.

## How To Obtain The Expected Values

Run this snippet in the current project's native framework/tooling to reproduce or inspect the assertion target.

```julia
using Test
include("docs/src/examples/literate/wang_magron_2021_example_6_2.jl")

observed = run_trace_case(broyden_tridiagonal_objective, TABLE_10_REF.n, 2; ts_algo=MMD())
@test observed.objective ≈ TABLE_10_REF.objective atol = 1e-4
@test TABLE_10_REF.paper_sparse_mb == 6
```

## Why These Values

The selected row identifier n = 20 and the paper reference values mb = 6 and objective = 0 are deterministic and are the stable literature targets echoed by the in-repo example constants. Reasoning: Table 10 reports those values directly, and docs/src/examples/literate/wang_magron_2021_example_6_2.jl binds the same row as TABLE_10_REF.

## Why Not Others

Table 10 also reports runtime, but runtime varies by machine and solver configuration, so it does not define correctness. Dense-column dashes are memory-limit outcomes, not deterministic assertion targets for this repo.
