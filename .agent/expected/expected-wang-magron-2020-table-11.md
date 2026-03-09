# Expected Value Extraction: `wang-magron-2020-table-11`

## Source

Section 6.2, Table 11 of Wang and Magron (2020): constrained trace minimization for the nc Broyden banded function over D, selected row n = 5.

## Problem Formulation

Section 6.2 studies constrained trace minimization for the noncommutative Broyden banded function over a semialgebraic set `D`. The objective family is the same banded polynomial defined on page 24, while pages 25-26 define
`D = {1 - X_1^2, ..., 1 - X_n^2, X_1 - 1/3, ..., X_n - 1/3}`.
Table 11 selects the sparse CS+TS row `n = 5` at minimum relaxation order `d_hat = 3`; the repository mirrors that case in `docs/src/examples/literate/wang_magron_2021_example_6_2.jl` with `TABLE_11_REF` and reproduces it via `run_trace_case(broyden_banded_objective, 5, 3; constrained = true, cs_algo = MF(), ts_algo = MMD())`. Reasoning: pages 25-26 define the feasible set, page 30 gives the selected table row, and the Literate example encodes the same constrained benchmark directly.

## How To Obtain The Expected Values

Run this snippet in the current project's native framework/tooling to reproduce or inspect the assertion target.

```julia
using Test
include("docs/src/examples/literate/wang_magron_2021_example_6_2.jl")

observed = run_trace_case(
    broyden_banded_objective,
    TABLE_11_REF.n,
    3;
    constrained=true,
    cs_algo=MF(),
    ts_algo=MMD(),
)
@test observed.objective ≈ TABLE_11_REF.objective atol = 1e-2
@test TABLE_11_REF.paper_sparse_mb == 19
```

## Why These Values

The selected row identifier n = 5 and the paper reference values mb = 19 and objective = 3.113 are deterministic and are the stable literature targets echoed by the in-repo example constants. Reasoning: Table 11 reports those values directly, and docs/src/examples/literate/wang_magron_2021_example_6_2.jl binds the same row as TABLE_11_REF.

## Why Not Others

Table 11 also reports runtime and comparison columns for CS and dense relaxations, but those values depend on solver environment or compare alternative formulations. The later regression target is the chosen CS+TS row, not every auxiliary column.
