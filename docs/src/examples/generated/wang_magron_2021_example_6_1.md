```@meta
EditURL = "../literate/wang_magron_2021_example_6_1.jl"
```

# Wang-Magron 2021 Example 6.1: Broyden-Banded Table 2 Row

This page reproduces one published row from Table 2 of Wang and Magron's
NCTSSOS paper [wangExploitingTermSparsity2021](@cite): the noncommutative
Broyden-banded benchmark at `n = 20`.

Broyden-type functions keep showing up in optimization papers for two reasons.
First, they isolate recognizable solver behavior instead of hiding it inside a
large application model. Second, reusing the same families keeps comparisons
across papers honest. The banded Broyden case is especially useful here
because each residual only talks to a short index neighborhood, so the example
directly tests whether a relaxation detects and exploits local sparsity.

## Setup

We solve the same benchmark twice:

- `dense`: no correlative sparsity, no term sparsity
- `sparse`: no correlative sparsity, but term sparsity via [`MMD`](@ref)

Table 2 reports the stable structural fields `mb` and `opt`. The paper's
notation table says printed `opt = 0` means the absolute objective is at most
`1e-4`, so that is the tolerance we check against.

````@example wang_magron_2021_example_6_1
using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)

const PAPER_ROW_N20 = (
    n=20,
    sparse_mb=15,
    sparse_opt_abs_le=1e-4,
    sparse_time_seconds=0.34,
    dense_mb=61,
    dense_opt_abs_le=1e-4,
    dense_time_seconds=1.42,
)
nothing #hide

flatten_sizes(sizes) = collect(Iterators.flatten(sizes))
max_block_size(result) = maximum(flatten_sizes(result.moment_matrix_sizes))

dense_config = SolverConfig(
    optimizer=SILENT_MOSEK,
    order=3,
    cs_algo=NoElimination(),
    ts_algo=NoElimination(),
)

sparse_config = SolverConfig(
    optimizer=SILENT_MOSEK,
    order=3,
    cs_algo=NoElimination(),
    ts_algo=MMD(),
)

function build_broyden_banded_problem(n::Int)
    registry, (X,) = create_noncommutative_variables([("X", 1:n)])
    poly_type = typeof(X[1] + X[min(n, 2)])
    objective = Float64(n) * one(poly_type)

    for i in 1:n
        neighbors = [j for j in max(1, i - 5):min(n, i + 1) if j != i]
        objective += 4.0 * X[i] + 4.0 * X[i]^2 + 10.0 * X[i]^3 + 20.0 * X[i]^4 + 25.0 * X[i]^6
        for j in neighbors
            objective += -2.0 * X[j] - 2.0 * X[j]^2 - 4.0 * X[i] * X[j] - 4.0 * X[i] * X[j]^2
            objective += -10.0 * X[i]^3 * X[j] - 10.0 * X[i]^3 * X[j]^2
        end
        for j in neighbors
            for k in neighbors
                objective += X[j] * X[k] + 2.0 * X[j]^2 * X[k] + X[j]^2 * X[k]^2
            end
        end
    end

    return polyopt(objective, registry)
end

function run_relaxation(pop, config)
    result = cs_nctssos(pop, config)
    return (
        mb=max_block_size(result),
        objective=result.objective,
    )
end
````

## Problem formulation

The noncommutative Broyden-banded objective is

```math
f = \sum_{i=1}^n g_i^\ast g_i,
\qquad
g_i = 2X_i + 5X_i^3 + 1 - \sum_{j \in J_i} (X_j + X_j^2),
```

with the banded neighbor set

```math
J_i = \{j : \max(1, i-5) \le j \le \min(n, i+1),\ j \ne i\}.
```

The structure matters. Each `g_i` only depends on `X_i` and a short band of
nearby generators, so a term-sparsity pass can keep the largest moment block
much smaller than the dense relaxation.

````@example wang_magron_2021_example_6_1
pop = build_broyden_banded_problem(PAPER_ROW_N20.n)
````

`pop` is a polynomial optimization problem with the paper's `n = 20` benchmark.

## Run the sparse and dense relaxations

We use relaxation order `3`, matching the extracted Table 2 reproduction note.

````@example wang_magron_2021_example_6_1
sparse_result = run_relaxation(pop, sparse_config)
dense_result = run_relaxation(pop, dense_config)

observed = (
    n=PAPER_ROW_N20.n,
    sparse_mb=sparse_result.mb,
    sparse_objective=sparse_result.objective,
    dense_mb=dense_result.mb,
    dense_objective=dense_result.objective,
)
````

`observed` holds the locally reproduced row on the stable fields.

````@example wang_magron_2021_example_6_1
observed
````

## Compare against the published row

The extracted Table 2 row for `n = 20` records:

- sparse: `mb = 15`, printed `opt = 0`, printed `time = 0.34`
- dense: `mb = 61`, printed `opt = 0`, printed `time = 1.42`

We compare the stable fields directly and treat the printed `0` objective as
`|opt| <= 1e-4`, exactly as the paper's notation table prescribes. We do not
assert runtime equality because elapsed time depends on hardware and solver
versions.

````@example wang_magron_2021_example_6_1
paper = PAPER_ROW_N20

checks = (
    sparse_mb_matches=observed.sparse_mb == paper.sparse_mb,
    dense_mb_matches=observed.dense_mb == paper.dense_mb,
    sparse_objective_matches=abs(observed.sparse_objective) <= paper.sparse_opt_abs_le,
    dense_objective_matches=abs(observed.dense_objective) <= paper.dense_opt_abs_le,
)
````

`checks` summarizes whether the local run matches the paper on stable fields.

````@example wang_magron_2021_example_6_1
checks

all(values(checks))
````

## Interpretation

The point of this benchmark is not that the optimum is exotic. The point is
that the sparse and dense relaxations certify the same near-zero minimum while
using very different moment blocks. On the published `n = 20` row, the sparse
hierarchy uses `mb = 15` while the dense hierarchy uses `mb = 61`. That gap is
exactly why Wang and Magron use this family to showcase term sparsity.

## Next steps

- [`cs_nctssos`](@ref): build the first relaxation
- [`cs_nctssos_higher`](@ref): continue from a sparse first step when needed
- [Sparsity Convergence](@ref sparsity-convergence): another view of how
  sparsity changes block structure across relaxation orders

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

