```@meta
EditURL = "../literate/sparsity_convergence.jl"
```

# [Stabilization vs. Exactness in the Sparse Hierarchy](@id sparsity-convergence)

!!! note "Prerequisites"
    This page assumes familiarity with the
    [CS-TSSOS hierarchy](@ref cs-tssos-hierarchy), where the bound
    ``\lambda^{\text{cs-ts}}_{d,k}`` is indexed by relaxation order ``d``
    and term-sparsity refinement order ``k``.

In the [CS-TSSOS hierarchy](@ref cs-tssos-hierarchy), increasing the
sparse order ``k`` at fixed ``d`` can only add edges to the term-sparsity
graph, weakly tightening the bound. At some ``k`` the graph stops
changing — its support and block structure become a fixed point. This is
**graph stabilization**.

It is tempting to read stabilization as "the hierarchy has converged."
But two independent properties hide behind that word:

| Property | Meaning | Type |
|:---|:---|:---|
| **Graph stabilization** | The term-sparsity graph at order ``k`` is a fixed point. | Structural |
| **Bound exactness** | ``\lambda^{\text{cs-ts}}_{d,k} = \lambda^{\text{cs}}_{d}`` (the sparse bound equals the dense bound at the same ``d``). | Numerical |

Graph stabilization does not imply bound exactness. Wang and Magron
exhibit an explicit counterexample in
[wangExploitingTermSparsity2020](@cite) (Example 3.4): the graph
stabilizes immediately, yet the sparse bound remains strictly loose.
They also show a problem where the two properties *happen to coincide*
(Example 3.8). We reproduce both below.

## Setup

We use a silent Mosek backend throughout. Two helpers keep the output compact:

- `term_summary` — extracts support sizes and block dimensions from one
  sparse-order step.
- `next_term_sparsities` — applies one additional support-extension step to a
  stored `SparsityResult`.

````julia
using NCTSSoS, MosekTools

const MOI = NCTSSoS.MOI
const SILENT_MOSEK = MOI.OptimizerWithAttributes(Mosek.Optimizer, MOI.Silent() => true)

function term_summary(term_sparsities)
    (
        support_sizes = [length(ts.term_sparse_graph_supp) for ts in term_sparsities],
        block_sizes = [length.(ts.block_bases) for ts in term_sparsities],
    )
end

function next_term_sparsities(sparsity, clique_idx, ts_algo)
    corr = sparsity.corr_sparsity
    current = sparsity.cliques_term_sparsities[clique_idx]
    activated_supp = reduce(
        NCTSSoS.sorted_union,
        [ts.term_sparse_graph_supp for ts in current]
    )
    cons_idx = corr.clq_cons[clique_idx]
    return NCTSSoS.term_sparsities(
        activated_supp,
        corr.cons[cons_idx],
        corr.clq_mom_mtx_bases[clique_idx],
        corr.clq_localizing_mtx_bases[clique_idx],
        ts_algo,
    )
end
````

## Stabilized graph, inexact bound (Example 3.4)

Consider the unconstrained noncommutative polynomial from
[wangExploitingTermSparsity2020](@cite), Example 3.4:

```math
f = X^2 - XY - YX + 3Y^2 - 2XYX + 2XY^2X - YZ - ZY + 6Z^2 + 9X^2Y + 9Z^2Y - 54ZYZ + 142ZY^2Z.
```

The `MMD()` chordal-extension heuristic is the closest built-in analogue of
the approximately-minimum strategy used in the paper. We compute two
successive sparse-order steps and compare their graph structure.

````julia
registry_34, (vars_34,) = create_noncommutative_variables([("X", 1:3)])
X_34, Y_34, Z_34 = vars_34

f_34 = X_34^2 - X_34 * Y_34 - Y_34 * X_34 + 3.0 * Y_34^2 -
    2.0 * X_34 * Y_34 * X_34 + 2.0 * X_34 * Y_34^2 * X_34 -
    Y_34 * Z_34 - Z_34 * Y_34 + 6.0 * Z_34^2 + 9.0 * X_34^2 * Y_34 +
    9.0 * Z_34^2 * Y_34 - 54.0 * Z_34 * Y_34 * Z_34 +
    142.0 * Z_34 * Y_34^2 * Z_34

pop_34 = polyopt(f_34, registry_34)
config_34 = SolverConfig(optimizer=SILENT_MOSEK, order=2, ts_algo=MMD())

sparsity_34_k1 = compute_sparsity(pop_34, config_34)
sparsity_34_k2 = next_term_sparsities(sparsity_34_k1, 1, MMD())

summary_34 = (
    algorithm = :MMD,
    k1 = term_summary(sparsity_34_k1.cliques_term_sparsities[1]),
    k2 = term_summary(sparsity_34_k2),
    stabilized = term_summary(sparsity_34_k1.cliques_term_sparsities[1]) ==
        term_summary(sparsity_34_k2),
    paper_lambda_1 = -0.00355,
    paper_lambda_min = 0.0,
)

@assert summary_34.stabilized

summary_34
````

````
(algorithm = :MMD, k1 = (support_sizes = [26], block_sizes = [[3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 1]]), k2 = (support_sizes = [26], block_sizes = [[3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 1]]), stabilized = true, paper_lambda_1 = -0.00355, paper_lambda_min = 0.0)
````

The support and block sizes at `k = 1` and `k = 2` are identical:
the graph is a fixed point. Yet the paper reports a sparse optimum of
``\lambda_1(f) \approx -0.00355`` versus a true minimum of
``\lambda_{\min}(f) = 0`` [wangExploitingTermSparsity2020](@cite).
Graph stabilization did not imply bound exactness.

## Stabilized graph, exact bound (Example 3.8)

The constrained problem from [wangExploitingTermSparsity2020](@cite),
Example 3.8, provides the positive counterpart:

```math
f = 2 - X^2 + XY^2X - Y^2, \qquad S = \{4 - X^2 - Y^2,\; XY + YX - 2\}.
```

Here the initial adjacency graph is already chordal, so no fill-in is
needed. We use `NoElimination()` and compare sparse orders `k = 1` and
`k = 2` via [`cs_nctssos`](@ref) and [`cs_nctssos_higher`](@ref).

````julia
registry_38, (vars_38,) = create_noncommutative_variables([("X", 1:2)])
X_38, Y_38 = vars_38

f_38 = 2.0 - X_38^2 + X_38 * Y_38^2 * X_38 - Y_38^2
g_ball_38 = 4.0 - X_38^2 - Y_38^2
g_link_38 = X_38 * Y_38 + Y_38 * X_38 - 2.0

pop_38 = polyopt(
    f_38,
    registry_38;
    ineq_constraints = [g_ball_38],
    eq_constraints = [g_link_38],
)
config_38 = SolverConfig(optimizer=SILENT_MOSEK, order=2, ts_algo=NoElimination())

result_38_k1 = cs_nctssos(pop_38, config_38)
result_38_k2 = cs_nctssos_higher(pop_38, result_38_k1, config_38)

summary_38 = (
    algorithm = :NoElimination,
    k1 = (
        objective = result_38_k1.objective,
        term_summary(result_38_k1.sparsity.cliques_term_sparsities[1])...,
    ),
    k2 = (
        objective = result_38_k2.objective,
        term_summary(result_38_k2.sparsity.cliques_term_sparsities[1])...,
    ),
    stabilized = term_summary(result_38_k1.sparsity.cliques_term_sparsities[1]) ==
        term_summary(result_38_k2.sparsity.cliques_term_sparsities[1]),
)

@assert summary_38.stabilized
@assert isapprox(summary_38.k1.objective, -1.0; atol=1e-6)
@assert isapprox(summary_38.k2.objective, -1.0; atol=1e-6)

summary_38
````

````
(algorithm = :NoElimination, k1 = (objective = -0.9999999212855436, support_sizes = [25, 17, 16], block_sizes = [[7], [3], [3]]), k2 = (objective = -0.9999999212855436, support_sizes = [25, 17, 16], block_sizes = [[7], [3], [3]]), stabilized = true)
````

Both sparse orders yield the same objective ``\approx -1.0``, matching the
dense SDP bound, and the graph is again a fixed point. This shows that
stabilization and exactness *can* coincide — but Example 3.4 above
demonstrates that they need not.

## Summary

The comparison table below collects both cases. Recall that *exactness*
here means the sparse SDP optimum equals the dense SDP optimum at the
same relaxation order ``d`` — not necessarily the true optimum of the
original problem.

Graph stabilization is a structural fixed point of the term-sparsity
refinement. It does not, by itself, certify that the sparse bound is
tight. Under additional algebraic conditions — detailed in
[wangExploitingTermSparsity2020](@cite) — stabilization does imply
exactness, but those conditions must be verified, not assumed.

````julia
comparison = [
    (
        example = "3.4",
        ts_algo = "MMD()",
        stabilized = string(summary_34.stabilized),
        outcome = "paper counterexample",
        objective_story = "paper: lambda_1(f) approx -0.00355 < 0 = lambda_min(f)",
    ),
    (
        example = "3.8",
        ts_algo = "NoElimination()",
        stabilized = string(summary_38.stabilized),
        outcome = "exact at the first sparse step",
        objective_story = "computed: k=1 = k=2 approx -1",
    ),
]

comparison
````

````
2-element Vector{@NamedTuple{example::String, ts_algo::String, stabilized::String, outcome::String, objective_story::String}}:
 (example = "3.4", ts_algo = "MMD()", stabilized = "true", outcome = "paper counterexample", objective_story = "paper: lambda_1(f) approx -0.00355 < 0 = lambda_min(f)")
 (example = "3.8", ts_algo = "NoElimination()", stabilized = "true", outcome = "exact at the first sparse step", objective_story = "computed: k=1 = k=2 approx -1")
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*
