# # [Wang & Magron (2021): Examples 3.4, 3.8, and 5.3](@id wang-magron-2021-paper-examples)
#
# Reproduce three small term-sparsity examples from Wang & Magron (2021)
# [wangExploitingTermSparsity2021](@cite):
#
# - Example 3.4: unconstrained eigenvalue optimization where the sparse order-2
#   relaxation loses tightness after chordal extension.
# - Example 3.8: constrained eigenvalue optimization where the sparse and dense
#   order-2 relaxations both return `-1`.
# - Example 5.3: the trace analogue of Example 3.8, again with value `-1`.
#
# This page uses the public APIs [`polyopt`](@ref), [`compute_sparsity`](@ref),
# and [`cs_nctssos`](@ref). Instead of reconstructing the paper figures edge by
# edge, we summarize the term-sparsity graph through the resulting moment-matrix
# block sizes.

# ## Setup
#
# `COSMO` is enough for these small reproductions, so the page runs without a
# commercial solver.

using NCTSSoS, COSMO, JuMP

const SOLVER = optimizer_with_attributes(
    COSMO.Optimizer,
    "verbose" => false,
    "eps_abs" => 1e-7,
    "eps_rel" => 1e-7,
);

flatten_sizes(sizes) = reduce(vcat, sizes)

function block_sizes(result)
    sort(flatten_sizes(result.moment_matrix_sizes))
end

function dense_order2_config()
    SolverConfig(
        optimizer=SOLVER,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=NoElimination()
    )
end

function sparse_order2_config()
    SolverConfig(
        optimizer=SOLVER,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=MMD()
    )
end

function build_example_3_4()
    reg, (x,) = create_noncommutative_variables([("x", 1:3)])
    f = 1.0 * x[1]^2 - x[1] * x[2] - x[2] * x[1] + 3.0 * x[2]^2 -
        2.0 * x[1] * x[2] * x[1] + 2.0 * x[1] * x[2]^2 * x[1] -
        x[2] * x[3] - x[3] * x[2] + 6.0 * x[3]^2 +
        9.0 * x[2]^2 * x[3] + 9.0 * x[3] * x[2]^2 -
        54.0 * x[3] * x[2] * x[3] + 142.0 * x[3] * x[2]^2 * x[3]
    return polyopt(f, reg)
end

function build_example_3_8()
    reg, (x,) = create_noncommutative_variables([("x", 1:2)])
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    g = 4.0 - x[1]^2 - x[2]^2
    h = x[1] * x[2] + x[2] * x[1] - 2.0
    return polyopt(f, reg; eq_constraints=[h], ineq_constraints=[g])
end

function build_example_5_3()
    reg, (x,) = create_noncommutative_variables([("x", 1:2)])
    f = 2.0 - x[1]^2 + x[1] * x[2]^2 * x[1] - x[2]^2
    objective = (1.0 * tr(f)) * one(typeof(x[1]))
    g = (1.0 * tr(4.0 - x[1]^2 - x[2]^2)) * one(typeof(x[1]))
    h = (1.0 * tr(x[1] * x[2] + x[2] * x[1] - 2.0)) * one(typeof(x[1]))
    return polyopt(objective, reg; eq_constraints=[h], ineq_constraints=[g])
end

nothing #hide

# ## Example 3.4: sparse order 2 is slightly looser
#
# The paper reports `lambda_1(f) ~= -0.00355`, while the dense optimum is `0`.

example_3_4 = build_example_3_4();
nothing #hide
# `example_3_4` is a `PolyOpt` for the unconstrained paper benchmark.

dense_3_4 = cs_nctssos(example_3_4, dense_order2_config(); dualize=false);
nothing #hide
# `dense_3_4` is a `PolyOptResult`; it approximates the dense order-2 bound.

sparse_3_4 = cs_nctssos(example_3_4, sparse_order2_config(); dualize=false);
nothing #hide
# `sparse_3_4` is the term-sparse order-2 result using `MMD()`.

@assert isapprox(dense_3_4.objective, 0.0; atol=1e-4)
@assert isapprox(sparse_3_4.objective, -0.00355; atol=1e-4)

(
    dense_objective=dense_3_4.objective,
    sparse_objective=sparse_3_4.objective,
    sparse_blocks=block_sizes(sparse_3_4),
)

# The package reproduces the same qualitative behavior as the paper:
# the dense relaxation stays at `0`, while the sparse relaxation drops to about
# `-3.55e-3` after the approximate chordal extension.

# ## Example 3.8: sparse and dense both reach `-1`
#
# This is the constrained eigenvalue example with the same order-2 relaxation as
# the paper's `lambda_ts^(2,1)(f, S) = lambda_min(f, S) = -1`.

example_3_8 = build_example_3_8();
nothing #hide
# `example_3_8` is a constrained `PolyOpt` with one equality and one inequality.

dense_3_8 = cs_nctssos(example_3_8, dense_order2_config(); dualize=false);
nothing #hide
# `dense_3_8` is the dense order-2 relaxation result.

sparse_3_8 = cs_nctssos(example_3_8, sparse_order2_config(); dualize=false);
nothing #hide
# `sparse_3_8` is the term-sparse order-2 relaxation result.

@assert isapprox(dense_3_8.objective, -1.0; atol=1e-6)
@assert isapprox(sparse_3_8.objective, -1.0; atol=1e-6)

(
    dense_objective=dense_3_8.objective,
    sparse_objective=sparse_3_8.objective,
    dense_blocks=block_sizes(dense_3_8),
    sparse_blocks=block_sizes(sparse_3_8),
)

# In this case the sparse decomposition is exact at order 2, so the dense and
# sparse objectives agree up to solver tolerance.

# ## Example 5.3: trace lift of Example 3.8
#
# The paper reuses the same `f` and `S` as Example 3.8 but minimizes the trace
# value `mu(f, S)` instead of the eigenvalue.
#
# In the current public API, the trace formulation is written by lifting both
# the objective and the constraints with [`tr`](@ref) before calling
# [`polyopt`](@ref). This is the package-level interpretation used here.

example_5_3 = build_example_5_3();
nothing #hide
# `example_5_3` is a trace-polynomial `PolyOpt` built from `tr(f)`, `tr(g)`, and `tr(h)`.

dense_5_3 = cs_nctssos(example_5_3, dense_order2_config(); dualize=false);
nothing #hide
# `dense_5_3` is the dense order-2 trace relaxation result.

sparse_5_3 = cs_nctssos(example_5_3, sparse_order2_config(); dualize=false);
nothing #hide
# `sparse_5_3` is the term-sparse order-2 trace relaxation result.

@assert isapprox(dense_5_3.objective, -1.0; atol=1e-6)
@assert isapprox(sparse_5_3.objective, -1.0; atol=1e-6)

(
    dense_objective=dense_5_3.objective,
    sparse_objective=sparse_5_3.objective,
    dense_blocks=block_sizes(dense_5_3),
    sparse_blocks=block_sizes(sparse_5_3),
)

# The trace lift stays exact at order 2 as well, matching the paper value
# `mu_ts^(2,1)(f, S) = mu_2(f, S) = -1`.

# ## Summary
#
# | Example | Objective type | Dense order-2 value | Sparse order-2 value |
# |:--|:--|--:|--:|
# | 3.4 | Eigenvalue | `0` | `-0.00355` |
# | 3.8 | Eigenvalue | `-1` | `-1` |
# | 5.3 | Trace | `-1` | `-1` |
#
# Next steps:
#
# - [Bell inequalities](@ref bell-inequalities) for a larger state-polynomial example.
# - [Tracial Polynomial Optimization](@ref tracial-polynomial-optimization) for more `tr(...)` workflows.
# - [`compute_sparsity`](@ref) and [`cs_nctssos`](@ref) for the full API surface.
