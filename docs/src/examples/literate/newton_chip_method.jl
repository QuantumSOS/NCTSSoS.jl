# # [Newton Chip Method](@id newton-chip-method)
#
# This page shows the Newton chip method for ordinary and tracial
# noncommutative optimization. `NCTSSoS.newton_chip_basis` reduces the rows of
# an unconstrained `Problem` over one `Free()` site. For a scalar trace
# objective it uses the Newton cyclic chip reduction. The background is in
# [burgdorfOptimizationPolynomialsNonCommuting2016](@cite) and
# [klep2022Optimization](@cite).
#
# **Prerequisites**: the [polynomial optimization API](@ref polynomial-optimization)
# and the [trace polynomial example](@ref tracial-polynomial-optimization).

# ## Setup
#
# The direct Mosek backend is silent by default. The cyclic reduction also
# needs a solver for its separating linear programs.

using NCTSSoS, Mosek
solver = MosekBackend();

pretty_basis(rows, registry) = [sprint(show, row; context=:registry => registry) for row in rows];

# ## Eigenvalue problems
#
# In the unconstrained eigenvalue setting, the certified polynomial is
# `f - lambda`. Keep the right chips of words `w` whose Hermitian square
# `w' * w` appears in its support
# [burgdorfOptimizationPolynomialsNonCommuting2016](@cite).

registry, (x,) = site_variables(Free(), 1, 2);
objective = 1.0 * x[1]^2 + 2.0 * x[2] * x[1] * x[1] * x[2] + 1.0;
problem = Problem(objective);
dense_basis = collect(NCTSSoS.basis(registry, 2));
chip_basis = NCTSSoS.newton_chip_basis(problem, 2);

pretty_basis(dense_basis, registry)

# The reduced basis keeps four rows:

pretty_basis(chip_basis, registry)

# ### Configure dense and custom relaxations
#
# `order=2` generates the full basis. `rows=chip_basis` supplies the reduced
# basis. Inspect the built rows with the public `NCTSSoS.block_rows` iterator.

dense = NCTSSoS.build(problem, Relaxation(order=2));
chip = NCTSSoS.build(problem, Relaxation(rows=chip_basis));
@assert Set(pretty_basis(only(NCTSSoS.block_rows(chip)), registry)) ==
        Set(pretty_basis((last(only(NCTSSoS.terms(p))) for p in chip_basis), registry))
(dense.sizes.rows, chip.sizes.rows)
@assert (dense.sizes.rows, chip.sizes.rows) == (7, 4) #src

# ### Solve both relaxations
#
# The objective is 1 plus two Hermitian squares. Setting both operators to
# zero attains 1, so the exact minimum is 1.

dense_result = solve(dense; solver);
chip_result = solve(chip; solver);
@assert isapprox(chip_result.bound, dense_result.bound; atol=1e-6)
@assert isapprox(dense_result.bound, 1; atol=1e-6) #src
@assert isapprox(chip_result.bound, 1; atol=1e-6) #src
(
    dense_result.bound,
    chip_result.bound,
    dense_result.sizes.moments,
    chip_result.sizes.moments,
)

# The ordinary workflow is:
# 1. build an unconstrained `Problem`,
# 2. call `NCTSSoS.newton_chip_basis`,
# 3. use `Relaxation(rows=...)`.

# ## Tracial problems
#
# Trace identifies cyclic words because `tr(AB) = tr(BA)`. For an
# unconstrained single-trace objective over one free site, the Newton cyclic
# chip reduction returns pure operator rows [klep2022Optimization](@cite).
# It does not cover products of traces, mixed state/operator objectives,
# algebraic relations such as `U^2 = I`, or constrained trace problems.

trace_registry, (y,) = site_variables(Free(), 1, 2; functional=NCTSSoS.RealTrace());
trace_objective = trace_symbol(1.0 * y[1]^2 + 1.0 * y[2] * y[1]^2 * y[2] + 1.0);
trace_problem = Problem(trace_objective);
dense_trace_basis = collect(NCTSSoS.basis(trace_registry, 2));
chip_trace_basis = NCTSSoS.newton_chip_basis(trace_problem, 2; solver);
(
    dense_size = length(dense_trace_basis),
    chip_words = pretty_basis(chip_trace_basis, trace_registry),
    chip_size = length(chip_trace_basis),
)

# ### Configure dense and custom relaxations
#
# The helper requires `solver=` for the trace case. If a numerical
# separating certificate is not accepted, it keeps the candidate row.
# The build checks that the custom rows cover the required moments.

dense_trace = NCTSSoS.build(trace_problem, Relaxation(order=2));
chip_trace = NCTSSoS.build(trace_problem, Relaxation(rows=chip_trace_basis));
@assert Set(pretty_basis(only(NCTSSoS.block_rows(chip_trace)), trace_registry)) ==
        Set(pretty_basis((last(only(NCTSSoS.terms(p))) for p in chip_trace_basis), trace_registry))
(dense_trace.sizes.rows, chip_trace.sizes.rows)
@assert (dense_trace.sizes.rows, chip_trace.sizes.rows) == (19, 4) #src

# ### Solve both relaxations
#
# Again, both nonconstant terms are traces of squares. Zero operators attain
# 1. The reduced and dense bounds must agree within `1e-6`.

dense_trace_result = solve(dense_trace; solver);
chip_trace_result = solve(chip_trace; solver);
@assert isapprox(chip_trace_result.bound, dense_trace_result.bound; atol=1e-6)
@assert isapprox(dense_trace_result.bound, 1; atol=1e-6) #src
@assert isapprox(chip_trace_result.bound, 1; atol=1e-6) #src
(
    dense_trace_result.bound,
    chip_trace_result.bound,
    dense_trace_result.sizes.moments,
    chip_trace_result.sizes.moments,
)

# The tracial workflow is:
# 1. build a scalar single-trace objective with `trace_symbol`,
# 2. call `NCTSSoS.newton_chip_basis(problem; solver)`,
# 3. use `Relaxation(rows=...)`.
#
# For general trace polynomials, or trace problems with relations or
# constraints, use the dense `Relaxation(order=...)` route instead.

# ## Next steps
#
# - [Trace Polynomial](@ref tracial-polynomial-optimization): larger tracial examples.
# - [Polynomial Optimization](@ref polynomial-optimization): modeling entry point.
