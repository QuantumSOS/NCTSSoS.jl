# test/problems/trace_polynomial/t2_broyden_tridiagonal_trace.jl
# T2: NC Broyden Tridiagonal (Trace)
#
# Source: Wang & Magron, "Exploiting Term Sparsity in Noncommutative
# Polynomial Optimization" (Table 10); benchmark metadata mirrored in
# `test-examples.typ` as T2.
#
# The literature source reports sparse optimum 0 and sparse max block size 6
# for this benchmark family. To keep the regression small and stable, we build
# the underlying ordinary polynomial, compute its Newton-chip basis, and inject
# that basis into the tracial relaxation through `moment_basis`. The solve path
# uses the dualized SOS formulation. The test only checks literature-backed
# quantities: the objective stays near 0 and the resulting term-sparse max
# block size does not exceed the published sparse block size.

using Test, TOML, NCTSSoS

const T2_OBJECTIVE_ATOL = 1e-5

function _t2_broyden_tridiagonal_trace_cases()
    path = joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "t2_broyden_tridiagonal_trace.toml")
    data = TOML.parsefile(path)
    haskey(data, "cases") || error("Missing key `cases` in T2 expectations TOML.")
    return data["cases"]
end

function broyden_tridiagonal_trace_problem(n::Int)
    n >= 2 || throw(ArgumentError("Broyden tridiagonal trace benchmark requires n >= 2."))

    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    poly_type = typeof(x[1] + x[2])
    one_poly = one(poly_type)

    first_residual = 3.0 * x[1] - 2.0 * x[1]^2 - 2.0 * x[2] + one_poly
    objective = adjoint(first_residual) * first_residual

    for i in 2:(n - 1)
        residual = 3.0 * x[i] - 2.0 * x[i]^2 - x[i - 1] - 2.0 * x[i + 1] + one_poly
        objective += adjoint(residual) * residual
    end

    last_residual = 3.0 * x[n] - 2.0 * x[n]^2 - x[n - 1] + one_poly
    objective += adjoint(last_residual) * last_residual

    poly_pop = polyopt(objective, reg)
    trace_pop = polyopt(tr(objective) * one(typeof(x[1])), reg)
    return poly_pop, trace_pop
end

function _t2_max_block_size(sparsity)
    return maximum(reduce(vcat, map(sparsity.cliques_term_sparsities) do ts
        length.(ts[1].block_bases)
    end))
end

@testset "T2 NC Broyden Tridiagonal (Trace)" begin
    cases = _t2_broyden_tridiagonal_trace_cases()
    @test [case["expected"]["n"] for case in cases] == [20, 40]

    @testset "n=$(case["expected"]["n"])" for case in cases
        expected = case["expected"]
        n = expected["n"]
        poly_pop, trace_pop = broyden_tridiagonal_trace_problem(n)

        sparsity_config = SolverConfig(
            optimizer=nothing,
            cs_algo=MF(),
            ts_algo=MMD(),
        )
        sparsity = compute_sparsity(trace_pop, sparsity_config)
        @test _t2_max_block_size(sparsity) == expected["sparse_mb"]

        solve_config = SolverConfig(
            optimizer=SOLVER,
            cs_algo=MF(),
            ts_algo=MMD(),
        )
        result = cs_nctssos(trace_pop, solve_config; dualize=true)

        @test result.objective ≈ expected["sparse_opt"] atol = T2_OBJECTIVE_ATOL
        @test maximum(reduce(vcat, result.moment_matrix_sizes)) <= expected["sparse_mb"]
    end
end
