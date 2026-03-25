# test/problems/trace_polynomial/t1_broyden_banded_trace.jl
# T1: NC Broyden Banded (Trace)
#
# Source: Wang & Magron, "Exploiting Term Sparsity in Noncommutative
# Polynomial Optimization" (Table 9); benchmark metadata mirrored in
# `test-examples.typ` as T1.
#
# This is a local-only standalone benchmark, not part of the curated CI suite in
# `test/runtests.jl`. The order-3 trace relaxation is too expensive for
# always-on CI, so the CI path keeps only the smaller T2 benchmark.
#
# Run manually with a local Mosek license:
#   make test-local
# or just this script:
#   make test-local-one SCRIPT=test/problems/trace_polynomial/t1_broyden_banded_trace.jl
# or
#   julia --project test/problems/trace_polynomial/t1_broyden_banded_trace.jl
#
# Benchmark path:
#   - objective: tr(f)
#   - relaxation order: 3  (deg(f) = 6)
#   - sparsity: CS = MF, TS = MMD
#   - no Newton-chip basis injection
#
# The reviewed literature row is n = 20. We also keep a larger local sweep point
# at n = 40 to exercise the same sparse trace hierarchy at a more useful size.

using Test, TOML, NCTSSoS, JuMP, MosekTools

const T1_SOLVER = optimizer_with_attributes(
    Mosek.Optimizer,
    "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
    "MSK_IPAR_LOG" => 0,
)
const T1_OBJECTIVE_ATOL = 1e-6
const T1_ORDER = 3

function _t1_broyden_banded_trace_cases()
    path = joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "t1_broyden_banded_trace.toml")
    data = TOML.parsefile(path)
    haskey(data, "cases") || error("Missing key `cases` in T1 expectations TOML.")
    return data["cases"]
end

function broyden_banded_trace_problem(n::Int)
    n >= 2 || throw(ArgumentError("Broyden banded trace benchmark requires n >= 2."))

    reg, (x,) = create_noncommutative_variables([("x", 1:n)])
    poly_type = typeof(x[1] + x[min(2, n)])
    objective = zero(poly_type)

    for i in 1:n
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)
        neighbor_sum = sum(x[j] + x[j]^2 for j in jset; init=zero(poly_type))
        residual = 2.0 * x[i] + 5.0 * x[i]^3 + 1.0 - neighbor_sum
        objective += adjoint(residual) * residual
    end

    return polyopt(tr(objective) * one(typeof(x[1])), reg)
end

function _t1_max_block_size(sparsity)
    return maximum(reduce(vcat, map(sparsity.cliques_term_sparsities) do ts
        length.(ts[1].block_bases)
    end))
end

@testset "T1 NC Broyden Banded (Trace)" begin
    cases = _t1_broyden_banded_trace_cases()
    @test [case["expected"]["n"] for case in cases] == [20, 40]

    @testset "n=$(case["expected"]["n"])" for case in cases
        expected = case["expected"]
        n = expected["n"]
        trace_pop = broyden_banded_trace_problem(n)

        sparsity_config = SolverConfig(
            optimizer=nothing,
            order=T1_ORDER,
            cs_algo=MF(),
            ts_algo=MMD(),
        )
        sparsity = compute_sparsity(trace_pop, sparsity_config)
        @test _t1_max_block_size(sparsity) == expected["sparse_mb"]

        solve_config = SolverConfig(
            optimizer=T1_SOLVER,
            order=T1_ORDER,
            cs_algo=MF(),
            ts_algo=MMD(),
        )
        result = cs_nctssos(trace_pop, solve_config; dualize=true)

        @test result.objective ≈ expected["sparse_opt"] atol = T1_OBJECTIVE_ATOL
        @test maximum(reduce(vcat, result.moment_matrix_sizes)) == expected["sparse_mb"]
    end
end
