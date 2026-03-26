# test/correlated_sparsity/csts_broyden_banded_polydisc.jl
# CS+TS Broyden banded on the NC polydisc
#
# Source: Wang & Magron, "Exploiting Term Sparsity in Noncommutative
# Polynomial Optimization" (Tables 7 and 11); benchmark metadata mirrored in
# `test-examples.typ` section 5.1.1 / `ex:csts-eig-bb` and section 5.1.3 /
# `ex:csts-tr-bb`.
#
# We keep only the smallest reviewed row `n = 5` in CI:
#   - eigenvalue CS+TS objective = 3.113, sparse mb = 11, dense mb = 156
#   - trace      CS+TS objective = 3.113, dense mb = 156
#
# For the trace benchmark we explicitly inject the embedded nc-word trace basis.
# Without that, the generic trace pipeline uses the broader state-word basis and
# the dense block size no longer matches the literature SOHS hierarchy.
#
# The current implementation contracts the sparse trace hierarchy slightly more
# than the paper's reported `mb = 19`: with the embedded operator basis we get
# `mb = 17` while preserving the literature objective. The fixture records both
# numbers so the test can protect the current implementation without hiding the
# paper value.

using Test, TOML, NCTSSoS

using NCTSSoS:
    MaxEntangled,
    NCStatePolynomial,
    NCStateWord,
    NonCommutativeAlgebra,
    StateWord,
    _embedded_trace_moment_basis,
    _uses_embedded_operator_basis,
    tr

if !@isdefined(SOLVER)
    include(joinpath(pkgdir(NCTSSoS), "test", "TestUtils.jl"))
end

const CSTS_BROYDEN_POLYDISC_EXPECTATIONS_PATH =
    joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "csts_broyden_banded_polydisc.toml")
const CSTS_BROYDEN_POLYDISC_OBJECTIVE_ATOL = 1e-3
const CSTS_BROYDEN_POLYDISC_ORDER = 3

function _csts_broyden_banded_polydisc_cases()
    data = TOML.parsefile(CSTS_BROYDEN_POLYDISC_EXPECTATIONS_PATH)
    haskey(data, "cases") || error("Missing key `cases` in CS+TS Broyden-polydisc expectations TOML.")
    return Dict(case["id"] => case["expected"] for case in data["cases"])
end

function _broyden_banded_polydisc_data(n::Int)
    n >= 2 || throw(ArgumentError("Broyden banded benchmark requires n >= 2."))

    reg, (x,) = create_noncommutative_variables([("X", 1:n)])
    poly_type = typeof(x[1] + x[min(2, n)])
    objective = zero(poly_type)

    for i in 1:n
        jset = setdiff(max(1, i - 5):min(n, i + 1), i)
        neighbor_sum = sum(x[j] + x[j]^2 for j in jset; init=zero(poly_type))
        residual = 2.0 * x[i] + 5.0 * x[i]^3 + 1.0 - neighbor_sum
        objective += adjoint(residual) * residual
    end

    constraints = vcat([1.0 - x[i]^2 for i in 1:n], [x[i] - 1.0 / 3 for i in 1:n])
    return reg, x, objective, constraints
end

function build_csts_broyden_banded_polydisc_problem(n::Int)
    reg, _, objective, constraints = _broyden_banded_polydisc_data(n)
    return polyopt(objective, reg; ineq_constraints=constraints)
end

function _lift_trace_constraint(poly, state_word_id)
    coeffs = [coef for (coef, _) in poly.terms]
    words = [NCStateWord(state_word_id, mono) for (_, mono) in poly.terms]
    return NCStatePolynomial(coeffs, words)
end

function build_csts_broyden_banded_trace_polydisc_problem(n::Int)
    reg, x, objective, constraints = _broyden_banded_polydisc_data(n)
    state_word_id = one(StateWord{MaxEntangled, NonCommutativeAlgebra, eltype(x[1].word)})
    trace_constraints = [_lift_trace_constraint(g, state_word_id) for g in constraints]
    trace_objective = tr(objective) * one(typeof(x[1]))
    return polyopt(trace_objective, reg; ineq_constraints=trace_constraints)
end

function _max_block_size(sparsity)
    block_sizes = Int[]
    for clique_ts in sparsity.cliques_term_sparsities
        for ts in clique_ts
            append!(block_sizes, length.(ts.block_bases))
        end
    end
    return maximum(block_sizes)
end

@testset "CS+TS Broyden banded on NC polydisc" begin
    cases = _csts_broyden_banded_polydisc_cases()

    @testset "Eigenvalue n=5" begin
        expected = cases["eigenvalue_n5"]
        pop = build_csts_broyden_banded_polydisc_problem(expected["n"])

        dense_structure = compute_sparsity(
            pop,
            SolverConfig(
                optimizer=nothing,
                order=CSTS_BROYDEN_POLYDISC_ORDER,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
            ),
        )
        @test _max_block_size(dense_structure) == expected["dense_mb"]

        sparse_structure = compute_sparsity(
            pop,
            SolverConfig(
                optimizer=nothing,
                order=CSTS_BROYDEN_POLYDISC_ORDER,
                cs_algo=MF(),
                ts_algo=MMD(),
            ),
        )
        @test _max_block_size(sparse_structure) == expected["sparse_mb"]

        result = cs_nctssos(
            pop,
            SolverConfig(
                optimizer=SOLVER,
                order=CSTS_BROYDEN_POLYDISC_ORDER,
                cs_algo=MF(),
                ts_algo=MMD(),
            );
            dualize=true,
        )

        @test result.objective ≈ expected["sparse_opt"] atol = CSTS_BROYDEN_POLYDISC_OBJECTIVE_ATOL
        @test maximum(flatten_sizes(result.moment_matrix_sizes)) == expected["sparse_mb"]
    end

    @testset "Trace n=5" begin
        expected = cases["trace_n5"]
        trace_pop = build_csts_broyden_banded_trace_polydisc_problem(expected["n"])

        trace_basis = _embedded_trace_moment_basis(trace_pop.registry, CSTS_BROYDEN_POLYDISC_ORDER)
        @test _uses_embedded_operator_basis(trace_basis)

        dense_structure = compute_sparsity(
            trace_pop,
            SolverConfig(
                optimizer=nothing,
                moment_basis=trace_basis,
                cs_algo=NoElimination(),
                ts_algo=NoElimination(),
            ),
        )
        @test _max_block_size(dense_structure) == expected["dense_mb"]

        sparse_structure = compute_sparsity(
            trace_pop,
            SolverConfig(
                optimizer=nothing,
                moment_basis=trace_basis,
                cs_algo=MF(),
                ts_algo=MMD(),
            ),
        )
        @test _max_block_size(sparse_structure) == expected["sparse_mb"]
        @test expected["sparse_mb"] <= expected["literature_sparse_mb"]

        result = cs_nctssos(
            trace_pop,
            SolverConfig(
                optimizer=SOLVER,
                moment_basis=trace_basis,
                cs_algo=MF(),
                ts_algo=MMD(),
            );
            dualize=true,
        )

        @test result.objective ≈ expected["sparse_opt"] atol = CSTS_BROYDEN_POLYDISC_OBJECTIVE_ATOL
        @test maximum(flatten_sizes(result.moment_matrix_sizes)) == expected["sparse_mb"]
        @test maximum(flatten_sizes(result.moment_matrix_sizes)) <= expected["literature_sparse_mb"]
    end
end
