# SOS (Sum of Squares) Dualization Tests
# Tests specific to SOS dualization components:
#   - Cαj coefficient extraction
#
# Note: Problem-specific SOS tests are in problems/ subdirectory.

using Test, NCTSSoS, JuMP
using SparseArrays
using NCTSSoS: get_Cαj

# SOLVER fallback for standalone/REPL execution
if !@isdefined(SOLVER)
    using MosekTools
    const SOLVER = optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_NUM_THREADS" => max(1, div(Sys.CPU_THREADS, 2)),
        "MSK_IPAR_LOG" => 0
    )
end

@testset "SOS Components" begin
    @testset "Cαj" begin
        # Test get_Cαj with new API: Vector{Monomial}, Matrix{Polynomial}
        registry, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Create a properly typed polynomial matrix
        poly_matrix = [
            1.0*x[1] 1.0*x[2]
            1.0*x[2] 1.0*x[1]+1.0*x[3]
        ]

        # Basis: the monomials themselves (sorted)
        basis_monomials = sort([one(x[1]), x[1], x[2], x[3]])

        C_α_js = get_Cαj(basis_monomials, poly_matrix)

        # Just verify structure - we have 5 non-zero entries
        @test length(C_α_js) == 5
        @test all(v -> v == 1.0, values(C_α_js))
    end

    @testset "Cαj complex" begin
        registry, (x,) = create_noncommutative_variables([("x", 1:2)])
        basis_polys = get_ncbasis(registry, 2)
        # For NC algebra, each basis poly is single-term; extract the monomial
        basis_monomials = [monomials(p)[1] for p in basis_polys]

        # Create properly typed polynomial matrix
        localizing_mtx = [
            1.0*x[1]-1.0*one(x[1])   1.0*x[2]^2-1.0*one(x[1])
            1.0*x[2]^2-1.0*one(x[1])   1.0*x[2]^3
        ]
        C_α_js = get_Cαj(basis_monomials, localizing_mtx)

        @test !isempty(C_α_js)
        @test all(v -> v == 1.0 || v == -1.0, values(C_α_js))
    end

    @testset "SOS dual objective uses affine constant coefficient" begin
        config = SolverConfig(optimizer=nothing, order=1, cs_algo=NoElimination(), ts_algo=NoElimination())
        n_symmetric_vars(n) = n * (n + 1) ÷ 2

        reg_real, (u,) = create_unipotent_variables([("u", 1:1)])
        real_pop = polyopt(1.0 * u[1], reg_real)
        real_sparsity = compute_sparsity(real_pop, config)
        real_mp = NCTSSoS.moment_relax(real_pop, real_sparsity.corr_sparsity, real_sparsity.cliques_term_sparsities)
        real_sos = NCTSSoS.sos_dualize(real_mp)
        real_expected_vars = sum(n_symmetric_vars(block.size) for block in real_mp.linear.psd_blocks_lin) +
            length(real_mp.linear.zero_constraints)
        @test num_variables(real_sos.model) == real_expected_vars

        reg_pauli, (σx, _, _) = create_pauli_variables(1:1)
        pauli_pop = polyopt(1.0 * σx[1], reg_pauli)
        pauli_sparsity = compute_sparsity(pauli_pop, config)
        pauli_mp = NCTSSoS.moment_relax(pauli_pop, pauli_sparsity.corr_sparsity, pauli_sparsity.cliques_term_sparsities)
        pauli_sos = NCTSSoS.sos_dualize(pauli_mp)
        pauli_expected_vars = sum(n_symmetric_vars(2 * block.size) for block in pauli_mp.linear.psd_blocks_lin) +
            length(pauli_mp.linear.zero_constraints)
        @test num_variables(pauli_sos.model) == pauli_expected_vars

        reg_state, (v,) = create_unipotent_variables([("v", 1:1)])
        state_objective = (1.0 * ς(v[1])) * one(typeof(v[1]))
        state_pop = polyopt(state_objective, reg_state)
        state_sparsity = compute_sparsity(state_pop, config)
        state_mp = NCTSSoS.moment_relax(state_pop, state_sparsity.corr_sparsity, state_sparsity.cliques_term_sparsities)
        state_sos = NCTSSoS.sos_dualize(state_mp)
        state_expected_vars = sum(n_symmetric_vars(size(mat, 1)) for (_cone, mat, _basis) in state_mp.constraints)
        @test num_variables(state_sos.model) == state_expected_vars
    end

    @testset "State SOS dualization rejects underspecified bases" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:1)])
        objective = tr(1.0 * x[1]) * one(typeof(x[1]))
        P = typeof(objective)
        M = typeof(only(monomials(objective)))
        basis = [one(M)]
        empty_constraints = Tuple{Symbol, Matrix{P}, Vector{M}}[]

        objective_mp = NCTSSoS.StateMomentProblem(objective, empty_constraints, basis, 0)
        @test_throws ArgumentError NCTSSoS.sos_dualize(objective_mp)

        identity_objective = one(P)
        constraint_matrix = reshape([objective], 1, 1)
        constraint_mp = NCTSSoS.StateMomentProblem(
            identity_objective,
            [(:PSD, constraint_matrix, basis)],
            basis,
            0,
        )
        @test_throws ArgumentError NCTSSoS.sos_dualize(constraint_mp)
    end
end
