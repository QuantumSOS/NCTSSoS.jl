# =============================================================================
# Fermionic Parity Superselection Rule Tests
# =============================================================================
# Tests for the fermionic parity superselection rule.
#
# Physical basis: In fermionic systems, parity superselection forbids observables
# with odd fermion number from having non-zero expectation values.
# Only operators with even parity (even number of creation/annihilation operators)
# can have non-zero physical expectation values.
#
# The fermionic operators satisfy the canonical anti-commutation relations (CAR):
#   {a_i, a†_j} = δ_{ij}
#   {a_i, a_j} = 0
#   {a†_i, a†_j} = 0
# =============================================================================

using NCTSSoS, Test
using JuMP
using NCTSSoS: variable_indices, simplify  # Disambiguate simplify from JuMP

# Load shared solver configuration
include("../setup.jl")

@testset "Fermionic Parity Superselection" begin

    @testset "has_even_parity helper function" begin
        # Create fermionic variables for testing
        registry, (a, a_dag) = create_fermionic_variables(1:3)

        # Test 1: Identity (empty word) has even parity (0 operators)
        M = typeof(a[1])
        m_identity = one(M)
        @test has_even_parity(m_identity) == true

        # Test 2: Single annihilation operator has odd parity (1 operator)
        @test has_even_parity(a[1]) == false

        # Test 3: Single creation operator has odd parity (1 operator)
        @test has_even_parity(a_dag[1]) == false

        # Test 4: Two operators (a_dag * a) have even parity (2 operators)
        # Before simplification, product has 2 operators
        pair_mono = a_dag[1] * a[1]
        @test has_even_parity(pair_mono) == true

        # Test 5: Three operators have odd parity
        # a_dag[1] * a_dag[2] * a[1] has 3 operators
        triple_mono = a_dag[1] * a_dag[2] * a[1]
        @test has_even_parity(triple_mono) == false

        # Test 6: Four operators have even parity
        # a_dag[1] * a_dag[2] * a[2] * a[1] has 4 operators
        quad_mono = a_dag[1] * a_dag[2] * a[2] * a[1]
        @test has_even_parity(quad_mono) == true
    end

    @testset "has_even_parity on simplified polynomials" begin
        # After simplification, polynomials may have multiple terms
        # All terms in even-parity expressions should still have even parity
        registry, (a, a_dag) = create_fermionic_variables(1:2)

        # Simplified even-parity polynomial: a_dag[1] * a[1] after normal ordering
        even_poly = simplify(a_dag[1] * a[1])
        for t in terms(even_poly)
            @test has_even_parity(t.monomial) == true
        end

        # Simplified polynomial: a[1] * a_dag[1] = delta - a_dag[1]*a[1]
        # Both identity (from delta) and number operator have even parity
        anticomm_poly = simplify(a[1] * a_dag[1])
        for t in terms(anticomm_poly)
            @test has_even_parity(t.monomial) == true
        end
    end

    @testset "Parity constraint injection in moment relaxation" begin
        # Create a simple fermionic problem and verify basis contains BOTH even and odd parity monomials
        # The parity superselection is enforced via constraints, not filtering
        registry, (a, a_dag) = create_fermionic_variables(1:2)

        # Simple even-parity Hamiltonian: H = a_1^dag a_1 + a_2^dag a_2
        ham = a_dag[1] * a[1] + a_dag[2] * a[2]

        # Create optimization problem
        pop = polyopt(ham, registry)

        # Use low order for testing with NoElimination to keep full basis
        solver_config = SolverConfig(optimizer=SOLVER, order=1, ts_algo=NoElimination())
        res = cs_nctssos(pop, solver_config)

        # The test passes if solution is obtained without errors
        # Ground state energy of number operators is >= 0
        @test res.objective >= -1e-6

        # Verify result type
        @test res isa NCTSSoS.PolyOptResult
    end

    @testset "Basis contains odd-parity monomials" begin
        # Verify that the basis includes odd-parity monomials (not filtered out)
        registry, (a, a_dag) = create_fermionic_variables(1:2)

        ham = a_dag[1] * a[1] + a_dag[2] * a[2]
        pop = polyopt(ham, registry)

        # Build moment problem to inspect basis
        order = 1
        corr_sparsity = NCTSSoS.correlative_sparsity(pop, order, NoElimination())

        cliques_objective = map(corr_sparsity.cliques) do clique_indices
            clique_set = Set(clique_indices)
            reduce(+, [
                issubset(variable_indices(mono), clique_set) ? coef * mono : zero(coef) * one(mono)
                for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))
            ])
        end

        initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
            NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
        end

        cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
            NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, NoElimination())
        end

        mp = NCTSSoS.moment_relax(pop, corr_sparsity, cliques_term_sparsities)

        # Check that basis contains BOTH even and odd parity monomials
        even_count = count(has_even_parity, mp.total_basis)
        odd_count = count(!has_even_parity, mp.total_basis)

        @test even_count > 0  # Has even-parity monomials
        @test odd_count > 0   # Has odd-parity monomials (NOT filtered out)
    end

    @testset "Zero constraints added for odd-parity entries" begin
        # Verify that Zero constraints are added for odd-parity moment entries
        registry, (a, a_dag) = create_fermionic_variables(1:2)

        ham = a_dag[1] * a[1] + a_dag[2] * a[2]
        pop = polyopt(ham, registry)

        order = 1
        corr_sparsity = NCTSSoS.correlative_sparsity(pop, order, NoElimination())

        cliques_objective = map(corr_sparsity.cliques) do clique_indices
            clique_set = Set(clique_indices)
            reduce(+, [
                issubset(variable_indices(mono), clique_set) ? coef * mono : zero(coef) * one(mono)
                for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))
            ])
        end

        initial_activated_supps = map(zip(cliques_objective, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases)) do (partial_obj, cons_idx, mom_mtx_base)
            NCTSSoS.init_activated_supp(partial_obj, corr_sparsity.cons[cons_idx], mom_mtx_base)
        end

        cliques_term_sparsities = map(zip(initial_activated_supps, corr_sparsity.clq_cons, corr_sparsity.clq_mom_mtx_bases, corr_sparsity.clq_localizing_mtx_bases)) do (init_act_supp, cons_idx, mom_mtx_bases, localizing_mtx_bases)
            NCTSSoS.term_sparsities(init_act_supp, corr_sparsity.cons[cons_idx], mom_mtx_bases, localizing_mtx_bases, NoElimination())
        end

        mp = NCTSSoS.moment_relax(pop, corr_sparsity, cliques_term_sparsities)

        # Count Zero constraints
        zero_count = count(c -> c[1] == :Zero, mp.constraints)

        # Should have Zero constraints for odd-parity entries
        @test zero_count > 0
    end

    @testset "Objective validation for odd parity" begin
        registry, (a, a_dag) = create_fermionic_variables(1:2)

        # Odd-parity objective: single annihilation operator (has 1 operator)
        # Need to create a Polynomial explicitly to pass to polyopt
        odd_poly = Polynomial([Term(1.0, a[1])])

        # This should throw an error because odd-parity operators have zero expectation
        @test_throws ErrorException polyopt(odd_poly, registry)

        # Even-parity objective should succeed
        even_poly = a_dag[1] * a[1] + a_dag[2] * a[2]
        pop = polyopt(even_poly, registry)
        @test pop isa NCTSSoS.PolyOpt
    end

    @testset "Pure odd-parity objective rejected" begin
        registry, (a, a_dag) = create_fermionic_variables(1:2)

        # Pure odd-parity: just a single creation operator
        odd_poly = Polynomial([Term(1.0, a_dag[1])])

        @test_throws ErrorException polyopt(odd_poly, registry)
    end

    @testset "Mixed parity objective rejected" begin
        registry, (a, a_dag) = create_fermionic_variables(1:2)

        # Mixed: has both even and odd parity terms
        # Even term: a_dag[1] * a[1] (2 operators)
        # Odd term: a[1] (1 operator)
        even_term = Term(1.0, a_dag[1] * a[1])
        odd_term = Term(1.0, a[1])
        mixed_poly = Polynomial([even_term, odd_term])

        # Should reject because of the odd-parity term
        @test_throws ErrorException polyopt(mixed_poly, registry)
    end

end
