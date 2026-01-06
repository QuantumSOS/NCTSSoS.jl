# NCTSSoS is loaded by parent runtests.jl
# Exported: simplify, create_*_variables, get_ncbasis, terms, monomials, has_even_parity
# Internal (not exported): encode_index, decode_site, normal_order_key, is_normal_ordered, find_first_out_of_order, combine_like_terms
using NCTSSoS:
    encode_index,
    decode_site,
    normal_order_key,
    is_normal_ordered,
    find_first_out_of_order,
    combine_like_terms

# Note: The new API uses AlgebraType dispatch for simplification instead of SimplifyAlgorithm
# Each algebra type (NonCommutativeAlgebra, PauliAlgebra, UnipotentAlgebra, etc.) has its own simplification rules

@testset "Simplification Interface" begin
    @testset "NonCommutative Simplification" begin
        # NonCommutativeAlgebra with encoded indices uses site-aware simplification
        # Operators on different sites commute (sorted by site)
        # Operators on same site preserve order

        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Simple case: multiplication now returns Monomial (word concatenation)
        result = x[1] * x[2]
        @test result isa Monomial
        @test degree(result) == 2

        # After simplify, we get a Monomial (NC algebra returns Monomial)
        simplified = simplify(result)
        @test simplified isa Monomial
        @test degree(simplified) == 2

        # Same variable twice
        result2 = x[1] * x[1]
        @test result2 isa Monomial
        @test degree(result2) == 2
    end

    @testset "Pauli Simplification" begin
        # Pauli algebra has specific simplification rules:
        # - Pauli matrices square to identity
        # - Products of different Pauli matrices produce phases

        reg, (σx, σy, σz) = create_pauli_variables(1:2)

        # Single Pauli operator
        @test degree(σx[1]) == 1
        @test degree(σy[1]) == 1
        @test degree(σz[1]) == 1

        # Pauli variables are monomials
        @test σx[1] isa Monomial{PauliAlgebra}
    end

    @testset "Projector Simplification" begin
        # Projector algebra: P^2 = P (idempotency)
        reg, (P,) = create_projector_variables([("P", 1:3)])

        @test P[1] isa Monomial{ProjectorAlgebra}
        @test degree(P[1]) == 1

        # Multiple projectors - multiplication returns Monomial
        result = P[1] * P[2]
        @test result isa Monomial
        @test degree(result) == 2

        # Equivalence test: P[2] * P[1]^2 * P[2] == P[2] * P[1] * P[2] (since P^2 = P)
        # Projector simplify returns Monomial
        lhs = simplify(P[2] * P[1] * P[1] * P[2])
        rhs = simplify(P[2] * P[1] * P[2])
        @test lhs isa Monomial
        @test rhs isa Monomial
        @test lhs == rhs
    end

    @testset "Unipotent Simplification" begin
        # Unipotent algebra: U^2 = I (squares to identity)
        reg, (U,) = create_unipotent_variables([("U", 1:3)])

        @test U[1] isa Monomial{UnipotentAlgebra}
        @test degree(U[1]) == 1

        # Multiplication of different unipotent variables - returns Monomial
        result = U[1] * U[2]
        @test result isa Monomial
        @test degree(result) == 2

        # Equivalence test: U[2] * U[1]^2 * U[2] == I (since U^2 = I)
        # U[2] * U[1]^2 * U[2] = U[2] * I * U[2] = U[2]^2 = I
        # Unipotent simplify returns Monomial
        m = simplify(U[2] * U[1] * U[1] * U[2])
        @test m isa Monomial
        @test isempty(m.word)  # Identity monomial has empty word
    end


    @testset "Term Structure" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        t = Term(2.0, m)

        @test t.coefficient == 2.0
        @test t.monomial == m

        # Term operations
        t_neg = -t
        @test t_neg.coefficient == -2.0

        t_scaled = 3.0 * t
        @test t_scaled.coefficient == 6.0
    end

    @testset "simplify returns Monomial" begin
        # Test that simplify returns a new Monomial
        # NonCommutative simplify returns Monomial
        m = Monomial{NonCommutativeAlgebra}(UInt8[2, 1])  # Will be sorted by site

        result = simplify(m)
        @test result isa Monomial
    end

    @testset "Adjoint and Simplification Interaction" begin
        # Test adjoint involution on directly created monomials
        m = Monomial{NonCommutativeAlgebra}(UInt8[5, 9])
        m_adj = adjoint(m)

        # adjoint is involutory: adjoint(adjoint(m)).word == m.word
        @test adjoint(adjoint(m)).word == m.word
        @test m' == m_adj  # Julia syntax shorthand
    end
end

@testset "FermionicAlgebra Simplification" begin
    # FermionicAlgebra uses Wick's theorem with contractions
    # CAR: {aᵢ, aⱼ†} = δᵢⱼ, {aᵢ, aⱼ} = 0, {aᵢ†, aⱼ†} = 0
    # Nilpotency: aᵢ² = 0, (aᵢ†)² = 0
    # Returns Polynomial due to delta corrections

    reg, (a, a_dag) = create_fermionic_variables(1:3)

    @testset "Basic Operations" begin
        # Empty word → identity PhysicsMonomial
        m_empty = Monomial{FermionicAlgebra}(Int32[])
        result = simplify(m_empty)
        @test result isa PhysicsMonomial
        @test length(result.coeffs) == 1
        @test result.coeffs[1] == 1
        @test isempty(result.monos[1].word)

        # Single annihilation a₁ → unchanged
        m_a1 = a[1]
        result_a1 = simplify(m_a1)
        @test length(result_a1.coeffs) == 1
        @test result_a1.coeffs[1] == 1
        @test result_a1.monos[1].word == [1]

        # Single creation a₁† → unchanged
        m_a1_dag = a_dag[1]
        result_a1_dag = simplify(m_a1_dag)
        @test length(result_a1_dag.coeffs) == 1
        @test result_a1_dag.coeffs[1] == 1
        @test result_a1_dag.monos[1].word == [-1]
    end

    @testset "Anticommutation (CAR)" begin
        # a₁ a₁† = 1 - a₁† a₁ (verify both terms)
        m = a[1] * a_dag[1]  # a₁ a₁†
        result = simplify(m)
        @test length(result.coeffs) == 2

        # Sort terms by degree (identity first, then normal-ordered)
        perm = sortperm(result.monos, by=m -> degree(m))
        sorted_coeffs = result.coeffs[perm]
        sorted_monos = result.monos[perm]

        # First term: identity (scalar contraction)
        @test sorted_coeffs[1] == 1
        @test isempty(sorted_monos[1].word)

        # Second term: -a₁† a₁ (normal-ordered with sign)
        @test sorted_coeffs[2] == -1
        @test sorted_monos[2].word == [-1, 1]

        # a₁† a₁ → unchanged (already normal)
        m_normal = a_dag[1] * a[1]
        result_normal = simplify(m_normal)
        @test length(result_normal.coeffs) == 1
        @test result_normal.coeffs[1] == 1
        @test result_normal.monos[1].word == [-1, 1]

        # a₁ a₂ → a₁ a₂ (annihilators already in order by mode)
        m_diff = a[1] * a[2]
        result_diff = simplify(m_diff)
        @test length(result_diff.coeffs) == 1
        @test result_diff.coeffs[1] == 1
        @test result_diff.monos[1].word == [1, 2]  # Sorted by mode

        # a₁† a₂† → a₁† a₂† (creators already in order)
        m_cr = a_dag[1] * a_dag[2]
        result_cr = simplify(m_cr)
        @test length(result_cr.coeffs) == 1
        @test result_cr.coeffs[1] == 1
        @test result_cr.monos[1].word == [-1, -2]
    end

    @testset "Nilpotency" begin
        # a₁ a₁ = 0
        m_a1a1 = a[1] * a[1]
        @test iszero(m_a1a1)
        result_a1a1 = simplify(m_a1a1)
        # Zero PhysicsMonomial has no coeffs or all zeros
        @test iszero(result_a1a1)

        # a₁† a₁† = 0
        m_dag_dag = a_dag[1] * a_dag[1]
        @test iszero(m_dag_dag)
        result_dag_dag = simplify(m_dag_dag)
        # Zero PhysicsMonomial has no coeffs or all zeros
        @test iszero(result_dag_dag)

        # Direct construction using PhysicsMonomial (Monomial validates normal-order)
        # a₁ a₁ = 0 (nilpotent - same mode annihilators)
        pm_nilp = PhysicsMonomial{FermionicAlgebra}(Int32[1, 1])
        @test iszero(pm_nilp)

        # Non-nilpotent case: a₁ a₁† a₁ a₁† ≠ 0 (alternating)
        m_alt = a[1] * a_dag[1] * a[1] * a_dag[1]
        @test !iszero(m_alt)
    end

    @testset "Surplus-based iszero (net flux >= 2)" begin
        # Use PhysicsMonomial for non-normal-ordered words (Monomial validates)
        # a₁ a₂ a₁ = -a₁ a₁ a₂ = 0 (surplus of 2 annihilations for mode 1)
        pm_cross_mode = PhysicsMonomial{FermionicAlgebra}(Int32[1, 2, 1])
        @test iszero(pm_cross_mode)

        # a₁† a₂† a₁† = 0 (surplus of 2 creations for mode 1)
        pm_cross_mode_dag = PhysicsMonomial{FermionicAlgebra}(Int32[-1, -2, -1])
        @test iszero(pm_cross_mode_dag)

        # a₁ a₁ a₁† has surplus 1 (2 ann - 1 cre = 1)
        # Actually: a₁ a₁ a₁† = (a₁ a₁) a₁† = 0 * a₁† = 0 due to nilpotency
        pm_surplus_1 = PhysicsMonomial{FermionicAlgebra}(Int32[1, 1, -1])
        # PhysicsMonomial auto-simplifies and detects nilpotency
        @test iszero(pm_surplus_1)

        # a₁† a₁ a₁† has surplus -1 (1 - 2 = -1), NOT zero
        # a₁† a₁ a₁† = a₁† (1 - a₁† a₁) = a₁† - a₁† a₁† a₁ = a₁† - 0 = a₁† ≠ 0
        pm_surplus_neg1 = PhysicsMonomial{FermionicAlgebra}(Int32[-1, 1, -1])
        @test !iszero(pm_surplus_neg1)
        # PhysicsMonomial already simplified - just verify it's a₁†
        @test length(pm_surplus_neg1.coeffs) == 1
        @test pm_surplus_neg1.coeffs[1] == 1
        @test pm_surplus_neg1.monos[1].word == [-1]  # Just a₁†

        # a₁ a₂ a₃ a₁ a₂ = 0 (modes 1 and 2 each have surplus 2)
        pm_multi_surplus = PhysicsMonomial{FermionicAlgebra}(Int32[1, 2, 3, 1, 2])
        @test iszero(pm_multi_surplus)

        # a₁ a₂ a₁† a₂† = non-zero (each mode has surplus 0)
        pm_balanced = PhysicsMonomial{FermionicAlgebra}(Int32[1, 2, -1, -2])
        @test !iszero(pm_balanced)
    end

    @testset "Multi-mode" begin
        # a₁ a₂† → -a₂† a₁ (different modes, need to swap an and cr)
        m_cross = a[1] * a_dag[2]
        result_cross = simplify(m_cross)
        @test length(result_cross.coeffs) == 1
        @test result_cross.coeffs[1] == -1  # Sign from anticommutation
        @test result_cross.monos[1].word == [-2, 1]

        # a₁ a₁† a₂ a₂† → 4 terms from two contractions
        m_two_mode = a[1] * a_dag[1] * a[2] * a_dag[2]
        result_two_mode = simplify(m_two_mode)
        @test length(result_two_mode.coeffs) == 4

        # Check total coefficient sum (should be 1 + 1 - 1 - 1 = 0 is wrong; it's more complex)
        # The correct expansion is: (1 - a₁† a₁)(1 - a₂† a₂) = 1 - a₁† a₁ - a₂† a₂ + a₁† a₁ a₂† a₂
        @test 1 in result_two_mode.coeffs  # Identity term
    end

    @testset "Complex Case" begin
        # a₁ a₂ a₃ a₃† a₂† a₁†
        m_complex = a[1] * a[2] * a[3] * a_dag[3] * a_dag[2] * a_dag[1]
        result_complex = simplify(m_complex)

        # Multiple contraction possibilities should produce multiple terms
        @test length(result_complex.coeffs) > 1

        # Should include fully contracted term (identity)
        identity_idx = findfirst(m -> isempty(m.word), result_complex.monos)
        @test !isnothing(identity_idx)
    end

    @testset "has_even_parity" begin
        # Identity (empty word) has 0 operators → even parity
        m_identity = Monomial{FermionicAlgebra}(Int32[])
        @test has_even_parity(m_identity) == true

        # Single annihilation a₁ → 1 operator → odd parity
        @test has_even_parity(a[1]) == false

        # Single creation a₁† → 1 operator → odd parity
        @test has_even_parity(a_dag[1]) == false

        # Two operators a₁ a₂ → even parity
        m_two_ann = a[1] * a[2]
        @test has_even_parity(m_two_ann) == true

        # Two operators a₁† a₂† → even parity
        m_two_cre = a_dag[1] * a_dag[2]
        @test has_even_parity(m_two_cre) == true

        # Number operator a₁† a₁ → 2 operators → even parity
        m_number = a_dag[1] * a[1]
        @test has_even_parity(m_number) == true

        # Three operators → odd parity
        m_three = a_dag[1] * a[1] * a[2]
        @test has_even_parity(m_three) == false

        # Four operators → even parity
        m_four = a_dag[1] * a_dag[2] * a[1] * a[2]
        @test has_even_parity(m_four) == true
    end
end

@testset "BosonicAlgebra Simplification" begin
    # BosonicAlgebra uses rook numbers on Ferrers boards
    # CCR: [cᵢ, cⱼ†] = δᵢⱼ, [cᵢ, cⱼ] = 0, [cᵢ†, cⱼ†] = 0
    # NOT nilpotent: cᵢ² ≠ 0
    # Returns Polynomial due to delta corrections

    reg, (c, c_dag) = create_bosonic_variables(1:3)

    @testset "Basic Operations" begin
        # Empty word → identity PhysicsMonomial
        m_empty = Monomial{BosonicAlgebra}(Int32[])
        result = simplify(m_empty)
        @test result isa PhysicsMonomial
        @test length(result.coeffs) == 1
        @test result.coeffs[1] == 1
        @test isempty(result.monos[1].word)

        # Single annihilation c₁ → unchanged
        m_c1 = c[1]
        result_c1 = simplify(m_c1)
        @test length(result_c1.coeffs) == 1
        @test result_c1.coeffs[1] == 1
        @test result_c1.monos[1].word == [1]

        # Single creation c₁† → unchanged
        m_c1_dag = c_dag[1]
        result_c1_dag = simplify(m_c1_dag)
        @test length(result_c1_dag.coeffs) == 1
        @test result_c1_dag.coeffs[1] == 1
        @test result_c1_dag.monos[1].word == [-1]
    end

    @testset "Commutation (CCR)" begin
        # c₁ c₁† = c₁† c₁ + 1 (verify both terms)
        m = c[1] * c_dag[1]  # c₁ c₁†
        result = simplify(m)
        @test length(result.coeffs) == 2

        # Sort terms by degree (identity first, then normal-ordered)
        perm = sortperm(result.monos, by=m -> degree(m))
        sorted_coeffs = result.coeffs[perm]
        sorted_monos = result.monos[perm]

        # First term: identity (delta correction)
        @test sorted_coeffs[1] == 1
        @test isempty(sorted_monos[1].word)

        # Second term: c₁† c₁ (normal-ordered, no sign for bosons)
        @test sorted_coeffs[2] == 1
        @test sorted_monos[2].word == [-1, 1]

        # c₁† c₁ → unchanged (already normal)
        m_normal = c_dag[1] * c[1]
        result_normal = simplify(m_normal)
        @test length(result_normal.coeffs) == 1
        @test result_normal.coeffs[1] == 1
        @test result_normal.monos[1].word == [-1, 1]

        # c₁ c₂ → c₁ c₂ (annihilators commute, just sort by mode)
        m_comm = c[1] * c[2]
        result_comm = simplify(m_comm)
        @test length(result_comm.coeffs) == 1
        @test result_comm.coeffs[1] == 1
        @test result_comm.monos[1].word == [1, 2]

        # c₁† c₂† → c₁† c₂† (creators commute)
        m_cr = c_dag[1] * c_dag[2]
        result_cr = simplify(m_cr)
        @test length(result_cr.coeffs) == 1
        @test result_cr.coeffs[1] == 1
        @test result_cr.monos[1].word == [-1, -2]
    end

    @testset "NOT Nilpotent" begin
        # c₁ c₁ ≠ 0 (bosons are NOT nilpotent)
        m_c1c1 = c[1] * c[1]
        # Note: iszero is only defined for FermionicAlgebra, not BosonicAlgebra
        result_c1c1 = simplify(m_c1c1)
        @test length(result_c1c1.coeffs) == 1
        @test result_c1c1.coeffs[1] == 1
        @test result_c1c1.monos[1].word == [1, 1]  # Stays as c₁ c₁

        # c₁† c₁† ≠ 0
        m_dag_dag = c_dag[1] * c_dag[1]
        result_dag_dag = simplify(m_dag_dag)
        @test length(result_dag_dag.coeffs) == 1
        @test result_dag_dag.coeffs[1] == 1
        @test result_dag_dag.monos[1].word == [-1, -1]
    end

    @testset "Multi-mode" begin
        # c₁ c₂† → c₂† c₁ (different modes, no delta)
        m_cross = c[1] * c_dag[2]
        result_cross = simplify(m_cross)
        @test length(result_cross.coeffs) == 1
        @test result_cross.coeffs[1] == 1
        @test result_cross.monos[1].word == [-2, 1]

        # c₁ c₂ c₁† c₂† → 4 terms
        # Expansion: (c₁ c₁† + 1)(c₂ c₂† + 1) - (direct but with commutations)
        # Normal form: 1 + c₁† c₁ + c₂† c₂ + c₁† c₂† c₁ c₂
        m_two_mode = c[1] * c[2] * c_dag[1] * c_dag[2]
        result_two_mode = simplify(m_two_mode)
        @test length(result_two_mode.coeffs) == 4

        # Check for identity term
        identity_idx = findfirst(m -> isempty(m.word), result_two_mode.monos)
        @test !isnothing(identity_idx)
        @test result_two_mode.coeffs[identity_idx] == 1

        # Check for fully normal-ordered term c₁† c₂† c₁ c₂
        full_term_idx = findfirst(m -> length(m.word) == 4, result_two_mode.monos)
        @test !isnothing(full_term_idx)
    end

    @testset "Rook Number Verification" begin
        # c c c† c† = c†² c² + 4 c† c + 2
        # This is a specific rook number identity
        m_rook = c[1] * c[1] * c_dag[1] * c_dag[1]
        result_rook = simplify(m_rook)

        # Should have 3 terms: c₁†² c₁², c₁† c₁, and identity
        @test length(result_rook.coeffs) == 3

        # Find and verify each term - build dict of degree -> (coeff, mono)
        coeffs_by_degree = Dict(degree(result_rook.monos[i]) => result_rook.coeffs[i]
                                 for i in 1:length(result_rook.coeffs))

        # Identity term (degree 0) with coefficient 2
        @test haskey(coeffs_by_degree, 0)
        @test coeffs_by_degree[0] == 2

        # c₁† c₁ term (degree 2) with coefficient 4
        @test haskey(coeffs_by_degree, 2)
        @test coeffs_by_degree[2] == 4

        # c₁†² c₁² term (degree 4) with coefficient 1
        @test haskey(coeffs_by_degree, 4)
        @test coeffs_by_degree[4] == 1
    end

    @testset "is_normal_ordered helper" begin
        # Empty word → normal ordered
        @test is_normal_ordered(Int32[]) == true

        # Single operator → normal ordered
        @test is_normal_ordered(Int32[1]) == true    # c₁
        @test is_normal_ordered(Int32[-1]) == true   # c₁†

        # Normal order: creators (negative) before annihilators (positive)
        # c₁† c₁ (normal)
        @test is_normal_ordered(Int32[-1, 1]) == true

        # c₁ c₁† (not normal - annihilator before creator)
        @test is_normal_ordered(Int32[1, -1]) == false

        # Multiple creators sorted by mode: c₁† c₂† (normal)
        @test is_normal_ordered(Int32[-1, -2]) == true

        # Multiple creators unsorted: c₂† c₁† (not normal)
        @test is_normal_ordered(Int32[-2, -1]) == false

        # Multiple annihilators sorted by mode: c₁ c₂ (normal)
        @test is_normal_ordered(Int32[1, 2]) == true

        # Multiple annihilators unsorted: c₂ c₁ (not normal)
        @test is_normal_ordered(Int32[2, 1]) == false

        # Full normal form: c₁† c₂† c₁ c₂
        @test is_normal_ordered(Int32[-1, -2, 1, 2]) == true

        # Not normal: c₁† c₁ c₂† c₂ (annihilator before second creator)
        @test is_normal_ordered(Int32[-1, 1, -2, 2]) == false
    end

    @testset "find_first_out_of_order helper" begin
        # Empty word → already normal (returns 0)
        @test find_first_out_of_order(Int32[]) == 0

        # Single operator → already normal (returns 0)
        @test find_first_out_of_order(Int32[1]) == 0

        # Normal order → returns 0
        @test find_first_out_of_order(Int32[-1, 1]) == 0
        @test find_first_out_of_order(Int32[-1, -2, 1, 2]) == 0

        # c₁ c₁† → position 1 is out of order (annihilator followed by creator)
        @test find_first_out_of_order(Int32[1, -1]) == 1

        # c₂† c₁† → position 1 is out of order (creator with higher mode before lower)
        @test find_first_out_of_order(Int32[-2, -1]) == 1

        # c₁† c₂ c₁ → position 2 is out of order (annihilator c₂ before c₁)
        @test find_first_out_of_order(Int32[-1, 2, 1]) == 2

        # c₁† c₁ c₂† c₂ → position 2 is out of order (annihilator before creator)
        @test find_first_out_of_order(Int32[-1, 1, -2, 2]) == 2
    end
end

@testset "Algebra Type Dispatch" begin
    @testset "Type Safety" begin
        # Monomials of different algebra types cannot be multiplied directly
        m_nc = Monomial{NonCommutativeAlgebra}([1])
        m_pauli = Monomial{PauliAlgebra}([1])

        @test typeof(m_nc) != typeof(m_pauli)
    end

    @testset "Polynomial with Algebra Types" begin
        m1 = Monomial{PauliAlgebra}([1])
        m2 = Monomial{PauliAlgebra}([2])

        p = Polynomial([Term(1.0 + 0.0im, m1), Term(2.0 + 0.0im, m2)])

        @test p isa Polynomial{PauliAlgebra}
        @test degree(p) == 1
    end
end

@testset "Immutable Monomial Simplification" begin
    # encode_index is imported at the top of this file
    # NOTE: SimpleAlgebra constructors now auto-canonicalize (sort by site, apply simplifications).
    # So the "original" is already in canonical form after construction.

    @testset "NonCommutative simplification returns new monomial" begin
        # Create monomial - constructor auto-sorts by site
        idx1_s1 = encode_index(UInt16, 1, 1)  # site 1
        idx1_s2 = encode_index(UInt16, 1, 2)  # site 2
        m_nc = Monomial{NonCommutativeAlgebra}(UInt16[idx1_s2, idx1_s1])

        result_nc = simplify(m_nc)

        # Verify different object returned (immutable - simplify always copies)
        @test result_nc !== m_nc

        # Verify result is sorted (site 1 before site 2)
        @test result_nc.word == [idx1_s1, idx1_s2]

        # Original is also sorted (constructor auto-canonicalizes)
        @test m_nc.word == [idx1_s1, idx1_s2]
    end

    @testset "Unipotent simplification with pair cancellation" begin
        # Create monomial - constructor auto-simplifies U² = I
        idx1_s1 = encode_index(UInt16, 1, 1)
        m_uni = Monomial{UnipotentAlgebra}(UInt16[idx1_s1, idx1_s1])

        result_uni = simplify(m_uni)

        # Verify different object returned (immutable - simplify always copies)
        @test result_uni !== m_uni

        # Both are empty (constructor applied U² = I simplification)
        @test isempty(result_uni.word)
        @test isempty(m_uni.word)
    end

    @testset "Unipotent simplification with site-based reordering" begin
        # UnipotentAlgebra: constructor auto-sorts by site
        idx1_s1 = encode_index(UInt16, 1, 1)  # site 1
        idx1_s2 = encode_index(UInt16, 1, 2)  # site 2
        m_uni2 = Monomial{UnipotentAlgebra}(UInt16[idx1_s2, idx1_s1])

        result_uni2 = simplify(m_uni2)

        # Verify different object returned (immutable - simplify always copies)
        @test result_uni2 !== m_uni2

        # Both are site-sorted (constructor auto-canonicalizes)
        @test result_uni2.word == [idx1_s1, idx1_s2]
        @test m_uni2.word == [idx1_s1, idx1_s2]
    end

    @testset "Projector simplification with duplicate removal" begin
        # Create monomial - constructor auto-simplifies P² = P
        idx1_s1 = encode_index(UInt16, 1, 1)
        m_proj = Monomial{ProjectorAlgebra}(UInt16[idx1_s1, idx1_s1, idx1_s1])

        result_proj = simplify(m_proj)

        # Verify different object returned (immutable - simplify always copies)
        @test result_proj !== m_proj

        # Both have single element (constructor applied P³ = P simplification)
        @test length(result_proj.word) == 1
        @test result_proj.word == [idx1_s1]
        @test length(m_proj.word) == 1
    end

    @testset "Projector simplification with reordering" begin
        # Create monomial - constructor auto-sorts by site
        idx1_s1 = encode_index(UInt16, 1, 1)  # site 1
        idx1_s2 = encode_index(UInt16, 1, 2)  # site 2
        m_proj2 = Monomial{ProjectorAlgebra}(UInt16[idx1_s2, idx1_s1])

        result_proj2 = simplify(m_proj2)

        # Verify different object returned (immutable - simplify always copies)
        @test result_proj2 !== m_proj2

        # Both are sorted (constructor auto-canonicalizes)
        @test result_proj2.word == [idx1_s1, idx1_s2]
        @test m_proj2.word == [idx1_s1, idx1_s2]
    end
end

@testset "Index Encoding: integration with simplification" begin
    # NCTSSOS Oracle: constraint reduction functions in src/utils.jl
    # - constraint_reduce_projector!: P² = P → delete one adjacent duplicate
    # - constraint_reduce_unipotent!: U² = I → delete both adjacent duplicates
    # See test/oracles/basis_counts.jl for runnable oracle code.

    @testset "ProjectorAlgebra encoding preserved through simplification" begin
        # NCTSSOS Oracle: constraint_reduce_projector!([1,1,2]) → [1,2]
        # Create encoded indices for projector operators on different sites
        p1_site1 = encode_index(UInt16, 1, 1)  # P₁ on site 1
        p1_site2 = encode_index(UInt16, 1, 2)  # P₁ on site 2

        # Create monomial P₁¹ P₁¹ P₁² (two P₁ on site 1, one on site 2)
        m = Monomial{ProjectorAlgebra}(UInt16[p1_site1, p1_site1, p1_site2])

        # Simplify: P² = P, so P₁¹ P₁¹ → P₁¹
        result = simplify(m)

        # Verify result preserves encoding (simplify returns Monomial for ProjectorAlgebra)
        @test result.word[1] == p1_site1  # First should be site 1
        @test result.word[2] == p1_site2  # Second should be site 2
        @test decode_site(result.word[1]) == 1
        @test decode_site(result.word[2]) == 2
    end

    @testset "UnipotentAlgebra encoding preserved through simplification" begin
        # NCTSSOS Oracle: constraint_reduce_unipotent!([1,1,2]) → [2]
        # Create encoded indices for unipotent operators on different sites
        u1_site1 = encode_index(UInt16, 1, 1)  # U₁ on site 1
        u1_site2 = encode_index(UInt16, 1, 2)  # U₁ on site 2

        # Create monomial U₁¹ U₁¹ U₁² (two U₁ on site 1, one on site 2)
        m = Monomial{UnipotentAlgebra}(UInt16[u1_site1, u1_site1, u1_site2])

        # Simplify: U² = I, so U₁¹ U₁¹ → identity (removed)
        result = simplify(m)

        # After U₁¹ U₁¹ cancels, only U₁² remains (simplify returns Monomial for UnipotentAlgebra)
        @test length(result.word) == 1
        @test result.word[1] == u1_site2
        @test decode_site(result.word[1]) == 2
    end

    @testset "NonCommutativeAlgebra site-based commutation" begin
        # NCTSSoS site-based commutation: operators on different sites commute.
        # Simplify performs stable sort by site (site 1 before site 2).
        # Within same site, order is preserved.
        #
        # NCTSSOS equivalent: _comm(word, partition, comm_var) in src/utils.jl
        # but NCTSSoS encodes site directly in each index.

        # Create encoded indices for operators on different sites
        x1_site1 = encode_index(UInt16, 1, 1)  # op 1 on site 1
        x2_site2 = encode_index(UInt16, 2, 2)  # op 2 on site 2
        x3_site1 = encode_index(UInt16, 3, 1)  # op 3 on site 1

        # Create monomial [op2@site2, op1@site1, op3@site1] (out of site order)
        m = Monomial{NonCommutativeAlgebra}(UInt16[x2_site2, x1_site1, x3_site1])

        # Simplify: stable sort by site → [op1@site1, op3@site1, op2@site2]
        # Site 1 elements [op1, op3] preserve relative order from input
        result = simplify(m)

        # Site 1 operators should come before site 2
        @test decode_site(result.word[1]) == 1
        @test decode_site(result.word[2]) == 1
        @test decode_site(result.word[3]) == 2

        # Verify the actual indices
        @test result.word[1] == x1_site1
        @test result.word[2] == x3_site1
        @test result.word[3] == x2_site2
    end
end

@testset "Mutable simplify! functions" begin
    # simplify! is exported; encode_index is imported at the top of this file

    @testset "NonCommutativeAlgebra simplify!" begin
        idx1_s1 = encode_index(UInt16, 1, 1)
        idx1_s2 = encode_index(UInt16, 1, 2)
        m = Monomial{NonCommutativeAlgebra}(UInt16[idx1_s2, idx1_s1])

        # simplify! mutates and returns the same monomial
        result = simplify!(m)

        @test result === m  # Same object
        @test m.word == [idx1_s1, idx1_s2]  # Mutated in place
    end

    @testset "ProjectorAlgebra simplify!" begin
        idx1_s1 = encode_index(UInt16, 1, 1)
        m = Monomial{ProjectorAlgebra}(UInt16[idx1_s1, idx1_s1, idx1_s1])

        # simplify! mutates and returns the same monomial
        result = simplify!(m)

        @test result === m  # Same object
        @test m.word == [idx1_s1]  # P³ = P, mutated in place
    end

    @testset "UnipotentAlgebra simplify!" begin
        idx1_s1 = encode_index(UInt16, 1, 1)
        m = Monomial{UnipotentAlgebra}(UInt16[idx1_s1, idx1_s1])

        # simplify! mutates and returns the same monomial
        result = simplify!(m)

        @test result === m  # Same object
        @test isempty(m.word)  # U² = I, mutated in place
    end

    @testset "simplify! vs simplify behavior" begin
        idx1_s1 = encode_index(UInt16, 1, 1)
        idx1_s2 = encode_index(UInt16, 1, 2)

        # NOTE: SimpleAlgebra constructors now auto-canonicalize.
        # Both simplify and simplify! operate on already-canonical monomials.

        # NonCommutative: constructor auto-sorts, so both start sorted
        m1 = Monomial{NonCommutativeAlgebra}(UInt16[idx1_s2, idx1_s1])
        m2 = Monomial{NonCommutativeAlgebra}(UInt16[idx1_s2, idx1_s1])

        result1 = simplify(m1)
        result2 = simplify!(m2)

        @test result1 !== m1  # simplify returns new object
        @test result2 === m2  # simplify! returns same object
        @test m1.word == [idx1_s1, idx1_s2]  # Already sorted by constructor
        @test m2.word == [idx1_s1, idx1_s2]  # simplify! on already-sorted is no-op

        # Projector: constructor auto-simplifies P² = P
        m3 = Monomial{ProjectorAlgebra}(UInt16[idx1_s1, idx1_s1])
        m4 = Monomial{ProjectorAlgebra}(UInt16[idx1_s1, idx1_s1])

        result3 = simplify(m3)
        result4 = simplify!(m4)

        @test result3 !== m3
        @test result4 === m4
        @test length(m3.word) == 1  # Already simplified by constructor
        @test length(m4.word) == 1  # simplify! on already-simplified is no-op

        # Unipotent: constructor auto-simplifies U² = I
        m5 = Monomial{UnipotentAlgebra}(UInt16[idx1_s1, idx1_s1])
        m6 = Monomial{UnipotentAlgebra}(UInt16[idx1_s1, idx1_s1])

        result5 = simplify(m5)
        result6 = simplify!(m6)

        @test result5 !== m5
        @test result6 === m6
        @test isempty(m5.word)  # Already simplified by constructor
        @test isempty(m6.word)  # simplify! on already-simplified is no-op
    end
end

@testset "Shared Normal Ordering Helpers (utils.jl)" begin
    @testset "normal_order_key" begin
        # Creation operators (negative) have key starting with 0
        @test normal_order_key(Int32(-1))[1] == 0  # c₁†
        @test normal_order_key(Int32(-3))[1] == 0  # c₃†

        # Annihilation operators (positive) have key starting with 1
        @test normal_order_key(Int32(1))[1] == 1   # c₁
        @test normal_order_key(Int32(2))[1] == 1   # c₂

        # Mode is the second element
        @test normal_order_key(Int32(-3))[2] == 3  # c₃† has mode 3
        @test normal_order_key(Int32(2))[2] == 2   # c₂ has mode 2

        # Sorting by normal_order_key gives normal order
        ops = Int32[2, -1, 1, -2]  # c₂ c₁† c₁ c₂†
        sorted_ops = sort(ops, by=normal_order_key)
        @test sorted_ops == Int32[-1, -2, 1, 2]  # c₁† c₂† c₁ c₂
    end

    @testset "combine_like_terms" begin
        # Test with FermionicAlgebra
        t1 = Term(1.0, Monomial{FermionicAlgebra}(Int32[-1, 1]))
        t2 = Term(2.0, Monomial{FermionicAlgebra}(Int32[-1, 1]))
        t3 = Term(1.0, Monomial{FermionicAlgebra}(Int32[]))

        result = combine_like_terms([t1, t2, t3])
        @test length(result) == 2

        # Find the combined term
        combined_term = findfirst(t -> t.monomial.word == Int32[-1, 1], result)
        @test !isnothing(combined_term)
        @test result[combined_term].coefficient == 3.0  # 1.0 + 2.0

        # Test cancellation
        t4 = Term(1.0, Monomial{FermionicAlgebra}(Int32[1]))
        t5 = Term(-1.0, Monomial{FermionicAlgebra}(Int32[1]))
        result_cancel = combine_like_terms([t4, t5])
        # Should return zero term when everything cancels
        @test length(result_cancel) == 1
        @test result_cancel[1].coefficient == 0.0

        # Test with BosonicAlgebra
        b1 = Term(2.0, Monomial{BosonicAlgebra}(Int32[-1, 1]))
        b2 = Term(3.0, Monomial{BosonicAlgebra}(Int32[-1, 1]))
        result_bos = combine_like_terms([b1, b2])
        @test length(result_bos) == 1
        @test result_bos[1].coefficient == 5.0

        # Test empty input
        empty_result = combine_like_terms(Term{Monomial{FermionicAlgebra,Int32},Float64}[])
        @test length(empty_result) == 1
        @test empty_result[1].coefficient == 0.0
    end
end

# =============================================================================
# NCTSSOS Oracle Comparison Tests
# These tests verify NCTSSoS simplification matches NCTSSOS behavior.
# Oracle reference: /Users/yushengzhao/projects/NCTSSOS/src/utils.jl
# =============================================================================

@testset "NCTSSOS Oracle: constraint_reduce! equivalence" begin
    @testset "UnipotentAlgebra matches NCTSSOS constraint_reduce!(unipotent)" begin
        # NCTSSOS constraint_reduce! with constraint="unipotent" removes consecutive
        # identical pairs with backtracking. NCTSSoS UnipotentAlgebra should match.

        # Reference: constraint_reduce!([1,1], constraint="unipotent") → []
        reg, (U,) = create_unipotent_variables([("U", 1:3)])

        # Case 1: Simple pair cancellation
        m1 = U[1] * U[1]
        s1 = simplify(m1)
        @test isempty(s1.word)  # [1,1] → []

        # Case 2: Multiple consecutive pairs
        m2 = U[1] * U[1] * U[2] * U[2]
        s2 = simplify(m2)
        @test isempty(s2.word)  # [1,1,2,2] → []

        # Case 3: Cascading cancellation
        m3 = U[1] * U[2] * U[2] * U[1]
        s3 = simplify(m3)
        @test isempty(s3.word)  # [1,2,2,1] → []

        # Case 4: Non-consecutive (no cancellation)
        m4 = U[1] * U[2] * U[1] * U[2]
        s4 = simplify(m4)
        @test length(s4.word) == 4  # [1,2,1,2] → [1,2,1,2]

        # Case 5: Mixed
        m5 = U[1] * U[1] * U[2] * U[3] * U[3] * U[2]
        s5 = simplify(m5)
        @test isempty(s5.word)  # [1,1,2,3,3,2] → [] (cascading)
    end

    @testset "ProjectorAlgebra matches NCTSSOS constraint_reduce!(projector)" begin
        # NCTSSOS constraint_reduce! with constraint≠"unipotent" removes consecutive
        # duplicates (P²=P idempotency).

        reg, (P,) = create_projector_variables([("P", 1:3)])

        # Case 1: Simple idempotency
        m1 = P[1] * P[1]
        s1 = simplify(m1)
        @test length(s1.word) == 1  # [1,1] → [1]
        @test s1.word == P[1].word

        # Case 2: Triple idempotency
        m2 = P[1] * P[1] * P[1]
        s2 = simplify(m2)
        @test length(s2.word) == 1  # [1,1,1] → [1]

        # Case 3: Mixed with idempotency
        m3 = P[1] * P[1] * P[2] * P[2] * P[1]
        s3 = simplify(m3)
        @test length(s3.word) == 3  # [1,1,2,2,1] → [1,2,1]

        # Case 4: Non-consecutive (no collapse)
        m4 = P[1] * P[2] * P[1]
        s4 = simplify(m4)
        @test length(s4.word) == 3  # [1,2,1] → [1,2,1]
    end
end

@testset "PauliAlgebra Algebraic Identities" begin
    reg, (σx, σy, σz) = create_pauli_variables(1:2)

    # Helper to get complex phase from PauliMonomial
    _phase(pm::PauliMonomial) = NCTSSoS.phase_to_complex(pm)

    @testset "Involution: σᵢ² = I" begin
        for σ in [σx[1], σy[1], σz[1], σx[2], σy[2], σz[2]]
            pm = simplify(σ * σ)
            @test isempty(pm.mono.word)
            @test _phase(pm) ≈ 1.0 + 0.0im
        end
    end

    @testset "Cyclic products: σₓσᵧ = iσz (same site)" begin
        # XY → iZ
        pm_xy = simplify(σx[1] * σy[1])
        @test _phase(pm_xy) ≈ im
        @test pm_xy.mono.word == σz[1].word

        # YZ → iX
        pm_yz = simplify(σy[1] * σz[1])
        @test _phase(pm_yz) ≈ im
        @test pm_yz.mono.word == σx[1].word

        # ZX → iY
        pm_zx = simplify(σz[1] * σx[1])
        @test _phase(pm_zx) ≈ im
        @test pm_zx.mono.word == σy[1].word
    end

    @testset "Anti-cyclic products: σᵧσₓ = -iσz (same site)" begin
        # YX → -iZ
        pm_yx = simplify(σy[1] * σx[1])
        @test _phase(pm_yx) ≈ -im
        @test pm_yx.mono.word == σz[1].word

        # ZY → -iX
        pm_zy = simplify(σz[1] * σy[1])
        @test _phase(pm_zy) ≈ -im
        @test pm_zy.mono.word == σx[1].word

        # XZ → -iY
        pm_xz = simplify(σx[1] * σz[1])
        @test _phase(pm_xz) ≈ -im
        @test pm_xz.mono.word == σy[1].word
    end

    @testset "Triple product: σₓσᵧσz = i·I" begin
        pm = simplify(σx[1] * σy[1] * σz[1])
        @test isempty(pm.mono.word)
        @test _phase(pm) ≈ im
    end

    @testset "Different sites commute" begin
        # σx₁ σy₂ should just be sorted by site
        pm = simplify(σx[2] * σy[1])
        @test _phase(pm) ≈ 1.0 + 0.0im
        @test length(pm.mono.word) == 2
    end
end

@testset "FermionicAlgebra CAR Identities" begin
    reg, (a, a_dag) = create_fermionic_variables(1:3)

    @testset "Nilpotency: aᵢ² = 0, (aᵢ†)² = 0" begin
        for i in 1:3
            @test iszero(a[i] * a[i])
            @test iszero(a_dag[i] * a_dag[i])
        end
    end

    @testset "Anticommutation: {aᵢ, aⱼ†} = δᵢⱼ" begin
        # Same mode: a₁ a₁† = 1 - a₁† a₁
        p = simplify(a[1] * a_dag[1])
        @test length(p.coeffs) == 2

        identity_idx = findfirst(m -> isempty(m.word), p.monos)
        @test !isnothing(identity_idx)
        @test p.coeffs[identity_idx] == 1

        # Different modes: a₁ a₂† = -a₂† a₁ (no delta)
        p_cross = simplify(a[1] * a_dag[2])
        @test length(p_cross.coeffs) == 1
        @test p_cross.coeffs[1] == -1  # Sign from anticommutation
    end

    @testset "Parity check" begin
        @test has_even_parity(Monomial{FermionicAlgebra}(Int32[]))  # Identity
        @test !has_even_parity(a[1])  # Single operator
        @test has_even_parity(a_dag[1] * a[1])  # Number operator
    end
end

@testset "BosonicAlgebra CCR Identities" begin
    reg, (c, c_dag) = create_bosonic_variables(1:3)

    @testset "Commutation: [cᵢ, cⱼ†] = δᵢⱼ" begin
        # Same mode: c₁ c₁† = c₁† c₁ + 1
        p = simplify(c[1] * c_dag[1])
        @test length(p.coeffs) == 2

        identity_idx = findfirst(m -> isempty(m.word), p.monos)
        @test !isnothing(identity_idx)
        @test p.coeffs[identity_idx] == 1

        normal_idx = findfirst(m -> !isempty(m.word), p.monos)
        @test !isnothing(normal_idx)
        @test p.coeffs[normal_idx] == 1  # No sign for bosons

        # Different modes: c₁ c₂† = c₂† c₁ (no delta)
        p_cross = simplify(c[1] * c_dag[2])
        @test length(p_cross.coeffs) == 1
        @test p_cross.coeffs[1] == 1  # No sign for bosons
    end

    @testset "NOT nilpotent: cᵢ² ≠ 0" begin
        p = simplify(c[1] * c[1])
        @test length(p.coeffs) == 1
        @test p.monos[1].word == Int8[1, 1]
    end

    @testset "Rook number identity: c c c† c† = 2 + 4c†c + c†²c²" begin
        p = simplify(c[1] * c[1] * c_dag[1] * c_dag[1])
        @test length(p.coeffs) == 3

        # Find terms by degree
        identity_idx = findfirst(m -> isempty(m.word), p.monos)
        @test !isnothing(identity_idx)
        @test p.coeffs[identity_idx] == 2

        deg2_idx = findfirst(m -> length(m.word) == 2, p.monos)
        @test !isnothing(deg2_idx)
        @test p.coeffs[deg2_idx] == 4

        deg4_idx = findfirst(m -> length(m.word) == 4, p.monos)
        @test !isnothing(deg4_idx)
        @test p.coeffs[deg4_idx] == 1
    end
end

@testset "Multi-site Commutation (Site-Based Simplification)" begin
    @testset "NonCommutativeAlgebra: different sites commute" begin
        reg, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 3:4)])

        # y₁ x₁ should become x₁ y₁ (site 1 before site 2)
        m = y[1] * x[1]
        s = simplify(m)
        @test decode_site(s.word[1]) == 1
        @test decode_site(s.word[2]) == 2
    end

    @testset "UnipotentAlgebra: different sites commute, then U²=I" begin
        reg, (U, V) = create_unipotent_variables([("U", 1:2), ("V", 3:4)])

        # V₁ U₁ V₁ U₁ should sort by site, then U pairs and V pairs cancel
        m = V[1] * U[1] * V[1] * U[1]
        s = simplify(m)
        @test isempty(s.word)  # After site sort: U₁U₁V₁V₁ → cancels
    end

    @testset "ProjectorAlgebra: different sites commute, then P²=P" begin
        reg, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 3:4)])

        # Q₁ P₁ Q₁ P₁ should sort by site: P₁P₁Q₁Q₁ → P₁Q₁
        m = Q[1] * P[1] * Q[1] * P[1]
        s = simplify(m)
        @test length(s.word) == 2  # P₁ Q₁
        @test decode_site(s.word[1]) == 1
        @test decode_site(s.word[2]) == 2
    end
end

# =============================================================================
# NCTSSOS Oracle: Word-Level Canonicalization
# =============================================================================
#
# Expected values generated from NCTSSOS oracle script (see end of file).
# Use signed Int16 to test generic algorithms without site-aware specializations.
#
# =============================================================================

# NCTSSOS _sym_canon: min(word, reverse(word))
const SYM_CANON_ORACLE = [
    (Int16[1, 2, 3], Int16[1, 2, 3]),
    (Int16[3, 2, 1], Int16[1, 2, 3]),
    (Int16[1, 2, 1], Int16[1, 2, 1]),
    (Int16[2, 1, 2], Int16[2, 1, 2]),
    (Int16[1, 3, 2], Int16[1, 3, 2]),
    (Int16[2, 3, 1], Int16[1, 3, 2]),
    (Int16[], Int16[]),
    (Int16[5], Int16[5]),
    (Int16[1, 2], Int16[1, 2]),
    (Int16[2, 1], Int16[1, 2]),
]

# NCTSSOS _cyclic_canon: lexicographically minimal rotation
const CYCLIC_CANON_ORACLE = [
    (Int16[1, 2, 3], Int16[1, 2, 3]),
    (Int16[2, 3, 1], Int16[1, 2, 3]),
    (Int16[3, 1, 2], Int16[1, 2, 3]),
    (Int16[3, 2, 1], Int16[1, 3, 2]),
    (Int16[2, 1, 3], Int16[1, 3, 2]),
    (Int16[1, 1, 2], Int16[1, 1, 2]),
    (Int16[2, 1, 1], Int16[1, 1, 2]),
    (Int16[], Int16[]),
    (Int16[5], Int16[5]),
    (Int16[1, 2], Int16[1, 2]),
    (Int16[2, 1], Int16[1, 2]),
]

# NCTSSOS min(_cyclic_canon(w), _cyclic_canon(reverse(w)))
const CYCLIC_SYM_CANON_ORACLE = [
    (Int16[1, 2, 3], Int16[1, 2, 3]),
    (Int16[3, 2, 1], Int16[1, 2, 3]),
    (Int16[2, 1, 3], Int16[1, 2, 3]),
    (Int16[1, 3, 2], Int16[1, 2, 3]),
    (Int16[1, 2, 1], Int16[1, 1, 2]),
    (Int16[2, 1, 2], Int16[1, 2, 2]),
    (Int16[1, 2, 3, 4], Int16[1, 2, 3, 4]),
    (Int16[4, 3, 2, 1], Int16[1, 2, 3, 4]),
]

# NCTSSOS constraint_reduce!(unipotent): consecutive pair removal with backtracking
const UNIPOTENT_REDUCE_ORACLE = [
    (Int16[1, 1], Int16[]),
    (Int16[1, 1, 2, 2], Int16[]),
    (Int16[1, 2, 2, 1], Int16[]),
    (Int16[1, 2, 1, 2], Int16[1, 2, 1, 2]),
    (Int16[1, 1, 2], Int16[2]),
    (Int16[1, 2, 2], Int16[1]),
    (Int16[1, 1, 1], Int16[1]),
    (Int16[1, 2, 3, 3, 2, 1], Int16[]),
    (Int16[], Int16[]),
    (Int16[1], Int16[1]),
    (Int16[1, 2, 3], Int16[1, 2, 3]),
]

# NCTSSOS constraint_reduce!(projector): consecutive duplicate removal (P²=P)
const PROJECTOR_REDUCE_ORACLE = [
    (Int16[1, 1], Int16[1]),
    (Int16[1, 1, 1], Int16[1]),
    (Int16[1, 2, 1], Int16[1, 2, 1]),
    (Int16[1, 1, 2, 2], Int16[1, 2]),
    (Int16[1, 1, 2, 2, 1], Int16[1, 2, 1]),
    (Int16[], Int16[]),
    (Int16[1], Int16[1]),
    (Int16[1, 2, 3], Int16[1, 2, 3]),
]

# NCTSSOS get_ncbasis counts: (n, d, count) where count = Σ_{k=0}^{d} n^k
const BASIS_COUNTS_ORACLE = [
    (1, 0, 1), (1, 1, 2), (1, 2, 3), (1, 3, 4),
    (2, 0, 1), (2, 1, 3), (2, 2, 7), (2, 3, 15),
    (3, 0, 1), (3, 1, 4), (3, 2, 13), (3, 3, 40),
    (4, 0, 1), (4, 1, 5), (4, 2, 21), (4, 3, 85),
]

@testset "NCTSSOS Oracle: symmetric_canon (word-level)" begin
    for (input, expected) in SYM_CANON_ORACLE
        @test symmetric_canon(copy(input)) == expected
    end
end

@testset "NCTSSOS Oracle: cyclic_canon (word-level)" begin
    for (input, expected) in CYCLIC_CANON_ORACLE
        @test cyclic_canon(copy(input)) == expected
    end
end

@testset "NCTSSOS Oracle: cyclic_symmetric_canon (word-level)" begin
    for (input, expected) in CYCLIC_SYM_CANON_ORACLE
        @test cyclic_symmetric_canon(copy(input)) == expected
    end
end

@testset "NCTSSOS Oracle: constraint_reduce! algorithm (word-level)" begin
    @testset "Unipotent reduction" begin
        for (input, expected) in UNIPOTENT_REDUCE_ORACLE
            # Apply NCTSSOS algorithm directly on word
            word = copy(input)
            i = 1
            while i < length(word)
                if word[i] == word[i+1]
                    deleteat!(word, i)
                    deleteat!(word, i)
                    i > 1 && (i -= 1)
                else
                    i += 1
                end
            end
            @test word == expected
        end
    end

    @testset "Projector reduction" begin
        for (input, expected) in PROJECTOR_REDUCE_ORACLE
            word = copy(input)
            i = 1
            while i < length(word)
                if word[i] == word[i+1]
                    deleteat!(word, i)
                else
                    i += 1
                end
            end
            @test word == expected
        end
    end
end

@testset "NCTSSOS Oracle: get_ncbasis counts" begin
    using NCTSSoS: decode_operator_id, get_ncbasis_deg

    for (n, d, expected_count) in BASIS_COUNTS_ORACLE
        reg, _ = create_noncommutative_variables([("x", 1:n)])
        basis = get_ncbasis(reg, d)
        @test length(basis) == expected_count
    end

    # Verify degree 0 returns identity
    reg, _ = create_noncommutative_variables([("x", 1:3)])
    basis = get_ncbasis(reg, 0)
    @test length(basis) == 1
    @test isone(basis[1])

    # Verify n=2, d=2 enumerates all 4 words
    reg2, _ = create_noncommutative_variables([("x", 1:2)])
    basis2 = get_ncbasis_deg(reg2, 2)
    @test length(basis2) == 4
    result_words = Set{Vector{UInt16}}()
    for poly in basis2
        for (_, mono) in poly
            push!(result_words, UInt16[decode_operator_id(idx) for idx in mono.word])
        end
    end
    @test result_words == Set([UInt16[1,1], UInt16[1,2], UInt16[2,1], UInt16[2,2]])
end

@testset "NCTSSOS Oracle: Site-aware simplification (two-site)" begin
    @testset "UnipotentAlgebra" begin
        reg, (U, V) = create_unipotent_variables([("U", 1:2), ("V", 1:2)])

        # Same site, same var → cancels
        @test isempty(simplify(U[1] * U[1]).word)

        # Same site, different vars → preserved
        @test length(simplify(U[2] * U[1]).word) == 2

        # Cross-site pairs cancel after site-sorting
        @test isempty(simplify(V[1] * U[1] * V[1] * U[1]).word)
    end

    @testset "ProjectorAlgebra" begin
        reg, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 1:2)])

        # Same site, same var → collapses
        @test length(simplify(P[1] * P[1]).word) == 1

        # Cross-site with idempotency
        s = simplify(Q[1] * P[1] * Q[1] * P[1])
        @test length(s.word) == 2
        @test decode_site(s.word[1]) == 1
        @test decode_site(s.word[2]) == 2
    end
end

# NCTSSOS Oracle: Two-site simplification (4 vars from 2 sites)
#
# NCTSSoS site-based commutation (simplify only, no canon):
#   - Different sites: sorted by site (site 1 before site 2)
#   - Same site: order PRESERVED (stable sort)
#
# NCTSSOS _comm(word, partition, comm_var):
#   - partition=2: vars 1,2 (site 1) placed before vars 3,4 (site 2)
#   - comm_var=2: only vars 1,2 sort among themselves (NOT 3,4)
#
# KEY DIFFERENCE:
#   NCTSSoS preserves within-site order (stable sort by site)
#   NCTSSOS _comm only sorts first comm_var variables

@testset "NCTSSOS Oracle: Two-site simplification (4 vars, 2 sites)" begin
    using NCTSSoS: decode_operator_id, decode_site

    # Create 4 variables on 2 sites:
    #   Site 1: x₁, x₂ (prefix group 1)
    #   Site 2: y₁, y₂ (prefix group 2)
    # Encoded indices are NOT 1,2,3,4 - they contain site info in the bits
    reg, ((x1, x2), (y1, y2)) = create_noncommutative_variables([("x", 1:2), ("y", 1:2)])

    # Verify site assignment
    @test decode_site(x1.word[1]) == 1
    @test decode_site(x2.word[1]) == 1
    @test decode_site(y1.word[1]) == 2
    @test decode_site(y2.word[1]) == 2

    # Test: y₁ x₁ → x₁ y₁ (site 2 moves after site 1)
    m1 = Monomial{NonCommutativeAlgebra}([y1.word[1], x1.word[1]])
    s1 = simplify(m1)
    @test s1.word == [x1.word[1], y1.word[1]]

    # Test: x₂ x₁ → x₂ x₁ (same site, order PRESERVED)
    m2 = Monomial{NonCommutativeAlgebra}([x2.word[1], x1.word[1]])
    s2 = simplify(m2)
    @test s2.word == [x2.word[1], x1.word[1]]  # NOT sorted to [x1, x2]!

    # Test: y₂ y₁ → y₂ y₁ (same site, order PRESERVED)
    m3 = Monomial{NonCommutativeAlgebra}([y2.word[1], y1.word[1]])
    s3 = simplify(m3)
    @test s3.word == [y2.word[1], y1.word[1]]  # NOT sorted to [y1, y2]!

    # Test: y₂ x₂ y₁ x₁ → x₂ x₁ y₂ y₁ (stable sort by site)
    # Site 1 elements [x₂, x₁] preserve relative order
    # Site 2 elements [y₂, y₁] preserve relative order
    m4 = Monomial{NonCommutativeAlgebra}([y2.word[1], x2.word[1], y1.word[1], x1.word[1]])
    s4 = simplify(m4)
    @test s4.word == [x2.word[1], x1.word[1], y2.word[1], y1.word[1]]

    # Test: y₁ x₁ y₂ x₂ → x₁ x₂ y₁ y₂ (stable sort by site)
    m5 = Monomial{NonCommutativeAlgebra}([y1.word[1], x1.word[1], y2.word[1], x2.word[1]])
    s5 = simplify(m5)
    @test s5.word == [x1.word[1], x2.word[1], y1.word[1], y2.word[1]]
end

@testset "NCTSSOS Oracle: Two-site symmetric_canon (4 vars, 2 sites)" begin
    using NCTSSoS: decode_operator_id, decode_site

    # symmetric_canon for Unsigned types:
    #   1. Sort word by site (stable)
    #   2. Sort reverse(word) by site (stable)
    #   3. Return min of the two (lexicographic on encoded indices)

    reg, ((x1, x2), (y1, y2)) = create_noncommutative_variables([("x", 1:2), ("y", 1:2)])

    # Test: x₂ x₁ y₂ y₁ vs reverse y₁ y₂ x₁ x₂
    # After site sort: [x₂,x₁,y₂,y₁] vs [x₁,x₂,y₁,y₂]
    # Min is lexicographically smaller (comparing encoded indices)
    m1 = Monomial{NonCommutativeAlgebra}([x2.word[1], x1.word[1], y2.word[1], y1.word[1]])
    s1 = simplify(m1)
    c1 = symmetric_canon(s1)
    # The reverse [y₁,y₂,x₁,x₂] sorted by site → [x₁,x₂,y₁,y₂]
    # Compare [x₂,x₁,y₂,y₁] vs [x₁,x₂,y₁,y₂] → min is [x₁,x₂,y₁,y₂] since x₁ < x₂
    @test c1.word == [x1.word[1], x2.word[1], y1.word[1], y2.word[1]]

    # Test: already canonical x₁ x₂ y₁ y₂
    m2 = Monomial{NonCommutativeAlgebra}([x1.word[1], x2.word[1], y1.word[1], y2.word[1]])
    s2 = simplify(m2)
    c2 = symmetric_canon(s2)
    @test c2.word == [x1.word[1], x2.word[1], y1.word[1], y2.word[1]]
end

# =============================================================================
# NCTSSOS Oracle Generation Script
# =============================================================================
#
# Run this in NCTSSOS to regenerate expected values:
#
# ```julia
# cd /Users/yushengzhao/projects/NCTSSOS && julia --project -e '
# function _sym_canon(a::Vector{UInt16})
#     i = 1
#     while i <= Int(ceil((length(a)-1)/2))
#         if a[i] < a[end+1-i]; return a
#         elseif a[i] > a[end+1-i]; return reverse(a)
#         else i += 1; end
#     end
#     return a
# end
#
# function _cyclic_canon(a::Vector{UInt16})
#     isempty(a) && return a
#     minimum([[a[i+1:end]; a[1:i]] for i=0:length(a)-1])
# end
#
# function constraint_reduce!(word; constraint="unipotent")
#     i = 1
#     while i < length(word)
#         if word[i] == word[i+1]
#             deleteat!(word, i)
#             if constraint == "unipotent"
#                 deleteat!(word, i)
#                 i > 1 && (i -= 1)
#             end
#         else i += 1; end
#     end
#     return word
# end
#
# function get_ncbasis(n, d)
#     basis = [UInt16[]]
#     for i = 1:d
#         for word in get_ncbasis(n, i-1)[2:end], idx in 1:n
#             push!(basis, [word; idx])
#         end
#         for idx in 1:n; push!(basis, UInt16[idx]); end
#     end
#     return basis
# end
#
# # Generate test data
# for w in [UInt16[1,2,3], UInt16[3,2,1], ...]
#     println("(Int16\$w, Int16\$(_sym_canon(copy(w)))),")
# end
# '
# ```
# =============================================================================

# =============================================================================
# NCTSSOS Oracle: Basis Generation with Constraint Reduction
# =============================================================================
#
# Oracle generation script for basis counts with constraint reduction:
#
# ```julia
# cd /Users/yushengzhao/projects/NCTSSOS && julia --project -e '
# function get_ncbasis(n, d; ind=Vector{UInt16}(1:n))
#     basis = [UInt16[]]
#     for i = 1:d
#         append!(basis, _get_ncbasis_deg(n, i, ind=ind))
#     end
#     return basis
# end
# function _get_ncbasis_deg(n, d; ind=Vector{UInt16}(1:n))
#     if d > 0
#         basis = Vector{UInt16}[]
#         for i = 1:n
#             temp = _get_ncbasis_deg(n, d-1, ind=ind)
#             push!.(temp, ind[i])
#             append!(basis, temp)
#         end
#         return basis
#     else
#         return [UInt16[]]
#     end
# end
# function constraint_reduce_projector!(word)
#     i = 1
#     while i < length(word)
#         if word[i] == word[i+1]; deleteat!(word, i)
#         else i += 1; end
#     end
#     return word
# end
# function constraint_reduce_unipotent!(word)
#     i = 1
#     while i < length(word)
#         if word[i] == word[i+1]
#             deleteat!(word, i); deleteat!(word, i)
#             i > 1 && (i -= 1)
#         else i += 1; end
#     end
#     return word
# end
# for n in 1:4, d in 0:3
#     basis = get_ncbasis(n, d)
#     proj = Set(constraint_reduce_projector!(copy(w)) for w in basis)
#     unip = Set(constraint_reduce_unipotent!(copy(w)) for w in basis)
#     println("(n=$n, d=$d): projector=$(length(proj)), unipotent=$(length(unip))")
# end
# '
# ```
#
# Output:
# (n=1, d=0): proj=1, unip=1
# (n=1, d=1): proj=2, unip=2
# (n=1, d=2): proj=2, unip=2
# (n=1, d=3): proj=2, unip=2
# (n=2, d=0): proj=1, unip=1
# (n=2, d=1): proj=3, unip=3
# (n=2, d=2): proj=5, unip=5
# (n=2, d=3): proj=7, unip=7
# (n=3, d=0): proj=1, unip=1
# (n=3, d=1): proj=4, unip=4
# (n=3, d=2): proj=10, unip=10
# (n=3, d=3): proj=22, unip=22
# (n=4, d=0): proj=1, unip=1
# (n=4, d=1): proj=5, unip=5
# (n=4, d=2): proj=17, unip=17
# (n=4, d=3): proj=53, unip=53
# =============================================================================

# NCTSSOS Oracle: unique word counts after constraint reduction
# (n, d, projector_count, unipotent_count)
# Generated from NCTSSOS oracle script (see above)
const CONSTRAINT_BASIS_COUNTS_ORACLE = [
    (1, 0, 1, 1), (1, 1, 2, 2), (1, 2, 2, 2), (1, 3, 2, 2),
    (2, 0, 1, 1), (2, 1, 3, 3), (2, 2, 5, 5), (2, 3, 7, 7),
    (3, 0, 1, 1), (3, 1, 4, 4), (3, 2, 10, 10), (3, 3, 22, 22),
    (4, 0, 1, 1), (4, 1, 5, 5), (4, 2, 17, 17), (4, 3, 53, 53),
]

@testset "NCTSSOS Oracle: Basis with constraint reduction (single-site)" begin
    using NCTSSoS: decode_operator_id

    @testset "ProjectorAlgebra (P²=P)" begin
        for (n, d, expected_proj, _) in CONSTRAINT_BASIS_COUNTS_ORACLE
            reg, _ = create_projector_variables([("P", 1:n)])
            basis = get_ncbasis(reg, d)
            # Extract unique operator-id words (ignore encoded site)
            words = Set([decode_operator_id.(m.word) for p in basis for m in monomials(p)])
            @test length(words) == expected_proj
        end
    end

    @testset "UnipotentAlgebra (U²=I)" begin
        for (n, d, _, expected_unip) in CONSTRAINT_BASIS_COUNTS_ORACLE
            reg, _ = create_unipotent_variables([("U", 1:n)])
            basis = get_ncbasis(reg, d)
            # Extract unique operator-id words (ignore encoded site)
            words = Set([decode_operator_id.(m.word) for p in basis for m in monomials(p)])
            @test length(words) == expected_unip
        end
    end
end

@testset "NCTSSOS Oracle: Multi-site commutation (NCTSSoS extension)" begin
    # NCTSSoS adds site-based commutation: operators on different sites commute.
    # This is NOT in NCTSSOS. Tests verify this NCTSSoS-specific behavior.

    using NCTSSoS: decode_site

    @testset "ProjectorAlgebra multi-site reduces basis size" begin
        # Single-site: 3 vars, degree 3 → 22 unique words (NCTSSOS oracle)
        reg_single, _ = create_projector_variables([("P", 1:3)])
        basis_single = get_ncbasis(reg_single, 3)
        monos_single = Set(m for p in basis_single for m in monomials(p))
        @test length(monos_single) == 22

        # Multi-site: V on site 1, W₁,W₂ on site 2 → 12 unique (site commutation)
        reg_multi, ((x,), (y, z)) = create_projector_variables([("V", 1:1), ("W", 1:2)])
        basis_multi = get_ncbasis(reg_multi, 3)
        monos_multi = Set(m for p in basis_multi for m in monomials(p))
        @test length(monos_multi) == 12  # Reduced by site-based commutation

        # Verify site-sorting: site 1 always before site 2
        for m in monos_multi
            if length(m.word) >= 2
                sites = decode_site.(m.word)
                # Check that site sequence is non-decreasing
                @test issorted(sites)
            end
        end
    end

    @testset "UnipotentAlgebra multi-site reduces basis size" begin
        # Single-site: 3 vars, degree 3 → 22 unique words (NCTSSOS oracle)
        reg_single, _ = create_unipotent_variables([("U", 1:3)])
        basis_single = get_ncbasis(reg_single, 3)
        monos_single = Set(m for p in basis_single for m in monomials(p))
        @test length(monos_single) == 22

        # Multi-site: V on site 1, W₁,W₂ on site 2 → 12 unique (site commutation)
        reg_multi, ((x,), (y, z)) = create_unipotent_variables([("V", 1:1), ("W", 1:2)])
        basis_multi = get_ncbasis(reg_multi, 3)
        monos_multi = Set(m for p in basis_multi for m in monomials(p))
        @test length(monos_multi) == 12  # Reduced by site-based commutation

        # Verify site-sorting: site 1 always before site 2
        for m in monos_multi
            if length(m.word) >= 2
                sites = decode_site.(m.word)
                @test issorted(sites)
            end
        end
    end
end

# =============================================================================
# Multi-site simplification tests
# =============================================================================

@testset "Multi-site simplification" begin
    # Setup: 2 variables per site, 2 sites
    # Site 1: a, b (a < b)
    # Site 2: c, d (c < d)
    # Raw word: [b, c, a, d] - out of order

    @testset "Pauli multi-site" begin
        # Create variables on 2 sites with 3 ops each (σx, σy, σz)
        reg, (σx, σy, σz) = create_pauli_variables(1:2)

        a, b = σx[1].word[1], σy[1].word[1]  # site 1, a < b
        c, d = σx[2].word[1], σy[2].word[1]  # site 2, c < d

        raw_word = [b, c, a, d]  # out of order across sites

        # Test canonical ordering: stable sort by site
        pm = PauliMonomial(raw_word)

        # Site 1 ops should come before site 2
        # Use Pauli site encoding: site = (idx-1) ÷ 3 + 1
        pauli_site(idx) = (idx - 1) ÷ 3 + 1
        sites = pauli_site.(pm.mono.word)
        @test issorted(sites)

        # Within each site, simplification occurs (σy σx = -i σz)
        # So we should have at most 1 op per site after simplification
        site_counts = Dict{Int,Int}()
        for s in sites
            site_counts[s] = get(site_counts, s, 0) + 1
        end
        for (_, count) in site_counts
            @test count <= 1
        end
    end

    @testset "Fermionic multi-mode" begin
        reg, (a_op,) = create_fermionic_variables(1:4)

        # Get raw indices
        a1 = a_op[1].word[1]  # mode 1 annihilation
        a2 = a_op[2].word[1]  # mode 2 annihilation

        # Different modes anticommute: a₁ a₂ = -a₂ a₁
        # PhysicsMonomial will normal-order these
        pm12 = PhysicsMonomial{FermionicAlgebra}(Int32[a1, a2])
        pm21 = PhysicsMonomial{FermionicAlgebra}(Int32[a2, a1])

        # Both should have one term in normal order
        @test length(pm12.coeffs) == 1
        @test length(pm21.coeffs) == 1

        # pm21 should have opposite sign from pm12
        # because swapping annihilation operators gives -1
        @test pm12.coeffs[1] == -pm21.coeffs[1]
    end

    @testset "Bosonic multi-mode" begin
        reg, (c_op,) = create_bosonic_variables(1:4)

        # Get raw indices
        c1 = c_op[1].word[1]  # mode 1 annihilation
        c2 = c_op[2].word[1]  # mode 2 annihilation

        # Different modes commute: c₁ c₂ = c₂ c₁
        pm12 = PhysicsMonomial{BosonicAlgebra}(Int32[c1, c2])
        pm21 = PhysicsMonomial{BosonicAlgebra}(Int32[c2, c1])

        # Both should be equal (modes commute)
        @test length(pm12.coeffs) == length(pm21.coeffs)
        @test length(pm12.coeffs) == 1
        @test pm12.coeffs[1] == pm21.coeffs[1]
    end

    @testset "NonCommutative multi-site" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:4)])

        a, b = x[1].word[1], x[2].word[1]  # different indices
        c, d = x[3].word[1], x[4].word[1]  # different indices

        raw_word = UInt16[b, c, a, d]  # out of order by site

        # Test: stable sort by site in constructor
        m = Monomial{NonCommutativeAlgebra}(raw_word)
        @test issorted(decode_site.(m.word))
    end
end

# =============================================================================
# Constructor simplification tests
# =============================================================================

@testset "Constructor simplification" begin
    @testset "Unipotent constructor simplification" begin
        reg, (u,) = create_unipotent_variables([("u", 1:2)])
        a, b = u[1].word[1], u[2].word[1]

        # u₁ u₁ should simplify to identity (U² = I)
        m = Monomial{UnipotentAlgebra}(UInt16[a, a])
        @test isone(m)

        # u₁ u₂ u₁ u₂ should simplify to identity
        m2 = Monomial{UnipotentAlgebra}(UInt16[a, b, a, b])
        @test isone(m2)

        # u₁ u₂ u₁ should simplify (cascade)
        m3 = Monomial{UnipotentAlgebra}(UInt16[a, b, a])
        @test m3.word == [b]
    end

    @testset "Projector constructor simplification" begin
        reg, (p,) = create_projector_variables([("p", 1:2)])
        a, b = p[1].word[1], p[2].word[1]

        # p₁ p₁ should simplify to p₁ (P² = P)
        m = Monomial{ProjectorAlgebra}(UInt16[a, a])
        @test m.word == [a]

        # p₁ p₁ p₁ should simplify to p₁
        m2 = Monomial{ProjectorAlgebra}(UInt16[a, a, a])
        @test m2.word == [a]

        # p₁ p₂ p₁ stays as-is after site-sorting (different sites don't simplify directly)
        # But order should be sorted by site
        m3 = Monomial{ProjectorAlgebra}(UInt16[a, b, a])
        @test issorted(decode_site.(m3.word))
    end

    @testset "Pauli constructor rejects non-canonical" begin
        reg, (σx, σy, σz) = create_pauli_variables(1:1)
        x, y = σx[1].word[1], σy[1].word[1]

        # σx₁ σy₁ is NOT canonical (should be iσz₁)
        @test_throws ArgumentError Monomial{PauliAlgebra}([x, y])

        # σx₁ alone IS canonical
        m = Monomial{PauliAlgebra}([x])
        @test m.word == [x]

        # Identity is canonical
        m_id = Monomial{PauliAlgebra}(Int[])
        @test isone(m_id)
    end

    @testset "Fermionic constructor rejects non-normal-ordered" begin
        reg, (a,) = create_fermionic_variables(1:1)
        ann = a[1].word[1]      # annihilation (positive)
        cre = -ann              # creation (negative)

        # a₁ a₁† is NOT normal-ordered (should be 1 - a₁†a₁)
        @test_throws ArgumentError Monomial{FermionicAlgebra}([ann, cre])

        # a₁†a₁ IS normal-ordered
        m = Monomial{FermionicAlgebra}([cre, ann])
        @test m.word == [cre, ann]

        # Identity is normal-ordered
        m_id = Monomial{FermionicAlgebra}(Int32[])
        @test isone(m_id)
    end

    @testset "Bosonic constructor rejects non-normal-ordered" begin
        reg, (c,) = create_bosonic_variables(1:1)
        ann = c[1].word[1]
        cre = -ann

        # c₁ c₁† is NOT normal-ordered
        @test_throws ArgumentError Monomial{BosonicAlgebra}([ann, cre])

        # c₁†c₁ IS normal-ordered
        m = Monomial{BosonicAlgebra}([cre, ann])
        @test m.word == [cre, ann]

        # Identity is normal-ordered
        m_id = Monomial{BosonicAlgebra}(Int32[])
        @test isone(m_id)
    end

    @testset "NonCommutative constructor auto-sorts" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])
        a, b, c = x[1].word[1], x[2].word[1], x[3].word[1]

        # Out of order by site
        raw = UInt16[c, a, b]
        m = Monomial{NonCommutativeAlgebra}(raw)

        # Should be sorted by site
        @test issorted(decode_site.(m.word))
    end
end
