# Note: FastPolynomials is loaded by setup.jl
using NCTSSoS.FastPolynomials:
    simplify,
    simplify!,
    create_noncommutative_variables,
    create_pauli_variables,
    create_projector_variables,
    create_unipotent_variables,
    create_fermionic_variables,
    create_bosonic_variables,
    get_ncbasis,
    terms,
    monomials

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

        # After simplify!, we get a Monomial (NC algebra returns Monomial)
        simplified = simplify!(result)
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

        # Basis simplification test: verify unique monomials after simplifying basis
        # For 3 projector variables up to degree 3, P_i^2 = P_i reduces consecutive duplicates
        reg3, ((x,), (y, z)) = create_projector_variables([("V", 1:1), ("W", 1:2)])
        basis = get_ncbasis(reg3, 3)

        # Extract unique monomials from the simplified basis
        basis_monos = Set(m for p in basis for m in monomials(p))

        # Expected unique monomials (identity + 3 degree-1 + 6 degree-2 + 2 valid degree-3 combos = 12 total)
        expected = Set([
            one(x),  # identity monomial
            x, y, z,
            x * y, x * z, y * z, z * y,
            x * y * z, x * z * y, y * z * y, z * y * z
        ])
        @test basis_monos == expected
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

        # Basis simplification test: verify unique monomials after simplifying basis
        # For 3 unipotent variables up to degree 3, U_i^2 = I reduces many words
        reg3, ((x,), (y, z)) = create_unipotent_variables([("V", 1:1), ("W", 1:2)])
        basis = get_ncbasis(reg3, 3)

        # Extract unique monomials from the simplified basis
        basis_monos = Set(m for p in basis for m in monomials(p))

        # Expected unique monomials (identity + 3 degree-1 + 6 degree-2 + 2 valid degree-3 combos = 12 total)
        # For unipotent: U_i^2 = I, so consecutive same-index pairs vanish
        expected = Set([
            one(x),  # identity monomial
            x, y, z,
            x * y, x * z, y * z, z * y,
            x * y * z, x * z * y, y * z * y, z * y * z
        ])
        @test basis_monos == expected
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

    @testset "simplify! Mutation" begin
        # Test that simplify! mutates the monomial
        # NonCommutative simplify! returns Monomial
        m = Monomial{NonCommutativeAlgebra}(UInt8[2, 1])  # Will be sorted by site

        result = simplify!(m)
        @test result isa Monomial
    end

    @testset "Star and Simplification Interaction" begin
        # Test star involution on directly created monomials (with proper hash)
        m = Monomial{NonCommutativeAlgebra}(UInt8[5, 9])
        m_star = star(m)

        # star is involutory: star(star(m)).word == m.word
        @test star(star(m)).word == m.word
    end
end

@testset "FermionicAlgebra Simplification" begin
    # FermionicAlgebra uses Wick's theorem with contractions
    # CAR: {aᵢ, aⱼ†} = δᵢⱼ, {aᵢ, aⱼ} = 0, {aᵢ†, aⱼ†} = 0
    # Nilpotency: aᵢ² = 0, (aᵢ†)² = 0
    # Returns Polynomial due to delta corrections

    reg, (a, a_dag) = create_fermionic_variables(1:3)

    @testset "Basic Operations" begin
        # Empty word → identity term
        m_empty = Monomial{FermionicAlgebra}(Int32[])
        result = simplify(m_empty)
        @test result isa Polynomial
        @test length(result.terms) == 1
        @test result.terms[1].coefficient == 1.0
        @test isempty(result.terms[1].monomial.word)

        # Single annihilation a₁ → unchanged
        m_a1 = a[1]
        result_a1 = simplify(m_a1)
        @test length(result_a1.terms) == 1
        @test result_a1.terms[1].coefficient == 1.0
        @test result_a1.terms[1].monomial.word == [1]

        # Single creation a₁† → unchanged
        m_a1_dag = a_dag[1]
        result_a1_dag = simplify(m_a1_dag)
        @test length(result_a1_dag.terms) == 1
        @test result_a1_dag.terms[1].coefficient == 1.0
        @test result_a1_dag.terms[1].monomial.word == [-1]
    end

    @testset "Anticommutation (CAR)" begin
        # a₁ a₁† = 1 - a₁† a₁ (verify both terms)
        m = a[1] * a_dag[1]  # a₁ a₁†
        result = simplify(m)
        @test length(result.terms) == 2

        # Sort terms by degree (identity first, then normal-ordered)
        terms_sorted = sort(result.terms, by=t -> degree(t.monomial))

        # First term: identity (scalar contraction)
        @test terms_sorted[1].coefficient == 1.0
        @test isempty(terms_sorted[1].monomial.word)

        # Second term: -a₁† a₁ (normal-ordered with sign)
        @test terms_sorted[2].coefficient == -1.0
        @test terms_sorted[2].monomial.word == [-1, 1]

        # a₁† a₁ → unchanged (already normal)
        m_normal = a_dag[1] * a[1]
        result_normal = simplify(m_normal)
        @test length(result_normal.terms) == 1
        @test result_normal.terms[1].coefficient == 1.0
        @test result_normal.terms[1].monomial.word == [-1, 1]

        # a₁ a₂ → a₁ a₂ (annihilators already in order by mode)
        m_diff = a[1] * a[2]
        result_diff = simplify(m_diff)
        @test length(result_diff.terms) == 1
        @test result_diff.terms[1].coefficient == 1.0
        @test result_diff.terms[1].monomial.word == [1, 2]  # Sorted by mode

        # a₁† a₂† → a₁† a₂† (creators already in order)
        m_cr = a_dag[1] * a_dag[2]
        result_cr = simplify(m_cr)
        @test length(result_cr.terms) == 1
        @test result_cr.terms[1].coefficient == 1.0
        @test result_cr.terms[1].monomial.word == [-1, -2]
    end

    @testset "Nilpotency" begin
        # a₁ a₁ = 0
        m_a1a1 = a[1] * a[1]
        @test iszero(m_a1a1)
        result_a1a1 = simplify(m_a1a1)
        # Zero polynomial has no terms (zero coefficients are filtered out)
        @test isempty(result_a1a1.terms) || (length(result_a1a1.terms) == 1 && iszero(result_a1a1.terms[1].coefficient))

        # a₁† a₁† = 0
        m_dag_dag = a_dag[1] * a_dag[1]
        @test iszero(m_dag_dag)
        result_dag_dag = simplify(m_dag_dag)
        # Zero polynomial has no terms (zero coefficients are filtered out)
        @test isempty(result_dag_dag.terms) || (length(result_dag_dag.terms) == 1 && iszero(result_dag_dag.terms[1].coefficient))

        # Direct monomial construction
        m_nilp = Monomial{FermionicAlgebra}(Int32[1, 1])
        @test iszero(m_nilp)

        # Non-nilpotent case: a₁ a₁† a₁ a₁† ≠ 0 (alternating)
        m_alt = a[1] * a_dag[1] * a[1] * a_dag[1]
        @test !iszero(m_alt)
    end

    @testset "Multi-mode" begin
        # a₁ a₂† → -a₂† a₁ (different modes, need to swap an and cr)
        m_cross = a[1] * a_dag[2]
        result_cross = simplify(m_cross)
        @test length(result_cross.terms) == 1
        @test result_cross.terms[1].coefficient == -1.0  # Sign from anticommutation
        @test result_cross.terms[1].monomial.word == [-2, 1]

        # a₁ a₁† a₂ a₂† → 4 terms from two contractions
        m_two_mode = a[1] * a_dag[1] * a[2] * a_dag[2]
        result_two_mode = simplify(m_two_mode)
        @test length(result_two_mode.terms) == 4

        # Check total coefficient sum (should be 1 + 1 - 1 - 1 = 0 is wrong; it's more complex)
        # The correct expansion is: (1 - a₁† a₁)(1 - a₂† a₂) = 1 - a₁† a₁ - a₂† a₂ + a₁† a₁ a₂† a₂
        coeffs = [t.coefficient for t in result_two_mode.terms]
        @test 1.0 in coeffs  # Identity term
    end

    @testset "Complex Case" begin
        # a₁ a₂ a₃ a₃† a₂† a₁†
        m_complex = a[1] * a[2] * a[3] * a_dag[3] * a_dag[2] * a_dag[1]
        result_complex = simplify(m_complex)

        # Multiple contraction possibilities should produce multiple terms
        @test length(result_complex.terms) > 1

        # Should include fully contracted term (identity)
        identity_term = findfirst(t -> isempty(t.monomial.word), result_complex.terms)
        @test !isnothing(identity_term)
    end
end

@testset "BosonicAlgebra Simplification" begin
    # BosonicAlgebra uses rook numbers on Ferrers boards
    # CCR: [cᵢ, cⱼ†] = δᵢⱼ, [cᵢ, cⱼ] = 0, [cᵢ†, cⱼ†] = 0
    # NOT nilpotent: cᵢ² ≠ 0
    # Returns Polynomial due to delta corrections

    reg, (c, c_dag) = create_bosonic_variables(1:3)

    @testset "Basic Operations" begin
        # Empty word → identity term
        m_empty = Monomial{BosonicAlgebra}(Int32[])
        result = simplify(m_empty)
        @test result isa Polynomial
        @test length(result.terms) == 1
        @test result.terms[1].coefficient == 1.0
        @test isempty(result.terms[1].monomial.word)

        # Single annihilation c₁ → unchanged
        m_c1 = c[1]
        result_c1 = simplify(m_c1)
        @test length(result_c1.terms) == 1
        @test result_c1.terms[1].coefficient == 1.0
        @test result_c1.terms[1].monomial.word == [1]

        # Single creation c₁† → unchanged
        m_c1_dag = c_dag[1]
        result_c1_dag = simplify(m_c1_dag)
        @test length(result_c1_dag.terms) == 1
        @test result_c1_dag.terms[1].coefficient == 1.0
        @test result_c1_dag.terms[1].monomial.word == [-1]
    end

    @testset "Commutation (CCR)" begin
        # c₁ c₁† = c₁† c₁ + 1 (verify both terms)
        m = c[1] * c_dag[1]  # c₁ c₁†
        result = simplify(m)
        @test length(result.terms) == 2

        # Sort terms by degree (identity first, then normal-ordered)
        terms_sorted = sort(result.terms, by=t -> degree(t.monomial))

        # First term: identity (delta correction)
        @test terms_sorted[1].coefficient == 1.0
        @test isempty(terms_sorted[1].monomial.word)

        # Second term: c₁† c₁ (normal-ordered, no sign for bosons)
        @test terms_sorted[2].coefficient == 1.0
        @test terms_sorted[2].monomial.word == [-1, 1]

        # c₁† c₁ → unchanged (already normal)
        m_normal = c_dag[1] * c[1]
        result_normal = simplify(m_normal)
        @test length(result_normal.terms) == 1
        @test result_normal.terms[1].coefficient == 1.0
        @test result_normal.terms[1].monomial.word == [-1, 1]

        # c₁ c₂ → c₁ c₂ (annihilators commute, just sort by mode)
        m_comm = c[1] * c[2]
        result_comm = simplify(m_comm)
        @test length(result_comm.terms) == 1
        @test result_comm.terms[1].coefficient == 1.0
        @test result_comm.terms[1].monomial.word == [1, 2]

        # c₁† c₂† → c₁† c₂† (creators commute)
        m_cr = c_dag[1] * c_dag[2]
        result_cr = simplify(m_cr)
        @test length(result_cr.terms) == 1
        @test result_cr.terms[1].coefficient == 1.0
        @test result_cr.terms[1].monomial.word == [-1, -2]
    end

    @testset "NOT Nilpotent" begin
        # c₁ c₁ ≠ 0 (bosons are NOT nilpotent)
        m_c1c1 = c[1] * c[1]
        # Note: iszero is only defined for FermionicAlgebra, not BosonicAlgebra
        result_c1c1 = simplify(m_c1c1)
        @test length(result_c1c1.terms) == 1
        @test result_c1c1.terms[1].coefficient == 1.0
        @test result_c1c1.terms[1].monomial.word == [1, 1]  # Stays as c₁ c₁

        # c₁† c₁† ≠ 0
        m_dag_dag = c_dag[1] * c_dag[1]
        result_dag_dag = simplify(m_dag_dag)
        @test length(result_dag_dag.terms) == 1
        @test result_dag_dag.terms[1].coefficient == 1.0
        @test result_dag_dag.terms[1].monomial.word == [-1, -1]
    end

    @testset "Multi-mode" begin
        # c₁ c₂† → c₂† c₁ (different modes, no delta)
        m_cross = c[1] * c_dag[2]
        result_cross = simplify(m_cross)
        @test length(result_cross.terms) == 1
        @test result_cross.terms[1].coefficient == 1.0
        @test result_cross.terms[1].monomial.word == [-2, 1]

        # c₁ c₂ c₁† c₂† → 4 terms
        # Expansion: (c₁ c₁† + 1)(c₂ c₂† + 1) - (direct but with commutations)
        # Normal form: 1 + c₁† c₁ + c₂† c₂ + c₁† c₂† c₁ c₂
        m_two_mode = c[1] * c[2] * c_dag[1] * c_dag[2]
        result_two_mode = simplify(m_two_mode)
        @test length(result_two_mode.terms) == 4

        # Check for identity term
        identity_terms = filter(t -> isempty(t.monomial.word), result_two_mode.terms)
        @test length(identity_terms) == 1
        @test identity_terms[1].coefficient == 1.0

        # Check for fully normal-ordered term c₁† c₂† c₁ c₂
        full_term = findfirst(t -> length(t.monomial.word) == 4, result_two_mode.terms)
        @test !isnothing(full_term)
    end

    @testset "Rook Number Verification" begin
        # c c c† c† = c†² c² + 4 c† c + 2
        # This is a specific rook number identity
        m_rook = c[1] * c[1] * c_dag[1] * c_dag[1]
        result_rook = simplify(m_rook)

        # Should have 3 terms: c₁†² c₁², c₁† c₁, and identity
        @test length(result_rook.terms) == 3

        # Find and verify each term
        terms_by_degree = Dict(degree(t.monomial) => t for t in result_rook.terms)

        # Identity term (degree 0) with coefficient 2
        @test haskey(terms_by_degree, 0)
        @test terms_by_degree[0].coefficient == 2.0

        # c₁† c₁ term (degree 2) with coefficient 4
        @test haskey(terms_by_degree, 2)
        @test terms_by_degree[2].coefficient == 4.0

        # c₁†² c₁² term (degree 4) with coefficient 1
        @test haskey(terms_by_degree, 4)
        @test terms_by_degree[4].coefficient == 1.0
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

@testset "Mutable Monomial Hash Consistency" begin
    using NCTSSoS.FastPolynomials: update_hash!, encode_index

    @testset "NonCommutative in-place mutation" begin
        # Create monomial with indices out of order by site
        idx1_s1 = encode_index(UInt16, 1, 1)  # site 1
        idx1_s2 = encode_index(UInt16, 1, 2)  # site 2
        m_nc = Monomial{NonCommutativeAlgebra}(UInt16[idx1_s2, idx1_s1])
        original_id = objectid(m_nc)

        result_nc = simplify!(m_nc)

        # Verify same object returned
        @test objectid(result_nc) == original_id

        # Verify hash is consistent with word
        @test m_nc.hash == hash(m_nc.word)

        # Verify word was sorted (site 1 before site 2)
        @test m_nc.word == [idx1_s1, idx1_s2]
    end

    @testset "Unipotent in-place mutation with pair cancellation" begin
        # Create monomial that will have pairs cancel: U[1] * U[1] = I
        idx1_s1 = encode_index(UInt16, 1, 1)
        m_uni = Monomial{UnipotentAlgebra}(UInt16[idx1_s1, idx1_s1])
        original_id = objectid(m_uni)

        result_uni = simplify!(m_uni)

        # Verify same object returned
        @test result_uni === m_uni
        @test objectid(result_uni) == original_id

        # Verify hash is consistent with word
        @test m_uni.hash == hash(m_uni.word)

        # Verify word is empty (U^2 = I = empty word)
        @test isempty(m_uni.word)
    end

    @testset "Unipotent in-place mutation with reordering" begin
        # Create monomial with different sites that will be reordered
        idx1_s1 = encode_index(UInt16, 1, 1)  # site 1
        idx1_s2 = encode_index(UInt16, 1, 2)  # site 2
        m_uni2 = Monomial{UnipotentAlgebra}(UInt16[idx1_s2, idx1_s1])
        original_id = objectid(m_uni2)

        result_uni2 = simplify!(m_uni2)

        # Verify same object returned
        @test result_uni2 === m_uni2

        # Verify hash is consistent with word
        @test m_uni2.hash == hash(m_uni2.word)

        # Verify word was sorted (site 1 before site 2)
        @test m_uni2.word == [idx1_s1, idx1_s2]
    end

    @testset "Projector in-place mutation with duplicate removal" begin
        # Create monomial with consecutive duplicates: P[1] * P[1] * P[1] = P[1]
        idx1_s1 = encode_index(UInt16, 1, 1)
        m_proj = Monomial{ProjectorAlgebra}(UInt16[idx1_s1, idx1_s1, idx1_s1])
        original_id = objectid(m_proj)

        result_proj = simplify!(m_proj)

        # Verify same object returned
        @test result_proj === m_proj
        @test objectid(result_proj) == original_id

        # Verify hash is consistent with word
        @test m_proj.hash == hash(m_proj.word)

        # Verify word has single element (P^3 = P)
        @test length(m_proj.word) == 1
        @test m_proj.word == [idx1_s1]
    end

    @testset "Projector in-place mutation with reordering" begin
        # Create monomial with different sites that will be reordered
        idx1_s1 = encode_index(UInt16, 1, 1)  # site 1
        idx1_s2 = encode_index(UInt16, 1, 2)  # site 2
        m_proj2 = Monomial{ProjectorAlgebra}(UInt16[idx1_s2, idx1_s1])
        original_id = objectid(m_proj2)

        result_proj2 = simplify!(m_proj2)

        # Verify same object returned
        @test result_proj2 === m_proj2

        # Verify hash is consistent with word
        @test m_proj2.hash == hash(m_proj2.word)

        # Verify word was sorted (site 1 before site 2)
        @test m_proj2.word == [idx1_s1, idx1_s2]
    end

    @testset "update_hash! function" begin
        # Test the update_hash! helper directly
        m = Monomial{NonCommutativeAlgebra}(UInt16[1, 2])
        old_hash = m.hash

        # Manually mutate the word
        m.word[1], m.word[2] = m.word[2], m.word[1]

        # Hash is now stale
        @test m.hash != hash(m.word)

        # Call update_hash!
        result = update_hash!(m)

        # Verify it returns the same object
        @test result === m

        # Hash is now correct
        @test m.hash == hash(m.word)
    end
end
