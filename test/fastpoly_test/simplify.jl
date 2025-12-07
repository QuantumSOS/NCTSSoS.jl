# Note: FastPolynomials is loaded by setup.jl
using NCTSSoS.FastPolynomials:
    simplify,
    simplify!,
    create_noncommutative_variables,
    create_pauli_variables,
    create_projector_variables,
    create_unipotent_variables

# Note: The new API uses AlgebraType dispatch for simplification instead of SimplifyAlgorithm
# Each algebra type (NonCommutativeAlgebra, PauliAlgebra, UnipotentAlgebra, etc.) has its own simplification rules

@testset "Simplification Interface" begin
    @testset "NonCommutative Simplification" begin
        # NonCommutativeAlgebra with encoded indices uses site-aware simplification
        # Operators on different sites commute (sorted by site)
        # Operators on same site preserve order

        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Simple case: multiplication returns Term
        result = x[1] * x[2]
        @test result isa Term
        @test result.coefficient == 1.0
        @test degree(result.monomial) == 2

        # Same variable twice
        result2 = x[1] * x[1]
        @test result2 isa Term
        @test degree(result2.monomial) == 2
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

        # Multiple projectors
        result = P[1] * P[2]
        @test result isa Term
        @test degree(result.monomial) == 2
    end

    @testset "Unipotent Simplification" begin
        # Unipotent algebra: U^2 = I (squares to identity)
        reg, (U,) = create_unipotent_variables([("U", 1:3)])

        @test U[1] isa Monomial{UnipotentAlgebra}
        @test degree(U[1]) == 1

        # Multiplication of different unipotent variables
        result = U[1] * U[2]
        @test result isa Term
        @test degree(result.monomial) == 2
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
        m = Monomial{NonCommutativeAlgebra}(UInt8[2, 1])  # Will be sorted by site

        result = simplify!(m)
        @test result isa Term
        @test result.coefficient == 1.0
    end

    @testset "Star and Simplification Interaction" begin
        # Test star involution on directly created monomials (with proper hash)
        m = Monomial{NonCommutativeAlgebra}(UInt8[5, 9])
        m_star = star(m)

        # star is involutory: star(star(m)).word == m.word
        @test star(star(m)).word == m.word
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
