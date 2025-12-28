# NCTSSoS is loaded by parent runtests.jl
# Exported: simplify, create_noncommutative_variables, get_ncbasis, coefficients, monomials, indices

# Note: Allocation tests for the new API
# The new implementation uses AlgebraType dispatch instead of SimplifyAlgorithm
# These tests verify that key operations remain allocation-free where expected

@testset "Monomial Allocations" begin
    @testset "Monomial Creation Cached" begin
        # After warmup, monomial operations should minimize allocations
        m1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        m2 = Monomial{NonCommutativeAlgebra}([1, 2, 3])

        # Equal monomials should have same hash (no re-computation needed)
        @test hash(m1) == hash(m2)
    end

    @testset "Degree Computation" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2, 3, 4, 5])

        # Warmup
        degree(m)

        # Degree is just length of word - should be allocation-free
        @test (@allocated degree(m)) == 0
    end

    @testset "isone Check" begin
        m_one = one(Monomial{NonCommutativeAlgebra,Int64})
        m_not_one = Monomial{NonCommutativeAlgebra}([1])

        # Warmup
        isone(m_one)
        isone(m_not_one)

        # isone checks empty word - should be allocation-free
        @test (@allocated isone(m_one)) == 0
        @test (@allocated isone(m_not_one)) == 0
    end
end

@testset "Term Allocations" begin
    @testset "iszero Check" begin
        m = Monomial{NonCommutativeAlgebra}([1])
        t = Term(0.0, m)
        t2 = Term(1.0, m)

        # Warmup
        iszero(t)
        iszero(t2)

        # iszero on Term just checks coefficient - should be allocation-free
        @test (@allocated iszero(t)) == 0
        @test (@allocated iszero(t2)) == 0
    end
end

@testset "Polynomial Accessor Allocations" begin
    m1 = Monomial{NonCommutativeAlgebra}([1])
    m2 = Monomial{NonCommutativeAlgebra}([2])

    p = Polynomial([Term(1.0, m1), Term(2.0, m2)])

    @testset "iszero and isone" begin
        p_zero = zero(typeof(p))
        p_one = one(typeof(p))

        # Warmup
        iszero(p)
        iszero(p_zero)
        isone(p_one)

        # These just check lengths/values
        @test (@allocated iszero(p)) == 0
        @test (@allocated iszero(p_zero)) == 0
    end

    @testset "Degree" begin
        # Warmup
        degree(p)

        # Degree iterates over terms - might have some allocation for iterator
        # but should be minimal
        alloc = @allocated degree(p)
        @test alloc <= 64  # Allow small allocation for iterator
    end
end

@testset "Basis Generation Consistency" begin
    # Test that repeated basis generation produces consistent results (registry-based API)
    reg, (x,) = create_noncommutative_variables([("x", 1:2)])
    basis1 = get_ncbasis(reg, 2)
    basis2 = get_ncbasis(reg, 2)

    @test length(basis1) == length(basis2)
    @test basis1 == basis2
end
