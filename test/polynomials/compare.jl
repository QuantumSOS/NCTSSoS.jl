# NCTSSoS is loaded by parent runtests.jl
# create_noncommutative_variables is exported from NCTSSoS

@testset "Comparison" begin
    @testset "Monomial Equality" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2, 3))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2, 3))
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 3, 2))

        @test m1 == m2
        @test m1 != m3
    end

    @testset "Monomial Ordering" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m4 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 1))
        m5 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 1))

        # Degree-first ordering
        @test isless(m1, m3)  # degree 1 < degree 2
        @test isless(m2, m3)  # degree 1 < degree 2

        # Same degree: lexicographic
        @test isless(m1, m2)  # [1] < [2]
        @test isless(m4, m3)  # [1,1] < [1,2]
        @test isless(m3, m5)  # [1,2] < [2,1]
    end

    @testset "Monomial Sorting" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))

        sorted_monos = sort([m1, m2, m3])
        @test sorted_monos == [m2, m3, m1]

        # In-set membership
        @test m1 in sorted_monos
        @test m2 in sorted_monos
    end

    @testset "Polynomial Equality" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        p1 = Polynomial([(1.0, m1), (2.0, m2)])
        p2 = Polynomial([(1.0, m1), (2.0, m2)])
        p3 = Polynomial([(1.00000001, m1), (2.0, m2)])

        @test p1 == p2
        @test p1 != p3  # Different coefficient
    end

    @testset "Polynomial Hash" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        p1 = Polynomial([(1.0, m1), (2.0, m2)])
        p2 = Polynomial([(1.0, m1), (2.0, m2)])
        p3 = Polynomial([(1.00000001, m1), (2.0, m2)])

        @test hash(p1) == hash(p2)

        # Unique should work correctly
        p_set = unique!([p1, p2, p1, p3])
        @test length(p_set) == 2
    end

    @testset "Monomial Hash" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m3 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 1))

        @test hash(m1) == hash(m2)
        @test hash(m1) != hash(m3)

        # Hash collisions handled correctly via equality
        @test m1 == m2
        @test m1 != m3
    end

    @testset "Identity Comparisons" begin
        m_identity = one(NormalMonomial{NonCommutativeAlgebra,UInt16})
        m_non_identity = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))

        @test isone(m_identity)
        @test !isone(m_non_identity)

        # Identity should be less than non-identity (degree 0 < degree 1)
        @test isless(m_identity, m_non_identity)
    end

    @testset "Cross-Algebra Type Comparison" begin
        # Monomials of different algebra types should not be equal
        # Use different sites for Pauli: 1,4 = σx on site 1, σx on site 2
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m_pauli = NormalMonomial{PauliAlgebra,Int}([1, 4])

        @test m_nc != m_pauli
    end

    @testset "Variable Registry Comparison" begin
        reg, (x,) = create_noncommutative_variables([("x", 1:5)])

        # Variables with same index should produce equal monomials
        m1 = x[1]
        m2 = x[1]
        m3 = x[2]

        @test m1 == m2
        @test m1 != m3
    end
end
