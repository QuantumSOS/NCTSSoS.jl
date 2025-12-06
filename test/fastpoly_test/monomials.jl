# Note: FastPolynomials is loaded by setup.jl
using .FastPolynomials: star, create_noncommutative_variables

@testset "Monomials" begin
    @testset "Creation" begin
        # Create monomials directly with word representation
        mono1 = Monomial{NonCommutativeAlgebra}([1, 2, 1, 3])
        @test mono1.word == [1, 2, 1, 3]
        @test degree(mono1) == 4

        # Empty monomial (identity)
        mono_identity = Monomial{NonCommutativeAlgebra}(Int[])
        @test isone(mono_identity)
        @test degree(mono_identity) == 0

        # Single element monomial
        mono_single = Monomial{NonCommutativeAlgebra}([5])
        @test mono_single.word == [5]
        @test degree(mono_single) == 1

        # Monomial with zeros should filter them out
        mono_with_zero = Monomial{NonCommutativeAlgebra}([1, 0, 2, 0, 3])
        @test mono_with_zero.word == [1, 2, 3]
        @test degree(mono_with_zero) == 3

        # one(Monomial) should return identity
        mono_one = one(Monomial{NonCommutativeAlgebra,Int64})
        @test isone(mono_one)
    end

    @testset "Degree" begin
        mono1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = Monomial{NonCommutativeAlgebra}(Int[])
        mono3 = Monomial{NonCommutativeAlgebra}([1, 1, 1, 1, 1])

        @test degree(mono1) == 3
        @test degree(mono2) == 0
        @test degree(mono3) == 5
    end

    @testset "Hash" begin
        mono1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = Monomial{NonCommutativeAlgebra}([1, 3, 3])
        mono3 = Monomial{NonCommutativeAlgebra}([1, 2, 3])

        @test hash(mono1) != hash(mono2)
        @test hash(mono1) == hash(mono3)

        mono_empty = Monomial{NonCommutativeAlgebra}(Int[])
        mono_zero_filtered = Monomial{NonCommutativeAlgebra}([0])
        @test hash(mono_empty) == hash(mono_zero_filtered)
    end

    @testset "Star Operation" begin
        # Star reverses word for self-adjoint algebras (unsigned types)
        mono1 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3])
        mono1_star = star(mono1)
        @test mono1_star.word == [3, 2, 1]

        # Empty monomial star should be empty
        mono_empty = Monomial{NonCommutativeAlgebra}(UInt8[])
        @test isone(star(mono_empty))

        # Single element star
        mono_single = Monomial{NonCommutativeAlgebra}(UInt8[5])
        @test star(mono_single).word == [5]

        # Star is involution: star(star(m)) == m
        mono2 = Monomial{NonCommutativeAlgebra}(UInt8[1, 2, 3, 4])
        @test star(star(mono2)) == mono2
    end

    @testset "Multiplication" begin
        # Note: Multiplication behavior depends on algebra type and index encoding
        # For basic NonCommutativeAlgebra with UInt indices, site-aware simplification applies
        reg, (x,) = create_noncommutative_variables([("x", 1:3)])

        # Same variable multiplication produces longer word
        result = x[1] * x[1]
        @test result isa Term
        @test degree(result.monomial) == 2

        # Different variables
        result2 = x[1] * x[2]
        @test degree(result2.monomial) == 2

        # Identity multiplication
        mono_id = one(Monomial{NonCommutativeAlgebra,UInt8})
        mono_x = Monomial{NonCommutativeAlgebra}(UInt8[1])
        result3 = mono_id * mono_x
        @test result3.monomial == mono_x

        result4 = mono_x * mono_id
        @test result4.monomial == mono_x
    end

    @testset "Comparison" begin
        mono1 = Monomial{NonCommutativeAlgebra}([1])
        mono2 = Monomial{NonCommutativeAlgebra}([1, 2])
        mono3 = Monomial{NonCommutativeAlgebra}([2])

        # Degree-first ordering
        @test isless(mono1, mono2)  # degree 1 < degree 2

        # Same degree: lexicographic
        @test isless(mono1, mono3)  # [1] < [2]

        # Sorting works
        monos = [mono2, mono3, mono1]
        sorted_monos = sort(monos)
        @test sorted_monos == [mono1, mono3, mono2]
    end

    @testset "Equality" begin
        mono1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        mono2 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        mono3 = Monomial{NonCommutativeAlgebra}([1, 3, 2])

        @test mono1 == mono2
        @test mono1 != mono3

        # Different algebra types are never equal
        mono_pauli = Monomial{PauliAlgebra}([1, 2, 3])
        mono_nc = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        @test mono_pauli != mono_nc
    end

    @testset "Algebra Type Preservation" begin
        # Monomials preserve their algebra type
        mono_pauli = Monomial{PauliAlgebra}([1, 2])
        mono_fermi = Monomial{FermionicAlgebra}(Int32[1, 2])
        mono_unipotent = Monomial{UnipotentAlgebra}([1, 2])

        @test mono_pauli isa Monomial{PauliAlgebra}
        @test mono_fermi isa Monomial{FermionicAlgebra}
        @test mono_unipotent isa Monomial{UnipotentAlgebra}
    end
end
