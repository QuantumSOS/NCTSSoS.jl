# NCTSSoS is loaded by parent runtests.jl
using Test, NCTSSoS
using NCTSSoS: cyclic_symmetric_canon

@testset "Canonicalization" begin
    @testset "symmetric_canon(Vector)" begin
        # Already canonical (ascending order)
        @test symmetric_canon([1, 2, 3]) == [1, 2, 3]

        # Needs reversal (descending becomes ascending)
        @test symmetric_canon([3, 2, 1]) == [1, 2, 3]

        # Palindrome (equal under reversal, returns copy)
        @test symmetric_canon([1, 2, 1]) == [1, 2, 1]
        @test symmetric_canon([1, 2, 3, 2, 1]) == [1, 2, 3, 2, 1]

        # Single element
        @test symmetric_canon([5]) == [5]

        # Empty vector
        @test symmetric_canon(Int[]) == Int[]

        # Equal elements (all same)
        @test symmetric_canon([2, 2, 2]) == [2, 2, 2]

        # First element decides
        @test symmetric_canon([1, 5, 5, 5]) == [1, 5, 5, 5]
        @test symmetric_canon([5, 5, 5, 1]) == [1, 5, 5, 5]

        # Middle element decides
        @test symmetric_canon([2, 1, 2]) == [2, 1, 2]
        @test symmetric_canon([2, 3, 2]) == [2, 3, 2]

        # Even length
        @test symmetric_canon([1, 2, 3, 4]) == [1, 2, 3, 4]
        @test symmetric_canon([4, 3, 2, 1]) == [1, 2, 3, 4]

        # Odd length
        @test symmetric_canon([1, 2, 3, 4, 5]) == [1, 2, 3, 4, 5]
        @test symmetric_canon([5, 4, 3, 2, 1]) == [1, 2, 3, 4, 5]
    end

    @testset "symmetric_canon(Monomial)" begin
        # Already canonical
        m1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        m1_canon = symmetric_canon(m1)
        @test m1_canon.word == [1, 2, 3]
        @test m1_canon isa Monomial{NonCommutativeAlgebra}

        # Needs reversal
        m2 = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        m2_canon = symmetric_canon(m2)
        @test m2_canon.word == [1, 2, 3]

        # Palindrome
        m3 = Monomial{NonCommutativeAlgebra}([1, 2, 1])
        m3_canon = symmetric_canon(m3)
        @test m3_canon.word == [1, 2, 1]

        # Single element
        m4 = Monomial{NonCommutativeAlgebra}([5])
        m4_canon = symmetric_canon(m4)
        @test m4_canon.word == [5]

        # Empty (identity)
        m5 = Monomial{NonCommutativeAlgebra}(Int[])
        m5_canon = symmetric_canon(m5)
        @test m5_canon.word == Int[]
        @test isone(m5_canon)

        # Test with PauliAlgebra (different algebra type)
        mp = Monomial{PauliAlgebra}(UInt16[3, 2, 1])
        mp_canon = symmetric_canon(mp)
        @test mp_canon.word == UInt16[1, 2, 3]
        @test mp_canon isa Monomial{PauliAlgebra}
    end

    @testset "cyclic_canon(Vector)" begin
        # Already canonical (smallest rotation)
        @test cyclic_canon([1, 2, 3]) == [1, 2, 3]

        # One rotation needed
        @test cyclic_canon([2, 3, 1]) == [1, 2, 3]

        # Two rotations needed
        @test cyclic_canon([3, 1, 2]) == [1, 2, 3]

        # All same elements
        @test cyclic_canon([2, 2, 2]) == [2, 2, 2]

        # Single element
        @test cyclic_canon([5]) == [5]

        # Empty vector
        @test cyclic_canon(Int[]) == Int[]

        # Test lexicographic ordering
        @test cyclic_canon([2, 1, 3]) == [1, 3, 2]
        @test cyclic_canon([3, 2, 1]) == [1, 3, 2]

        # Test with repeated patterns
        @test cyclic_canon([2, 1, 2, 1]) == [1, 2, 1, 2]

        # Multiple equal minimums (return first found)
        @test cyclic_canon([1, 2, 1, 2]) == [1, 2, 1, 2]

        # Longer sequence
        @test cyclic_canon([5, 1, 2, 3, 4]) == [1, 2, 3, 4, 5]
        @test cyclic_canon([4, 5, 1, 2, 3]) == [1, 2, 3, 4, 5]
    end

    @testset "cyclic_canon(Monomial)" begin
        # Already canonical
        m1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        m1_canon = cyclic_canon(m1)
        @test m1_canon.word == [1, 2, 3]
        @test m1_canon isa Monomial{NonCommutativeAlgebra}

        # One rotation
        m2 = Monomial{NonCommutativeAlgebra}([2, 3, 1])
        m2_canon = cyclic_canon(m2)
        @test m2_canon.word == [1, 2, 3]

        # Two rotations
        m3 = Monomial{NonCommutativeAlgebra}([3, 1, 2])
        m3_canon = cyclic_canon(m3)
        @test m3_canon.word == [1, 2, 3]

        # Single element
        m4 = Monomial{NonCommutativeAlgebra}([5])
        m4_canon = cyclic_canon(m4)
        @test m4_canon.word == [5]

        # Empty (identity)
        m5 = Monomial{NonCommutativeAlgebra}(Int[])
        m5_canon = cyclic_canon(m5)
        @test m5_canon.word == Int[]
        @test isone(m5_canon)

        # Test with PauliAlgebra
        mp = Monomial{PauliAlgebra}(UInt16[2, 3, 1])
        mp_canon = cyclic_canon(mp)
        @test mp_canon.word == UInt16[1, 2, 3]
        @test mp_canon isa Monomial{PauliAlgebra}
    end

    @testset "cyclic_symmetric_canon(Vector)" begin
        # Simple case: both directions give same result
        @test cyclic_symmetric_canon([3, 2, 1]) == [1, 2, 3]

        # Case from docstring: [3, 1, 4, 2] -> [1, 3, 2, 4]
        # Original: [3, 1, 4, 2]
        #   cyclic_canon([3, 1, 4, 2]) rotations: [3,1,4,2], [1,4,2,3], [4,2,3,1], [2,3,1,4] -> min = [1,4,2,3]
        # Reverse: [2, 4, 1, 3]
        #   cyclic_canon([2, 4, 1, 3]) rotations: [2,4,1,3], [4,1,3,2], [1,3,2,4], [3,2,4,1] -> min = [1,3,2,4]
        # min([1,4,2,3], [1,3,2,4]) = [1,3,2,4]
        @test cyclic_symmetric_canon([3, 1, 4, 2]) == [1, 3, 2, 4]

        # Already canonical
        @test cyclic_symmetric_canon([1, 2, 3]) == [1, 2, 3]

        # Reverse gives smaller result
        @test cyclic_symmetric_canon([2, 1]) == [1, 2]

        # Empty
        @test cyclic_symmetric_canon(Int[]) == Int[]

        # Single element
        @test cyclic_symmetric_canon([5]) == [5]

        # All same
        @test cyclic_symmetric_canon([2, 2, 2]) == [2, 2, 2]

        # Palindromic pattern (cyclic rotation gives [1, 1, 2])
        # [1, 2, 1] rotations: [1,2,1], [2,1,1], [1,1,2] -> min = [1,1,2]
        # reverse([1,2,1]) = [1,2,1], same rotations -> min = [1,1,2]
        @test cyclic_symmetric_canon([1, 2, 1]) == [1, 1, 2]

        # Test case: [3, 2, 4, 1]
        # cyclic_canon([3,2,4,1]) = [1,3,2,4]
        # cyclic_canon(reverse([3,2,4,1])) = cyclic_canon([1,4,2,3]) = [1,4,2,3]
        # min([1,3,2,4], [1,4,2,3]) = [1,3,2,4]
        @test cyclic_symmetric_canon([3, 2, 4, 1]) == [1, 3, 2, 4]
    end

    @testset "cyclic_symmetric_canon(Monomial)" begin
        # Simple case
        m1 = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        m1_canon = cyclic_symmetric_canon(m1)
        @test m1_canon.word == [1, 2, 3]
        @test m1_canon isa Monomial{NonCommutativeAlgebra}

        # Docstring example
        m2 = Monomial{NonCommutativeAlgebra}([3, 1, 4, 2])
        m2_canon = cyclic_symmetric_canon(m2)
        @test m2_canon.word == [1, 3, 2, 4]

        # Already canonical
        m3 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        m3_canon = cyclic_symmetric_canon(m3)
        @test m3_canon.word == [1, 2, 3]

        # Empty (identity)
        m4 = Monomial{NonCommutativeAlgebra}(Int[])
        m4_canon = cyclic_symmetric_canon(m4)
        @test m4_canon.word == Int[]
        @test isone(m4_canon)

        # Test with PauliAlgebra
        mp = Monomial{PauliAlgebra}(UInt16[3, 2, 1])
        mp_canon = cyclic_symmetric_canon(mp)
        @test mp_canon.word == UInt16[1, 2, 3]
        @test mp_canon isa Monomial{PauliAlgebra}
    end

    @testset "canonicalize(Monomial) - symmetric mode" begin
        # Default mode (cyclic=false) uses symmetric_canon
        m1 = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        m1_canon = canonicalize(m1)
        @test m1_canon.word == [1, 2, 3]
        @test m1_canon isa Monomial{NonCommutativeAlgebra}

        # Explicit cyclic=false
        m2 = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        m2_canon = canonicalize(m2; cyclic=false)
        @test m2_canon.word == [1, 2, 3]

        # Test difference between symmetric and cyclic
        # [2, 3, 1]: symmetric_canon([2,3,1]) = [1,3,2] (reverse is smaller)
        m3 = Monomial{NonCommutativeAlgebra}([2, 3, 1])
        m3_sym = canonicalize(m3; cyclic=false)
        @test m3_sym.word == [1, 3, 2]

        # Empty
        m4 = Monomial{NonCommutativeAlgebra}(Int[])
        m4_canon = canonicalize(m4)
        @test isone(m4_canon)
    end

    @testset "canonicalize(Monomial) - cyclic mode" begin
        # cyclic=true uses cyclic_symmetric_canon
        m1 = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        m1_canon = canonicalize(m1; cyclic=true)
        @test m1_canon.word == [1, 2, 3]

        # Test docstring example from cyclic_symmetric_canon
        m2 = Monomial{NonCommutativeAlgebra}([3, 1, 4, 2])
        m2_canon = canonicalize(m2; cyclic=true)
        @test m2_canon.word == [1, 3, 2, 4]

        # Compare cyclic vs symmetric mode
        # [2, 3, 1]:
        #   symmetric: reverse([2,3,1]) = [1,3,2] (smaller) -> [1,3,2]
        #   cyclic: cyclic_canon([2,3,1]) = [1,2,3], cyclic_canon([1,3,2]) = [1,3,2], min = [1,2,3]
        m3 = Monomial{NonCommutativeAlgebra}([2, 3, 1])
        m3_sym = canonicalize(m3; cyclic=false)
        m3_cyc = canonicalize(m3; cyclic=true)
        @test m3_sym.word == [1, 3, 2]
        @test m3_cyc.word == [1, 2, 3]

        # Empty
        m4 = Monomial{NonCommutativeAlgebra}(Int[])
        m4_canon = canonicalize(m4; cyclic=true)
        @test isone(m4_canon)
    end

    @testset "canonicalize(Polynomial) - term combining" begin
        # Test case from docstring: [3,2,1] and [1,2,3] combine
        m1 = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        m2 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        p = Polynomial([Term(1.0, m1), Term(2.0, m2)])

        p_canon = canonicalize(p)

        # Both canonicalize to [1, 2, 3], coefficients should sum
        @test length(terms(p_canon)) == 1
        @test coefficients(p_canon)[1] == 3.0
        @test monomials(p_canon)[1].word == [1, 2, 3]

        # Test with complex coefficients
        p_complex = Polynomial([Term(1.0+1.0im, m1), Term(2.0+0.0im, m2)])
        p_complex_canon = canonicalize(p_complex)
        @test length(terms(p_complex_canon)) == 1
        @test coefficients(p_complex_canon)[1] == 3.0 + 1.0im
    end

    @testset "canonicalize(Polynomial) - term cancellation" begin
        # Terms that cancel to zero
        m1 = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        m2 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        p = Polynomial([Term(2.0, m1), Term(-2.0, m2)])

        p_canon = canonicalize(p)

        # Should result in zero polynomial
        @test iszero(p_canon)
        @test isempty(terms(p_canon))
    end

    @testset "canonicalize(Polynomial) - empty polynomial" begin
        # Empty polynomial
        p_empty = zero(Polynomial{NonCommutativeAlgebra,Int64,Float64})
        p_empty_canon = canonicalize(p_empty)
        @test iszero(p_empty_canon)
        @test isempty(terms(p_empty_canon))
    end

    @testset "canonicalize(Polynomial) - single term" begin
        # Single term
        m = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        p = Polynomial([Term(5.0, m)])

        p_canon = canonicalize(p)

        @test length(terms(p_canon)) == 1
        @test coefficients(p_canon)[1] == 5.0
        @test monomials(p_canon)[1].word == [1, 2, 3]
    end

    @testset "canonicalize(Polynomial) - already canonical" begin
        # Already canonical polynomial
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([1, 3])
        p = Polynomial([Term(1.0, m1), Term(2.0, m2)])

        p_canon = canonicalize(p)

        # Should have same structure (no combining)
        @test length(terms(p_canon)) == 2
        @test monomials(p_canon)[1].word == [1, 2]
        @test monomials(p_canon)[2].word == [1, 3]
    end

    @testset "canonicalize(Polynomial) - cyclic mode" begin
        # Test cyclic mode with rotations
        # [2, 3, 1] rotates to [1, 2, 3]
        # [3, 1, 2] rotates to [1, 2, 3]
        m1 = Monomial{NonCommutativeAlgebra}([2, 3, 1])
        m2 = Monomial{NonCommutativeAlgebra}([3, 1, 2])
        p = Polynomial([Term(1.5, m1), Term(2.5, m2)])

        p_canon = canonicalize(p; cyclic=true)

        # Both should canonicalize to [1, 2, 3] and combine
        @test length(terms(p_canon)) == 1
        @test coefficients(p_canon)[1] == 4.0
        @test monomials(p_canon)[1].word == [1, 2, 3]
    end

    @testset "canonicalize(Polynomial) - multiple term combining" begin
        # Three terms that all canonicalize to the same monomial
        m1 = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        m2 = Monomial{NonCommutativeAlgebra}([3, 2, 1])  # symmetric -> [1,2,3]
        m3 = Monomial{NonCommutativeAlgebra}([1, 2, 3])  # already canonical

        p = Polynomial([Term(1.0, m1), Term(2.0, m2), Term(3.0, m3)])
        p_canon = canonicalize(p)

        @test length(terms(p_canon)) == 1
        @test coefficients(p_canon)[1] == 6.0
        @test monomials(p_canon)[1].word == [1, 2, 3]
    end

    @testset "canonicalize(Polynomial) - partial combining" begin
        # Some terms combine, others don't
        m1 = Monomial{NonCommutativeAlgebra}([3, 2, 1])  # -> [1,2,3]
        m2 = Monomial{NonCommutativeAlgebra}([1, 2, 3])  # -> [1,2,3]
        m3 = Monomial{NonCommutativeAlgebra}([4, 5, 6])  # -> [4,5,6] (distinct)

        p = Polynomial([Term(1.0, m1), Term(2.0, m2), Term(3.0, m3)])
        p_canon = canonicalize(p)

        @test length(terms(p_canon)) == 2
        # Terms should be sorted, so [1,2,3] comes before [4,5,6]
        @test monomials(p_canon)[1].word == [1, 2, 3]
        @test coefficients(p_canon)[1] == 3.0
        @test monomials(p_canon)[2].word == [4, 5, 6]
        @test coefficients(p_canon)[2] == 3.0
    end

    @testset "canonicalize - type preservation" begin
        # Test that algebra types are preserved throughout
        m_pauli = Monomial{PauliAlgebra}(UInt16[3, 2, 1])
        m_nc = Monomial{NonCommutativeAlgebra}([3, 2, 1])

        @test symmetric_canon(m_pauli) isa Monomial{PauliAlgebra}
        @test cyclic_canon(m_pauli) isa Monomial{PauliAlgebra}
        @test cyclic_symmetric_canon(m_pauli) isa Monomial{PauliAlgebra}
        @test canonicalize(m_pauli) isa Monomial{PauliAlgebra}

        @test symmetric_canon(m_nc) isa Monomial{NonCommutativeAlgebra}
        @test cyclic_canon(m_nc) isa Monomial{NonCommutativeAlgebra}
        @test cyclic_symmetric_canon(m_nc) isa Monomial{NonCommutativeAlgebra}
        @test canonicalize(m_nc) isa Monomial{NonCommutativeAlgebra}
    end
end
