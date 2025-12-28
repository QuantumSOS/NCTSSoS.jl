# NCTSSoS is loaded by parent runtests.jl
# Exported: StateWord, StatePolynomial, ς, tr
# Internal (not exported): NCStateWord, Arbitrary, MaxEntangled, expval, neat_dot
using NCTSSoS: NCStateWord, Arbitrary, MaxEntangled, expval, neat_dot

@testset "StateWord{MaxEntangled} (Tracial)" begin
    @testset "Creation and Equality" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([3])

        sw = StateWord{MaxEntangled}([m1, m2])
        @test sw isa StateWord{MaxEntangled}
        @test length(sw.state_monos) == 2

        # StateWords are sorted by monomial ordering
        sw_sorted = StateWord{MaxEntangled}([m2, m1])
        @test sw == sw_sorted  # Order doesn't matter, they get sorted

        # tr() convenience constructor
        sw_tr = tr(m1)
        @test sw_tr isa StateWord{MaxEntangled}
    end

    @testset "Cyclic Canonicalization for Trace" begin
        # For MaxEntangled (trace) states, cyclic permutations should be equivalent
        # tr(ABC) = tr(BCA) = tr(CAB)
        m_abc = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        m_bca = Monomial{NonCommutativeAlgebra}([2, 3, 1])
        m_cab = Monomial{NonCommutativeAlgebra}([3, 1, 2])

        sw_abc = tr(m_abc)
        sw_bca = tr(m_bca)
        sw_cab = tr(m_cab)

        @test sw_abc == sw_bca
        @test sw_bca == sw_cab
        @test sw_abc == sw_cab
        @test hash(sw_abc) == hash(sw_bca) == hash(sw_cab)

        # Also test with reversal (tr(ABC) = tr(C†B†A†))
        # For NonCommutativeAlgebra, adjoint reverses the word
        m_cba = Monomial{NonCommutativeAlgebra}([3, 2, 1])
        sw_cba = tr(m_cba)
        @test sw_abc == sw_cba

        # Test with longer words
        m_12345 = Monomial{NonCommutativeAlgebra}([1, 2, 3, 4, 5])
        m_23451 = Monomial{NonCommutativeAlgebra}([2, 3, 4, 5, 1])
        m_34512 = Monomial{NonCommutativeAlgebra}([3, 4, 5, 1, 2])

        @test tr(m_12345) == tr(m_23451)
        @test tr(m_12345) == tr(m_34512)
    end

    @testset "Degree" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])  # degree 2
        m2 = Monomial{NonCommutativeAlgebra}([3])     # degree 1

        sw = StateWord{MaxEntangled}([m1, m2])
        @test degree(sw) == 3  # sum of degrees
    end

    @testset "Variables" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([2, 3])

        sw = StateWord{MaxEntangled}([m1, m2])
        vars = variables(sw)

        @test 1 in vars
        @test 2 in vars
        @test 3 in vars
    end

    @testset "One and IsOne" begin
        sw_one = one(StateWord{MaxEntangled,NonCommutativeAlgebra,Int64})
        @test isone(sw_one)

        m = Monomial{NonCommutativeAlgebra}([1])
        sw_not_one = StateWord{MaxEntangled}([m])
        @test !isone(sw_not_one)
    end

    @testset "Multiplication" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        sw1 = StateWord{MaxEntangled}([m1])
        sw2 = StateWord{MaxEntangled}([m2])

        # Multiplication returns StateWord (commutative, no phase)
        result = sw1 * sw2
        @test result isa StateWord
        @test length(result.state_monos) == 2
    end

    @testset "Comparison" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([1, 2])

        sw1 = StateWord{MaxEntangled}([m1])
        sw2 = StateWord{MaxEntangled}([m2])

        # Degree-first ordering
        @test isless(sw1, sw2)

        # Sorting works
        sorted_sws = sort([sw2, sw1])
        @test sorted_sws[1] == sw1
    end

    @testset "Hash and Uniqueness" begin
        m1 = Monomial{NonCommutativeAlgebra}([1, 2])
        m2 = Monomial{NonCommutativeAlgebra}([3])

        sw1 = StateWord{MaxEntangled}([m1, m2])
        sw2 = StateWord{MaxEntangled}([m2, m1])  # Same content, different order

        @test sw1 == sw2
        @test hash(sw1) == hash(sw2)

        @test length(unique([sw1, sw2, sw1])) == 1
    end
end

@testset "StateWord{Arbitrary}" begin
    @testset "Creation with ς" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])

        sw = ς(m)
        @test sw isa StateWord{Arbitrary}
        @test length(sw.state_monos) == 1
    end

    @testset "Multiplication of Arbitrary StateWords" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([2])

        sw1 = ς(m1)
        sw2 = ς(m2)

        # Multiplication returns StateWord (commutative, no phase)
        result = sw1 * sw2
        @test result isa StateWord
        @test length(result.state_monos) == 2
    end

    @testset "Involution (not Cyclic) Canonicalization for Arbitrary" begin
        # For Arbitrary states, only involution canonicalization is applied
        # <M> = <M†> but <ABC> ≠ <BCA> in general

        m_abc = Monomial{NonCommutativeAlgebra}([1, 2, 3])
        m_bca = Monomial{NonCommutativeAlgebra}([2, 3, 1])
        # adjoint([1, 2, 3]) = [-3, -2, -1] for NonCommutativeAlgebra
        m_abc_adj = Monomial{NonCommutativeAlgebra}([-3, -2, -1])

        sw_abc = ς(m_abc)
        sw_bca = ς(m_bca)
        sw_abc_adj = ς(m_abc_adj)

        # Cyclic permutations should NOT be equal for Arbitrary states
        @test sw_abc != sw_bca

        # But m and adjoint(m) should be equal (involution symmetry)
        @test sw_abc == sw_abc_adj
    end
end

@testset "NCStateWord{MaxEntangled}" begin
    @testset "Creation" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])
        m_nc = Monomial{NonCommutativeAlgebra}([3])

        sw = StateWord{MaxEntangled}([m_state])
        ncsw = NCStateWord(sw, m_nc)

        @test ncsw isa NCStateWord{MaxEntangled}
        @test ncsw.sw == sw
        @test ncsw.nc_word == m_nc
    end

    @testset "Degree" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])  # degree 2
        m_nc = Monomial{NonCommutativeAlgebra}([3, 4])      # degree 2

        sw = StateWord{MaxEntangled}([m_state])
        ncsw = NCStateWord(sw, m_nc)

        @test degree(ncsw) == 4  # 2 + 2
    end

    @testset "Variables" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])
        m_nc = Monomial{NonCommutativeAlgebra}([3, 4])

        sw = StateWord{MaxEntangled}([m_state])
        ncsw = NCStateWord(sw, m_nc)

        vars = variables(ncsw)
        @test all(i in vars for i in 1:4)
    end

    @testset "Multiplication" begin
        # Use UInt8 monomials - multiplication is only defined for Unsigned types
        m1 = Monomial{NonCommutativeAlgebra}(UInt8[1])
        m2 = Monomial{NonCommutativeAlgebra}(UInt8[2])

        sw1 = StateWord{MaxEntangled}([m1])
        sw2 = StateWord{MaxEntangled}([m2])

        ncsw1 = NCStateWord(sw1, m1)
        ncsw2 = NCStateWord(sw2, m2)

        result = ncsw1 * ncsw2
        @test result isa NCStateWord
    end

    @testset "Adjoint" begin
        m_state = Monomial{NonCommutativeAlgebra}(UInt8[1, 2])
        m_nc = Monomial{NonCommutativeAlgebra}(UInt8[3, 4])

        sw = StateWord{MaxEntangled}([m_state])
        ncsw = NCStateWord(sw, m_nc)

        ncsw_adj = adjoint(ncsw)
        @test ncsw_adj isa NCStateWord
    end

    @testset "neat_dot" begin
        # Use UInt8 monomials for multiplication support
        m = Monomial{NonCommutativeAlgebra}(UInt8[1])
        sw = StateWord{MaxEntangled}([m])
        ncsw = NCStateWord(sw, m)

        nd = neat_dot(ncsw, ncsw)
        @test nd == adjoint(ncsw) * ncsw
    end

    @testset "expval" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])
        m_nc = Monomial{NonCommutativeAlgebra}([3])

        sw = StateWord{MaxEntangled}([m_state])
        ncsw = NCStateWord(sw, m_nc)

        ev = expval(ncsw)
        @test ev isa StateWord{MaxEntangled}
        @test length(ev.state_monos) == 2  # original + nc_word
    end

    @testset "One and IsOne" begin
        ncsw_one = one(NCStateWord{MaxEntangled,NonCommutativeAlgebra,Int64})
        @test isone(ncsw_one)
    end
end

@testset "NCStateWord{Arbitrary}" begin
    @testset "Creation" begin
        m_state = Monomial{NonCommutativeAlgebra}([1, 2])
        m_nc = Monomial{NonCommutativeAlgebra}([3])

        sw = ς(m_state)
        ncsw = NCStateWord(sw, m_nc)

        @test ncsw isa NCStateWord{Arbitrary}
    end

    @testset "Comparison" begin
        m1 = Monomial{NonCommutativeAlgebra}([1])
        m2 = Monomial{NonCommutativeAlgebra}([1, 2])

        sw1 = ς(m1)
        sw2 = ς(m2)

        ncsw1 = NCStateWord(sw1, m1)
        ncsw2 = NCStateWord(sw2, m2)

        @test isless(ncsw1, ncsw2)  # degree 2 < degree 4
    end

    @testset "Hash and Uniqueness" begin
        m = Monomial{NonCommutativeAlgebra}([1, 2])
        sw = ς(m)
        nc_word = Monomial{NonCommutativeAlgebra}([3])

        ncsw1 = NCStateWord(sw, nc_word)
        ncsw2 = NCStateWord(sw, nc_word)

        @test ncsw1 == ncsw2
        @test hash(ncsw1) == hash(ncsw2)
    end
end
