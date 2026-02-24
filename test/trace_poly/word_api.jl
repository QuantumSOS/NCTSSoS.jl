# NCTSSoS is loaded by parent runtests.jl
# Exported: StateSymbol, StateWord, StatePolynomial, ς, tr
# Internal (not exported): NCStateWord, Arbitrary, MaxEntangled, expval, neat_dot
using NCTSSoS: NCStateWord, Arbitrary, MaxEntangled, expval, neat_dot

_sym(::Type{ST}, m::NormalMonomial{A,T}) where {ST<:NCTSSoS.StateType,A<:NCTSSoS.AlgebraType,T<:Integer} =
    StateSymbol{ST}(m)

_sw(::Type{ST}, monos::Vector{NormalMonomial{A,T}}) where {ST<:NCTSSoS.StateType,A<:NCTSSoS.AlgebraType,T<:Integer} =
    StateWord{ST,A,T}([_sym(ST, m) for m in monos])

@testset "StateWord{MaxEntangled} (Tracial)" begin
    @testset "Creation and Equality" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        sw = _sw(MaxEntangled, [m1, m2])
        @test sw isa StateWord{MaxEntangled}
        @test length(sw.state_syms) == 2

        # StateWords are sorted by monomial ordering
        sw_sorted = _sw(MaxEntangled, [m2, m1])
        @test sw == sw_sorted  # Order doesn't matter, they get sorted

        # Single expectation symbol
        sym_tr = _sym(MaxEntangled, m1)
        @test sym_tr isa StateSymbol{MaxEntangled}
    end

    @testset "Cyclic Canonicalization for Trace" begin
        # For MaxEntangled (trace) states, cyclic permutations should be equivalent
        # tr(ABC) = tr(BCA) = tr(CAB)
        m_abc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2, 3))
        m_bca = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 3, 1))
        m_cab = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 1, 2))

        sym_abc = _sym(MaxEntangled, m_abc)
        sym_bca = _sym(MaxEntangled, m_bca)
        sym_cab = _sym(MaxEntangled, m_cab)

        @test sym_abc == sym_bca
        @test sym_bca == sym_cab
        @test sym_abc == sym_cab
        @test hash(sym_abc) == hash(sym_bca) == hash(sym_cab)

        # Also test with reversal (tr(ABC) = tr(C†B†A†))
        # For NonCommutativeAlgebra, adjoint reverses the word
        m_cba = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 2, 1))
        sym_cba = _sym(MaxEntangled, m_cba)
        @test sym_abc == sym_cba

        # Test with longer words
        m_12345 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2, 3, 4, 5))
        m_23451 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 3, 4, 5, 1))
        m_34512 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 4, 5, 1, 2))

        @test _sym(MaxEntangled, m_12345) == _sym(MaxEntangled, m_23451)
        @test _sym(MaxEntangled, m_12345) == _sym(MaxEntangled, m_34512)
    end

    @testset "Degree" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))  # degree 2
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))     # degree 1

        sw = _sw(MaxEntangled, [m1, m2])
        @test degree(sw) == 3  # sum of degrees
    end

    @testset "Variables" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 3))

        sw = _sw(MaxEntangled, [m1, m2])
        vars = variables(sw)

        @test nc_idx(1) in vars
        @test nc_idx(2) in vars
        @test nc_idx(3) in vars
    end

    @testset "One and IsOne" begin
        sw_one = one(StateWord{MaxEntangled,NonCommutativeAlgebra,UInt16})
        @test isone(sw_one)

        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        sw_not_one = _sw(MaxEntangled, [m])
        @test !isone(sw_not_one)
    end

    @testset "Multiplication" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))

        sw1 = _sw(MaxEntangled, [m1])
        sw2 = _sw(MaxEntangled, [m2])

        # Multiplication returns StateWord (commutative, no phase)
        result = sw1 * sw2
        @test result isa StateWord
        @test length(result.state_syms) == 2
    end

    @testset "Comparison" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))

        sw1 = _sw(MaxEntangled, [m1])
        sw2 = _sw(MaxEntangled, [m2])

        # Degree-first ordering
        @test isless(sw1, sw2)

        # Sorting works
        sorted_sws = sort([sw2, sw1])
        @test sorted_sws[1] == sw1
    end

    @testset "Hash and Uniqueness" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        sw1 = _sw(MaxEntangled, [m1, m2])
        sw2 = _sw(MaxEntangled, [m2, m1])  # Same content, different order

        @test sw1 == sw2
        @test hash(sw1) == hash(sw2)

        @test length(unique([sw1, sw2, sw1])) == 1
    end
end

@testset "expval" begin
    @testset "NormalMonomial -> StateSymbol{Arbitrary}" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        ev = expval(m)
        @test ev isa StateSymbol{Arbitrary}
        @test ev == StateSymbol{Arbitrary}(m)
    end

    @testset "TwistedGroupAlgebra -> StatePolynomial" begin
        _, (σx, σy, _) = create_pauli_variables(1:1)
        p = σx[1] * σy[1]  # i*σz₁
        ev = expval(p)

        expected = StatePolynomial(coefficients(p), [_sw(Arbitrary, [mono]) for (_, mono) in p.terms])

        @test ev isa StatePolynomial
        @test ev == expected
    end

    @testset "PBWAlgebra -> StatePolynomial" begin
        _, (a, a_dag) = create_fermionic_variables(1:1)
        p = a[1] * a_dag[1]  # 1 - a₁†a₁
        ev = expval(p)

        expected = StatePolynomial(coefficients(p), [_sw(Arbitrary, [mono]) for (_, mono) in p.terms])

        @test ev isa StatePolynomial
        @test ev == expected
    end
end

@testset "StateWord{Arbitrary}" begin
    @testset "Creation with ς" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))

        sym = _sym(Arbitrary, m)
        @test sym isa StateSymbol{Arbitrary}
    end

    @testset "Multiplication of Arbitrary StateWords" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))

        sym1 = _sym(Arbitrary, m1)
        sym2 = _sym(Arbitrary, m2)

        # Multiplication returns StateWord (commutative, no phase)
        result = sym1 * sym2
        @test result isa StateWord
        @test length(result.state_syms) == 2
    end

    @testset "Involution (not Cyclic) Canonicalization for Arbitrary" begin
        # For Arbitrary states (supported real monoid algebras), symmetric canonicalization is used.
        # Cyclic permutations are NOT identified, but reversal is.

        # Odd degree (3): involution NOT applied
        m_abc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2, 3))
        m_bca = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 3, 1))
        m_abc_adj = adjoint(m_abc)

        sym_abc = _sym(Arbitrary, m_abc)
        sym_bca = _sym(Arbitrary, m_bca)
        sym_abc_adj = _sym(Arbitrary, m_abc_adj)

        # Cyclic permutations should NOT be equal for Arbitrary states
        @test sym_abc != sym_bca

        # Reversal is identified under symmetric canonicalization
        @test sym_abc == sym_abc_adj

        # Even degree (2): involution IS applied, so m == adjoint(m)
        # For NonCommutativeAlgebra: adjoint([1, 2]) = [2, 1] (reverse, no negation)
        m_ab = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m_ab_adj = adjoint(m_ab)  # [2, 1]
        sym_ab = _sym(Arbitrary, m_ab)
        sym_ab_adj = _sym(Arbitrary, m_ab_adj)
        @test sym_ab == sym_ab_adj
    end
end

@testset "NCStateWord{MaxEntangled}" begin
    @testset "Creation" begin
        m_state = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        sw = _sw(MaxEntangled, [m_state])
        ncsw = NCStateWord(sw, m_nc)

        @test ncsw isa NCStateWord{MaxEntangled}
        @test ncsw.sw == sw
        @test ncsw.nc_word == m_nc
    end

    @testset "Degree" begin
        m_state = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))  # degree 2
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 4))      # degree 2

        sw = _sw(MaxEntangled, [m_state])
        ncsw = NCStateWord(sw, m_nc)

        @test degree(ncsw) == 4  # 2 + 2
    end

    @testset "Variables" begin
        m_state = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 4))

        sw = _sw(MaxEntangled, [m_state])
        ncsw = NCStateWord(sw, m_nc)

        vars = variables(ncsw)
        @test all(nc_idx(i) in vars for i in 1:4)
    end

    @testset "Multiplication" begin
        # Multiplication is only defined for unsigned (site-encoded) indices.
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2))

        sw1 = _sw(MaxEntangled, [m1])
        sw2 = _sw(MaxEntangled, [m2])

        ncsw1 = NCStateWord(sw1, m1)
        ncsw2 = NCStateWord(sw2, m2)

        result = ncsw1 * ncsw2
        @test result isa NCStateWord
    end

    @testset "Adjoint" begin
        m_state = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 4))

        sw = _sw(MaxEntangled, [m_state])
        ncsw = NCStateWord(sw, m_nc)

        ncsw_adj = adjoint(ncsw)
        @test ncsw_adj isa NCStateWord
    end

    @testset "neat_dot" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        sw = _sw(MaxEntangled, [m])
        ncsw = NCStateWord(sw, m)

        nd = neat_dot(ncsw, ncsw)
        @test nd == adjoint(ncsw) * ncsw
    end

    @testset "expval" begin
        m_state = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        sw = _sw(MaxEntangled, [m_state])
        ncsw = NCStateWord(sw, m_nc)

        ev = expval(ncsw)
        @test ev isa StateWord{MaxEntangled}
        @test length(ev.state_syms) == 2  # original + nc_word
    end

    @testset "One and IsOne" begin
        ncsw_one = one(NCStateWord{MaxEntangled,NonCommutativeAlgebra,UInt16})
        @test isone(ncsw_one)
    end
end

@testset "NCStateWord{Arbitrary}" begin
    @testset "Creation" begin
        m_state = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        sw = _sw(Arbitrary, [m_state])
        ncsw = NCStateWord(sw, m_nc)

        @test ncsw isa NCStateWord{Arbitrary}
    end

    @testset "Comparison" begin
        m1 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1))
        m2 = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))

        sw1 = _sw(Arbitrary, [m1])
        sw2 = _sw(Arbitrary, [m2])

        ncsw1 = NCStateWord(sw1, m1)
        ncsw2 = NCStateWord(sw2, m2)

        @test isless(ncsw1, ncsw2)  # degree 2 < degree 4
    end

    @testset "Hash and Uniqueness" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        sw = _sw(Arbitrary, [m])
        m_nc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3))

        ncsw1 = NCStateWord(sw, m_nc)
        ncsw2 = NCStateWord(sw, m_nc)

        @test ncsw1 == ncsw2
        @test hash(ncsw1) == hash(ncsw2)
    end
end

@testset "StateSymbol Canonicalization" begin
    @testset "StateSymbol{Arbitrary} - Symmetric Canon (real monoid algebras)" begin
        # Symmetric canonicalization identifies reversal for all degrees.
        m_odd = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2, 3))
        m_odd_adj = adjoint(m_odd)
        sym_odd1 = StateSymbol{Arbitrary}(m_odd)
        sym_odd2 = StateSymbol{Arbitrary}(m_odd_adj)
        @test sym_odd1 == sym_odd2

        # Even degree (2): involution IS applied, so <M> = <M†>
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        m_adj = adjoint(m)
        
        sym1 = StateSymbol{Arbitrary}(m)
        sym2 = StateSymbol{Arbitrary}(m_adj)
        
        @test sym1 == sym2
        @test hash(sym1) == hash(sym2)
        
        # Both should canonicalize to the same monomial (min of m, adjoint(m))
        @test sym1.mono == sym2.mono
    end

    @testset "StateSymbol{Arbitrary} - No canon (complex / PBW algebras)" begin
        _, (a, a_dag) = create_fermionic_variables(1:1)
        @test StateSymbol{Arbitrary}(a[1]) != StateSymbol{Arbitrary}(a_dag[1])
    end
    
    @testset "StateSymbol{MaxEntangled} - Cyclic Symmetric Canon" begin
        # For MaxEntangled states: tr(ABC) = tr(BCA) = tr(CAB) = tr(CBA) = ...
        m_abc = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2, 3))
        m_bca = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(2, 3, 1))
        m_cab = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 1, 2))
        m_cba = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(3, 2, 1))  # reverse (adjoint for NC)
        
        sym_abc = StateSymbol{MaxEntangled}(m_abc)
        sym_bca = StateSymbol{MaxEntangled}(m_bca)
        sym_cab = StateSymbol{MaxEntangled}(m_cab)
        sym_cba = StateSymbol{MaxEntangled}(m_cba)
        
        # All cyclic rotations and reversal should be equivalent
        @test sym_abc == sym_bca
        @test sym_bca == sym_cab
        @test sym_abc == sym_cba
        
        # All should have the same canonical monomial
        @test sym_abc.mono == sym_bca.mono
        @test sym_bca.mono == sym_cab.mono
        @test sym_cab.mono == sym_cba.mono
    end
    
    @testset "StateSymbol - Basic Properties" begin
        m = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2))
        sym = StateSymbol{Arbitrary}(m)
        
        # degree
        @test degree(sym) == 2
        
        # variables
        vars = variables(sym)
        @test nc_idx(1) in vars
        @test nc_idx(2) in vars
        
        # one and isone
        sym_one = one(StateSymbol{Arbitrary,NonCommutativeAlgebra,UInt16})
        @test isone(sym_one)
        @test !isone(sym)
        
        # adjoint returns self (due to canonicalization)
        @test adjoint(sym) === sym
        
        # isless (degree-first)
        m_longer = NormalMonomial{NonCommutativeAlgebra,UInt16}(nc_word(1, 2, 3))
        sym_longer = StateSymbol{Arbitrary}(m_longer)
        @test isless(sym, sym_longer)
    end
end
