# NCTSSoS is loaded by parent runtests.jl
using Test, NCTSSoS
using NCTSSoS: encode_index, cyclic_symmetric_canon

const T = UInt16

# Encoded based on operator id & site
o1s1 = encode_index(T, 1, 1)
o2s1 = encode_index(T, 2, 1)
o3s2 = encode_index(T, 3, 2)
o4s2 = encode_index(T, 4, 2)

@testset "Symmetric Canonicalization" begin
    @testset "NoncommutativeAlgebra" begin
        # m.word is already site-sorted: [site 1..., site 2...]
        m = NormalMonomial{NonCommutativeAlgebra,T}(T[o2s1, o1s1, o4s2, o3s2])

        # reverse(m.word) = [o3s2, o4s2, o1s1, o2s1]
        # simplify(NonCommutative): stable sort by site → [o1s1, o2s1, o3s2, o4s2]
        # compare [o2s1,...] vs [o1s1,...] → pick reversed representative
        c = symmetric_canon(m)
        @test c == T[o1s1, o2s1, o3s2, o4s2]
        @test c !== m.word
    end

    @testset "ProjectorAlgebra" begin
        m = NormalMonomial{ProjectorAlgebra,T}(T[o2s1, o1s1, o3s2])
        # reverse = [o3s2,o1s1,o2s1] → site-sort → [o1s1,o2s1,o3s2]
        c = symmetric_canon(m)
        @test c == T[o1s1, o2s1, o3s2]
        @test c !== m.word
    end

    @testset "UnipotentAlgebra" begin
        m = NormalMonomial{UnipotentAlgebra,T}(T[o2s1, o1s1, o3s2])
        # reverse = [o3s2,o1s1,o2s1] → site-sort → [o1s1,o2s1,o3s2]
        c = symmetric_canon(m)
        @test c == T[o1s1, o2s1, o3s2]
        @test c !== m.word
    end

    @testset "PauliAlgebra" begin
        # Pauli uses (idx-1)÷3+1 as site; pick one operator per site.
        m = NormalMonomial{PauliAlgebra,Int}(Int[2, 4, 7])
        c = symmetric_canon(m)
        @test c == m.word
        @test c !== m.word
    end

    @testset "FermionicAlgebra" begin
        # Normal ordered: creators (negative) first, descending by mode; annihilators ascending.
        m = NormalMonomial{FermionicAlgebra,Int32}(Int32[-2, -1, 1, 3])
        c = symmetric_canon(m)
        @test c == m.word
        @test c !== m.word
    end

    @testset "BosonicAlgebra" begin
        m = NormalMonomial{BosonicAlgebra,Int32}(Int32[-2, -2, -1, 1, 1, 3])
        c = symmetric_canon(m)
        @test c == m.word
        @test c !== m.word
    end
end

@testset "Cyclic Canonicalization" begin
    @testset "NoncommutativeAlgebra" begin
        m = NormalMonomial{NonCommutativeAlgebra,T}(T[o2s1, o1s1, o3s2])
        # Rotations (then site-simplify):
        #   [o1s1,o3s2,o2s1] → [o1s1,o2s1,o3s2]
        #   [o3s2,o2s1,o1s1] → [o2s1,o1s1,o3s2]
        # Min is [o1s1,o2s1,o3s2].
        c = cyclic_canon(m)
        @test c == T[o1s1, o2s1, o3s2]
        @test c !== m.word
    end

    @testset "ProjectorAlgebra" begin
        # Chosen so a rotation creates adjacent duplicates after site-sorting, then P²=P shortens.
        m = NormalMonomial{ProjectorAlgebra,T}(T[o2s1, o1s1, o2s1, o3s2])
        # Rotation by 1: [o1s1,o2s1,o3s2,o2s1] → site-sort → [o1s1,o2s1,o2s1,o3s2] → [o1s1,o2s1,o3s2]
        # This is lexicographically smaller than the original (starts with o1s1).
        c = cyclic_canon(m)
        @test c == T[o1s1, o2s1, o3s2]
        @test c !== m.word
    end

    @testset "UnipotentAlgebra" begin
        # Chosen so a rotation creates adjacent duplicates after site-sorting, then U²=I cancels the pair.
        m = NormalMonomial{UnipotentAlgebra,T}(T[o2s1, o1s1, o2s1, o3s2])
        # Rotation by 1: [o1s1,o2s1,o3s2,o2s1] → site-sort → [o1s1,o2s1,o2s1,o3s2] → cancel o2s1,o2s1 → [o1s1,o3s2]
        c = cyclic_canon(m)
        @test c == T[o1s1, o3s2]
        @test c !== m.word

        # Trace + involution must cyclically reduce U_a w U_a to w even when
        # the shorter representative is lexicographically larger. This is the
        # PHBB17 trace-polynomial count bug: tr(U₁U₂U₁) = tr(U₂), not a new
        # degree-3 trace symbol.
        boundary_reducible = NormalMonomial{UnipotentAlgebra,T}(T[o1s1, o2s1, o1s1])
        @test cyclic_canon(boundary_reducible) == T[o2s1]
    end

    @testset "Two commuting sites: $A" for A in (NonCommutativeAlgebra, UnipotentAlgebra, ProjectorAlgebra)
        # tr(u v) = tr(u' v') when u' is a rotation of u on site 1 and v' is a
        # rotation of v on site 2. One rotation of the whole word moves the
        # letters of one site only, so both sites need their own rotation.
        o5s1 = encode_index(T, 5, 1)
        o6s2 = encode_index(T, 6, 2)
        canon(w) = cyclic_canon(NormalMonomial{A,T}(NCTSSoS.simplify!(A, copy(w))))
        @test canon(T[o2s1, o5s1, o1s1, o4s2, o3s2, o6s2]) == canon(T[o1s1, o2s1, o5s1, o4s2, o3s2, o6s2])
        words(letters, n) = n == 0 ? [T[]] : [[w; l] for w in words(letters, n - 1) for l in letters]
        rotations(w) = isempty(w) ? [w] : [circshift(w, k) for k in 0:length(w)-1]
        pairs = [(u, v) for nu in 1:3 for nv in 1:3 for u in words([o1s1, o2s1, o5s1], nu) for v in words([o3s2, o4s2, o6s2], nv)]
        @test all(pairs) do (u, v)
            forms = unique!([canon([ru; rv]) for ru in rotations(u) for rv in rotations(v)])
            length(forms) == 1 && canon(only(forms)) == only(forms)
        end
    end

    @testset "PauliAlgebra" begin
        m = NormalMonomial{PauliAlgebra,Int}(Int[1, 4, 7])
        c = cyclic_canon(m)
        @test c == m.word
        @test c !== m.word
    end

    @testset "FermionicAlgebra" begin
        m = NormalMonomial{FermionicAlgebra,Int32}(Int32[-2, -1, 1, 3])
        c = cyclic_canon(m)
        @test c == m.word
        @test c !== m.word
    end

    @testset "BosonicAlgebra" begin
        m = NormalMonomial{BosonicAlgebra,Int32}(Int32[-2, -2, -1, 1, 1, 3])
        c = cyclic_canon(m)
        @test c == m.word
        @test c !== m.word
    end
end

@testset "Cyclic Symmetric Canonicalization" begin
    @testset "NoncommutativeAlgebra" begin
        m = NormalMonomial{NonCommutativeAlgebra,T}(T[o2s1, o1s1, o4s2, o3s2])
        # cyclic_canon(m) yields [o1s1,o2s1,o4s2,o3s2]
        # adjoint(m) = reverse + simplify → [o1s1,o2s1,o3s2,o4s2]
        # cyclic_canon(adjoint(m)) yields [o1s1,o2s1,o3s2,o4s2]
        # min picks [o1s1,o2s1,o3s2,o4s2].
        c = cyclic_symmetric_canon(m)
        @test c == T[o1s1, o2s1, o3s2, o4s2]
        @test c !== m.word
    end

    @testset "ProjectorAlgebra" begin
        m = NormalMonomial{ProjectorAlgebra,T}(T[o2s1, o1s1, o4s2, o3s2])
        c = cyclic_symmetric_canon(m)
        @test c == T[o1s1, o2s1, o3s2, o4s2]
        @test c !== m.word
    end

    @testset "UnipotentAlgebra" begin
        m = NormalMonomial{UnipotentAlgebra,T}(T[o2s1, o1s1, o4s2, o3s2])
        c = cyclic_symmetric_canon(m)
        @test c == T[o1s1, o2s1, o3s2, o4s2]
        @test c !== m.word
    end

    @testset "PauliAlgebra" begin
        m = NormalMonomial{PauliAlgebra,Int}(Int[2, 4, 7])
        c = cyclic_symmetric_canon(m)
        @test c == m.word
        @test c !== m.word
    end

    @testset "FermionicAlgebra" begin
        m = NormalMonomial{FermionicAlgebra,Int32}(Int32[-2, -1, 1, 3])
        c = cyclic_symmetric_canon(m)
        @test c == m.word
        @test c !== m.word
    end

    @testset "BosonicAlgebra" begin
        m = NormalMonomial{BosonicAlgebra,Int32}(Int32[-2, -2, -1, 1, 1, 3])
        c = cyclic_symmetric_canon(m)
        @test c == m.word
        @test c !== m.word
    end
end
