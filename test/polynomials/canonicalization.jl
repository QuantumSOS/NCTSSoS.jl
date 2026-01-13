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
