# NCTSSoS is loaded by parent runtests.jl
using Test, NCTSSoS
import LinearAlgebra

function _kron_all(mats::AbstractVector{<:AbstractMatrix})
    isempty(mats) && return Matrix{ComplexF64}(LinearAlgebra.I, 1, 1)
    return reduce(LinearAlgebra.kron, mats)
end

# -----------------------------------------------------------------------------
# Pauli oracle (exact 2^n representation)
# -----------------------------------------------------------------------------

const _PAULI_X = ComplexF64[0 1; 1 0]
const _PAULI_Y = ComplexF64[0 -im; im 0]
const _PAULI_Z = ComplexF64[1 0; 0 -1]
const _PAULI_I = Matrix{ComplexF64}(LinearAlgebra.I, 2, 2)

@inline _pauli_site_oracle(idx::Integer) = (Int(idx) - 1) ÷ 3 + 1
@inline _pauli_type_oracle(idx::Integer) = (Int(idx) - 1) % 3  # 0=X, 1=Y, 2=Z

function _pauli_op_oracle(idx::Integer, nsites::Int)
    site = _pauli_site_oracle(idx)
    type = _pauli_type_oracle(idx)
    local_op = type == 0 ? _PAULI_X : type == 1 ? _PAULI_Y : _PAULI_Z

    mats = [_PAULI_I for _ in 1:nsites]
    mats[site] = local_op
    return _kron_all(mats)
end

function _pauli_word_oracle(word::AbstractVector{<:Integer}, nsites::Int)
    dim = 2^nsites
    acc = Matrix{ComplexF64}(LinearAlgebra.I, dim, dim)
    for idx in word
        acc = acc * _pauli_op_oracle(idx, nsites)
    end
    return acc
end

@testset "Matrix Oracles" begin
    @testset "PauliAlgebra: simplify matches matrix model" begin
        # 2 sites -> 4×4 matrices
        nsites = 2
        idxs = 1:(3 * nsites)

        # Exhaustive low-degree verification (fast: dim=4)
        for a in idxs
            for b in idxs
                for c in idxs
                    word = [a, b, c]
                    m = Monomial{PauliAlgebra}(word)
                    t = simplify(m)

                    lhs = _pauli_word_oracle(word, nsites)
                    rhs = ComplexF64(t.coefficient) * _pauli_word_oracle(t.monomial.word, nsites)
                    @test isapprox(lhs, rhs; atol=1e-12, rtol=0)
                end
            end
        end
    end
end

# -----------------------------------------------------------------------------
# Fermionic oracle (exact Jordan–Wigner representation on (C^2)^⊗n)
# -----------------------------------------------------------------------------

const _JW_I = _PAULI_I
const _JW_Z = _PAULI_Z
const _JW_SP = ComplexF64[0 1; 0 0]  # σ+
const _JW_SM = ComplexF64[0 0; 1 0]  # σ-

function _jw_op(mode::Int, kind::Symbol, nmodes::Int)
    mats = Vector{Matrix{ComplexF64}}(undef, nmodes)
    for i in 1:nmodes
        if i < mode
            mats[i] = _JW_Z
        elseif i == mode
            mats[i] = kind === :annihilate ? _JW_SM : _JW_SP
        else
            mats[i] = _JW_I
        end
    end
    return _kron_all(mats)
end

function _fermion_word_oracle(word::AbstractVector{<:Integer}, nmodes::Int)
    dim = 2^nmodes
    acc = Matrix{ComplexF64}(LinearAlgebra.I, dim, dim)
    for op in word
        mode = abs(Int(op))
        kind = op > 0 ? :annihilate : :create
        acc = acc * _jw_op(mode, kind, nmodes)
    end
    return acc
end

function _fermion_poly_oracle(p::Polynomial{FermionicAlgebra,T,C}, nmodes::Int) where {T<:Integer,C<:Number}
    dim = 2^nmodes
    acc = zeros(ComplexF64, dim, dim)
    for t in p.terms
        acc .+= ComplexF64(t.coefficient) * _fermion_word_oracle(t.monomial.word, nmodes)
    end
    return acc
end

@testset "Matrix Oracles" begin
    @testset "FermionicAlgebra: simplify matches Jordan–Wigner model" begin
        nmodes = 3

        # Nilpotency: a₁a₁ = 0, a₁†a₁† = 0
        for word in (Int32[1, 1], Int32[-1, -1])
            m = Monomial{FermionicAlgebra}(word)
            p = simplify(m)
            lhs = _fermion_word_oracle(word, nmodes)
            rhs = _fermion_poly_oracle(p, nmodes)
            @test isapprox(lhs, rhs; atol=1e-12, rtol=0)
        end

        # CAR: a₁ a₁† = 1 - a₁† a₁
        m = Monomial{FermionicAlgebra}(Int32[1, -1])
        p = simplify(m)
        lhs = _fermion_word_oracle(m.word, nmodes)
        rhs = _fermion_poly_oracle(p, nmodes)
        @test isapprox(lhs, rhs; atol=1e-12, rtol=0)

        # Cross-mode: a₁ a₂† = -a₂† a₁
        m = Monomial{FermionicAlgebra}(Int32[1, -2])
        p = simplify(m)
        lhs = _fermion_word_oracle(m.word, nmodes)
        rhs = _fermion_poly_oracle(p, nmodes)
        @test isapprox(lhs, rhs; atol=1e-12, rtol=0)
    end
end
