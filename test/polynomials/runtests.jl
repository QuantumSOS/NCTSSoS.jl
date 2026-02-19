using Test, NCTSSoS

# Shared helpers for NonCommutativeAlgebra tests (site-encoded unsigned indices).
const NC_INDEX_T = UInt16
nc_idx(::Type{T}, op_id::Integer, site::Integer=1) where {T<:Unsigned} =
    NCTSSoS.encode_index(T, op_id, site)
nc_idx(op_id::Integer, site::Integer=1) = nc_idx(NC_INDEX_T, op_id, site)

nc_word(::Type{T}, ids::Integer...; site::Integer=1) where {T<:Unsigned} =
    T[nc_idx(T, i, site) for i in ids]
nc_word(ids::Integer...; site::Integer=1) = nc_word(NC_INDEX_T, ids...; site=site)

@testset "Polynomials" begin
    include("algebra_types.jl")
    include("variables.jl")
    include("monomials.jl")
    include("composed_monomial.jl")
    include("polynomial.jl")
    include("arithmetic.jl")
    include("compare.jl")
    include("canonicalization.jl")
    include("simplify.jl")
    include("matrix_oracles.jl")
    include("basis.jl")
    include("term_sparsity_structural.jl")
    include("state_word.jl")
    include("statepolynomial.jl")
    include("state_basis.jl")
    include("utils.jl")
end
