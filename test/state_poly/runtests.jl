using Test, NCTSSoS

# Shared helpers for NonCommutativeAlgebra tests (site-encoded unsigned indices).
if !isdefined(@__MODULE__, :NC_INDEX_T)
    const NC_INDEX_T = UInt16
end

if !isdefined(@__MODULE__, :nc_idx)
    nc_idx(::Type{T}, op_id::Integer, site::Integer=1) where {T<:Unsigned} =
        NCTSSoS.encode_index(T, op_id, site)
    nc_idx(op_id::Integer, site::Integer=1) = nc_idx(NC_INDEX_T, op_id, site)
end

if !isdefined(@__MODULE__, :nc_word)
    nc_word(::Type{T}, ids::Integer...; site::Integer=1) where {T<:Unsigned} =
        T[nc_idx(T, i, site) for i in ids]
    nc_word(ids::Integer...; site::Integer=1) = nc_word(NC_INDEX_T, ids...; site=site)
end

@testset "State Polynomials" begin
    include("state_word.jl")
    include("statepolynomial.jl")
    include("state_basis.jl")
    include("polyopt_constructor.jl")
    include("integration_chsh.jl")
    include("coverage_edges.jl")
end
