# Lightweight BPSDP-native SDP serialization helpers.
#
# This is intentionally not a general interchange format. It freezes the exact
# sparse operator data that BPSDP.jl consumes: A, b, c, and the cone block layout.
# Use it for fast tuning/preconditioner loops after NCTSSoS has built and lowered
# a PQG-V2RDM instance once.
module BPSDPNativeIO

using BPSDP
using Serialization
using SparseArrays

const FORMAT_NAME = "bpsdp-native-sparse"
const FORMAT_VERSION = 1

function _operator_sparse_matrix(operator)
    hasfield(typeof(operator), :A) || throw(ArgumentError(
        "BPSDP-native export requires a SparseOperator-like object with field `A`, got $(typeof(operator))",
    ))
    A = getfield(operator, :A)
    A isa SparseMatrixCSC || throw(ArgumentError("operator.A must be a SparseMatrixCSC, got $(typeof(A))"))
    return A
end

function block_kinds(problem::BPSDP.Problem)
    return String[
        block isa BPSDP.HermitianPSDBlock ? "HPSD" :
        block isa BPSDP.PSDBlock ? "PSD" :
        throw(ArgumentError("unsupported BPSDP block type $(typeof(block))"))
        for block in problem.blocks
    ]
end

function block_specs(problem::BPSDP.Problem)
    dims = Int.(problem.layout.dims)
    return collect(zip(block_kinds(problem), dims))
end

function payload(problem::BPSDP.Problem; metadata = Dict{String,Any}())
    A = _operator_sparse_matrix(problem.operator)
    return Dict{Symbol,Any}(
        :format => FORMAT_NAME,
        :format_version => FORMAT_VERSION,
        :storage_type => string(eltype(problem.c)),
        :block_kinds => block_kinds(problem),
        :block_dims => Int.(problem.layout.dims),
        :A => copy(A),
        :b => copy(problem.b),
        :c => copy(problem.c),
        :metadata => metadata,
    )
end

function write_native(path::AbstractString, problem::BPSDP.Problem; metadata = Dict{String,Any}())
    mkpath(dirname(path))
    data = payload(problem; metadata = metadata)
    open(path, "w") do io
        serialize(io, data)
    end
    return path
end

function _block_from_kind(kind, n::Integer)
    k = Symbol(kind)
    if k in (:HPSD, :HermitianPSDBlock, :hermitian)
        return BPSDP.HermitianPSDBlock(n)
    elseif k in (:PSD, :PSDBlock, :real)
        return BPSDP.PSDBlock(n)
    else
        throw(ArgumentError("unsupported block kind $(repr(kind))"))
    end
end

function read_native(path::AbstractString)
    data = open(deserialize, path)
    get(data, :format, nothing) == FORMAT_NAME || throw(ArgumentError("not a $FORMAT_NAME payload: $path"))
    get(data, :format_version, nothing) == FORMAT_VERSION || throw(ArgumentError(
        "unsupported $FORMAT_NAME version $(get(data, :format_version, nothing)); expected $FORMAT_VERSION",
    ))

    kinds = data[:block_kinds]
    dims = data[:block_dims]
    length(kinds) == length(dims) || throw(ArgumentError("block_kinds/block_dims length mismatch"))
    blocks = [_block_from_kind(kind, n) for (kind, n) in zip(kinds, dims)]

    A = data[:A]
    A isa SparseMatrixCSC || throw(ArgumentError("payload A must be SparseMatrixCSC, got $(typeof(A))"))
    problem = BPSDP.Problem(blocks, data[:c], data[:b], BPSDP.SparseOperator(A))
    return problem, data
end

function _dimensions(block_kinds, block_dims, A::SparseMatrixCSC, c)
    return Dict{String,Any}(
        "num_blocks" => length(block_dims),
        "block_kinds" => block_kinds,
        "block_dims" => block_dims,
        "primal_dimension" => size(A, 2),
        "dual_dimension" => size(A, 1),
        "A_nnz" => nnz(A),
        "A_eltype" => string(eltype(A)),
        "c_eltype" => string(eltype(c)),
        "sqrt_dual_dimension" => sqrt(float(size(A, 1))),
    )
end

function dimensions(problem::BPSDP.Problem)
    A = _operator_sparse_matrix(problem.operator)
    return _dimensions(block_kinds(problem), Int.(problem.layout.dims), A, problem.c)
end

function dimensions(data::AbstractDict)
    A = data[:A]
    return _dimensions(data[:block_kinds], data[:block_dims], A, data[:c])
end

end # module
