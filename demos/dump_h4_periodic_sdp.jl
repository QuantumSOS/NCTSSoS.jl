#!/usr/bin/env julia
# Dump the H4/Nk=2 periodic V2RDM moment relaxation as standalone SDP data.
#
# This script includes h4_periodic_moment_sos.jl to formulate the symbolic
# NCTSSoS MomentProblem once, then writes solver-native artifacts that can be
# loaded without NCTSSoS:
#
#   * libsdp: real-lifted SDPA sparse format (.dat-s), use procedure=maximize.
#   * BPSDP.jl: serialized SparseMatrixCSC callback problem plus a small loader.
#
# Usage from this repository root:
#   julia --project=. demos/dump_h4_periodic_sdp.jl
#   julia --project=. demos/dump_h4_periodic_sdp.jl --format=bpsdp
#   julia --project=. demos/dump_h4_periodic_sdp.jl --integrals=output/...

using Dates
using Printf
using Serialization
using SparseArrays

include(joinpath(@__DIR__, "h4_periodic_moment_sos.jl"))

const DEFAULT_LIBSDP_OUT = "/Users/exaclior/QuantumSOS/libsdp/examples/python_interface/h4_periodic_moment_sos.dat-s"
const DEFAULT_BPSDP_OUT = "/Users/exaclior/QuantumSOS/BPSDP.jl/examples/h4_periodic_moment_sos_native.jls"
const DEFAULT_BPSDP_LOADER = "/Users/exaclior/QuantumSOS/BPSDP.jl/examples/h4_periodic_moment_sos_native_loader.jl"
const DEFAULT_LIBSDP_RUNNER = "/Users/exaclior/QuantumSOS/libsdp/examples/python_interface/solve_h4_periodic_moment_sos.py"

struct DumpOptions
    formulation_args::Vector{String}
    format::Symbol
    libsdp_out::String
    libsdp_runner::String
    bpsdp_out::String
    bpsdp_loader::String
    metadata_out::String
    atol::Float64
end

mutable struct SparseBuilder{T}
    I::Vector{Int}
    J::Vector{Int}
    V::Vector{T}
    atol::Float64
end

SparseBuilder(::Type{T}; atol::Real = 1e-14) where {T} =
    SparseBuilder{T}(Int[], Int[], T[], Float64(atol))

function add_entry!(builder::SparseBuilder{T}, row::Integer, col::Integer, value) where {T}
    abs(value) <= builder.atol && return nothing
    push!(builder.I, Int(row))
    push!(builder.J, Int(col))
    push!(builder.V, T(value))
    return nothing
end

function parse_dump_options(argv)
    format = :all
    libsdp_out = DEFAULT_LIBSDP_OUT
    libsdp_runner = DEFAULT_LIBSDP_RUNNER
    bpsdp_out = DEFAULT_BPSDP_OUT
    bpsdp_loader = DEFAULT_BPSDP_LOADER
    metadata_out = ""
    atol = 1e-14
    formulation_args = String[]

    for arg in argv
        if startswith(arg, "--format=")
            value = split(arg, "=", limit = 2)[2]
            if value in ("all", "both")
                format = :all
            elseif value == "libsdp"
                format = :libsdp
            elseif value == "bpsdp"
                format = :bpsdp
            else
                throw(ArgumentError("unknown --format=$value; expected all, libsdp, or bpsdp"))
            end
        elseif startswith(arg, "--libsdp-out=")
            libsdp_out = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--libsdp-runner=")
            libsdp_runner = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--bpsdp-out=")
            bpsdp_out = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--bpsdp-loader=")
            bpsdp_loader = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--metadata-out=")
            metadata_out = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--atol=")
            atol = parse(Float64, split(arg, "=", limit = 2)[2])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. demos/dump_h4_periodic_sdp.jl [dump options] [formulation options]\n")
            println("Dump options:")
            println("  --format=all|libsdp|bpsdp   artifacts to write; default all")
            println("  --libsdp-out=PATH           SDPA sparse file for libsdp")
            println("  --libsdp-runner=PATH        small Python loader/runner for libsdp")
            println("  --bpsdp-out=PATH            serialized BPSDP.jl native sparse callback data")
            println("  --bpsdp-loader=PATH         Julia loader for the BPSDP.jl native data")
            println("  --metadata-out=PATH         metadata TOML-ish text; default next to primary artifact")
            println("  --atol=FLOAT                drop coefficients with abs(value) <= atol; default 1e-14")
            println()
            println("Formulation options are forwarded to h4_periodic_moment_sos.jl:")
            println("  --integrals=PATH --blocking=momentum|spin|none")
            exit(0)
        else
            push!(formulation_args, arg)
        end
    end

    if isempty(metadata_out)
        primary = format == :libsdp ? libsdp_out : bpsdp_out
        metadata_out = primary * ".metadata.txt"
    end

    return DumpOptions(formulation_args, format, libsdp_out, libsdp_runner,
        bpsdp_out, bpsdp_loader, metadata_out, atol)
end

function mkparent(path::AbstractString)
    dir = dirname(path)
    !isempty(dir) && mkpath(dir)
    return path
end

function symmetric_basis_data(mp)
    symmetric_basis = NCTSSoS._sorted_symmetric_basis(mp.total_basis)
    sorted_unsym_basis = sort(mp.total_basis)
    basis_to_sym_idx = Dict(
        m => searchsortedfirst(symmetric_basis, NCTSSoS.symmetric_canon(m))
        for m in mp.total_basis
    )
    return symmetric_basis, sorted_unsym_basis, basis_to_sym_idx
end

function objective_rhs(mp, symmetric_basis)
    n_basis = length(symmetric_basis)
    rhs = zeros(Float64, 2 * n_basis)
    for (coef, mono) in mp.objective.terms
        canon_mono = NCTSSoS.symmetric_canon(mono)
        idx = searchsortedfirst(symmetric_basis, canon_mono)
        if idx <= n_basis && symmetric_basis[idx] == canon_mono
            rhs[idx] -= real(coef)
            rhs[n_basis + idx] -= imag(coef)
        end
    end
    return rhs
end

function block_offsets(block_dims::AbstractVector{<:Integer})
    offsets = Vector{Int}(undef, length(block_dims))
    offset = 1
    for (i, dim) in pairs(block_dims)
        offsets[i] = offset
        offset += Int(dim)^2
    end
    return offsets, offset - 1
end

@inline function flat_index(offsets, block_dims, block::Integer, row::Integer, col::Integer)
    n = Int(block_dims[block])
    return offsets[block] + (Int(col) - 1) * n + (Int(row) - 1)
end

function flat_to_block_coordinate(col::Integer, block_ends::Vector{Int}, block_offsets::Vector{Int}, block_dims::Vector{Int})
    block = searchsortedfirst(block_ends, Int(col))
    local_index = Int(col) - block_offsets[block]
    n = block_dims[block]
    row = mod(local_index, n) + 1
    column = div(local_index, n) + 1
    return block, row, column
end

function coefficient_maps(mp)
    hpsd_constraints = Int[]
    zero_constraints = Int[]
    for (i, (cone, _)) in pairs(mp.constraints)
        if cone == :HPSD
            push!(hpsd_constraints, i)
        elseif cone == :Zero
            push!(zero_constraints, i)
        else
            throw(ArgumentError("unsupported cone $cone in H4 SDP export"))
        end
    end
    return hpsd_constraints, zero_constraints
end

function compact_rows(A::SparseMatrixCSC{T,Int}, rhs::Vector{Float64}; atol::Real) where {T}
    SparseArrays.droptol!(A, atol)
    dropzeros!(A)
    row_has_nnz = falses(size(A, 1))
    for row in rowvals(A)
        row_has_nnz[row] = true
    end

    keep = row_has_nnz .| (abs.(rhs) .> atol)
    bad = findall((.!row_has_nnz) .& (abs.(rhs) .> atol))
    isempty(bad) || throw(ArgumentError("inconsistent empty SDP equality rows: $(bad[1:min(end, 5)])"))

    row_map = zeros(Int, length(keep))
    kept_rhs = Float64[]
    for i in eachindex(keep)
        if keep[i]
            row_map[i] = length(kept_rhs) + 1
            push!(kept_rhs, rhs[i])
        end
    end

    I, J, V = findnz(A)
    I2 = Int[]
    J2 = Int[]
    V2 = T[]
    sizehint!(I2, length(I))
    sizehint!(J2, length(J))
    sizehint!(V2, length(V))
    for k in eachindex(I)
        mapped = row_map[I[k]]
        if mapped != 0
            push!(I2, mapped)
            push!(J2, J[k])
            push!(V2, V[k])
        end
    end
    A2 = sparse(I2, J2, V2, length(kept_rhs), size(A, 2))
    SparseArrays.droptol!(A2, atol)
    dropzeros!(A2)
    return A2, kept_rhs, count(!, keep)
end

function add_native_hermitian_functional!(builder, row_index, offsets, block_dims, block, r, c, d)
    # Encode real(d * X[r,c]) as real(conj(A) ⋅ X) with a Hermitian A block.
    # BPSDP's complex dual residual lives in full block storage, so emitting
    # non-Hermitian row matrices is not harmless: it creates fake skew residuals.
    if r == c
        add_entry!(builder, row_index, flat_index(offsets, block_dims, block, r, c), real(d) + 0.0im)
    else
        add_entry!(builder, row_index, flat_index(offsets, block_dims, block, r, c), conj(d) / 2)
        add_entry!(builder, row_index, flat_index(offsets, block_dims, block, c, r), d / 2)
    end
    return nothing
end

function add_native_complex_coeff!(builder, row_re, row_im, offsets, block_dims, block, row, col, coef)
    # SOS rows contain -real(coef * X[row,col]) and -imag(coef * X[row,col]).
    add_native_hermitian_functional!(builder, row_re, offsets, block_dims, block, row, col, -coef)
    add_native_hermitian_functional!(builder, row_im, offsets, block_dims, block, row, col, im * coef)
    return nothing
end

function assemble_bpsdp_native(mp; atol::Real = 1e-14)
    symmetric_basis, sorted_unsym_basis, basis_to_sym_idx = symmetric_basis_data(mp)
    n_basis = length(symmetric_basis)
    rhs = objective_rhs(mp, symmetric_basis)
    hpsd_constraints, zero_constraints = coefficient_maps(mp)

    block_dims = Int[]
    block_roles = NamedTuple[]
    hpsd_block = Dict{Int,Int}()
    zero_blocks = Dict{Int,Tuple{Int,Int}}()

    for ci in hpsd_constraints
        dim = size(mp.constraints[ci][2], 1)
        push!(block_dims, dim)
        hpsd_block[ci] = length(block_dims)
        push!(block_roles, (; role = "hpsd", source_constraint = ci, sign = 1, dim))
    end

    push!(block_dims, 1)
    bplus_block = length(block_dims)
    push!(block_roles, (; role = "objective_bound_plus", source_constraint = 0, sign = 1, dim = 1))
    push!(block_dims, 1)
    bminus_block = length(block_dims)
    push!(block_roles, (; role = "objective_bound_minus", source_constraint = 0, sign = -1, dim = 1))

    for ci in zero_constraints
        dim = size(mp.constraints[ci][2], 1)
        push!(block_dims, dim)
        plus_block = length(block_dims)
        push!(block_roles, (; role = "zero_free_plus", source_constraint = ci, sign = 1, dim))
        push!(block_dims, dim)
        minus_block = length(block_dims)
        push!(block_roles, (; role = "zero_free_minus", source_constraint = ci, sign = -1, dim))
        zero_blocks[ci] = (plus_block, minus_block)
    end

    offsets, ncols = block_offsets(block_dims)
    builder = SparseBuilder(ComplexF64; atol)

    c = zeros(ComplexF64, ncols) # low-level BPSDP minimizes; this is -max(b⁺ - b⁻)
    bplus_col = flat_index(offsets, block_dims, bplus_block, 1, 1)
    bminus_col = flat_index(offsets, block_dims, bminus_block, 1, 1)
    c[bplus_col] = -1.0 + 0.0im
    c[bminus_col] = 1.0 + 0.0im
    add_entry!(builder, 1, bplus_col, -1.0 + 0.0im)
    add_entry!(builder, 1, bminus_col, 1.0 + 0.0im)

    for (ci, (cone, mat)) in pairs(mp.constraints)
        Cαjs = NCTSSoS.get_Cαj(sorted_unsym_basis, mat)
        if cone == :HPSD
            block = hpsd_block[ci]
            for ((basis_idx, row, col), coef) in Cαjs
                unsym_mono = sorted_unsym_basis[basis_idx]
                sym_idx = basis_to_sym_idx[unsym_mono]
                add_native_complex_coeff!(builder, sym_idx, n_basis + sym_idx,
                    offsets, block_dims, block, row, col, coef)
            end
        elseif cone == :Zero
            plus_block, minus_block = zero_blocks[ci]
            for ((basis_idx, row, col), coef) in Cαjs
                unsym_mono = sorted_unsym_basis[basis_idx]
                sym_idx = basis_to_sym_idx[unsym_mono]
                add_native_complex_coeff!(builder, sym_idx, n_basis + sym_idx,
                    offsets, block_dims, plus_block, row, col, coef)
                add_native_complex_coeff!(builder, sym_idx, n_basis + sym_idx,
                    offsets, block_dims, minus_block, row, col, -coef)
            end
        end
    end

    A = sparse(builder.I, builder.J, builder.V, 2 * n_basis, ncols)
    A, b, dropped_rows = compact_rows(A, rhs; atol)
    return (; A, b, c, block_dims, block_roles, n_basis, dropped_rows,
        nrows_before = 2 * n_basis, hpsd_constraints, zero_constraints)
end

function add_real_symmetric_coeff!(builder, row, offsets, block_dims, block, p, q, value)
    if p == q
        add_entry!(builder, row, flat_index(offsets, block_dims, block, p, q), value)
    else
        half = value / 2
        add_entry!(builder, row, flat_index(offsets, block_dims, block, p, q), half)
        add_entry!(builder, row, flat_index(offsets, block_dims, block, q, p), half)
    end
    return nothing
end

function add_real_lift_coeff!(builder, row_re, row_im, offsets, block_dims, block, dim, row, col, coef, sign)
    c_re = real(coef) * sign
    c_im = imag(coef) * sign

    # Same adjoint as NCTSSoS._sos_dualize_hermitian:
    # X1 = Z[row,col] + Z[dim+row,dim+col]
    # X2 = Z[dim+row,col] - Z[row,dim+col]
    add_real_symmetric_coeff!(builder, row_re, offsets, block_dims, block, row, col, -c_re)
    add_real_symmetric_coeff!(builder, row_re, offsets, block_dims, block, dim + row, dim + col, -c_re)
    add_real_symmetric_coeff!(builder, row_re, offsets, block_dims, block, dim + row, col, c_im)
    add_real_symmetric_coeff!(builder, row_re, offsets, block_dims, block, row, dim + col, -c_im)

    add_real_symmetric_coeff!(builder, row_im, offsets, block_dims, block, row, col, -c_im)
    add_real_symmetric_coeff!(builder, row_im, offsets, block_dims, block, dim + row, dim + col, -c_im)
    add_real_symmetric_coeff!(builder, row_im, offsets, block_dims, block, dim + row, col, -c_re)
    add_real_symmetric_coeff!(builder, row_im, offsets, block_dims, block, row, dim + col, c_re)
    return nothing
end

function assemble_sdpa_real_lift(mp; atol::Real = 1e-14)
    symmetric_basis, sorted_unsym_basis, basis_to_sym_idx = symmetric_basis_data(mp)
    n_basis = length(symmetric_basis)
    rhs = objective_rhs(mp, symmetric_basis)
    hpsd_constraints, zero_constraints = coefficient_maps(mp)

    block_dims = Int[]
    block_roles = NamedTuple[]
    hpsd_block = Dict{Int,Int}()
    zero_blocks = Dict{Int,Tuple{Int,Int}}()

    for ci in hpsd_constraints
        dim = size(mp.constraints[ci][2], 1)
        push!(block_dims, 2 * dim)
        hpsd_block[ci] = length(block_dims)
        push!(block_roles, (; role = "real_lift_hpsd", source_constraint = ci, sign = 1, dim = 2 * dim, source_dim = dim))
    end

    push!(block_dims, 1)
    bplus_block = length(block_dims)
    push!(block_roles, (; role = "objective_bound_plus", source_constraint = 0, sign = 1, dim = 1, source_dim = 1))
    push!(block_dims, 1)
    bminus_block = length(block_dims)
    push!(block_roles, (; role = "objective_bound_minus", source_constraint = 0, sign = -1, dim = 1, source_dim = 1))

    for ci in zero_constraints
        source_dim = size(mp.constraints[ci][2], 1)
        dim = 2 * source_dim
        push!(block_dims, dim)
        plus_block = length(block_dims)
        push!(block_roles, (; role = "real_lift_zero_free_plus", source_constraint = ci, sign = 1, dim, source_dim))
        push!(block_dims, dim)
        minus_block = length(block_dims)
        push!(block_roles, (; role = "real_lift_zero_free_minus", source_constraint = ci, sign = -1, dim, source_dim))
        zero_blocks[ci] = (plus_block, minus_block)
    end

    offsets, ncols = block_offsets(block_dims)
    builder = SparseBuilder(Float64; atol)

    bplus_col = flat_index(offsets, block_dims, bplus_block, 1, 1)
    bminus_col = flat_index(offsets, block_dims, bminus_block, 1, 1)
    add_entry!(builder, 1, bplus_col, -1.0)
    add_entry!(builder, 1, bminus_col, 1.0)

    for (ci, (cone, mat)) in pairs(mp.constraints)
        Cαjs = NCTSSoS.get_Cαj(sorted_unsym_basis, mat)
        source_dim = size(mat, 1)
        if cone == :HPSD
            block = hpsd_block[ci]
            for ((basis_idx, row, col), coef) in Cαjs
                unsym_mono = sorted_unsym_basis[basis_idx]
                sym_idx = basis_to_sym_idx[unsym_mono]
                add_real_lift_coeff!(builder, sym_idx, n_basis + sym_idx,
                    offsets, block_dims, block, source_dim, row, col, coef, 1)
            end
        elseif cone == :Zero
            plus_block, minus_block = zero_blocks[ci]
            for ((basis_idx, row, col), coef) in Cαjs
                unsym_mono = sorted_unsym_basis[basis_idx]
                sym_idx = basis_to_sym_idx[unsym_mono]
                add_real_lift_coeff!(builder, sym_idx, n_basis + sym_idx,
                    offsets, block_dims, plus_block, source_dim, row, col, coef, 1)
                add_real_lift_coeff!(builder, sym_idx, n_basis + sym_idx,
                    offsets, block_dims, minus_block, source_dim, row, col, coef, -1)
            end
        end
    end

    A = sparse(builder.I, builder.J, builder.V, 2 * n_basis, ncols)
    A, b, dropped_rows = compact_rows(A, rhs; atol)
    return (; A, b, block_dims, block_roles, n_basis, dropped_rows,
        nrows_before = 2 * n_basis, hpsd_constraints, zero_constraints)
end

function write_sdpa_sparse(path::AbstractString, assembly; atol::Real = 1e-14)
    mkparent(path)
    A = assembly.A
    b = assembly.b
    block_dims = Vector{Int}(assembly.block_dims)
    block_offsets_vec, _ = block_offsets(block_dims)
    block_ends = block_offsets_vec .+ block_dims .^ 2 .- 1

    AT = SparseMatrixCSC(transpose(A))

    open(path, "w") do io
        println(io, "* H4/Nk=2 periodic V2RDM SOS dual dumped from NCTSSoS.jl")
        println(io, "* Real-lifted SDPA sparse format. Use libsdp/BPSDP procedure=maximize.")
        @printf(io, "%10d\n", length(b))
        @printf(io, "%10d\n", length(block_dims))
        for dim in block_dims
            @printf(io, "%10d ", dim)
        end
        write(io, "\n")
        for value in b
            @printf(io, "% .16e ", value)
        end
        write(io, "\n")

        # Objective row F0: maximize b⁺ - b⁻.  Blocks are recorded in roles.
        for (block, role) in pairs(assembly.block_roles)
            if role.role == "objective_bound_plus"
                @printf(io, "%5d %5d %5d %5d % .16e\n", 0, block, 1, 1, 1.0)
            elseif role.role == "objective_bound_minus"
                @printf(io, "%5d %5d %5d %5d % .16e\n", 0, block, 1, 1, -1.0)
            end
        end

        for matrix_index in 1:size(A, 1)
            for ptr in AT.colptr[matrix_index]:(AT.colptr[matrix_index + 1] - 1)
                flat_col = AT.rowval[ptr]
                value = AT.nzval[ptr]
                abs(value) <= atol && continue
                block, row, col = flat_to_block_coordinate(flat_col, block_ends, block_offsets_vec, block_dims)
                row <= col || continue
                @printf(io, "%5d %5d %5d %5d % .16e\n", matrix_index, block, row, col, value)
            end
        end
    end
    return path
end

function write_bpsdp_payload(path::AbstractString, assembly, metadata)
    mkparent(path)
    payload = Dict{Symbol,Any}(
        :format_version => 1,
        :created_at => string(Dates.now()),
        :description => "H4/Nk=2 periodic V2RDM SOS dual from NCTSSoS.jl; BPSDP.jl low-level ComplexBPSDPProblem data",
        :A => assembly.A,
        :b => assembly.b,
        :c => assembly.c,
        :block_dims => assembly.block_dims,
        :block_roles => assembly.block_roles,
        :objective_convention => "problem minimizes c⋅x = -(SOS lower bound); report -state.objective_primal",
        :metadata => metadata,
    )
    open(path, "w") do io
        serialize(io, payload)
    end
    return path
end

function write_bpsdp_loader(path::AbstractString, default_data::AbstractString)
    mkparent(path)
    rel = basename(default_data)
    content = """
# Loader for the H4/Nk=2 periodic V2RDM SDP dumped from NCTSSoS.jl.
# This file intentionally depends only on BPSDP.jl plus Julia stdlibs.

module H4PeriodicMomentSOSBPSDP

using BPSDP
using LinearAlgebra
using Serialization
using SparseArrays

const DEFAULT_DATA = joinpath(@__DIR__, $(repr(rel)))

function load_h4_periodic_data(path::AbstractString = DEFAULT_DATA)
    return open(deserialize, path)
end

function h4_periodic_problem(data = load_h4_periodic_data())
    A = SparseMatrixCSC{ComplexF64,Int}(data[:A])
    A_conj = conj.(A)
    b = Float64.(data[:b])
    c = ComplexF64.(data[:c])
    block_dims = Int.(data[:block_dims])
    dual_work = zeros(ComplexF64, size(A, 1))

    Au! = function (out, x)
        mul!(dual_work, A_conj, x)
        @. out = real(dual_work)
        return nothing
    end

    ATu! = function (out, y)
        mul!(out, transpose(A), y)
        return nothing
    end

    return BPSDP.ComplexBPSDPProblem(block_dims, c, b, Au!, ATu!)
end

function solve_h4_periodic(; path::AbstractString = DEFAULT_DATA, options::BPSDP.BPSDPOptions = BPSDP.BPSDPOptions())
    data = load_h4_periodic_data(path)
    problem = h4_periodic_problem(data)
    state = BPSDP.solve(problem; options)
    return (; state, lower_bound = -state.objective_primal, problem, data)
end

end # module
"""
    write(path, content)
    return path
end

function write_libsdp_runner(path::AbstractString, sdpa_path::AbstractString)
    mkparent(path)
    rel = basename(sdpa_path)
    content = """
#!/usr/bin/env python3
# Load and solve the H4/Nk=2 periodic V2RDM SDP dump with libsdp.
# The SDPA file is a maximization-form SOS dual; libsdp needs procedure=maximize.

import os
import numpy as np
from libsdp import sdp_options
from libsdp.sdpa_file_io import read_sdpa_problem
from libsdp.sdp_helper import sdp_solver


def main():
    filename = os.path.join(os.path.dirname(__file__), $(repr(rel)))
    b, matrices, block_dim = read_sdpa_problem(filename)

    options = sdp_options()
    options.sdp_algorithm = "bpsdp"
    options.procedure = "maximize"
    options.guess_type = "zero"

    sdp = sdp_solver(options, matrices, block_dim)
    x = sdp.solve(b, 50000)
    c = np.array(sdp.get_c())
    lower_bound = -float(np.dot(c, x))
    print(f"SOS lower bound: {lower_bound:.12f}")


if __name__ == "__main__":
    main()
"""
    write(path, content)
    try
        chmod(path, 0o755)
    catch
    end
    return path
end

function write_metadata(path::AbstractString, data, options::Options, dump_options::DumpOptions, assemblies)
    mkparent(path)
    stats = constraint_stats(data.moment_problem)
    open(path, "w") do io
        println(io, "# H4 periodic SDP dump metadata")
        println(io, "created_at = \"$(Dates.now())\"")
        println(io, "integrals = \"$(options.integrals_path)\"")
        println(io, "blocking = \"$(options.blocking)\"")
        println(io, "formulation = \"1D+PQG, total N, trace constraints\"")
        println(io, "particle_sector = \"total\"")
        println(io, "trace_constraints = true")
        println(io, "include_one_d = true")
        println(io, "zero_objective = false")
        println(io, "atol = $(dump_options.atol)")
        println(io, "moment_hpsd_sizes = $(stats.hpsd_sizes)")
        println(io, "moment_zero_constraints = $(length(stats.zero_sizes))")
        println(io, "canonical_monomials = $(stats.total_canonical_moments)")
        println(io, "direct_real_moment_variables = $(stats.direct_real_moment_variables)")
        for (name, assembly) in assemblies
            println(io)
            println(io, "[$name]")
            println(io, "rows = $(length(assembly.b))")
            println(io, "rows_before_clean = $(assembly.nrows_before)")
            println(io, "dropped_zero_rows = $(assembly.dropped_rows)")
            println(io, "cols = $(size(assembly.A, 2))")
            println(io, "nnz = $(nnz(assembly.A))")
            println(io, "blocks = $(length(assembly.block_dims))")
            println(io, "max_block_dim = $(maximum(assembly.block_dims))")
        end
    end
    return path
end

function main(argv = ARGS)
    dump_options = parse_dump_options(argv)
    formulation_options = parse_options(dump_options.formulation_args)

    @printf("%-48s %s\n", "max RSS at start", rss_string())
    build_seconds = @elapsed data = build_h4_pqg_moment_problem(formulation_options)
    print_summary(data, formulation_options, build_seconds)

    assemblies = Pair{String,Any}[]
    metadata = Dict{Symbol,Any}(
        :integrals_path => formulation_options.integrals_path,
        :blocking => string(formulation_options.blocking),
        :formulation => "1D+PQG, total N, trace constraints",
        :particle_sector => "total",
        :trace_constraints => true,
        :include_one_d => true,
        :zero_objective => false,
    )

    if dump_options.format in (:all, :bpsdp)
        @printf("%-48s %s\n", "assembling BPSDP.jl native data", dump_options.bpsdp_out)
        seconds = @elapsed bpsdp_assembly = assemble_bpsdp_native(data.moment_problem; atol = dump_options.atol)
        @printf("%-48s %.3f s\n", "BPSDP assembly walltime", seconds)
        @printf("%-48s rows=%d cols=%d nnz=%d blocks=%d\n", "BPSDP assembly size",
            length(bpsdp_assembly.b), size(bpsdp_assembly.A, 2), nnz(bpsdp_assembly.A), length(bpsdp_assembly.block_dims))
        write_bpsdp_payload(dump_options.bpsdp_out, bpsdp_assembly, metadata)
        write_bpsdp_loader(dump_options.bpsdp_loader, dump_options.bpsdp_out)
        push!(assemblies, "bpsdp_native" => bpsdp_assembly)
        @printf("%-48s %s\n", "wrote BPSDP data", dump_options.bpsdp_out)
        @printf("%-48s %s\n", "wrote BPSDP loader", dump_options.bpsdp_loader)
        GC.gc()
    end

    if dump_options.format in (:all, :libsdp)
        @printf("%-48s %s\n", "assembling libsdp SDPA data", dump_options.libsdp_out)
        seconds = @elapsed sdpa_assembly = assemble_sdpa_real_lift(data.moment_problem; atol = dump_options.atol)
        @printf("%-48s %.3f s\n", "SDPA assembly walltime", seconds)
        @printf("%-48s rows=%d cols=%d nnz=%d blocks=%d\n", "SDPA assembly size",
            length(sdpa_assembly.b), size(sdpa_assembly.A, 2), nnz(sdpa_assembly.A), length(sdpa_assembly.block_dims))
        seconds = @elapsed write_sdpa_sparse(dump_options.libsdp_out, sdpa_assembly; atol = dump_options.atol)
        write_libsdp_runner(dump_options.libsdp_runner, dump_options.libsdp_out)
        push!(assemblies, "libsdp_sdpa" => sdpa_assembly)
        @printf("%-48s %.3f s\n", "SDPA write walltime", seconds)
        @printf("%-48s %s\n", "wrote libsdp SDPA", dump_options.libsdp_out)
        @printf("%-48s %s\n", "wrote libsdp runner", dump_options.libsdp_runner)
        GC.gc()
    end

    write_metadata(dump_options.metadata_out, data, formulation_options, dump_options, assemblies)
    @printf("%-48s %s\n", "wrote metadata", dump_options.metadata_out)
    @printf("%-48s %s\n", "max RSS at end", rss_string())
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
