#!/usr/bin/env julia
# Solve a libsdp-complex SDPA-sparse `.dat-c` file through BPSDP.jl's public
# primitive API: Problem(blocks, c, b, SparseOperator(A)).
#
# The parser hermitianizes every sparse functional.  That preserves
# Re(<F,X>) for Hermitian X and avoids fake skew-Hermitian dual residuals in
# BPSDP's full-storage complex blocks.

using BPSDP
using LinearAlgebra
using Printf
using SparseArrays

try
    @eval using CUDA
catch err
    @warn "CUDA.jl is not available; CUDA backend will be rejected" exception = (err, catch_backtrace())
end

const ATOL = 1e-14

struct StopSolve <: Exception
    reason::Symbol
end

mutable struct LastMetrics
    outer::Int
    inner::Int
    primal::Float64
    dual::Float64
    mu::Float64
    primal_error::Float64
    dual_error::Float64
    walltime::Float64
end

LastMetrics() = LastMetrics(0, 0, NaN, NaN, NaN, Inf, Inf, 0.0)

function next_data_line(io)
    for line in eachline(io)
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "*") || startswith(stripped, "\"")) && continue
        return stripped
    end
    error("unexpected EOF while parsing SDPA-c file")
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

function add_hermitian_functional!(I, J, V, matrix_row::Int, offsets, block_dims, block, row, col, value::ComplexF64)
    abs(value) <= ATOL && return nothing
    if row == col
        diag_value = ComplexF64(real(value), 0.0)
        abs(diag_value) <= ATOL && return nothing
        push!(I, matrix_row)
        push!(J, flat_index(offsets, block_dims, block, row, col))
        push!(V, diag_value)
    else
        half = value / 2
        push!(I, matrix_row)
        push!(J, flat_index(offsets, block_dims, block, row, col))
        push!(V, half)
        push!(I, matrix_row)
        push!(J, flat_index(offsets, block_dims, block, col, row))
        push!(V, conj(half))
    end
    return nothing
end

function add_hermitian_objective!(c, offsets, block_dims, block, row, col, value::ComplexF64)
    abs(value) <= ATOL && return nothing
    if row == col
        c[flat_index(offsets, block_dims, block, row, col)] += ComplexF64(real(value), 0.0)
    else
        half = value / 2
        c[flat_index(offsets, block_dims, block, row, col)] += half
        c[flat_index(offsets, block_dims, block, col, row)] += conj(half)
    end
    return nothing
end

function parse_dat_c_as_bpsdp(path::AbstractString)
    open(path, "r") do io
        m = parse(Int, next_data_line(io))
        nblocks = parse(Int, next_data_line(io))
        block_dims = parse.(Int, split(next_data_line(io)))
        length(block_dims) == nblocks || error("block dimension count mismatch")
        b = parse.(Float64, split(next_data_line(io)))
        length(b) == m || error("RHS length mismatch")

        offsets, ncols = block_offsets(block_dims)
        c = zeros(ComplexF64, ncols)
        I = Int[]
        J = Int[]
        V = ComplexF64[]

        for line in eachline(io)
            stripped = strip(line)
            (isempty(stripped) || startswith(stripped, "*") || startswith(stripped, "\"")) && continue
            tok = split(stripped)
            length(tok) >= 6 || error("bad sparse entry: $stripped")
            k = parse(Int, tok[1])
            block = parse(Int, tok[2])
            row = parse(Int, tok[3])
            col = parse(Int, tok[4])
            value = ComplexF64(parse(Float64, tok[5]), parse(Float64, tok[6]))
            if k == 0
                add_hermitian_objective!(c, offsets, block_dims, block, row, col, value)
            else
                add_hermitian_functional!(I, J, V, k, offsets, block_dims, block, row, col, value)
            end
        end

        A = sparse(I, J, V, m, ncols)
        dropzeros!(A)
        return (; A, b, c, block_dims)
    end
end

function parse_args(argv)
    opts = Dict{String,String}()
    flags = Set{String}()
    for arg in argv
        if startswith(arg, "--") && occursin("=", arg)
            key, value = split(arg[3:end], "=", limit = 2)
            opts[String(key)] = String(value)
        elseif startswith(arg, "--")
            push!(flags, String(arg[3:end]))
        else
            error("unexpected positional argument: $arg")
        end
    end
    haskey(opts, "problem") || error("--problem=PATH is required")
    return (; problem = opts["problem"],
        backend = Symbol(get(opts, "backend", "cpu")),
        max_iter = parse(Int, get(opts, "max-iter", "2000")),
        cg_max_iter = parse(Int, get(opts, "cg-max-iter", "1000")),
        max_seconds = parse(Float64, get(opts, "max-seconds", "900")),
        monitor_period = parse(Int, get(opts, "monitor-period", "25")),
        plateau_window = parse(Int, get(opts, "plateau-window", "250")),
        plateau_abs = parse(Float64, get(opts, "plateau-abs", "1e-6")),
        objective_tol = parse(Float64, get(opts, "objective-tol", "1e-5")),
        error_tol = parse(Float64, get(opts, "error-tol", "1e-5")),
        cg_tol = parse(Float64, get(opts, "cg-tol", "1e-8")),
        penalty_parameter = parse(Float64, get(opts, "penalty", "0.1")),
        mu_update_frequency = parse(Int, get(opts, "mu-update-frequency", "250")),
        energy_shift = parse(Float64, get(opts, "energy-shift", "NaN")),
        threads = parse(Int, get(opts, "threads", string(Threads.nthreads()))),
        summary = get(opts, "summary", ""),
        fail_fast = "fail-fast" in flags)
end

function backend_value(backend::Symbol)
    if backend === :cpu
        return BPSDP.CPUBackend()
    elseif backend === :cuda || backend === :gpu
        return BPSDP.CUDABackend()
    else
        throw(ArgumentError("unknown backend $backend; use cpu or cuda"))
    end
end

function cuda_functional_string()
    if isdefined(Main, :CUDA)
        try
            return string(CUDA.functional())
        catch err
            return "error: $err"
        end
    end
    return "CUDA.jl not loaded"
end

function cuda_device_string()
    if isdefined(Main, :CUDA)
        try
            CUDA.functional() || return "none"
            return string(CUDA.device())
        catch err
            return "error: $err"
        end
    end
    return "none"
end

function write_summary(path::AbstractString, rows::Vector{Pair{String,String}})
    isempty(path) && return nothing
    mkpath(dirname(path))
    open(path, "w") do io
        for (key, value) in rows
            println(io, key, " = ", value)
        end
    end
    return nothing
end

function main(argv = ARGS)
    args = parse_args(argv)
    args.backend in (:cpu, :cuda, :gpu) || error("--backend must be cpu or cuda")

    println("== BPSDP.jl primitive dat-c solve ==")
    @printf("%-36s %s\n", "problem", args.problem)
    @printf("%-36s %s\n", "backend", string(args.backend))
    @printf("%-36s %d\n", "Julia threads", Threads.nthreads())
    @printf("%-36s %s\n", "CUDA functional", cuda_functional_string())
    @printf("%-36s %s\n", "CUDA device", cuda_device_string())
    flush(stdout)

    parse_seconds = @elapsed data = parse_dat_c_as_bpsdp(args.problem)
    @printf("%-36s %.3f s\n", "parse/assemble walltime", parse_seconds)
    @printf("%-36s %d\n", "blocks", length(data.block_dims))
    @printf("%-36s %d / %d / %d\n", "block min/max/sum", minimum(data.block_dims), maximum(data.block_dims), sum(data.block_dims))
    @printf("%-36s %d\n", "primal full-storage variables", length(data.c))
    @printf("%-36s %d\n", "constraints", length(data.b))
    @printf("%-36s %d\n", "A nnz", nnz(data.A))
    @printf("%-36s %d\n", "objective nnz", count(!iszero, data.c))
    flush(stdout)

    problem = BPSDP.Problem(
        BPSDP.HermitianPSDBlock.(Int.(data.block_dims)),
        ComplexF64.(data.c),
        Float64.(data.b),
        BPSDP.SparseOperator(SparseMatrixCSC{ComplexF64,Int}(data.A)),
    )

    start = time()
    last = LastMetrics()
    objective_history = Float64[]
    status = :not_started
    stop_reason = :none
    solution = nothing

    monitor = function (_print_level, outer, inner, primal, dual, mu, primal_error, dual_error)
        wall = time() - start
        last.outer = outer
        last.inner += inner
        last.primal = primal
        last.dual = dual
        last.mu = mu
        last.primal_error = primal_error
        last.dual_error = dual_error
        last.walltime = wall
        push!(objective_history, primal)
        if outer % args.monitor_period == 0 || outer <= 5
            @printf("%8d %6d %14.6f %14.6f %12.6e %9.3f %12.5e %12.5e %9.1f\n",
                outer, inner, primal, dual, abs(primal - dual), mu, primal_error, dual_error, wall)
            flush(stdout)
        end
        if wall >= args.max_seconds
            throw(StopSolve(:timeout))
        end
        if args.plateau_window > 0 && length(objective_history) > args.plateau_window
            old = objective_history[end - args.plateau_window]
            if isfinite(old) && isfinite(primal) && abs(primal - old) <= args.plateau_abs
                throw(StopSolve(:plateau))
            end
        end
        return nothing
    end

    settings = BPSDP.Settings(
        max_iter = args.max_iter,
        cg_max_iter = args.cg_max_iter,
        mu_update_frequency = args.mu_update_frequency,
        penalty_parameter = args.penalty_parameter,
        cg_convergence = args.cg_tol,
        dynamic_cg_convergence = true,
        sdp_objective_convergence_abs = args.objective_tol,
        sdp_error_convergence_abs = args.error_tol,
        guess_type = :zero,
        print_level = 0,
        progress_monitor = monitor,
        backend = backend_value(args.backend),
    )

    println()
    println("   outer  inner            c.x            b.y    |c.x-b.y|        mu     ||Ax-b||  ||ATy-c+z||      wall")
    flush(stdout)

    solve_seconds = @elapsed begin
        try
            solution = BPSDP.solve(problem; settings)
            status = solution.status
            stop_reason = status
        catch err
            if err isa StopSolve
                status = :stopped
                stop_reason = err.reason
            else
                status = :error
                stop_reason = Symbol(nameof(typeof(err)))
                if args.fail_fast
                    rethrow()
                end
                println()
                @printf("ERROR: %s\n", sprint(showerror, err))
            end
        end
    end

    println()
    println("== BPSDP.jl primitive results ==")
    @printf("%-36s %.3f s\n", "solve walltime", solve_seconds)
    @printf("%-36s %s\n", "backend", string(args.backend))
    @printf("%-36s %s\n", "status", string(status))
    @printf("%-36s %s\n", "stop_reason", string(stop_reason))

    primal_objective = solution !== nothing ? solution.primal_objective : last.primal
    dual_objective = solution !== nothing ? solution.dual_objective : last.dual

    if solution !== nothing
        @printf("%-36s %d\n", "outer_iterations", solution.iterations)
        @printf("%-36s %d\n", "inner_iterations", solution.inner_iterations)
        @printf("%-36s %.12f\n", "primal objective", primal_objective)
        @printf("%-36s %.12f\n", "dual objective", dual_objective)
        @printf("%-36s %.6e\n", "objective gap", solution.objective_gap)
        @printf("%-36s %.6e\n", "primal error", solution.primal_error)
        @printf("%-36s %.6e\n", "dual error", solution.dual_error)
    else
        @printf("%-36s %d\n", "outer_iterations", last.outer)
        @printf("%-36s %d\n", "inner_iterations", last.inner)
        @printf("%-36s %.12f\n", "last primal objective", primal_objective)
        @printf("%-36s %.12f\n", "last dual objective", dual_objective)
        @printf("%-36s %.6e\n", "last objective gap", abs(last.primal - last.dual))
        @printf("%-36s %.6e\n", "last primal error", last.primal_error)
        @printf("%-36s %.6e\n", "last dual error", last.dual_error)
    end
    if isfinite(args.energy_shift)
        @printf("%-36s %.12f\n", "recovered total primal", primal_objective + args.energy_shift)
        @printf("%-36s %.12f\n", "recovered total dual", dual_objective + args.energy_shift)
    end

    rows = Pair{String,String}[
        "problem" => args.problem,
        "backend" => string(args.backend),
        "status" => string(status),
        "stop_reason" => string(stop_reason),
        "parse_seconds" => @sprintf("%.6f", parse_seconds),
        "solve_seconds" => @sprintf("%.6f", solve_seconds),
        "blocks" => string(length(data.block_dims)),
        "block_min" => string(minimum(data.block_dims)),
        "block_max" => string(maximum(data.block_dims)),
        "block_sum" => string(sum(data.block_dims)),
        "primal_dim" => string(length(data.c)),
        "constraints" => string(length(data.b)),
        "A_nnz" => string(nnz(data.A)),
    ]
    if isfinite(args.energy_shift)
        append!(rows, [
            "energy_shift" => @sprintf("%.12f", args.energy_shift),
            "recovered_total_primal" => @sprintf("%.12f", primal_objective + args.energy_shift),
            "recovered_total_dual" => @sprintf("%.12f", dual_objective + args.energy_shift),
        ])
    end
    if solution !== nothing
        append!(rows, [
            "outer_iterations" => string(solution.iterations),
            "inner_iterations" => string(solution.inner_iterations),
            "primal_objective" => @sprintf("%.12f", solution.primal_objective),
            "dual_objective" => @sprintf("%.12f", solution.dual_objective),
            "objective_gap" => @sprintf("%.6e", solution.objective_gap),
            "primal_error" => @sprintf("%.6e", solution.primal_error),
            "dual_error" => @sprintf("%.6e", solution.dual_error),
        ])
    else
        append!(rows, [
            "outer_iterations" => string(last.outer),
            "inner_iterations" => string(last.inner),
            "primal_objective" => @sprintf("%.12f", last.primal),
            "dual_objective" => @sprintf("%.12f", last.dual),
            "objective_gap" => @sprintf("%.6e", abs(last.primal - last.dual)),
            "primal_error" => @sprintf("%.6e", last.primal_error),
            "dual_error" => @sprintf("%.6e", last.dual_error),
        ])
    end
    write_summary(args.summary, rows)
    flush(stdout)
    return status == :error ? 2 : 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main())
end
