#!/usr/bin/env julia
# Phase-2 H2/Nk=2 A*A' diagnostic sweep.
#
# Consumes block operators dumped by demos/h2_periodic_nk2_moment_sos.jl and
# writes the per-block decision table required by PHASE2_PLAN.md.

using LinearAlgebra
using SparseArrays
using Printf
using NPZ
using JSON3
using Dates

const DEFAULT_INPUT_DIR = normpath(joinpath(@__DIR__, "..", "output", "phase2", "h2_nk2"))
const DEFAULT_RANK_TOLS = [1e-10, 1e-12, 1e-14]

struct Options
    input_dir::String
    rank_tols::Vector{Float64}
    ram_bytes::Float64
end

function parse_rank_tols(value::AbstractString)
    return [parse(Float64, strip(tok)) for tok in split(value, ",")]
end

function parse_options(argv)
    input_dir = DEFAULT_INPUT_DIR
    rank_tols = copy(DEFAULT_RANK_TOLS)
    ram_bytes = Float64(Sys.total_memory())

    for arg in argv
        if startswith(arg, "--input-dir=") || startswith(arg, "--output-dir=")
            input_dir = split(arg, "=", limit = 2)[2]
        elseif startswith(arg, "--rank-tol=")
            rank_tols = parse_rank_tols(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--ram-gib=")
            ram_bytes = parse(Float64, split(arg, "=", limit = 2)[2]) * 1024.0^3
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=<env-with-NPZ-JSON3> probes/h2_nk2_aastar_diagnostic.jl [options]\n")
            println("Options:")
            println("  --input-dir=PATH   default output/phase2/h2_nk2")
            println("  --output-dir=PATH  alias for --input-dir")
            println("  --rank-tol=a,b,c   default 1e-10,1e-12,1e-14")
            println("  --ram-gib=GB       working-RAM threshold for direct factor decision")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end
    return Options(input_dir, rank_tols, ram_bytes)
end

function read_sparse_npz(path::AbstractString)
    d = NPZ.npzread(path)
    rows = Int.(d["row"]) .+ 1
    cols = Int.(d["col"]) .+ 1
    vals = Float64.(d["data"])
    shape = Int.(d["shape"])
    return sparse(rows, cols, vals, shape[1], shape[2])
end

function write_sparse_npz(path::AbstractString, A::SparseMatrixCSC)
    mkpath(dirname(path))
    I, J, V = findnz(A)
    NPZ.npzwrite(path, Dict(
        "row" => Int64.(I .- 1),
        "col" => Int64.(J .- 1),
        "data" => Float64.(V),
        "shape" => Int64[size(A, 1), size(A, 2)],
        "index_base" => Int64[0],
    ))
end

function write_json(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, obj)
        write(io, '\n')
    end
end

finite_or_string(x::Real) = isfinite(Float64(x)) ? Float64(x) : string(Float64(x))
tol_key(tol::Real) = @sprintf("%.0e", tol)

function structural_redundancy(A::SparseMatrixCSC)
    M = Matrix(A)
    signatures = Dict{String,Int}()
    zero_rows = 0
    for r in axes(M, 1)
        parts = String[]
        for c in axes(M, 2)
            v = M[r, c]
            abs(v) <= 1e-13 && continue
            push!(parts, string(c, ":", @sprintf("%.12e", v)))
        end
        if isempty(parts)
            zero_rows += 1
        else
            sig = join(parts, ",")
            signatures[sig] = get(signatures, sig, 0) + 1
        end
    end
    nonzero_rows = size(A, 1) - zero_rows
    duplicate_rows = nonzero_rows - length(signatures)
    return (; zero_rows, duplicate_rows, predicted_kernel_dim = zero_rows + duplicate_rows)
end

function spectrum_diagnostic(A::SparseMatrixCSC, rank_tols::Vector{Float64})
    m = size(A, 1)
    n = size(A, 2)
    if m == 0
        return Dict(
            "rows" => 0, "cols" => n, "sigma1" => 0.0,
            "rank_by_tol" => Dict{String,Any}(),
            "singular_values" => Float64[],
            "threshold_mode" => "relative_to_sigma1",
        ), zeros(Float64, 0, 0)
    end

    G = Matrix(A * transpose(A))
    G = 0.5 .* (G .+ transpose(G))
    F = eigen(Symmetric(G))
    vals_asc = max.(F.values, 0.0)
    order = sortperm(vals_asc; rev = true)
    sigma = vals_asc[order]
    vecs_desc = F.vectors[:, order]
    sigma1 = isempty(sigma) ? 0.0 : sigma[1]

    ranks = Dict{String,Any}()
    for tol in rank_tols
        threshold = tol * sigma1
        r = count(>(threshold), sigma)
        sigma_r = r > 0 ? sigma[r] : 0.0
        sigma_next = r < length(sigma) ? sigma[r + 1] : 0.0
        gap = sigma_next > 0 ? sigma_r / sigma_next : Inf
        kappa = sigma_r > 0 ? sigma1 / sigma_r : Inf
        ranks[tol_key(tol)] = Dict(
            "threshold" => threshold,
            "rank" => r,
            "sigma_rank" => finite_or_string(sigma_r),
            "sigma_next" => finite_or_string(sigma_next),
            "kernel_dim" => m - r,
            "kappa_range" => finite_or_string(kappa),
            "spectral_gap" => finite_or_string(gap),
            "estimated_cg_iterations_1e-8" => finite_or_string(isfinite(kappa) ? 9.0 * sqrt(kappa) : Inf),
        )
    end

    default_tol = 1e-12 in rank_tols ? 1e-12 : rank_tols[min(2, length(rank_tols))]
    kernel_cols = findall(<=(default_tol * sigma1), sigma)
    kernel_basis = isempty(kernel_cols) ? zeros(Float64, m, 0) : vecs_desc[:, kernel_cols]

    return Dict(
        "rows" => m,
        "cols" => n,
        "sigma1" => finite_or_string(sigma1),
        "rank_by_tol" => ranks,
        "singular_values" => sigma,
        "threshold_mode" => "relative_to_sigma1",
    ), kernel_basis
end

function sparsity_diagnostic(A::SparseMatrixCSC)
    G = sparse(A * transpose(A))
    m = size(G, 1)
    nnz_G = nnz(G)
    density = m == 0 ? 0.0 : nnz_G / (m * m)
    if m == 0
        return Dict(
            "aat_rows" => 0, "aat_nnz" => 0, "aat_density" => 0.0,
            "amd_nnz_L" => 0, "amd_memory_bytes" => 0,
            "metis_available" => false,
        )
    end

    delta = max(1.0, maximum(abs, G.nzval; init = 0.0)) * 1e-12
    shifted = G + delta * sparse(I, m, m)
    amd_nnz_L = try
        F = cholesky(Symmetric(shifted); check = false)
        nnz(sparse(F.L))
    catch err
        @warn "AMD Cholesky fill measurement failed" err
        -1
    end

    return Dict(
        "aat_rows" => m,
        "aat_nnz" => nnz_G,
        "aat_density" => density,
        "regularization_delta" => delta,
        "amd_nnz_L" => amd_nnz_L,
        "amd_memory_bytes" => amd_nnz_L < 0 ? -1 : 8 * amd_nnz_L,
        "metis_available" => false,
        "metis_note" => "METIS not required/available in this Phase-2 sweep; AMD measured with CHOLMOD.",
    )
end

function rank_record(spec::Dict, tol::Real)
    return spec["rank_by_tol"][tol_key(tol)]
end

function asfloat(x)
    x isa Real && return Float64(x)
    x isa AbstractString && return x == "Inf" ? Inf : parse(Float64, x)
    return Float64(x)
end

function classify(m::Int, rank::Int, predicted_kernel::Int, gap::Float64, kappa::Float64, amd_bytes::Real, ram_bytes::Real)
    m == 0 && return "Direct Cholesky"
    empirical_kernel = m - rank
    excess_kernel = max(empirical_kernel - predicted_kernel, 0)
    rank_fraction = rank / m

    if excess_kernel > 0.10 * m || rank_fraction < 0.5
        return "Facial reduction needed"
    elseif gap > 10 && amd_bytes >= 0 && amd_bytes < ram_bytes
        return "Direct Cholesky"
    elseif gap > 10
        return "Tikhonov+CG"
    elseif isfinite(kappa) && sqrt(kappa) <= 1e3
        return "Mazziotti-CG"
    else
        return "Punt to BM"
    end
end

# -----------------------------------------------------------------------------
# Tiny pure-Julia PNG writer. This avoids dragging a plotting stack into the
# server temp environment just to draw small singular-value histograms.
# -----------------------------------------------------------------------------

be32(x::Integer) = UInt8[(UInt32(x) >> 24) & 0xff, (UInt32(x) >> 16) & 0xff, (UInt32(x) >> 8) & 0xff, UInt32(x) & 0xff]

function crc32(bytes::Vector{UInt8})
    crc = UInt32(0xffffffff)
    for b in bytes
        crc ⊻= UInt32(b)
        for _ in 1:8
            if (crc & UInt32(1)) == UInt32(1)
                crc = (crc >> 1) ⊻ UInt32(0xedb88320)
            else
                crc >>= 1
            end
        end
    end
    return ~crc
end

function adler32(bytes::Vector{UInt8})
    a = UInt32(1)
    b = UInt32(0)
    modv = UInt32(65521)
    for x in bytes
        a = (a + UInt32(x)) % modv
        b = (b + a) % modv
    end
    return (b << 16) | a
end

function png_chunk(io, typ::AbstractString, data::Vector{UInt8})
    typebytes = Vector{UInt8}(codeunits(typ))
    write(io, be32(length(data)))
    write(io, typebytes)
    write(io, data)
    write(io, be32(crc32(vcat(typebytes, data))))
end

function zlib_store(raw::Vector{UInt8})
    out = UInt8[0x78, 0x01]
    offset = 1
    while offset <= length(raw)
        len = min(65535, length(raw) - offset + 1)
        final = offset + len > length(raw)
        push!(out, final ? 0x01 : 0x00)
        l = UInt16(len)
        nl = ~l
        push!(out, UInt8(l & 0xff), UInt8((l >> 8) & 0xff), UInt8(nl & 0xff), UInt8((nl >> 8) & 0xff))
        append!(out, raw[offset:(offset + len - 1)])
        offset += len
    end
    write_adler = adler32(raw)
    append!(out, be32(write_adler))
    return out
end

function write_png_rgb(path::AbstractString, img::Array{UInt8,3})
    h, w, c = size(img)
    c == 3 || error("RGB image must have three channels")
    raw = UInt8[]
    sizehint!(raw, h * (1 + 3w))
    for y in 1:h
        push!(raw, 0x00)
        for x in 1:w, ch in 1:3
            push!(raw, img[y, x, ch])
        end
    end
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a])
        ihdr = UInt8[]
        append!(ihdr, be32(w)); append!(ihdr, be32(h))
        append!(ihdr, UInt8[0x08, 0x02, 0x00, 0x00, 0x00])
        png_chunk(io, "IHDR", ihdr)
        png_chunk(io, "IDAT", zlib_store(raw))
        png_chunk(io, "IEND", UInt8[])
    end
end

function setpixel!(img, x, y, color)
    h, w, _ = size(img)
    1 <= x <= w && 1 <= y <= h || return nothing
    img[y, x, 1] = color[1]
    img[y, x, 2] = color[2]
    img[y, x, 3] = color[3]
    return nothing
end

function fillrect!(img, x1, y1, x2, y2, color)
    h, w, _ = size(img)
    xa, xb = max(1, min(x1, x2)), min(w, max(x1, x2))
    ya, yb = max(1, min(y1, y2)), min(h, max(y1, y2))
    for y in ya:yb, x in xa:xb
        setpixel!(img, x, y, color)
    end
end

function write_histogram_png(path::AbstractString, sigma::Vector{Float64})
    w, h = 720, 420
    img = fill(UInt8(255), h, w, 3)
    black = UInt8[0x20, 0x20, 0x20]
    blue = UInt8[0x2b, 0x6c, 0xb0]
    gray = UInt8[0xdd, 0xdd, 0xdd]
    left, right, top, bottom = 60, w - 25, 25, h - 55
    fillrect!(img, left, bottom, right, bottom + 1, black)
    fillrect!(img, left, top, left + 1, bottom, black)
    for k in 0:4
        y = round(Int, bottom - k * (bottom - top) / 4)
        fillrect!(img, left, y, right, y, gray)
    end
    if !isempty(sigma)
        positive = [s for s in sigma if s > 0]
        if !isempty(positive)
            logs = log10.(positive)
            lo, hi = minimum(logs), maximum(logs)
            hi == lo && (hi = lo + 1)
            nbins = 40
            counts = zeros(Int, nbins)
            for v in logs
                b = clamp(fld(Int(floor((v - lo) / (hi - lo) * nbins)), 1) + 1, 1, nbins)
                counts[b] += 1
            end
            maxc = max(maximum(counts), 1)
            plotw = right - left - 8
            barw = max(1, plotw ÷ nbins)
            for b in 1:nbins
                x1 = left + 4 + (b - 1) * plotw ÷ nbins
                x2 = min(right - 2, x1 + barw - 1)
                y1 = round(Int, bottom - counts[b] / maxc * (bottom - top - 8))
                fillrect!(img, x1, y1, x2, bottom - 1, blue)
            end
        end
    end
    write_png_rgb(path, img)
end

function analyze_matrix(A::SparseMatrixCSC, rank_tols::Vector{Float64})
    spec, kernel_basis = spectrum_diagnostic(A, rank_tols)
    sparse_stats = sparsity_diagnostic(A)
    predicted = structural_redundancy(A)
    rec12 = rank_record(spec, 1e-12)
    m = size(A, 1)
    rank12 = Int(rec12["rank"])
    kernel12 = Int(rec12["kernel_dim"])
    gap12 = asfloat(rec12["spectral_gap"])
    kappa12 = asfloat(rec12["kappa_range"])
    cg12 = asfloat(rec12["estimated_cg_iterations_1e-8"])
    decision = classify(m, rank12, predicted.predicted_kernel_dim, gap12, kappa12,
        sparse_stats["amd_memory_bytes"], 0.0) # caller rewrites with RAM below
    return (; spec, kernel_basis, sparse_stats, predicted, rank12, kernel12, gap12, kappa12, cg12, decision)
end

function block_dirs(input_dir::AbstractString)
    root = joinpath(input_dir, "blocks")
    isdir(root) || error("missing blocks directory: $root")
    return sort!([joinpath(root, name) for name in readdir(root) if isdir(joinpath(root, name))])
end

function decision_with_ram(result, ram_bytes)
    return classify(result.spec["rows"], result.rank12, result.predicted.predicted_kernel_dim,
        result.gap12, result.kappa12, result.sparse_stats["amd_memory_bytes"], ram_bytes)
end

function format_float(x; digits=3)
    if x isa AbstractString
        return x
    end
    xf = Float64(x)
    if isinf(xf)
        return "Inf"
    elseif isnan(xf)
        return "NaN"
    elseif xf == 0
        return "0"
    elseif abs(xf) < 1e-3 || abs(xf) > 1e4
        return @sprintf("%.*e", digits, xf)
    else
        return @sprintf("%.*f", digits, xf)
    end
end

function try_json(path)
    isfile(path) || return nothing
    try
        return JSON3.read(read(path, String))
    catch err
        @warn "could not parse JSON" path err
        return nothing
    end
end

function object_get(obj, keys...)
    cur = obj
    for key in keys
        cur === nothing && return nothing
        try
            cur = cur[key]
        catch
            try
                cur = getproperty(cur, Symbol(key))
            catch
                return nothing
            end
        end
    end
    return cur
end

function write_summary(input_dir::AbstractString, rows::Vector{Dict{String,Any}}, ram_bytes::Real)
    path = joinpath(input_dir, "summary.md")
    solve = try_json(joinpath(input_dir, "solve", "objective.json"))

    counts = Dict{String,Int}()
    for row in rows
        decision = row["decision"]
        counts[decision] = get(counts, decision, 0) + 1
    end

    open(path, "w") do io
        println(io, "# Phase 2 H₂/Nk=2 A A* diagnostic summary")
        println(io)
        println(io, "Generated: `", Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"), "`")
        println(io, "Input dir: `", input_dir, "`")
        println(io, "Working RAM threshold for direct Cholesky: ", format_float(ram_bytes / 1024.0^3), " GiB")
        println(io)

        if solve !== nothing
            println(io, "## BPSDP solve check")
            println(io)
            println(io, "| quantity | value |")
            println(io, "|---|---:|")
            for (label, keys) in [
                ("termination", ("termination_status",)),
                ("objective active Ha", ("objective_active_Ha",)),
                ("objective shifted total Ha", ("objective_total_shifted_Ha",)),
                ("HF active Ha", ("hf_active_reconstructed_Ha",)),
                ("primal feasibility linf", ("primal_feasibility_linf",)),
                ("primal feasibility l2", ("primal_feasibility_l2",)),
                ("BPSDP outer iterations", ("bpsdp_state", "outer_iterations")),
                ("BPSDP inner iterations", ("bpsdp_state", "inner_iterations")),
                ("BPSDP objective gap", ("bpsdp_state", "objective_gap")),
            ]
                value = object_get(solve, keys...)
                value === nothing && continue
                println(io, "| ", label, " | `", value, "` |")
            end
            println(io)
        end

        println(io, "## Decision counts")
        println(io)
        println(io, "| decision | count |")
        println(io, "|---|---:|")
        for decision in sort!(collect(keys(counts)))
            println(io, "| ", decision, " | ", counts[decision], " |")
        end
        println(io)

        eq_rows = [row for row in rows if row["operator"] == "A_eq"]
        full_rows = [row for row in rows if row["operator"] == "A_full"]
        eq_facial = count(row -> row["decision"] == "Facial reduction needed", eq_rows)
        full_facial = count(row -> row["decision"] == "Facial reduction needed", full_rows)
        full_direct = count(row -> row["decision"] == "Direct Cholesky", full_rows)
        termination = solve === nothing ? nothing : object_get(solve, "termination_status")

        println(io, "## Instance-level decision")
        println(io)
        println(io, "- `A_eq`: ", eq_facial, "/", length(eq_rows), " blocks classified as `Facial reduction needed`.")
        println(io, "- `A_full`: ", full_direct, "/", length(full_rows), " blocks classified as `Direct Cholesky`; ", full_facial, "/", length(full_rows), " still require facial reduction/presolve by the Phase-2 rule.")
        if termination !== nothing
            println(io, "- BPSDP PSD-block solve status: `", termination, "`.")
        end
        println(io, "- Verdict: row/nullspace elimination or facial-reduction-style presolve is a blocking workstream before treating H₄/Nk=2 as a plain BPSDP/Mazziotti-CG run. Direct Cholesky remains viable only for the cleaned `A_full` subblocks that survive that presolve.")
        raw_status = solve === nothing ? nothing : object_get(solve, "raw_status")
        if string(raw_status) == "cg_failure"
            println(io, "- Causality note: the BPSDP `cg_failure` is consistent with the measured singular/equality-kernel structure, but this table alone does not prove facial reduction is the only cause; scaling and implementation checks remain follow-up work.")
        else
            println(io, "- Causality note: no immediate BPSDP `cg_failure` was observed; any remaining singularity is residual algebraic/scaling structure rather than orphan-variable fallout.")
        end
        println(io)

        println(io, "## Per-block decision table")
        println(io)
        println(io, "| block | A | m | n | nnz(A) | nnz(AA*) | σ₁ | rank@1e-12 | ker | pred ker | excess | κ_range | gap | AMD nnz(L) | CG est | decision |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|")
        for row in rows
            println(io, "| ", row["block"], " | ", row["operator"], " | ",
                row["m"], " | ", row["n"], " | ", row["nnz_A"], " | ", row["nnz_AAt"], " | ",
                format_float(row["sigma1"]), " | ", row["rank12"], " | ", row["kernel12"], " | ",
                row["predicted_kernel"], " | ", row["excess_kernel"], " | ", format_float(row["kappa12"]), " | ",
                format_float(row["gap12"]), " | ", row["amd_nnz_L"], " | ", format_float(row["cg12"]), " | ",
                row["decision"], " |")
        end
        println(io)

        println(io, "## Moment-matrix rank at the BPSDP optimum")
        println(io)
        termination = solve === nothing ? nothing : object_get(solve, "termination_status")
        if termination !== nothing && string(termination) != "OPTIMAL"
            println(io, "No optimum was obtained (`termination_status = ", termination, "`), so primal moment-rank validation is unavailable for this run.")
            println(io)
        end
        println(io, "| block | rank@1e-10 | rank@1e-12 | rank@1e-14 | min positive eig@1e-12 |")
        println(io, "|---|---:|---:|---:|---:|")
        for dir in block_dirs(input_dir)
            label = basename(dir)
            y = try_json(joinpath(dir, "Y_eigvals.json"))
            y === nothing && continue
            r10 = object_get(y, "rank_by_tol", "1e-10", "rank")
            r12 = object_get(y, "rank_by_tol", "1e-12", "rank")
            r14 = object_get(y, "rank_by_tol", "1e-14", "rank")
            min12 = object_get(y, "rank_by_tol", "1e-12", "smallest_nonzero_eigenvalue")
            println(io, "| ", label, " | ", r10, " | ", r12, " | ", r14, " | ", min12, " |")
        end
        println(io)

        println(io, "## Notes")
        println(io)
        println(io, "- `pred ker` is the construction-level predicted redundancy: exact zero rows plus exact duplicate rows after local block restriction.")
        println(io, "- `excess = empirical kernel@1e-12 - pred ker`; this is the number to chase for missed algebraic/facial-reduction structure.")
        println(io, "- METIS fill was not used; AMD fill is CHOLMOD's ordering on `AA* + δI`.")
        println(io, "- Histogram PNGs live under `plots/sigma_histogram_<block>_<operator>.png`.")
    end
    return path
end

function main(argv = ARGS)
    options = parse_options(argv)
    mkpath(joinpath(options.input_dir, "plots"))

    rows = Dict{String,Any}[]
    for dir in block_dirs(options.input_dir)
        label = basename(dir)
        spectra_out = Dict{String,Any}()
        sparsity_out = Dict{String,Any}()
        kernels = Dict{String,Any}()

        for opname in ("A_full", "A_eq")
            path = joinpath(dir, opname * ".npz")
            isfile(path) || continue
            A = read_sparse_npz(path)
            result = analyze_matrix(A, options.rank_tols)
            decision = decision_with_ram(result, options.ram_bytes)
            spectra_out[opname] = result.spec
            sparsity_out[opname] = result.sparse_stats
            kernels[opname] = result.kernel_basis
            write_histogram_png(joinpath(options.input_dir, "plots", "sigma_histogram_$(label)_$(opname).png"), Float64.(result.spec["singular_values"]))

            empirical_kernel = result.kernel12
            predicted_kernel = result.predicted.predicted_kernel_dim
            push!(rows, Dict{String,Any}(
                "block" => label,
                "operator" => opname,
                "m" => size(A, 1),
                "n" => size(A, 2),
                "nnz_A" => nnz(A),
                "nnz_AAt" => result.sparse_stats["aat_nnz"],
                "sigma1" => result.spec["sigma1"],
                "rank12" => result.rank12,
                "kernel12" => result.kernel12,
                "predicted_kernel" => predicted_kernel,
                "excess_kernel" => max(empirical_kernel - predicted_kernel, 0),
                "kappa12" => result.kappa12,
                "gap12" => result.gap12,
                "amd_nnz_L" => result.sparse_stats["amd_nnz_L"],
                "cg12" => result.cg12,
                "decision" => decision,
            ))
        end

        write_json(joinpath(dir, "spectrum.json"), spectra_out)
        write_json(joinpath(dir, "sparsity.json"), sparsity_out)
        NPZ.npzwrite(joinpath(dir, "kernel_basis.npz"), kernels)
    end

    summary = write_summary(options.input_dir, rows, options.ram_bytes)
    println("wrote ", summary)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
