#!/usr/bin/env julia
# Fast structured H4 periodic PQG V2RDM exporter.
#
# This bypasses the generic Polynomial/moment_relax construction and emits the
# same libsdp/BPSDP `.dat-c` problem from closed-form PQG index algebra.

using LinearAlgebra
using NCTSSoS
using Printf

include(joinpath(@__DIR__, "..", "demos", "SDPALibsdpExport.jl"))
using .SDPALibsdpExport: export_libsdp

struct Options
    nk::Int
    norb::Int
    nelec_per_cell::Int
    integrals_path::Union{Nothing,String}
    blocking::Symbol
    include_one_d::Bool
    spin_resolved_trace::Bool
    singlet_s2::Bool
    outdir::String
    basename::String
    sense::Symbol
    orphans_per_block::Int
end

function parse_options(argv)
    nk = 2
    norb = 4
    nelec_per_cell = 4
    integrals_path = nothing
    blocking = :momentum
    include_one_d = false
    spin_resolved_trace = true
    singlet_s2 = true
    outdir = normpath(joinpath(@__DIR__, "..", "output", "h4_nk_structured"))
    basename = ""
    sense = :min
    orphans_per_block = 32

    for arg in argv
        if startswith(arg, "--nk=")
            nk = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--norb=")
            norb = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--nelec-per-cell=")
            nelec_per_cell = parse(Int, split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--integrals=")
            value = String(split(arg, "=", limit = 2)[2])
            integrals_path = isempty(value) || value == "mock" ? nothing : value
        elseif startswith(arg, "--blocking=")
            value = split(arg, "=", limit = 2)[2]
            if value in ("momentum", "k")
                blocking = :momentum
            elseif value == "spin"
                blocking = :spin
            elseif value == "none"
                blocking = :none
            else
                throw(ArgumentError("unknown --blocking=$value; expected momentum, spin, or none"))
            end
        elseif arg == "--include-1d"
            include_one_d = true
        elseif arg == "--no-1d"
            include_one_d = false
        elseif arg in ("--paper-spin", "--spin-singlet")
            spin_resolved_trace = true
            singlet_s2 = true
        elseif arg == "--spin-resolved-trace"
            spin_resolved_trace = true
        elseif arg == "--no-spin-resolved-trace"
            spin_resolved_trace = false
        elseif arg == "--singlet-s2"
            singlet_s2 = true
        elseif arg == "--no-singlet-s2"
            singlet_s2 = false
        elseif startswith(arg, "--outdir=")
            outdir = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--basename=")
            basename = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--sense=")
            sense = Symbol(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--orphans-per-block=")
            orphans_per_block = parse(Int, split(arg, "=", limit = 2)[2])
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=. probes/h4_nk_periodic_structured_export.jl [options]\n")
            println("Core options:")
            println("  --nk=N                         k-point count; default 2")
            println("  --norb=N                       active spatial orbitals per k; default 4")
            println("  --nelec-per-cell=N             active electrons per cell; default 4")
            println("  --integrals=PATH|mock          PySCF .txt or .npz integral dump; default mock")
            println("  --blocking=momentum|spin|none  PQG block grouping; default momentum")
            println("  --include-1d | --no-1d         include/drop 1D PSD blocks; default drop")
            println("  --spin-resolved-trace | --no-spin-resolved-trace")
            println("  --singlet-s2 | --no-singlet-s2")
            println("Export options:")
            println("  --outdir=PATH                  output directory; default output/h4_nk_structured")
            println("  --basename=STEM                output stem; default h4_nk<N>_<source>_structured_pqg")
            println("  --sense=min|max                default min")
            println("  --orphans-per-block=N          default 32")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    nk >= 1 || throw(ArgumentError("--nk must be positive"))
    norb >= 1 || throw(ArgumentError("--norb must be positive"))
    nelec_per_cell >= 1 || throw(ArgumentError("--nelec-per-cell must be positive"))
    nelec_per_cell <= 2 * norb || throw(ArgumentError("--nelec-per-cell cannot exceed 2*--norb spin capacity"))
    iseven(nelec_per_cell) || throw(ArgumentError("--nelec-per-cell must be even for the closed-shell active-HF summary and spin-singlet constraints"))
    sense in (:min, :max) || throw(ArgumentError("--sense must be min or max"))
    orphans_per_block >= 1 || throw(ArgumentError("--orphans-per-block must be >= 1"))

    if isempty(basename)
        source = integrals_path === nothing ? "mock" : splitext(Base.basename(String(integrals_path)))[1]
        basename = "h4_nk$(nk)_$(source)_structured_pqg"
    end

    return Options(nk, norb, nelec_per_cell, integrals_path, blocking,
        include_one_d, spin_resolved_trace, singlet_s2, outdir, basename,
        sense, orphans_per_block)
end

bytes_to_gib(bytes::Integer) = bytes / 1024.0^3
rss_string() = @sprintf("%.3f GiB", bytes_to_gib(Sys.maxrss()))

@inline complex_entry(fields, re_idx::Int, im_idx::Int) =
    complex(parse(Float64, fields[re_idx]), parse(Float64, fields[im_idx]))

function load_integrals_txt(path::AbstractString; nk::Int, norb::Int)
    h1e = Dict(k => zeros(ComplexF64, norb, norb) for k in 0:(nk - 1))
    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()

    for (line_no, raw_line) in enumerate(eachline(path))
        line = strip(raw_line)
        (isempty(line) || startswith(line, '#')) && continue

        fields = split(line)
        tag = first(fields)
        if tag == "h1e"
            length(fields) == 6 || error("malformed h1e line $line_no in $(repr(path))")
            k = parse(Int, fields[2])
            0 <= k < nk || error("h1e line $line_no has k=$k outside 0:$(nk - 1)")
            p = parse(Int, fields[3]) + 1
            q = parse(Int, fields[4]) + 1
            1 <= p <= norb && 1 <= q <= norb || error("h1e line $line_no has orbital outside 0:$(norb - 1)")
            try
                h1e[k][p, q] = complex_entry(fields, 5, 6)
            catch err
                error("failed parsing h1e numeric value at line $line_no in $(repr(path)): $(err)")
            end
        elseif tag == "eri"
            length(fields) == 11 || error("malformed eri line $line_no in $(repr(path))")
            k1 = parse(Int, fields[2])
            k2 = parse(Int, fields[3])
            k3 = parse(Int, fields[4])
            k4 = parse(Int, fields[5])
            all(k -> 0 <= k < nk, (k1, k2, k3, k4)) ||
                error("eri line $line_no has k outside 0:$(nk - 1)")
            p = parse(Int, fields[6]) + 1
            r = parse(Int, fields[7]) + 1
            q = parse(Int, fields[8]) + 1
            s = parse(Int, fields[9]) + 1
            all(i -> 1 <= i <= norb, (p, r, q, s)) ||
                error("eri line $line_no has orbital outside 0:$(norb - 1)")
            block = get!(eri, (k1, k2, k3, k4)) do
                zeros(ComplexF64, norb, norb, norb, norb)
            end
            try
                block[p, r, q, s] = complex_entry(fields, 10, 11)
            catch err
                error("failed parsing eri numeric value at line $line_no in $(repr(path)): $(err)")
            end
        else
            error("unexpected integral tag $(repr(tag)) at line $line_no in $(repr(path))")
        end
    end

    return h1e, eri
end

function load_integrals_npz(path::AbstractString; nk::Int, norb::Int)
    try
        @eval import NPZ
    catch err
        error("Loading .npz integrals requires NPZ.jl in the active environment. Use the .txt dump or add NPZ. Original error: $(err)")
    end

    data = NPZ.npzread(path)
    h1e = Dict{Int,Matrix{ComplexF64}}()
    for k in 0:(nk - 1)
        key = "h1e_k$(k)"
        haskey(data, key) || error("npz file $(repr(path)) is missing $key")
        h1e[k] = ComplexF64.(data[key][1:norb, 1:norb])
    end

    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()
    for k1 in 0:(nk - 1), k2 in 0:(nk - 1), k3 in 0:(nk - 1), k4 in 0:(nk - 1)
        key = "eri_$(k1)_$(k2)_$(k3)_$(k4)"
        haskey(data, key) || continue
        eri[(k1, k2, k3, k4)] = ComplexF64.(data[key][1:norb, 1:norb, 1:norb, 1:norb])
    end

    return h1e, eri
end

function load_integrals(path::AbstractString; nk::Int, norb::Int)
    ext = lowercase(splitext(path)[2])
    if ext == ".npz"
        return load_integrals_npz(path; nk, norb)
    else
        return load_integrals_txt(path; nk, norb)
    end
end

function mock_h4_integrals(; nk::Int, norb::Int)
    h1e = Dict{Int,Matrix{ComplexF64}}()
    for k in 0:(nk - 1)
        θ = 2π * k / nk
        H = zeros(ComplexF64, norb, norb)
        for p in 1:norb
            H[p, p] = -1.30 - 0.11p + 0.035cos(θ + 0.37p)
        end
        for p in 1:(norb - 1)
            amp = -0.055 / (1 + 0.2p)
            val = amp * cis(θ * (p % 2 == 0 ? -1 : 1))
            H[p, p + 1] = val
            H[p + 1, p] = conj(val)
        end
        h1e[k] = H
    end

    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()
    for k1 in 0:(nk - 1), k2 in 0:(nk - 1), k3 in 0:(nk - 1)
        k4 = mod(k1 + k2 - k3, nk)
        V = zeros(ComplexF64, norb, norb, norb, norb)
        transfer = min(mod(k1 - k3, nk), mod(k3 - k1, nk))
        qscale = 1.0 / (1.0 + transfer)
        for p in 1:norb, q in 1:norb
            direct = 0.18 * qscale / (1.0 + abs(p - q))
            exchange = 0.045 * qscale / (1.0 + abs(p - q))
            V[p, p, q, q] += direct
            V[p, q, q, p] -= exchange
        end
        if any(!iszero, V)
            eri[(k1, k2, k3, k4)] = V
        end
    end
    return h1e, eri
end

function active_hf_energy(h1e, eri; nk::Int, nelec_per_cell::Int)
    nocc = nelec_per_cell ÷ 2
    e1 = 0.0 + 0.0im
    for k in 0:(nk - 1), i in 1:nocc
        e1 += 2.0 * h1e[k][i, i] / nk
    end

    e2 = 0.0 + 0.0im
    for ((k1, k2, k3, k4), Vraw) in eri
        V = Vraw / nk^2
        for i in 1:nocc, j in 1:nocc
            if k1 == k3 && k2 == k4
                e2 += 2.0 * V[i, i, j, j]
            end
            if k1 == k4 && k2 == k3
                e2 -= V[i, j, j, i]
            end
        end
    end
    return real(e1 + e2)
end

function print_summary(linear, options::Options, hf_active::Real, build_seconds::Real)
    block_dims = [block.size for block in linear.psd_blocks_lin]
    println("== H4/Nk=$(options.nk) structured periodic PQG V2RDM ==")
    @printf("%-44s %s\n", "integrals", options.integrals_path === nothing ? "Julia deterministic mock" : options.integrals_path)
    @printf("%-44s %s\n", "blocking", string(options.blocking))
    @printf("%-44s %s\n", "extra ¹D PSD block", options.include_one_d ? "included" : "disabled")
    @printf("%-44s %s\n", "spin-resolved D traces", string(options.spin_resolved_trace))
    @printf("%-44s %s\n", "singlet S² constraint", string(options.singlet_s2))
    println()
    @printf("%-44s %d\n", "k-points", options.nk)
    @printf("%-44s %d\n", "active spatial orbitals per k", options.norb)
    @printf("%-44s %d\n", "spin-orbital modes", 2 * options.nk * options.norb)
    @printf("%-44s %d\n", "active electrons total", options.nk * options.nelec_per_cell)
    @printf("%-44s %.12f Ha\n", "active HF energy", hf_active)
    println()
    @printf("%-44s %d\n", "HPSD blocks", length(block_dims))
    @printf("%-44s %d\n", "largest HPSD block", maximum(block_dims))
    @printf("%-44s %d\n", "Σ block n²", sum(n -> n * n, block_dims; init = 0))
    @printf("%-44s %d\n", "canonical moments", length(linear.moments))
    @printf("%-44s %d\n", "HPSD pivots", length(linear.pivots))
    @printf("%-44s %d\n", "free/orphan moment keys", length(linear.free_keys))
    @printf("%-44s %d\n", "zero scalar rows", length(linear.zero_constraints))
    @printf("%-44s %d\n", "objective terms", length(linear.objective_lin))
    @printf("%-44s %.3f s\n", "structured build walltime", build_seconds)
    @printf("%-44s %s\n", "max RSS after build", rss_string())
    println()
    return nothing
end

function main(argv = ARGS)
    options = parse_options(argv)
    @printf("%-44s %s\n", "max RSS at start", rss_string())

    h1e, eri = options.integrals_path === nothing ?
        mock_h4_integrals(; nk = options.nk, norb = options.norb) :
        load_integrals(options.integrals_path; nk = options.nk, norb = options.norb)
    hf_active = active_hf_energy(h1e, eri; nk = options.nk, nelec_per_cell = options.nelec_per_cell)

    build_seconds = @elapsed linear = build_pqg_moment_data(
        h1e,
        eri;
        nk = options.nk,
        norb = options.norb,
        nelec_per_cell = options.nelec_per_cell,
        blocking = options.blocking,
        spin_resolved_trace = options.spin_resolved_trace,
        singlet_s2 = options.singlet_s2,
        include_one_d = options.include_one_d,
    )

    print_summary(linear, options, hf_active, build_seconds)

    @printf("%-44s %s\n", "writing libsdp export to", options.outdir)
    @printf("%-44s %s\n", "basename", options.basename)
    @printf("%-44s %s\n", "sense", string(options.sense))
    @printf("%-44s %d\n", "orphans per aux block", options.orphans_per_block)
    flush(stdout)

    export_seconds = @elapsed summary = export_libsdp(
        linear;
        outdir = options.outdir,
        basename = options.basename,
        sense = options.sense,
        orphans_per_block = options.orphans_per_block,
        legacy_zero_split = true,
    )

    dat_path = joinpath(options.outdir, options.basename * ".dat-c")
    meta_path = joinpath(options.outdir, options.basename * "_meta.txt")
    println()
    println("== libsdp export summary ==")
    @printf("%-44s %s\n", "SDPA file", dat_path)
    @printf("%-44s %s\n", "metadata", meta_path)
    @printf("%-44s %d (HPSD=%d, aux=%d)\n", "PSD blocks", summary.n_blocks, summary.n_hpsd_blocks, summary.n_aux_blocks)
    @printf("%-44s %d / %d\n", "min/max block dim", minimum(summary.block_dims), maximum(summary.block_dims))
    @printf("%-44s %d\n", "unique canonical moments", summary.n_unique_moments)
    @printf("%-44s %d\n", "equality constraints", summary.n_constraints)
    @printf("%-44s %d\n", "objective sparse entries", summary.n_objective_entries)
    @printf("%-44s %d\n", "constraint sparse entries", summary.n_constraint_entries)
    @printf("%-44s %.3f s\n", "export walltime", export_seconds)
    @printf("%-44s %s\n", "max RSS after export", rss_string())
    return summary
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
