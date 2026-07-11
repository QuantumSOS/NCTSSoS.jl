#!/usr/bin/env julia

# Structural benchmark for the sparse degree-4 Pauli chain basis used in
# Wang et al. arXiv:2604.01555v1, Table "Maximal block sizes".
#
# This does not call an SDP solver. It measures the symbolic basis and symmetry
# block construction for the paper's hard 1D target:
#   N=100, d=4, r=1 -> sparse basis size 12_001, target max block size 31.
#
# Usage:
#   NCTS_PERF_NS=8,12 julia --project=. --startup-file=no perf/pauli_sparse_chain_d4_blocks.jl
#   NCTS_PERF_ALLOW_LARGE=true NCTS_PERF_TARGET_ONLY=true NCTS_PERF_NS=100 NCTS_PERF_ESTIMATED_WALL_SECONDS=<estimate-before-launch> NCTS_PERF_ESTIMATED_RSS_GIB=<estimate-before-launch> julia --project=. --startup-file=no perf/pauli_sparse_chain_d4_blocks.jl
#
# Safety:
#   By default this refuses N > 13.  For intentional structural benchmarks,
#   set NCTS_PERF_ALLOW_LARGE=true only with explicit wall/RSS estimates and
#   safe load/memory telemetry.  Target-only mode is checked too; use the same
#   pressure gate for large analytic shape reports.

Base.include(@__MODULE__, joinpath(@__DIR__, "shared_load_guard.jl"))

function _sparse_chain_preimport_large_run_pressure_guard()
    _ncts_load_guard_parse_bool("NCTS_PERF_ALLOW_LARGE", false) ||
        return nothing
    ns_label = strip(get(ENV, "NCTS_PERF_NS", "8,12"))
    _ncts_check_large_run_pressure_guard(
        env_prefix="NCTS_PERF",
        label="sparse Pauli chain d4 Ns=$ns_label",
    )
    return nothing
end

_sparse_chain_preimport_large_run_pressure_guard()

using Dates
using Printf
using JuMP
using NCTSSoS

function _fmt_bytes(bytes::Integer)
    units = ("B", "KiB", "MiB", "GiB", "TiB")
    value = Float64(bytes)
    unit = 1
    while value >= 1024 && unit < length(units)
        value /= 1024
        unit += 1
    end
    return @sprintf("%.2f %s", value, units[unit])
end

function _timed(label::AbstractString, f)
    GC.gc(true)
    stats = @timed f()
    @printf("| `%s` | %.6f | %s | %.6f |\n", label, stats.time, _fmt_bytes(stats.bytes), stats.gctime)
    flush(stdout)
    return stats.value
end

function _parse_bool_env(name::AbstractString, default::Bool)
    raw = lowercase(strip(get(ENV, name, string(default))))
    raw in ("true", "1", "yes", "y") && return true
    raw in ("false", "0", "no", "n") && return false
    throw(ArgumentError("$name must be a boolean value, got $(repr(raw))."))
end

function _parse_int_list_env(name::AbstractString, default::AbstractString)
    values = Int[]
    for raw in split(get(ENV, name, default), ",")
        item = strip(raw)
        isempty(item) && continue
        push!(values, parse(Int, item))
    end
    isempty(values) && throw(ArgumentError("$name did not contain any integer sizes."))
    return values
end

function _check_size_guard(ns::AbstractVector{Int})
    allow_large = _parse_bool_env("NCTS_PERF_ALLOW_LARGE", false)
    max_n = parse(Int, get(ENV, "NCTS_PERF_MAX_N", "13"))
    if allow_large
        _ncts_check_large_run_pressure_guard(
            env_prefix="NCTS_PERF",
            label="sparse Pauli chain d4 Ns=$(join(ns, ","))",
        )
        return nothing
    end
    for n in ns
        n <= max_n || throw(ArgumentError(
            "Refusing N=$n because NCTS_PERF_ALLOW_LARGE=false and " *
            "NCTS_PERF_MAX_N=$max_n. Set NCTS_PERF_ALLOW_LARGE=true " *
            "only for intentional large structural runs."
        ))
    end
    return nothing
end

_block_variable_count(sizes) = sum(n * (n + 1) ÷ 2 for n in sizes)

function _print_target_only(N::Integer; degree::Integer=4)
    target = pauli_translation_structural_targets(N, degree)

    println("\n## N = $N, d = $degree")
    println("\n#### Analytic translation-DFT target")
    println("\n- construction performed: `false`")
    println("- basis size: `$(target.basis_size)`")
    println("- orbit basis size: `$(target.orbit_basis_size)`")
    println("- PSD blocks: `$(target.n_blocks)`")
    println("- largest logical block: `$(target.logical_max_block)`")
    println("- largest solver-facing block: `$(target.psd_max_block)`")
    println("- solver-facing symmetric entries: `$(target.psd_symmetric_entries)`")
    println("- solver-facing symmetric bytes: `$(_fmt_bytes(target.psd_symmetric_bytes))`")
    println("- product-cache hit rate: `$(@sprintf("%.6f", target.product_cache_hit_rate))`")
    println("- max closed contiguous RDM k: `$(target.max_contiguous_rdm_k)`")
    println("- max linear state-opt width: `$(target.max_linear_state_opt_width)`")
    println("- max PSD state-opt width: `$(target.max_psd_state_opt_width)`")
    println("- solve supported: `$(target.solve_supported)`")
    println("- solve blocker: `$(target.solve_blocker)`")
    println("- solve blocker reason: `$(target.solve_blocker_reason)`")
    println("- solve unsupported block features: `$(target.solve_unsupported_block_features)`")
    println("- solve unsupported zero-row features: `$(target.solve_unsupported_zero_features)`")
    println("- logical block histogram: `$(target.logical_block_histogram)`")
    println("- solver-facing block histogram: `$(target.psd_block_histogram)`")
    println("- axis-orbit size histogram: `$(target.axis_orbit_size_histogram)`")
end

function _run_one(N::Integer; degree::Integer=4)
    println("\n## N = $N, d = $degree")
    println("\n| phase | wall time (s) | allocated | GC time (s) |")
    println("|:--|--:|--:|--:|")

    registry = _timed("create Pauli registry", () -> first(create_pauli_variables(1:Int(N))))
    basis = _timed("pauli_contiguous_chain_basis", () -> pauli_contiguous_chain_basis(registry, Int(degree)))
    T = eltype(basis[1].word)

    group = _timed("build translation/reflection/sign group", () -> begin
        translation = pauli_site_permutation([2:Int(N); 1])
        reflection = pauli_site_permutation(reverse(1:Int(N)))
        sign = pauli_sign_symmetry(Int(N); integer_type=T)
        CliffordSymmetryGroup([translation, reflection, sign]; nqubits=Int(N), integer_type=T)
    end)

    charge_groups = _timed("charge/spatial/sign block decomposition", () -> NCTSSoS._pauli_charge_transform_groups(
        basis,
        PauliChargeSectorSpec(nqubits=Int(N), max_degree=Int(degree)),
        group,
    ))

    block_sizes = sort([size(block.row_basis, 1) for group in charge_groups for block in group])
    target = nothing
    if 2 * Int(degree) < Int(N)
        target = pauli_translation_structural_targets(N, degree)
    end
    basis_target = target === nothing ? "n/a for N ≤ 2d" : string(target.basis_size)
    max_target = target === nothing ? "n/a for N ≤ 2d" : string(target.logical_max_block)

    println("\n#### Summary")
    println("\n- basis size: `$(length(basis))` (paper sparse target: `$basis_target`)")
    println("- finite group order: `$(length(group))`")
    println("- charge sectors: `$(length(charge_groups))`")
    println("- PSD blocks: `$(length(block_sizes))`")
    println("- largest PSD block: `$(maximum(block_sizes))` (translation-DFT logical target: `$max_target`)")
    println("- reduced PSD scalar variables: `$(_block_variable_count(block_sizes))`")
    println("- largest 20 block sizes: `$(block_sizes[max(1, length(block_sizes)-19):end])`")
end

function main()
    ns = _parse_int_list_env("NCTS_PERF_NS", "8,12")
    degree = parse(Int, get(ENV, "NCTS_PERF_DEGREE", "4"))
    target_only = _parse_bool_env("NCTS_PERF_TARGET_ONLY", false)
    _check_size_guard(ns)

    println("# Sparse Pauli chain degree-$degree structural benchmark")
    println("\n- generated: `$(Dates.now())`")
    println("- Julia: `$(VERSION)`")
    println("- threads: `$(Threads.nthreads())`")
    println("- CPU: `$(Sys.cpu_info()[1].model)`")
    println("- solver calls: none")
    println("- target_only: `$target_only`")

    if target_only
        for N in ns
            _print_target_only(N; degree)
        end
        return nothing
    end

    # Warm a tiny case so the reported sizes are not dominated by first-use JIT.
    _run_one(4; degree=min(degree, 3))
    for N in ns
        _run_one(N; degree)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
