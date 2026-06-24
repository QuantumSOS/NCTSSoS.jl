#!/usr/bin/env julia

# Structural benchmark for the sparse degree-4 Pauli chain basis used in
# Wang et al. arXiv:2604.01555v1, Table "Maximal block sizes".
#
# This does not call an SDP solver. It measures the symbolic basis and symmetry
# block construction for the paper's hard 1D target:
#   N=100, d=4, r=1 -> sparse basis size 12_001, target max block size 31.
#
# Usage:
#   NCTS_PERF_NS=16,32,100 julia --project=. --startup-file=no perf/pauli_sparse_chain_d4_blocks.jl

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

function _paper_sparse_basis_size(N::Integer, d::Integer)
    # The closed form assumes no periodic-window collisions, i.e. N > d.
    Int(N) > Int(d) || return nothing
    return 1 + 3 * Int(N) * (3^Int(d) - 1) ÷ 2
end

function _paper_symmetry_max_block(d::Integer)
    d = Int(d)
    return isodd(d) ? (3^(d + 1) - 1) ÷ 8 : (3^(d + 1) + 5) ÷ 8
end

_block_variable_count(sizes) = sum(n * (n + 1) ÷ 2 for n in sizes)

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
    expected_basis = _paper_sparse_basis_size(N, degree)
    expected_max = _paper_symmetry_max_block(degree)
    basis_target = isnothing(expected_basis) ? "n/a for N ≤ d" : string(expected_basis)

    println("\n#### Summary")
    println("\n- basis size: `$(length(basis))` (paper sparse target: `$basis_target`)")
    println("- finite group order: `$(length(group))`")
    println("- charge sectors: `$(length(charge_groups))`")
    println("- PSD blocks: `$(length(block_sizes))`")
    println("- largest PSD block: `$(maximum(block_sizes))` (paper symmetry target: `$expected_max`)")
    println("- reduced PSD scalar variables: `$(_block_variable_count(block_sizes))`")
    println("- largest 20 block sizes: `$(block_sizes[max(1, length(block_sizes)-19):end])`")
end

function main()
    ns = parse.(Int, split(get(ENV, "NCTS_PERF_NS", "16,32,100"), ","))
    degree = parse(Int, get(ENV, "NCTS_PERF_DEGREE", "4"))

    println("# Sparse Pauli chain degree-$degree structural benchmark")
    println("\n- generated: `$(Dates.now())`")
    println("- Julia: `$(VERSION)`")
    println("- threads: `$(Threads.nthreads())`")
    println("- CPU: `$(Sys.cpu_info()[1].model)`")
    println("- solver calls: none")

    # Warm a tiny case so the reported sizes are not dominated by first-use JIT.
    _run_one(4; degree=min(degree, 4))
    for N in ns
        _run_one(N; degree)
    end
end

main()
