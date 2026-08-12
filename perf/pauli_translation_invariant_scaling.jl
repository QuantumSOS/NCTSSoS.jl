#!/usr/bin/env julia

# Structural scaling benchmark for the specialized translation-invariant
# (momentum/DFT) Pauli-chain relaxation against the paper target of
# Wang et al. arXiv:2604.01555v1:
#   N=100, d=4, r=1 -> sparse basis size 12_001, max real PSD block <= 31.
#
# This does not call an SDP solver; it measures symbolic model construction.
#
# Usage:
#   NCTS_PERF_NS=16,32,64,100 julia --project=. --startup-file=no perf/pauli_translation_invariant_scaling.jl
#   NCTS_PERF_LEGACY=1 additionally measures reflection=false, conjugate_symmetry=false.

using Printf
using NCTSSoS

function _run_one(n::Int, order::Int; reflection::Bool, conjugate_symmetry::Bool)
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops; coupling=1 / 4, periodic=true), registry)
    GC.gc(true)
    stats = @timed pauli_translation_invariant_moment_relaxation(
        pop, ops, order; reflection, conjugate_symmetry
    )
    _, report = stats.value
    @printf(
        "| %d | %d | %s | %d | %d | %d | %d | %d | %.3f | %.3f |\n",
        n, order, reflection && conjugate_symmetry ? "mirror+conj" : "legacy",
        report.basis_size, report.orbit_basis_size,
        length(report.psd_block_sizes), maximum(report.psd_block_sizes),
        report.n_unique_moment_matrix_elements, stats.time, stats.bytes / 2^30,
    )
    flush(stdout)
    return nothing
end

ns = [parse(Int, s) for s in split(get(ENV, "NCTS_PERF_NS", "16,32,64,100"), ",")]
order = parse(Int, get(ENV, "NCTS_PERF_ORDER", "4"))
legacy = get(ENV, "NCTS_PERF_LEGACY", "0") == "1"

println("| N | d | rules | basis | orbit reps | blocks | max block | moments | time (s) | alloc (GiB) |")
println("|--:|--:|:--|--:|--:|--:|--:|--:|--:|--:|")
_run_one(4, min(order, 2); reflection=true, conjugate_symmetry=true)  # warmup
for n in ns
    legacy && _run_one(n, order; reflection=false, conjugate_symmetry=false)
    _run_one(n, order; reflection=true, conjugate_symmetry=true)
end
