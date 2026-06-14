#!/usr/bin/env julia

using Printf
using NCTSSoS
using ManiDSDP

const MOI = NCTSSoS.MOI
const CHAIN_N = parse(Int, get(ENV, "HEISENBERG_N", "16"))
const RELAX_ORDER = parse(Int, get(ENV, "HEISENBERG_ORDER", "3"))

function heisenberg_periodic_pop(n::Integer)
    registry, (σx, σy, σz) = create_pauli_variables(1:n)
    H = sum(ComplexF64(0.25) * op[i] * op[mod1(i + 1, n)]
            for op in (σx, σy, σz) for i in 1:n)
    return polyopt(H, registry)
end

function ts_algorithm(name::Symbol)
    name === :none && return NoElimination()
    name === :mmd && return MMD()
    name === :mf && return MF()
    error("unsupported ts=$name")
end

println("Heisenberg sparsity probe N=", CHAIN_N, " order=", RELAX_ORDER)
flush(stdout)
t_pop = @elapsed pop = heisenberg_periodic_pop(CHAIN_N)
@printf("pop %.3fs\n", t_pop)
flush(stdout)

for ts in (:none, :mmd, :mf)
    GC.gc()
    optimizer = MOI.OptimizerWithAttributes(ManiDSDP.Optimizer, MOI.RawOptimizerAttribute("rank") => 1, MOI.Silent() => true)
    cfg = SolverConfig(optimizer=optimizer, order=RELAX_ORDER, ts_algo=ts_algorithm(ts))
    t = @elapsed sparsity = compute_sparsity(pop, cfg)
    sizes = NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities)
    total_dim = sum(sum(block) for block in sizes)
    max_block = maximum(vcat(sizes...))
    @printf("ts=%s compute=%.3fs block_sizes=%s total_dim=%d max_block=%d\n",
            String(ts), t, sprint(show, sizes), total_dim, max_block)
    flush(stdout)
end
