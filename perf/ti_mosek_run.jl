#!/usr/bin/env julia

# Mosek-backed solve of the mirror-adapted translation-invariant (DFT) XXX
# Heisenberg relaxation. Runs N=16 order 4 as a sanity anchor (COSMO
# reference: -7.1443433630), then the paper-target N=100 order 4.
#
# Measured on 32 Mosek threads (Xeon, 128 cores):
#   N=16  order=4 bound=-7.1443484721  OPTIMAL, 23.4 s
#   N=100 order=4 bound=-44.4239941277 OPTIMAL, 1640.5 s, peak RSS ~270 GiB
#
# Env knobs: NCTS_MOSEK_THREADS (default cpu/4), NCTS_MOSEK_LOG (default 0).
# Self-bootstraps a solver environment in ~/.nctssos_ti_mosek_env with this
# checkout dev'ed, so it can run on a bare host via easy-ssh.

using Pkg
const REPO = dirname(@__DIR__)
Pkg.activate(joinpath(homedir(), ".nctssos_ti_mosek_env"))
Pkg.develop(path=REPO)
Pkg.add(["MosekTools", "JuMP"])
Pkg.instantiate()

using NCTSSoS, JuMP, MosekTools, Printf, Dates
using NCTSSoS: pauli_translation_invariant_nctssos, heisenberg_chain_hamiltonian,
    create_pauli_variables, polyopt

function mosek_optimizer()
    threads = parse(Int, get(ENV, "NCTS_MOSEK_THREADS", string(max(1, div(Sys.CPU_THREADS, 4)))))
    return optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_LOG" => parse(Int, get(ENV, "NCTS_MOSEK_LOG", "0")),
        "MSK_IPAR_NUM_THREADS" => threads,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => 1e-8,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => 1e-8,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => 1e-8,
    )
end

function run_instance(n::Int, order::Int)
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
    GC.gc(true)
    t0 = time()
    res = pauli_translation_invariant_nctssos(pop, ops, order, mosek_optimizer(); dualize=true)
    wall = time() - t0
    @printf("RESULT N=%d order=%d bound=%.10f per_site=%.10f status=%s wall=%.1fs max_block=%d n_blocks=%d unique_moments=%d\n",
        n, order, res.objective, res.objective / n,
        string(termination_status(res.model)), wall,
        maximum(res.report.psd_block_sizes), length(res.report.psd_block_sizes),
        res.report.n_unique_moment_matrix_elements)
    flush(stdout)
    return res.objective
end

println("host=$(gethostname()) julia=$(VERSION) started=$(Dates.now())")
println("threads=$(get(ENV, "NCTS_MOSEK_THREADS", "auto(cpu/4)")) cpu_threads=$(Sys.CPU_THREADS)")
flush(stdout)

b16 = run_instance(16, 4)
@printf("SANITY N=16 |mosek - cosmo| = %.2e (cosmo ref -7.1443433630)\n", abs(b16 - (-7.1443433630)))
flush(stdout)

run_instance(100, 4)
println("done=$(Dates.now())")
