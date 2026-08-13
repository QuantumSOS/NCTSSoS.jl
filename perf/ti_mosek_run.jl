#!/usr/bin/env julia

# Mosek-backed solve of the mirror-adapted translation-invariant (DFT) XXX
# Heisenberg relaxation. Runs N=16 order 4 as a sanity anchor (COSMO
# reference: -7.1443433630), then the paper-target N=100 order 4.
#
# Measured with the default moment formulation (dualize=false), 16 Mosek
# threads, and fresh processes:
#   N=32  bound=-14.2306029477 OPTIMAL,  48.1 s, peak RSS 2.69 GiB (+1.98 GiB)
#   N=48  bound=-21.3325283611 OPTIMAL,  69.8 s, peak RSS 4.08 GiB (+3.37 GiB)
#   N=64  bound=-28.4368184925 OPTIMAL, 110.0 s, peak RSS 5.31 GiB (+4.60 GiB)
#   N=100 bound=-44.4239799523 OPTIMAL, 237.5 s, peak RSS 8.98 GiB (+8.27 GiB)
# Previous dual SOS run (dualize=true, 32 threads): N=16 bound=-7.1443484721
# in 23.4 s; N=100 bound=-44.4239941277 in 1640.5 s with ~270 GiB peak RSS.
#
# Env knobs: NCTS_MOSEK_THREADS (default cpu/4), NCTS_MOSEK_LOG (default 0),
# NCTS_TI_DUALIZE (default false; set true to reproduce the dual SOS form).
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

function run_instance(n::Int, order::Int; dualize::Bool)
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
    GC.gc(true)
    t0 = time()
    res = pauli_translation_invariant_nctssos(pop, ops, order, mosek_optimizer(); dualize)
    wall = time() - t0
    @printf("RESULT N=%d order=%d bound=%.10f per_site=%.10f status=%s wall=%.1fs max_block=%d n_blocks=%d unique_moments=%d\n",
        n, order, res.objective, res.objective / n,
        string(termination_status(res.model)), wall,
        maximum(res.report.psd_block_sizes), length(res.report.psd_block_sizes),
        res.report.n_unique_moment_matrix_elements)
    flush(stdout)
    return res.objective
end

dualize = parse(Bool, get(ENV, "NCTS_TI_DUALIZE", "false"))
println("host=$(gethostname()) julia=$(VERSION) started=$(Dates.now())")
println("threads=$(get(ENV, "NCTS_MOSEK_THREADS", "auto(cpu/4)")) cpu_threads=$(Sys.CPU_THREADS) dualize=$dualize")
flush(stdout)

b16 = run_instance(16, 4; dualize)
@printf("SANITY N=16 |mosek - cosmo| = %.2e (cosmo ref -7.1443433630)\n", abs(b16 - (-7.1443433630)))
flush(stdout)

run_instance(100, 4; dualize)
println("done=$(Dates.now())")
