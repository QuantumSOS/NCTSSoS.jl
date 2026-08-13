#!/usr/bin/env julia

# Mosek-backed solve of the mirror-adapted translation-invariant (DFT) XXX
# Heisenberg relaxation. Defaults to N=16 order 4 as a sanity anchor (COSMO
# reference: -7.1443433630), then the paper-target N=100 order 4; NCTS_TI_NS
# selects any comma-separated set of chain lengths.
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
# NCTS_TI_DUALIZE (default false; set true to reproduce the dual SOS form),
# NCTS_TI_NS (default 16,100), NCTS_TI_RDM (default 0), NCTS_TI_STATE
# (default none), NCTS_TI_STATE_RANGE (default 5), and NCTS_TI_AXIS_SYMMETRY
# (default false).  The practical arXiv:2604.01555 reproduction in the 31 GiB
# orb uses NS=10,20,30,100, RDM=9, STATE=linear_psd, STATE_RANGE=10, and axis
# symmetry; the reference k=10 RDM exceeds the orb memory limit.
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

function run_instance(
    n::Int,
    order::Int;
    dualize::Bool,
    rdm_level::Int,
    state_optimality::Symbol,
    state_optimality_range::Int,
    axis_permutation_symmetry::Bool,
)
    registry, ops = create_pauli_variables(1:n)
    pop = polyopt(heisenberg_chain_hamiltonian(ops), registry)
    GC.gc(true)
    t0 = time()
    res = pauli_translation_invariant_nctssos(
        pop,
        ops,
        order,
        mosek_optimizer();
        dualize,
        rdm_levels=rdm_level == 0 ? Int[] : [rdm_level],
        state_optimality,
        state_optimality_range,
        axis_permutation_symmetry,
    )
    wall = time() - t0
    dual_bound = objective_bound(res.model)
    term_status = termination_status(res.model)
    primal_point = primal_status(res.model)
    dual_point = dual_status(res.model)
    @printf("RESULT N=%d order=%d rdm=%d state=%s state_range=%d axis_symmetry=%s objective=%.10f per_site=%.10f objective_bound=%.10f bound_per_site=%.10f status=%s primal=%s dual=%s wall=%.1fs max_block=%d n_blocks=%d unique_moments=%d\n",
        n, order, rdm_level, string(state_optimality), state_optimality_range,
        string(axis_permutation_symmetry), res.objective, res.objective / n,
        dual_bound, dual_bound / n, string(term_status), string(primal_point),
        string(dual_point), wall,
        maximum(res.report.psd_block_sizes), length(res.report.psd_block_sizes),
        res.report.n_unique_moment_matrix_elements)
    flush(stdout)
    isfinite(res.objective) && isfinite(dual_bound) || error(
        "Mosek returned a non-finite objective for N=$n (termination status $term_status)."
    )
    primal_point == MOI.FEASIBLE_POINT || error(
        "Mosek did not return a primal feasible point for N=$n: $primal_point."
    )
    dual_point == MOI.FEASIBLE_POINT || error(
        "Mosek did not return a dual feasible point for N=$n: $dual_point."
    )
    return res.objective
end

dualize = parse(Bool, get(ENV, "NCTS_TI_DUALIZE", "false"))
ns = parse.(Int, split(get(ENV, "NCTS_TI_NS", "16,100"), ','))
rdm_level = parse(Int, get(ENV, "NCTS_TI_RDM", "0"))
state_optimality = Symbol(get(ENV, "NCTS_TI_STATE", "none"))
state_optimality_range = parse(Int, get(ENV, "NCTS_TI_STATE_RANGE", "5"))
axis_permutation_symmetry = parse(Bool, get(ENV, "NCTS_TI_AXIS_SYMMETRY", "false"))
println("host=$(gethostname()) julia=$(VERSION) started=$(Dates.now())")
println("threads=$(get(ENV, "NCTS_MOSEK_THREADS", "auto(cpu/4)")) cpu_threads=$(Sys.CPU_THREADS) dualize=$dualize ns=$ns rdm=$rdm_level state=$state_optimality state_range=$state_optimality_range axis_symmetry=$axis_permutation_symmetry")
flush(stdout)

for n in ns
    bound = run_instance(
        n,
        4;
        dualize,
        rdm_level,
        state_optimality,
        state_optimality_range,
        axis_permutation_symmetry,
    )
    if n == 16 && rdm_level == 0 && state_optimality == :none && !axis_permutation_symmetry
        delta = abs(bound - (-7.1443433630))
        @printf("SANITY N=16 |mosek - cosmo| = %.2e (cosmo ref -7.1443433630)\n", delta)
        flush(stdout)
        delta <= 1e-4 || error("N=16 Mosek sanity bound differs from the COSMO reference by $delta.")
    end
end
println("done=$(Dates.now())")
