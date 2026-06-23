#!/usr/bin/env julia

# Mosek-backed 1D periodic Heisenberg scaling benchmark.
#
# This intentionally decomposes the NCTSSoS pipeline into phases instead of
# calling cs_nctssos directly, so setup/model-build/solver costs are visible.
#
# Usage examples:
#   NCTS_PERF_MODE=order2_singlet NCTS_PERF_NS=8,10,12 julia --project=. --startup-file=no perf/heisenberg_mosek_scaling.jl
#   NCTS_PERF_MODE=sparse_d4 NCTS_PERF_NS=8,10,12,16 julia --project=. --startup-file=no perf/heisenberg_mosek_scaling.jl
#
# Modes:
#   order2_singlet  automatic order-2 basis, translation/reflection, charge sectors, singlet constraints
#   sparse_d4       contiguous degree-4 basis, translation/reflection/sign, charge sectors

using Dates
using Printf
using JuMP
using NCTSSoS
using MosekTools

const MOI = JuMP.MOI

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

function _proc_status_kb(field::AbstractString)
    path = "/proc/self/status"
    isfile(path) || return nothing
    prefix = field * ":"
    for line in eachline(path)
        startswith(line, prefix) || continue
        parts = split(line)
        length(parts) >= 2 || return nothing
        return parse(Int, parts[2])
    end
    return nothing
end

function _rss_string()
    rss = _proc_status_kb("VmRSS")
    hwm = _proc_status_kb("VmHWM")
    isnothing(rss) && isnothing(hwm) && return "n/a"
    rss_s = isnothing(rss) ? "n/a" : _fmt_bytes(rss * 1024)
    hwm_s = isnothing(hwm) ? "n/a" : _fmt_bytes(hwm * 1024)
    return "rss=$rss_s, hwm=$hwm_s"
end

function _timed(label::AbstractString, f)
    GC.gc(true)
    stats = @timed f()
    @printf("| `%s` | %.6f | %s | %.6f | %s |\n", label, stats.time, _fmt_bytes(stats.bytes), stats.gctime, _rss_string())
    flush(stdout)
    return stats.value, stats
end

function _safe_show(x)
    try
        return string(x)
    catch err
        return "<show failed: $(typeof(err))>"
    end
end

function _safe_int(f, default=-1)
    try
        return Int(f())
    catch
        return default
    end
end

function _block_variable_count(sizes)
    return sum(Int(n) * (Int(n) + 1) ÷ 2 for n in sizes)
end

function _largest_tail(sizes; k::Int=20)
    isempty(sizes) && return Int[]
    sorted = sort(Int.(sizes))
    return sorted[max(1, length(sorted) - k + 1):end]
end

function _make_solver()
    # Keep logging off; use half the machine by default so the benchmark does not
    # turn a shared box into soup. Override with NCTS_MOSEK_THREADS if needed.
    default_threads = max(1, Sys.CPU_THREADS ÷ 2)
    threads = parse(Int, get(ENV, "NCTS_MOSEK_THREADS", string(default_threads)))
    return optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_LOG" => 0,
        "MSK_IPAR_NUM_THREADS" => threads,
    )
end

function _heisenberg_pop(N::Integer)
    registry, (σx, σy, σz) = create_pauli_variables(1:Int(N))
    H = sum(
        ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, Int(N))]
        for op in (σx, σy, σz) for i in 1:Int(N)
    )
    return (; registry, H, pop=polyopt(H, registry))
end

function _translation_reflection(N::Integer)
    n = Int(N)
    return (
        pauli_site_permutation([2:n; 1]),
        pauli_site_permutation(reverse(1:n)),
    )
end

function _build_order2_singlet_case(N::Integer; offblock_check::Symbol, check_invariance::Bool)
    base = _heisenberg_pop(N)
    translation, reflection = _translation_reflection(N)
    symmetry = SymmetrySpec(
        [translation, reflection];
        pauli_charge=PauliChargeSectorSpec(nqubits=Int(N), max_degree=2),
        pauli_singlet=PauliSingletConstraintSpec(nqubits=Int(N), max_degree=2),
        check_invariance=check_invariance,
        offblock_check=offblock_check,
    )
    cfg = SolverConfig(
        optimizer=_make_solver(),
        order=2,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry,
    )
    return (; N=Int(N), mode=:order2_singlet, degree=2, base..., moment_basis=nothing, symmetry, cfg)
end

function _build_sparse_d4_case(N::Integer; degree::Integer, offblock_check::Symbol, check_invariance::Bool)
    base = _heisenberg_pop(N)
    basis = pauli_contiguous_chain_basis(base.registry, Int(degree))
    T = eltype(basis[1].word)
    translation, reflection = _translation_reflection(N)
    sign = pauli_sign_symmetry(Int(N); integer_type=T)
    symmetry = SymmetrySpec(
        [translation, reflection, sign];
        pauli_charge=PauliChargeSectorSpec(nqubits=Int(N), max_degree=Int(degree)),
        check_invariance=check_invariance,
        offblock_check=offblock_check,
    )
    cfg = SolverConfig(
        optimizer=_make_solver(),
        moment_basis=basis,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry,
    )
    return (; N=Int(N), mode=:sparse_d4, degree=Int(degree), base..., moment_basis=basis, symmetry, cfg)
end

function _build_case(N::Integer, mode::Symbol; degree::Integer, offblock_check::Symbol, check_invariance::Bool)
    if mode === :order2_singlet
        return _build_order2_singlet_case(N; offblock_check, check_invariance)
    elseif mode === :sparse_d4
        return _build_sparse_d4_case(N; degree, offblock_check, check_invariance)
    else
        throw(ArgumentError("unknown NCTS_PERF_MODE=$(repr(mode)); use order2_singlet or sparse_d4"))
    end
end

function _model_counts(model)
    nvars = _safe_int(() -> JuMP.num_variables(model))
    ncons = _safe_int(() -> JuMP.num_constraints(model; count_variable_in_set_constraints=true))
    return (; nvars, ncons)
end

function _objective_value(model)
    try
        return JuMP.objective_value(model)
    catch err
        @warn "objective_value unavailable" exception=(err, catch_backtrace())
        return NaN
    end
end

function _run_one(N::Integer, mode::Symbol; degree::Integer, offblock_check::Symbol, check_invariance::Bool)
    println("\n## N = $N")
    println("\n| phase | wall time (s) | allocated | GC time (s) | process memory |")
    println("|:--|--:|--:|--:|:--|")

    case, setup_stats = _timed("setup Heisenberg case + basis/symmetry", () -> _build_case(
        N,
        mode;
        degree,
        offblock_check,
        check_invariance,
    ))

    sparsity, sparsity_stats = _timed("compute_sparsity", () -> compute_sparsity(case.pop, case.cfg))
    pre_basis = only(sparsity.corr_sparsity.clq_mom_mtx_bases)
    pre_basis_size = length(pre_basis)

    mp_report, relax_stats = _timed("moment_relax_symmetric", () -> NCTSSoS.moment_relax_symmetric(
        case.pop,
        sparsity.corr_sparsity,
        sparsity.cliques_term_sparsities,
        case.symmetry,
    ))
    mp, report = mp_report

    sos, model_stats = _timed("sos_dualize (JuMP/MOI model build)", () -> NCTSSoS.sos_dualize(mp))

    _, setopt_stats = _timed("set Mosek optimizer", () -> begin
        JuMP.set_optimizer(sos.model, case.cfg.optimizer)
        nothing
    end)

    _, solve_stats = _timed("Mosek optimize!", () -> begin
        JuMP.optimize!(sos.model)
        nothing
    end)

    status = JuMP.termination_status(sos.model)
    primal_status = JuMP.primal_status(sos.model)
    dual_status = JuMP.dual_status(sos.model)
    objective = _objective_value(sos.model)
    solve_time = try
        JuMP.solve_time(sos.model)
    catch
        solve_stats.time
    end
    counts = _model_counts(sos.model)
    block_sizes = Int.(report.psd_block_sizes)
    pre_block_sizes = length.(first(first(sparsity.cliques_term_sparsities)).block_bases)
    pre_psd_vars = _block_variable_count(pre_block_sizes)
    reduced_psd_vars = _block_variable_count(block_sizes)
    total_wall = setup_stats.time + sparsity_stats.time + relax_stats.time + model_stats.time + setopt_stats.time + solve_stats.time

    println("\n#### Summary")
    println("\n- mode: `$(case.mode)`")
    println("- degree/order: `$(case.degree)`")
    println("- pre-symmetry moment basis size: `$pre_basis_size`")
    if !isnothing(case.moment_basis)
        println("- supplied contiguous basis size: `$(length(case.moment_basis))`")
    end
    println("- pre-symmetry dense PSD scalar variables estimate: `$pre_psd_vars`")
    println("- symmetry group order: `$(report.group_order)`")
    println("- PSD blocks: `$(length(block_sizes))`")
    println("- largest PSD block: `$(isempty(block_sizes) ? 0 : maximum(block_sizes))`")
    println("- reduced PSD scalar variables: `$reduced_psd_vars`")
    println("- largest 20 PSD block sizes: `$(_largest_tail(block_sizes))`")
    println("- invariant moment count: `$(report.invariant_moment_count)`")
    println("- symbolic constraints: `$(length(mp.constraints))`")
    println("- linear moments: `$(length(mp.linear.moments))`")
    println("- unique moment-matrix elements: `$(mp.n_unique_moment_matrix_elements)`")
    println("- JuMP variables: `$(counts.nvars)`")
    println("- JuMP constraints: `$(counts.ncons)`")
    println("- termination status: `$(status)`")
    println("- primal status: `$(primal_status)`")
    println("- dual status: `$(dual_status)`")
    println("- objective: `$objective`")
    println("- objective per spin: `$(objective / Int(N))`")
    println("- JuMP-reported solve time: `$solve_time` seconds")
    println("- measured total wall time: `$total_wall` seconds")
    println("- process memory after solve: `$(_rss_string())`")

    return (; N=Int(N), mode=case.mode, objective, objective_per_spin=objective / Int(N), status, total_wall,
        setup_time=setup_stats.time, sparsity_time=sparsity_stats.time, relax_time=relax_stats.time,
        model_time=model_stats.time, solve_time=solve_stats.time, jump_solve_time=solve_time,
        pre_basis_size, psd_blocks=length(block_sizes), largest_block=isempty(block_sizes) ? 0 : maximum(block_sizes),
        reduced_psd_vars, jump_variables=counts.nvars, jump_constraints=counts.ncons)
end

function _print_rollup(rows)
    isempty(rows) && return nothing
    println("\n## Rollup")
    println("\n| N | total wall (s) | setup | sparsity | symmetry/relax | model build | Mosek wall | Mosek reported | objective/spin | largest block | PSD vars | JuMP vars | status |")
    println("|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|")
    for r in rows
        @printf(
            "| %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.12g | %d | %d | %d | `%s` |\n",
            r.N,
            r.total_wall,
            r.setup_time,
            r.sparsity_time,
            r.relax_time,
            r.model_time,
            r.solve_time,
            r.jump_solve_time,
            r.objective_per_spin,
            r.largest_block,
            r.reduced_psd_vars,
            r.jump_variables,
            _safe_show(r.status),
        )
    end
end

function main()
    ns = parse.(Int, split(get(ENV, "NCTS_PERF_NS", "8,10,12,16"), ","))
    mode = Symbol(get(ENV, "NCTS_PERF_MODE", "order2_singlet"))
    degree = parse(Int, get(ENV, "NCTS_PERF_DEGREE", "4"))
    offblock_check = Symbol(get(ENV, "NCTS_PERF_OFFBLOCK_CHECK", "off"))
    check_invariance = parse(Bool, get(ENV, "NCTS_PERF_CHECK_INVARIANCE", "true"))
    reasonable_wall_s = parse(Float64, get(ENV, "NCTS_REASONABLE_WALL_S", string(30 * 60)))
    reasonable_mem_gib = parse(Float64, get(ENV, "NCTS_REASONABLE_MEM_GIB", "96"))

    println("# Mosek-backed periodic Heisenberg chain scaling")
    println("\n- generated: `$(Dates.now())`")
    println("- Julia: `$(VERSION)`")
    println("- threads: `$(Threads.nthreads())`")
    println("- CPU: `$(Sys.cpu_info()[1].model)`")
    println("- mode: `$mode`")
    println("- N list: `$(join(ns, ","))`")
    println("- offblock_check: `$offblock_check`")
    println("- check_invariance: `$check_invariance`")
    println("- reasonable wall budget: `$reasonable_wall_s` seconds")
    println("- reasonable memory budget: `$reasonable_mem_gib` GiB")
    println("- initial process memory: `$(_rss_string())`")
    println("\nReasonable means setup + JuMP/MOI build + Mosek solve completes within the wall/memory budget above on this host. Runs past that are not engineering wins; they are just heat.")

    rows = NamedTuple[]
    for N in ns
        try
            row = _run_one(N, mode; degree, offblock_check, check_invariance)
            push!(rows, row)
            if row.total_wall > reasonable_wall_s
                println("\nStopping after N=$N: measured wall time exceeded the reasonable budget.")
                break
            end
        catch err
            println("\n## N = $N failed")
            println("\n- error type: `$(typeof(err))`")
            println("- error: `$(_safe_show(err))`")
            Base.show_backtrace(stdout, catch_backtrace())
            println()
            break
        end
        flush(stdout)
    end
    _print_rollup(rows)
end

main()
