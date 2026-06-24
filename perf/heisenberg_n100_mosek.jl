#!/usr/bin/env julia

# Mosek-backed order-2 (degree-4) XXX Heisenberg chain run with all currently
# implemented Pauli reductions: U(1) charge sectors, spatial dihedral symmetry,
# and SU(2)-singlet moment equalities.

using Dates
using Printf
using TOML
using JuMP
using NCTSSoS
using MosekTools

function _readchomp_or(cmd::Cmd, fallback::AbstractString)
    try
        return readchomp(cmd)
    catch
        return fallback
    end
end

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
    started = Dates.now()
    stats = @timed f()
    finished = Dates.now()
    @printf("%-28s %10.3f s  %12s allocated  %8.3f s gc\n",
        label, stats.time, _fmt_bytes(stats.bytes), stats.gctime)
    flush(stdout)
    return (value=stats.value, time=stats.time, bytes=stats.bytes,
        gctime=stats.gctime, started=started, finished=finished)
end

_block_variable_count(sizes) = sum(n * (n + 1) ÷ 2 for n in sizes; init=0)

function _heisenberg_case(N::Integer)
    registry, (sx, sy, sz) = create_pauli_variables(1:N)
    H = sum(
        ComplexF64(1 / 4) * op[i] * op[mod1(i + 1, N)]
        for op in (sx, sy, sz) for i in 1:N
    )
    pop = polyopt(H, registry)

    translation = pauli_site_permutation([2:N; 1])
    reflection = pauli_site_permutation(reverse(1:N))
    symmetry = SymmetrySpec(
        [translation, reflection];
        pauli_charge=PauliChargeSectorSpec(nqubits=Int(N), max_degree=2),
        pauli_singlet=PauliSingletConstraintSpec(nqubits=Int(N), max_degree=2),
        check_invariance=parse(Bool, get(ENV, "NCTS_CHECK_INVARIANCE", "true")),
        offblock_check=Symbol(get(ENV, "NCTS_OFFBLOCK_CHECK", "off")),
    )
    cfg = SolverConfig(
        optimizer=nothing,
        order=2,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry,
    )
    return (; pop, symmetry, cfg)
end

function _mosek_optimizer()
    threads = parse(Int, get(ENV, "NCTS_MOSEK_THREADS", string(max(1, div(Sys.CPU_THREADS, 2)))))
    return optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_IPAR_LOG" => parse(Int, get(ENV, "NCTS_MOSEK_LOG", "0")),
        "MSK_IPAR_NUM_THREADS" => threads,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => parse(Float64, get(ENV, "NCTS_MOSEK_PFEAS", "1e-8")),
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => parse(Float64, get(ENV, "NCTS_MOSEK_DFEAS", "1e-8")),
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => parse(Float64, get(ENV, "NCTS_MOSEK_REL_GAP", "1e-8")),
    )
end

function _write_report(path::AbstractString, data::Dict)
    open(path, "w") do io
        TOML.print(io, data; sorted=true)
    end
end

function _write_markdown(path::AbstractString, data::Dict)
    open(path, "w") do io
        println(io, "# Heisenberg N=$(data["N"]) sparse-degree-$(data["polynomial_degree"]) Mosek symmetry run")
        println(io)
        println(io, "- generated: `$(data["generated"])`")
        println(io, "- git branch: `$(data["git_branch"])`")
        println(io, "- git commit: `$(data["git_commit"])`")
        println(io, "- git diff stat: `$(data["git_diff_stat"])`")
        println(io, "- relaxation order: `$(data["relaxation_order"])`")
        println(io, "- sparse polynomial degree: `$(data["polynomial_degree"])`")
        println(io, "- Julia: `$(data["julia_version"])`")
        println(io, "- threads: `$(data["threads"])`")
        println(io, "- CPU: `$(data["cpu"])`")
        println(io, "- run directory: `$(data["run_dir"])`")
        println(io, "- solver: `Mosek`")
        println(io, "- SDP form: `$(data["sdp_form"])`")
        println(io)
        println(io, "## Result")
        println(io)
        println(io, "- objective: `$(data["objective"])`")
        println(io, "- objective per site: `$(data["objective_per_site"])`")
        println(io, "- solver status: `$(data["solver_status"])`")
        println(io, "- primal status: `$(data["primal_status"])`")
        println(io, "- dual status: `$(data["dual_status"])`")
        println(io)
        println(io, "## Symmetry")
        println(io)
        println(io, "- ingredients: `$(join(data["symmetry_ingredients"], ", "))`")
        println(io, "- dense half-basis dimension: `$(data["dense_half_basis_dim"])`")
        println(io, "- dense full moment-basis estimate: `$(data["dense_full_moment_basis_dim"])`")
        println(io, "- dense PSD scalar variables: `$(data["dense_psd_scalar_variables"])`")
        println(io, "- spatial group order: `$(data["group_order"])`")
        println(io, "- PSD blocks: `$(data["psd_blocks"])`")
        println(io, "- largest PSD block: `$(data["largest_psd_block"])`")
        println(io, "- reduced PSD scalar variables: `$(data["reduced_psd_scalar_variables"])`")
        println(io, "- block provenance: `$(join(data["block_provenance"], ", "))`")
        println(io, "- invariant moment count: `$(data["invariant_moment_count"])`")
        println(io, "- symbolic constraints: `$(data["symbolic_constraints"])`")
        println(io, "- zero scalar constraints: `$(data["zero_scalar_constraints"])`")
        println(io, "- JuMP variables: `$(data["jump_variables"])`")
        println(io, "- JuMP constraints: `$(data["jump_constraints"])`")
        println(io)
        println(io, "## Timings")
        println(io)
        println(io, "| phase | seconds | allocated bytes | GC seconds |")
        println(io, "|:--|--:|--:|--:|")
        for phase in data["timings"]
            println(io, "| $(phase["phase"]) | $(phase["seconds"]) | $(phase["bytes"]) | $(phase["gc_seconds"]) |")
        end
    end
end

function main()
    N = parse(Int, get(ENV, "NCTS_HEISENBERG_N", "100"))
    outdir = get(ENV, "NCTS_RESULTS_DIR", "results")
    mkpath(outdir)

    println("# Heisenberg N=$N sparse-degree-4/order-2 Pauli symmetry + Mosek run")
    println("started: $(Dates.now())")
    println("threads: $(Threads.nthreads())")
    println("cpu: $(Sys.cpu_info()[1].model)")
    flush(stdout)

    case_t = _timed("build case", () -> _heisenberg_case(N))
    sparsity_t = _timed("compute_sparsity", () -> compute_sparsity(case_t.value.pop, case_t.value.cfg))
    basis = only(sparsity_t.value.corr_sparsity.clq_mom_mtx_bases)
    dense_dim = length(basis)

    relax_t = _timed("moment_relax_symmetric", () -> NCTSSoS.moment_relax_symmetric(
        case_t.value.pop,
        sparsity_t.value.corr_sparsity,
        sparsity_t.value.cliques_term_sparsities,
        case_t.value.symmetry,
    ))
    mp, report = relax_t.value

    solve_t = _timed("solve_sdp Mosek", () -> NCTSSoS.solve_sdp(mp, _mosek_optimizer(); dualize=true))
    solve = solve_t.value

    sizes = report.psd_block_sizes
    data = Dict{String,Any}(
        "N" => N,
        "relaxation_order" => 2,
        "polynomial_degree" => 4,
        "sdp_form" => "dualized JuMP conic SDP",
        "symmetry_ingredients" => [
            "Pauli U(1) charge sectors through degree 2",
            "dihedral spatial translation/reflection group",
            "order-2 SU(2)-singlet moment equalities",
            "lazy charge-orbit representative blocks for large sectors",
        ],
        "generated" => string(Dates.now()),
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
        "cpu" => Sys.cpu_info()[1].model,
        "run_dir" => pwd(),
        "git_branch" => _readchomp_or(`git branch --show-current`, "unknown"),
        "git_commit" => _readchomp_or(`git rev-parse --short HEAD`, "unknown"),
        "git_diff_stat" => _readchomp_or(`git diff --shortstat`, ""),
        "objective" => solve.objective,
        "objective_per_site" => solve.objective / N,
        "solver_status" => string(solve.status),
        "primal_status" => string(primal_status(solve.model)),
        "dual_status" => string(dual_status(solve.model)),
        "dense_half_basis_dim" => dense_dim,
        "dense_full_moment_basis_dim" => report.basis_full_size,
        "dense_psd_scalar_variables" => dense_dim * (dense_dim + 1) ÷ 2,
        "group_order" => report.group_order,
        "psd_blocks" => length(sizes),
        "largest_psd_block" => maximum(sizes),
        "reduced_psd_scalar_variables" => _block_variable_count(sizes),
        "block_provenance" => string.(report.block_provenance),
        "invariant_moment_count" => report.invariant_moment_count,
        "symbolic_constraints" => length(mp.constraints),
        "zero_scalar_constraints" => length(mp.linear.zero_constraints),
        "jump_variables" => JuMP.num_variables(solve.model),
        "jump_constraints" => JuMP.num_constraints(solve.model; count_variable_in_set_constraints=true),
        "timings" => [
            Dict("phase" => "build case", "seconds" => case_t.time, "bytes" => case_t.bytes, "gc_seconds" => case_t.gctime),
            Dict("phase" => "compute_sparsity", "seconds" => sparsity_t.time, "bytes" => sparsity_t.bytes, "gc_seconds" => sparsity_t.gctime),
            Dict("phase" => "moment_relax_symmetric", "seconds" => relax_t.time, "bytes" => relax_t.bytes, "gc_seconds" => relax_t.gctime),
            Dict("phase" => "solve_sdp Mosek", "seconds" => solve_t.time, "bytes" => solve_t.bytes, "gc_seconds" => solve_t.gctime),
        ],
    )

    stem = joinpath(outdir, "heisenberg_n$(N)_order2_pauli_symmetry_mosek")
    _write_report(stem * ".toml", data)
    _write_markdown(stem * ".md", data)

    println("finished: $(Dates.now())")
    println("objective: $(solve.objective)")
    println("objective per site: $(solve.objective / N)")
    println("status: $(solve.status)")
    println("report: $(stem).md")
end

main()
