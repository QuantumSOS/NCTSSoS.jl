#!/usr/bin/env julia
# Shared H2/H4 periodic V2RDM solver-evidence driver. Manual HAI artifact, not a test.
using NCTSSoS
using JuMP
const MOI = JuMP.MOI
using LinearAlgebra
using Printf
using Dates
using JSON3
import BPSDP
import COSMO
import IterativeSolvers
import LinearMaps
import SCS
Base.@kwdef mutable struct Options
    system::String = ""
    nk::Int = 0
    norb::Int = 0
    nelec_per_cell::Int = 0
    integrals::String = ""
    meta::String = ""
    output_dir::String = ""
    solvers::Vector{Symbol} = [:cosmo, :scs, :bpsdp]
    time_limit::Float64 = 3600.0
    repeats::Int = 1
    force::Bool = false
    include_one_d::Bool = false
    cosmo_max_iter::Int = 200_000
    cosmo_eps_abs::Float64 = 1e-5
    cosmo_eps_rel::Float64 = 1e-5
    cosmo_tol_constant::Float64 = 1.0
    cosmo_tol_exponent::Float64 = 1.5
    scs_max_iters::Int = 200_000
    scs_eps_abs::Float64 = 1e-5
    scs_eps_rel::Float64 = 1e-5
    bpsdp_max_iter::Int = 50_000
    bpsdp_cg_max_iter::Int = 10_000
    bpsdp_obj_tol::Float64 = 1e-5
    bpsdp_err_tol::Float64 = 1e-5
    # Keep BPSDP stopping absolute for fair raw ∞-norm tuning; server HEAD defaults
    # relative tolerances to 1e-5, which otherwise stops above the requested raw tol.
    bpsdp_obj_tol_rel::Float64 = 0.0
    bpsdp_err_tol_rel::Float64 = 0.0
    # cg_convergence: 1e-8 matches BPSDP.jl's recommended config with dynamic_cg off.
    # Old value 1e-10 was tighter than needed and slowed inner CG without helping outer convergence.
    bpsdp_cg_tol::Float64 = 1e-8
    # mu_update_frequency: 25 matches BPSDP.jl test/reference_cpp.jl. Old value 500 left the
    # penalty parameter stale between updates and stalled outer convergence.
    bpsdp_mu_update_frequency::Int = 25
    bpsdp_penalty_parameter::Float64 = 0.1
    # dynamic_cg_convergence: BPSDP.jl reference disables it; adaptive CG tol at early outer
    # iters can over-loosen and cause inaccurate steps.
    bpsdp_dynamic_cg::Bool = false
    bpsdp_dependent_rows::Symbol = :keep
    bpsdp_dependent_row_tol::Float64 = 1e-12
    solver_threads::Int = parse(Int, get(ENV, "SOLVER_THREADS", "8"))
end
function parse_solvers(s::AbstractString)
    out = Symbol[]
    for raw in split(s, ',')
        t = lowercase(strip(raw)); isempty(t) && continue
        if t == "all"
            append!(out, [:cosmo, :scs, :bpsdp])
        elseif t in ("cosmo", "cosmo-cg", "cosmo_cg")
            push!(out, :cosmo)
        elseif t in ("scs", "scs-indirect", "scs_indirect")
            push!(out, :scs)
        elseif t == "bpsdp"
            push!(out, :bpsdp)
        else
            throw(ArgumentError("unknown solver $(repr(t)); expected cosmo, scs, bpsdp, or all"))
        end
    end
    isempty(out) && throw(ArgumentError("--solvers names no solvers"))
    return unique(out)
end
function parse_options(argv)
    o = Options()
    for arg in argv
        val() = split(arg, "=", limit = 2)[2]
        if startswith(arg, "--system="); o.system = lowercase(val())
        elseif startswith(arg, "--nk="); o.nk = parse(Int, val())
        elseif startswith(arg, "--norb="); o.norb = parse(Int, val())
        elseif startswith(arg, "--nelec-per-cell="); o.nelec_per_cell = parse(Int, val())
        elseif startswith(arg, "--integrals="); o.integrals = val()
        elseif startswith(arg, "--meta="); o.meta = val()
        elseif startswith(arg, "--output-dir="); o.output_dir = val()
        elseif startswith(arg, "--solvers="); o.solvers = parse_solvers(val())
        elseif startswith(arg, "--time-limit="); o.time_limit = parse(Float64, val())
        elseif startswith(arg, "--repeats="); o.repeats = parse(Int, val())
        elseif startswith(arg, "--include-one-d="); o.include_one_d = (lowercase(val()) in ("true","1","yes"))
        elseif startswith(arg, "--cosmo-max-iter="); o.cosmo_max_iter = parse(Int, val())
        elseif startswith(arg, "--cosmo-eps-abs="); o.cosmo_eps_abs = parse(Float64, val())
        elseif startswith(arg, "--cosmo-eps-rel="); o.cosmo_eps_rel = parse(Float64, val())
        elseif startswith(arg, "--cosmo-tol-constant="); o.cosmo_tol_constant = parse(Float64, val())
        elseif startswith(arg, "--cosmo-tol-exponent="); o.cosmo_tol_exponent = parse(Float64, val())
        elseif startswith(arg, "--scs-max-iters="); o.scs_max_iters = parse(Int, val())
        elseif startswith(arg, "--scs-eps-abs="); o.scs_eps_abs = parse(Float64, val())
        elseif startswith(arg, "--scs-eps-rel="); o.scs_eps_rel = parse(Float64, val())
        elseif startswith(arg, "--bpsdp-max-iter="); o.bpsdp_max_iter = parse(Int, val())
        elseif startswith(arg, "--bpsdp-cg-max-iter="); o.bpsdp_cg_max_iter = parse(Int, val())
        elseif startswith(arg, "--bpsdp-obj-tol="); o.bpsdp_obj_tol = parse(Float64, val())
        elseif startswith(arg, "--bpsdp-err-tol="); o.bpsdp_err_tol = parse(Float64, val())
        elseif startswith(arg, "--bpsdp-obj-tol-rel="); o.bpsdp_obj_tol_rel = parse(Float64, val())
        elseif startswith(arg, "--bpsdp-err-tol-rel="); o.bpsdp_err_tol_rel = parse(Float64, val())
        elseif startswith(arg, "--bpsdp-cg-tol="); o.bpsdp_cg_tol = parse(Float64, val())
        elseif startswith(arg, "--bpsdp-mu-update-frequency="); o.bpsdp_mu_update_frequency = parse(Int, val())
        elseif startswith(arg, "--bpsdp-penalty-parameter="); o.bpsdp_penalty_parameter = parse(Float64, val())
        elseif startswith(arg, "--bpsdp-dynamic-cg="); o.bpsdp_dynamic_cg = (lowercase(val()) in ("true","1","yes"))
        elseif startswith(arg, "--bpsdp-dependent-rows="); o.bpsdp_dependent_rows = Symbol(lowercase(val()))
        elseif startswith(arg, "--bpsdp-dependent-row-tol="); o.bpsdp_dependent_row_tol = parse(Float64, val())
        elseif startswith(arg, "--solver-threads="); o.solver_threads = parse(Int, val())
        elseif arg == "--force"; o.force = true
        elseif arg in ("-h", "--help"); print_help(); exit(0)
        else; throw(ArgumentError("unknown argument: $arg"))
        end
    end
    o.system in ("h2", "h4") || throw(ArgumentError("--system must be h2 or h4"))
    o.nk > 0 && o.norb > 0 && o.nelec_per_cell > 0 || throw(ArgumentError("--nk, --norb, --nelec-per-cell must be positive"))
    for (name, path) in (("--integrals", o.integrals), ("--meta", o.meta), ("--output-dir", o.output_dir))
        isempty(path) && throw(ArgumentError("$name is required"))
    end
    o.repeats >= 1 || throw(ArgumentError("--repeats must be positive"))
    return o
end
function print_help()
    println("Usage: julia --project=. driver.jl --system=h2|h4 --nk=N --norb=N --nelec-per-cell=N --integrals=PATH --meta=PATH --output-dir=PATH [--include-one-d=true|false] [--solvers=cosmo,scs,bpsdp]")
end
@inline complex_entry(f, r, i) = complex(parse(Float64, f[r]), parse(Float64, f[i]))
function load_integrals_txt(path::AbstractString; norb::Int)
    h1e = Dict{Int,Matrix{ComplexF64}}()
    eri = Dict{NTuple{4,Int},Array{ComplexF64,4}}()
    for (line_no, raw) in enumerate(eachline(path))
        line = strip(raw); (isempty(line) || startswith(line, '#')) && continue
        f = split(line); tag = first(f)
        if tag == "h1e"
            length(f) == 6 || error("malformed h1e line $line_no in $path")
            k = parse(Int, f[2]); p = parse(Int, f[3]) + 1; q = parse(Int, f[4]) + 1
            block = get!(h1e, k) do; zeros(ComplexF64, norb, norb); end
            block[p, q] = complex_entry(f, 5, 6)
        elseif tag == "eri"
            length(f) == 11 || error("malformed eri line $line_no in $path")
            key = ntuple(i -> parse(Int, f[i + 1]), 4)
            p = parse(Int, f[6]) + 1; r = parse(Int, f[7]) + 1; q = parse(Int, f[8]) + 1; s = parse(Int, f[9]) + 1
            block = get!(eri, key) do; zeros(ComplexF64, norb, norb, norb, norb); end
            block[p, r, q, s] = complex_entry(f, 10, 11)
        else
            error("unexpected integral tag $(repr(tag)) at line $line_no in $path")
        end
    end
    return h1e, eri
end
function lookup_mono(L, key)
    v = get(L.key_to_monomial, key, nothing)
    v !== nothing && return v
    for (k, mono) in L.key_to_monomial
        NCTSSoS.key_isequal(k, key) && return mono
    end
    error("missing monomial for key $(repr(key))")
end
function poly_from_form(form, L)
    M = typeof(lookup_mono(L, L.identity))
    terms = Tuple{ComplexF64,M}[]
    for (key, coef) in form
        push!(terms, (ComplexF64(coef), lookup_mono(L, key)))
    end
    return Polynomial(terms)
end
function moment_problem_from_linear(L)
    M = typeof(lookup_mono(L, L.identity))
    T = eltype(L.identity)
    P = Polynomial{FermionicAlgebra,T,ComplexF64}
    constraints = Tuple{Symbol,Matrix{P}}[]  # PSD placeholders; cached L carries entries and zero rows.
    z = zero(P)
    for block in L.psd_blocks_lin
        push!(constraints, (block.meta.cone, fill(z, block.size, block.size)))
    end
    total_basis = M[lookup_mono(L, key) for key in L.moments]
    return NCTSSoS._moment_problem_with_linear(
        FermionicAlgebra, T, M, P, poly_from_form(L.objective_lin, L),
        constraints, total_basis, length(L.moments), L)
end
function build_problem(opts::Options)
    h1e, eri = load_integrals_txt(opts.integrals; norb = opts.norb)
    build_seconds = @elapsed linear = build_pqg_moment_data(h1e, eri;
        nk = opts.nk, norb = opts.norb, nelec_per_cell = opts.nelec_per_cell,
        blocking = :momentum, spin_resolved_trace = true, singlet_s2 = true,
        include_one_d = opts.include_one_d)
    NCTSSoS.assert_moment_linear_data_invariants(linear)
    mp = moment_problem_from_linear(linear)
    stats = Dict{String,Any}(
        "include_one_d" => opts.include_one_d,
        "hpsd_block_sizes" => [b.size for b in linear.psd_blocks_lin if b.meta.cone == :HPSD],
        "total_canonical_moments" => length(linear.moments),
        "zero_scalar_constraints" => length(linear.zero_constraints),
        "n_orphans" => length(NCTSSoS.orphan_keys(mp)),
    )
    @info "moment-problem build" n_orphans=stats["n_orphans"] hpsd_block_sizes=stats["hpsd_block_sizes"] total_canonical_moments=stats["total_canonical_moments"]
    return mp, stats, build_seconds
end
has_field(x, s::Symbol) = s in fieldnames(typeof(x))
field_or_nothing(x, s::Symbol) = x === nothing || !has_field(x, s) ? nothing : getfield(x, s)
first_field(x, names::Symbol...) = (for n in names; v = field_or_nothing(x, n); v === nothing || return v; end; nothing)
finite(x) = x === nothing ? nothing : (x isa Real ? (isfinite(Float64(x)) ? Float64(x) : string(x)) : x)
function set_solver!(model, solver::Symbol, opts::Options, cg_iters::Vector{Int}, bpsdp_progress::Vector{Any} = Any[])
    if solver == :cosmo
        isdefined(COSMO, :CGIndirectKKTSolver) || error("COSMO.CGIndirectKKTSolver unavailable; install IterativeSolvers and LinearMaps")
        kkt = COSMO.with_options(COSMO.CGIndirectKKTSolver; tol_constant = opts.cosmo_tol_constant, tol_exponent = opts.cosmo_tol_exponent)
        JuMP.set_optimizer(model, JuMP.optimizer_with_attributes(COSMO.Optimizer,
            "kkt_solver" => kkt, "max_iter" => opts.cosmo_max_iter,
            "eps_abs" => opts.cosmo_eps_abs, "eps_rel" => opts.cosmo_eps_rel,
            "time_limit" => opts.time_limit, "verbose" => false,
            "verbose_timing" => false, "decompose" => false, "complete_dual" => false))
        JuMP.set_silent(model)
    elseif solver == :scs
        SCS.is_available(SCS.IndirectSolver) || error("SCS.IndirectSolver missing in this build")
        JuMP.set_optimizer(model, JuMP.optimizer_with_attributes(SCS.Optimizer,
            "linear_solver" => SCS.IndirectSolver, "max_iters" => opts.scs_max_iters,
            "eps_abs" => opts.scs_eps_abs, "eps_rel" => opts.scs_eps_rel,
            "time_limit_secs" => opts.time_limit, "verbose" => 0))
    elseif solver == :bpsdp
        monitor = (_lvl, outer, inner, objective_primal, objective_dual, mu, primal_error, dual_error) -> begin
            push!(cg_iters, Int(inner))
            if outer <= 10 || outer % 100 == 0
                push!(bpsdp_progress, Dict(
                    "outer" => Int(outer),
                    "inner" => Int(inner),
                    "mu" => finite(mu),
                    "primal_error" => finite(primal_error),
                    "dual_error" => finite(dual_error),
                    "objective_primal" => finite(objective_primal),
                    "objective_dual" => finite(objective_dual),
                ))
            end
            return nothing
        end
        JuMP.set_optimizer(model, () -> BPSDP.Optimizer(
            max_iter = opts.bpsdp_max_iter, cg_max_iter = opts.bpsdp_cg_max_iter,
            mu_update_frequency = opts.bpsdp_mu_update_frequency,
            penalty_parameter = opts.bpsdp_penalty_parameter,
            cg_convergence = opts.bpsdp_cg_tol,
            dynamic_cg_convergence = opts.bpsdp_dynamic_cg,
            sdp_objective_convergence = opts.bpsdp_obj_tol,
            sdp_objective_convergence_rel = opts.bpsdp_obj_tol_rel,
            sdp_error_convergence = opts.bpsdp_err_tol,
            sdp_error_convergence_rel = opts.bpsdp_err_tol_rel, guess_type = :zero,
            print_level = 1, progress_monitor = monitor,
            dependent_rows = opts.bpsdp_dependent_rows,
            dependent_row_tol = opts.bpsdp_dependent_row_tol,
            projection_threads = opts.solver_threads))
    else
        error("unknown solver $solver")
    end
end
raw_optimizer(model) = try JuMP.unsafe_backend(model) catch; nothing end
safe_status(model) = try string(JuMP.termination_status(model)) catch e; "unavailable: $(sprint(showerror, e))" end
safe_raw_status(model) = try string(MOI.get(model, MOI.RawStatusString())) catch e; "unavailable: $(sprint(showerror, e))" end
safe_obj(model) = try JuMP.objective_value(model) catch e; "unavailable: $(sprint(showerror, e))" end
safe_solve_time(model) = try finite(MOI.get(model, MOI.SolveTimeSec())) catch; nothing end
function bpsdp_state_summary(raw)
    st = field_or_nothing(raw, :state); st === nothing && return nothing
    return Dict("outer_iterations" => first_field(st, :outer_iterations, :iterations),
        "primal_error" => finite(field_or_nothing(st, :primal_error)),
        "dual_error" => finite(field_or_nothing(st, :dual_error)),
        "objective_primal" => finite(first_field(st, :objective_primal, :primal_objective)),
        "objective_dual" => finite(first_field(st, :objective_dual, :dual_objective)),
        "termination_reason" => string(first_field(st, :termination_reason, :status)))
end
function scs_state_summary(raw)
    sol = field_or_nothing(raw, :sol); sol === nothing && return nothing
    po = field_or_nothing(sol, :objective_value); do_ = field_or_nothing(sol, :dual_objective_value)
    info = field_or_nothing(sol, :info)  # SCS sol.info exposes residuals on newer SCS.jl builds
    return Dict("iterations" => field_or_nothing(sol, :iterations),
        "solve_time_sec_reported" => finite(field_or_nothing(sol, :solve_time_sec)),
        "primal_objective" => finite(po), "dual_objective" => finite(do_),
        "objective_gap" => (po isa Real && do_ isa Real ? finite(abs(po - do_)) : nothing),
        "raw_status" => string(field_or_nothing(sol, :raw_status)),
        "ret_val" => field_or_nothing(sol, :ret_val),
        "res_pri" => finite(first_field(sol, :res_pri, :res_primal)),
        "res_dual" => finite(first_field(sol, :res_dual)),
        "res_infeas" => finite(field_or_nothing(sol, :res_infeas)),
        "res_unbdd_a" => finite(field_or_nothing(sol, :res_unbdd_a)),
        "res_unbdd_p" => finite(field_or_nothing(sol, :res_unbdd_p)),
        "info_res_pri" => finite(field_or_nothing(info, :res_pri)),
        "info_res_dual" => finite(field_or_nothing(info, :res_dual)),
        "info_gap" => finite(field_or_nothing(info, :gap)))
end
function cosmo_state_summary(raw)
    res = field_or_nothing(raw, :results); res === nothing && return nothing
    info = field_or_nothing(res, :info); times = field_or_nothing(res, :times)
    inner = field_or_nothing(raw, :inner); kkt = field_or_nothing(inner, :kkt_solver)
    indirect = field_or_nothing(kkt, :indirect_kktsolver); mult = field_or_nothing(indirect, :multiplications)
    return Dict("status" => string(field_or_nothing(res, :status)),
        "iterations" => field_or_nothing(res, :iter),
        "r_prim" => finite(field_or_nothing(info, :r_prim)),
        "r_dual" => finite(field_or_nothing(info, :r_dual)),
        "solver_time" => finite(field_or_nothing(times, :solver_time)),
        "kkt_matvecs_total" => mult === nothing ? nothing : sum(mult),
        "kkt_matvecs_max" => mult === nothing || isempty(mult) ? nothing : maximum(mult))
end
function truncate_log_file!(path; max_bytes = 262_144)
    isfile(path) || return Dict("path" => path, "exists" => false, "truncated" => false)
    bytes = filesize(path); bytes <= max_bytes && return Dict("path" => path, "exists" => true, "bytes" => bytes, "truncated" => false)
    head_n = max_bytes ÷ 2; tail_n = max_bytes - head_n
    head = tail = UInt8[]
    open(path, "r") do io
        head = read(io, head_n); seek(io, max(Int(bytes) - tail_n, 0)); tail = read(io)
    end
    open(path, "w") do io
        write(io, head); write(io, "\n... [log truncated from $bytes to $max_bytes bytes] ...\n"); write(io, tail)
    end
    return Dict("path" => path, "exists" => true, "bytes" => filesize(path), "original_bytes" => bytes, "truncated" => true)
end
function optimize_with_redirected_stdio!(model, out_path, err_path)
    mkpath(dirname(out_path)); opt_err = nothing
    open(out_path, "w") do out; open(err_path, "w") do err
        try
            redirect_stdout(out) do; redirect_stderr(err) do; JuMP.optimize!(model); end; end
        catch e
            opt_err = sprint(showerror, e)
        end
    end; end
    return opt_err, Dict("stdout" => truncate_log_file!(out_path), "stderr" => truncate_log_file!(err_path))
end
function json_ready(x)
    x === nothing && return nothing; x isa Missing && return nothing
    (x isa AbstractString || x isa Bool || x isa Integer) && return x
    x isa AbstractFloat && return isfinite(x) ? x : string(x)
    x isa Real && return finite(x); x isa Symbol && return string(x); x isa DateTime && return string(x)
    x isa Pair && return string(x.first) => json_ready(x.second)
    x isa NamedTuple && return Dict(String(k) => json_ready(v) for (k, v) in pairs(x))
    x isa AbstractDict && return Dict(string(k) => json_ready(v) for (k, v) in x)
    x isa Tuple && return Any[json_ready(v) for v in x]
    x isa AbstractVector && return Any[json_ready(v) for v in x]
    return string(x)
end
function write_json(path, obj)
    mkpath(dirname(path)); open(path, "w") do io; JSON3.write(io, json_ready(obj)); write(io, '\n'); end
end
read_meta(path) = JSON3.read(read(path, String), Dict{String,Any})
rss_string() = @sprintf("%.3f GiB", Sys.maxrss() / 1024.0^3)
git_rev() = try readchomp(pipeline(`git rev-parse --short HEAD`, stderr = devnull)) catch; nothing end
common_runtime_json(opts::Union{Options,Nothing} = nothing) = Dict(
    "hostname" => gethostname(), "julia_version" => string(VERSION),
    "julia_threads" => Threads.nthreads(),
    "julia_num_threads_env" => get(ENV, "JULIA_NUM_THREADS", nothing),
    "blas_threads" => BLAS.get_num_threads(),
    "openblas_num_threads_env" => get(ENV, "OPENBLAS_NUM_THREADS", nothing),
    "omp_num_threads_env" => get(ENV, "OMP_NUM_THREADS", nothing),
    "solver_threads_env" => get(ENV, "SOLVER_THREADS", nothing),
    "solver_threads_target" => opts === nothing ? nothing : opts.solver_threads,
    "cpu_threads" => Sys.CPU_THREADS, "max_rss" => rss_string(), "git_rev" => git_rev())
function jump_model_size(model)
    sz = Dict{String,Any}(
        "num_variables" => JuMP.num_variables(model),
        "num_constraints_total" => JuMP.num_constraints(model; count_variable_in_set_constraints = true),
        "num_constraints_no_var" => JuMP.num_constraints(model; count_variable_in_set_constraints = false),
    )
    by_type = Dict{String,Int}()
    for (F, S) in JuMP.list_of_constraint_types(model)
        by_type[string("(", F, ", ", S, ")")] = JuMP.num_constraints(model, F, S)
    end
    sz["by_constraint_type"] = by_type
    return sz
end
solver_label(s::Symbol) = string(s)
function add_objectives!(r, meta)
    active = get(r, "active_objective_Ha", nothing); shift = get(meta, "energy_shift", nothing); ehf = get(meta, "ehf", nothing)
    r["energy_shift_Ha"] = shift; r["ehf_Ha"] = ehf
    if active isa Real && shift isa Real
        total = active + shift; r["total_per_cell_Ha"] = total
        r["gap_to_hf_Ha"] = ehf isa Real ? total - ehf : nothing
    else
        r["total_per_cell_Ha"] = nothing; r["gap_to_hf_Ha"] = nothing
    end
end
function solve_one(mp, problem_stats, build_seconds, meta, opts::Options, solver::Symbol, result_path::String, repeat::Int)
    started = Dates.now(); cg = Int[]; bpsdp_progress = Any[]
    r = Dict{String,Any}("system" => opts.system, "nk" => opts.nk, "solver" => solver_label(solver),
        "repeat" => repeat, "norb" => opts.norb, "nelec_per_cell" => opts.nelec_per_cell,
        "integrals" => opts.integrals, "meta" => opts.meta, "started_at" => string(started),
        "runtime" => common_runtime_json(opts), "problem" => problem_stats,
        "build_wall_seconds" => build_seconds,
        "pyscf" => Dict("ehf_Ha" => get(meta, "ehf", nothing), "energy_shift_Ha" => get(meta, "energy_shift", nothing), "active_hf_Ha" => get(meta, "active_hf", nothing), "basis" => get(meta, "basis", nothing), "density_fit" => get(meta, "density_fit", nothing)))
    model = nothing; lowering_error = nothing
    r["lowering_wall_seconds"] = @elapsed try
        model, _ = build_jump_model(mp; formulation = :psd_blocks, representation = :complex)
    catch e
        lowering_error = sprint(showerror, e)
    end
    r["lowering_error"] = lowering_error
    if model !== nothing
        try
            r["sdp_size"] = jump_model_size(model)
        catch e
            r["sdp_size_error"] = sprint(showerror, e)
        end
    end
    if lowering_error === nothing
        seterr = nothing
        try set_solver!(model, solver, opts, cg, bpsdp_progress) catch e; seterr = sprint(showerror, e) end
        r["set_solver_error"] = seterr
        if seterr === nothing
            solver_dir = dirname(result_path)
            opt_err = nothing; stdio = nothing
            r["solve_wall_seconds"] = @elapsed opt_err, stdio = optimize_with_redirected_stdio!(model, joinpath(solver_dir, "stdout.log"), joinpath(solver_dir, "stderr.log"))
            raw = raw_optimizer(model); r["termination_status"] = safe_status(model); r["raw_status"] = safe_raw_status(model)
            r["solve_time_sec_reported"] = safe_solve_time(model); r["active_objective_Ha"] = safe_obj(model)
            r["optimize_error"] = opt_err; r["optimize_stdio"] = stdio
            if solver == :bpsdp
                st = bpsdp_state_summary(raw); r["bpsdp_state"] = st; r["inner_iterations"] = sum(cg); r["iterations"] = st === nothing ? nothing : get(st, "outer_iterations", nothing); r["r_prim"] = st === nothing ? nothing : get(st, "primal_error", nothing); r["r_dual"] = st === nothing ? nothing : get(st, "dual_error", nothing); r["cg_iterations_per_outer"] = cg; r["bpsdp_progress_samples"] = bpsdp_progress
                # BPSDP MOI bridge returns no solution at index 1 when termination is ITERATION_LIMIT,
                # so JuMP.objective_value throws. Recover the primal objective from BPSDP's own state.
                if !(r["active_objective_Ha"] isa Real) && st !== nothing
                    op = get(st, "objective_primal", nothing)
                    if op isa Real
                        r["active_objective_recovered_from"] = "bpsdp_state.objective_primal"
                        r["active_objective_Ha"] = op
                    end
                end
            elseif solver == :scs
                st = scs_state_summary(raw); r["scs_state"] = st; r["iterations"] = st === nothing ? nothing : get(st, "iterations", nothing); r["inner_iterations"] = r["iterations"]
            else
                st = cosmo_state_summary(raw); r["cosmo_state"] = st; r["iterations"] = st === nothing ? nothing : get(st, "iterations", nothing); r["inner_iterations"] = st === nothing ? nothing : get(st, "kkt_matvecs_total", nothing); r["r_prim"] = st === nothing ? nothing : get(st, "r_prim", nothing); r["r_dual"] = st === nothing ? nothing : get(st, "r_dual", nothing)
            end
            add_objectives!(r, meta)
        end
    end
    add_objectives!(r, meta)
    r["finished_at"] = string(Dates.now()); r["runtime_after"] = common_runtime_json(opts)
    write_json(result_path, r)
    return r
end
function validate_output_dir!(path)
    abs = abspath(path)
    abs in Set([abspath(pwd()), homedir(), dirname(homedir()), "/"]) && throw(ArgumentError("refusing dangerous --output-dir=$path"))
    mkpath(path)
end
function read_result(path)
    JSON3.read(read(path, String), Dict{String,Any})
end
function fmt(x)
    x === nothing && return ""
    x isa Real && return @sprintf("%.8g", x)
    return replace(string(x), "\n" => " ", "|" => "¦")
end
function write_markdown_summary(system_root::String, system::String)
    rows = Any[]
    for nkdir in sort(filter(p -> isdir(joinpath(system_root, p)) && startswith(p, "nk"), readdir(system_root)))
        nk = parse(Int, replace(nkdir, "nk" => ""))
        for solver in ("cosmo", "scs", "bpsdp")
            path = joinpath(system_root, nkdir, solver, "result.json")
            isfile(path) && push!(rows, (nk, solver, read_result(path)))
        end
    end
    path = joinpath(system_root, "summary.md")
    open(path, "w") do io
        println(io, "# ", uppercase(system), " periodic V2RDM solver evidence\n")
        println(io, "Lowering: `build_jump_model(mp; formulation=:psd_blocks, representation=:complex)`. Default orphan handling; no fallback.\n")
        println(io, "| Nk | solver | status | iter | inner iter | solve s | active Ha | shift Ha | total Ha | HF Ha | gap to HF Ha | r_prim | r_dual |")
        println(io, "|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for (nk, solver, r) in rows
            println(io, "| ", nk, " | `", solver, "` | `", fmt(get(r, "termination_status", get(r, "lowering_error", ""))), "` | ",
                fmt(get(r, "iterations", nothing)), " | ", fmt(get(r, "inner_iterations", nothing)), " | ", fmt(get(r, "solve_wall_seconds", nothing)), " | ",
                fmt(get(r, "active_objective_Ha", nothing)), " | ", fmt(get(r, "energy_shift_Ha", nothing)), " | ", fmt(get(r, "total_per_cell_Ha", nothing)), " | ",
                fmt(get(r, "ehf_Ha", nothing)), " | ", fmt(get(r, "gap_to_hf_Ha", nothing)), " | ", fmt(get(r, "r_prim", nothing)), " | ", fmt(get(r, "r_dual", nothing)), " |")
        end
    end
end
function pin_threads!(opts::Options)
    opts.solver_threads >= 1 || throw(ArgumentError("--solver-threads must be >= 1"))
    BLAS.set_num_threads(opts.solver_threads)
    actual = BLAS.get_num_threads()
    actual == opts.solver_threads || @warn "BLAS.set_num_threads did not stick" requested=opts.solver_threads actual=actual
    @info "pinned solver threads" target=opts.solver_threads blas=actual julia=Threads.nthreads() openblas_env=get(ENV, "OPENBLAS_NUM_THREADS", nothing) omp_env=get(ENV, "OMP_NUM_THREADS", nothing)
    return nothing
end
function main(argv = ARGS)
    opts = parse_options(argv); validate_output_dir!(opts.output_dir)
    pin_threads!(opts)
    meta = read_meta(opts.meta)
    mp, stats, build_seconds = build_problem(opts)
    for repeat in 1:opts.repeats, solver in opts.solvers
        leaf = opts.repeats == 1 ? solver_label(solver) : joinpath(solver_label(solver), "repeat_$repeat")
        result_path = joinpath(opts.output_dir, leaf, "result.json")
        if isfile(result_path) && !opts.force
            @info "skip existing result" result_path
        else
            @printf("== %s Nk=%d %s repeat=%d ==\n", uppercase(opts.system), opts.nk, solver_label(solver), repeat); flush(stdout)
            r = solve_one(mp, stats, build_seconds, meta, opts, solver, result_path, repeat)
            @printf("%-6s status=%s active=%s total=%s gap=%s\n", solver_label(solver), fmt(get(r, "termination_status", get(r, "lowering_error", ""))), fmt(get(r, "active_objective_Ha", nothing)), fmt(get(r, "total_per_cell_Ha", nothing)), fmt(get(r, "gap_to_hf_Ha", nothing)))
        end
        write_markdown_summary(dirname(opts.output_dir), opts.system)
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
