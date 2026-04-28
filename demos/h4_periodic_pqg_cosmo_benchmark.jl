# demos/h4_periodic_pqg_cosmo_benchmark.jl
#
# Purpose
# -------
# Build the H4 periodic Nk=2 PQG / V2RDM-style basis, compare the
# single-fat-block and paper-style (ΔN, K) block routes, add the explicit
# particle-number and trace equalities through `moment_eq_constraints`, and try
# to solve the selected route with COSMO.
#
# Default route
# -------------
# `--route=paper` is the only sensible default. The single-fat-block route is
# still useful as a comparison point, but it is structurally much larger and can
# stall before even finishing `moment_relax`.
#
# Usage
# -----
#   julia --project demos/h4_periodic_pqg_cosmo_benchmark.jl
#   julia --project demos/h4_periodic_pqg_cosmo_benchmark.jl \
#       --route=paper demos/results/h4_periodic_pqg_cosmo_benchmark.md
#
# Flags
# -----
#   --route=paper|fat      Choose the manual 6-block route or the single PQG block.
#   --direct=true|false    Direct COSMO assemble/optimize vs JuMP primal model (default: true).
#   --max-iter=N           COSMO max iterations (default: 10000).
#   --eps=TOL              COSMO eps_abs / eps_rel (default: 1e-4).
#   --decompose=true|false COSMO chordal decomposition switch (default: false).

using SparseArrays
using NCTSSoS

include(joinpath(@__DIR__, "H4PeriodicPQGHelpers.jl"))
include(joinpath(@__DIR__, "results", "h4_periodic_probe_common.jl"))

const PQG = H4PeriodicPQGHelpers

function parse_flags(args)
    route = :paper
    direct = true
    max_iter = 10_000
    eps = 1e-4
    decompose = false
    rest = String[]

    for arg in args
        if startswith(arg, "--route=")
            route_str = lowercase(split(arg, "="; limit = 2)[2])
            route = route_str == "paper" ? :paper :
                    route_str == "fat" ? :fat :
                    error("Unknown route $(route_str). Use --route=paper or --route=fat.")
        elseif startswith(arg, "--direct=")
            v = lowercase(split(arg, "="; limit = 2)[2])
            direct = v in ("1", "true", "yes")
        elseif startswith(arg, "--max-iter=")
            max_iter = parse(Int, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--eps=")
            eps = parse(Float64, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--decompose=")
            v = lowercase(split(arg, "="; limit = 2)[2])
            decompose = v in ("1", "true", "yes")
        elseif startswith(arg, "--")
            error("Unknown flag: $(arg)")
        else
            push!(rest, arg)
        end
    end

    return (; route, direct, max_iter, eps, decompose, positional = rest)
end

const CLI = parse_flags(ARGS)

function static_block_stats(block_sizes)
    return (
        n_psd_blocks = length(block_sizes),
        max_psd_block = isempty(block_sizes) ? 0 : maximum(block_sizes),
        sum_sq = sum(d^2 for d in block_sizes),
        sum_tri_slots = sum(d * (d + 1) ÷ 2 for d in block_sizes),
        psd_block_sizes = block_sizes,
    )
end

function route_static_stats(data, route::Symbol)
    sizes = sort(length.(PQG.route_block_bases(data, route)); rev = true)
    return static_block_stats(sizes)
end

function cosmo_optimizer(; max_iter::Int, eps::Float64, decompose::Bool)
    return optimizer_with_attributes(
        COSMO.Optimizer,
        MOI.Silent() => false,
        "verbose" => true,
        "verbose_timing" => true,
        "max_iter" => max_iter,
        "eps_abs" => eps,
        "eps_rel" => eps,
        "eps_prim_inf" => 1e-5,
        "eps_dual_inf" => 1e-5,
        "check_termination" => 200,
        "decompose" => decompose,
    )
end

function cosmo_settings(; max_iter::Int, eps::Float64, decompose::Bool)
    return COSMO.Settings(
        verbose = true,
        verbose_timing = true,
        max_iter = max_iter,
        eps_abs = eps,
        eps_rel = eps,
        eps_prim_inf = 1e-5,
        eps_dual_inf = 1e-5,
        check_termination = 200,
        decompose = decompose,
    )
end

function evaluate_constraint_residuals(data, monomap)
    n_up_val = real(PQG.evaluate_polynomial(data.number_ops.n_up, monomap))
    n_dn_val = real(PQG.evaluate_polynomial(data.number_ops.n_dn, monomap))
    d_trace_val = real(PQG.evaluate_polynomial(data.trace_ops.D_trace, monomap))
    q_trace_val = real(PQG.evaluate_polynomial(data.trace_ops.Q_trace, monomap))
    g_trace_val = real(PQG.evaluate_polynomial(data.trace_ops.G_trace, monomap))

    return (
        n_up = n_up_val,
        n_dn = n_dn_val,
        d_trace = d_trace_val,
        q_trace = q_trace_val,
        g_trace = g_trace_val,
        n_up_residual = n_up_val - 4.0,
        n_dn_residual = n_dn_val - 4.0,
        d_trace_residual = d_trace_val - data.trace_targets.D,
        q_trace_residual = q_trace_val - data.trace_targets.Q,
        g_trace_residual = g_trace_val - data.trace_targets.G,
    )
end

function markdown_report(; data, route, static_route_stats, static_other_stats,
                           route_mp_stats, route_moment_relax_time,
                           direct_mode, build_stats, solve_stats,
                           energy = nothing, residuals = nothing,
                           term_status = nothing, primal_status = nothing,
                           dual_status = nothing)
    io = IOBuffer()
    println(io, "# H4 periodic Nk=2 — PQG basis + COSMO benchmark")
    println(io)
    println(io, "Selected route: `", PQG.route_label(route), "`.")
    println(io)
    println(io, "## Shared formulation")
    println(io)
    println(io, "- explicit PQG basis size: `", length(data.basis_data.pqg_basis), "`")
    println(io, "- explicit equality constraints via `moment_eq_constraints`:")
    println(io, "  - `N_up = 4`")
    println(io, "  - `N_dn = 4`")
    println(io, "  - `Tr(²D) = ", data.trace_targets.D, "`")
    println(io, "  - `Tr(²Q) = ", data.trace_targets.Q, "`")
    println(io, "  - `Tr(²G) = ", data.trace_targets.G, "`")
    println(io)
    println(io, "The CAR-derived PQG linear maps are **implicit** in this formulation because `²D`, `²Q`, and `²G` are not independent variables here; they are different slices of one shared canonical moment map.")
    println(io)
    println(io, "## Static block comparison")
    println(io)
    println(io, "| Route | PSD blocks | largest block | Σ dim² | upper-tri slots |")
    println(io, "| :--- | ---: | ---: | ---: | ---: |")
    println(io, "| single PQG block | $(static_other_stats.n_psd_blocks) | $(static_other_stats.max_psd_block) | $(static_other_stats.sum_sq) | $(static_other_stats.sum_tri_slots) |")
    println(io, "| paper-style (ΔN, K) blocks | $(static_route_stats.n_psd_blocks) | $(static_route_stats.max_psd_block) | $(static_route_stats.sum_sq) | $(static_route_stats.sum_tri_slots) |")
    println(io)
    println(io, "Paper-style sector rows: `", [row.dim for row in data.basis_data.sector_rows], "`")
    println(io)
    println(io, "## Selected route: symbolic relaxation")
    println(io)
    println(io, "| Quantity | Value |")
    println(io, "| :--- | ---: |")
    println(io, "| route | `", PQG.route_label(route), "` |")
    println(io, "| moment_relax wall (s) | $(round(route_moment_relax_time.time; digits = 2)) |")
    println(io, "| moment_relax alloc | $(format_bytes(route_moment_relax_time.bytes)) |")
    println(io, "| PSD blocks | $(route_mp_stats.n_psd_blocks) |")
    println(io, "| largest PSD block | $(route_mp_stats.max_psd_block) |")
    println(io, "| Σ dim² | $(route_mp_stats.sum_sq) |")
    println(io, "| upper-tri PSD slots | $(route_mp_stats.sum_tri_slots) |")
    println(io, "| unique moments | $(route_mp_stats.n_unique_moments) |")
    println(io, "| total_basis length | $(route_mp_stats.total_basis_len) |")
    println(io, "| zero scalar equalities | $(route_mp_stats.n_zero_scalar_eqs) |")
    println(io)

    if direct_mode
        println(io, "## Direct COSMO path")
        println(io)
        println(io, "| Quantity | Value |")
        println(io, "| :--- | ---: |")
        println(io, "| direct variables | $(build_stats.n_vars) |")
        println(io, "| direct rows | $(build_stats.total_rows) |")
        println(io, "| direct A nnz | $(build_stats.total_nnz) |")
        println(io, "| direct constraint blocks | $(build_stats.n_constraint_blocks) |")
        println(io, "| COSMO assemble wall (s) | $(round(build_stats.assemble_time; digits = 2)) |")
        println(io, "| COSMO optimize wall (s) | $(round(solve_stats.solve_time; digits = 2)) |")
    else
        println(io, "## JuMP primal path")
        println(io)
        println(io, "| Quantity | Value |")
        println(io, "| :--- | ---: |")
        println(io, "| JuMP basis length | $(build_stats.n_basis) |")
        println(io, "| JuMP variables | $(build_stats.n_variables) |")
        println(io, "| JuMP constraints | $(build_stats.n_constraints) |")
        println(io, "| JuMP build wall (s) | $(round(build_stats.build_time; digits = 2)) |")
        println(io, "| COSMO optimize wall (s) | $(round(solve_stats.solve_time; digits = 2)) |")
    end
    println(io)

    if !isnothing(term_status)
        println(io, "## Solver termination")
        println(io)
        println(io, "- `termination_status`: `", term_status, "`")
        println(io, "- `primal_status`: `", primal_status, "`")
        println(io, "- `dual_status`: `", dual_status, "`")
        println(io)
    end

    if !isnothing(energy)
        figure_v2rdm = data.asset.preflight["figure_v2rdm_nk2"]
        println(io, "## Energy")
        println(io)
        println(io, "| Quantity | Value |")
        println(io, "| :--- | ---: |")
        println(io, "| active-space objective | $(round(energy.objective; digits = 8)) |")
        println(io, "| HF constant shift | $(round(energy.hf_shift; digits = 8)) |")
        println(io, "| total energy | $(round(energy.total; digits = 8)) |")
        println(io, "| paper V2RDM figure | $(round(figure_v2rdm; digits = 3)) |")
        println(io)
    end

    if !isnothing(residuals)
        println(io, "## Equality residuals")
        println(io, "| Quantity | Value | residual |")
        println(io, "| :--- | ---: | ---: |")
        println(io, "| N_up | $(round(residuals.n_up; digits = 8)) | $(residuals.n_up_residual) |")
        println(io, "| N_dn | $(round(residuals.n_dn; digits = 8)) | $(residuals.n_dn_residual) |")
        println(io, "| Tr(²D) | $(round(residuals.d_trace; digits = 8)) | $(residuals.d_trace_residual) |")
        println(io, "| Tr(²Q) | $(round(residuals.q_trace; digits = 8)) | $(residuals.q_trace_residual) |")
        println(io, "| Tr(²G) | $(round(residuals.g_trace; digits = 8)) | $(residuals.g_trace_residual) |")
        println(io)
    end

    return String(take!(io))
end

flushln("="^78)
flushln("H4 periodic Nk=2 — PQG basis + COSMO benchmark")
flushln("="^78)
flushln("Selected route       : ", PQG.route_label(CLI.route))
flushln("COSMO path           : ", CLI.direct ? "direct assemble/optimize" : "JuMP primal")
flushln(@sprintf("COSMO settings       : max_iter = %d, eps_abs/rel = %.1e, decompose = %s",
                 CLI.max_iter, CLI.eps, CLI.decompose))

shared_time = @timed PQG.build_h4_pqg_problem_data()
data = shared_time.value
flushln(@sprintf("[stage] shared build        : %.2f s  |  alloc %s",
                 shared_time.time, format_bytes(shared_time.bytes)))
flushln("Shared basis size      : ", length(data.basis_data.pqg_basis))
flushln("Trace targets          : Tr(²D) = ", data.trace_targets.D,
        ", Tr(²Q) = ", data.trace_targets.Q,
        ", Tr(²G) = ", data.trace_targets.G)
flushln("The CAR-derived PQG maps are implicit in the shared moment map; only particle number + traces are injected explicitly.")

fat_stats = route_static_stats(data, :fat)
paper_stats = route_static_stats(data, :paper)
flushln("\n-- static block comparison --")
flushln(@sprintf("  single PQG block            : blocks = %d | max = %d | Σ dim² = %d | tri slots = %d",
                 fat_stats.n_psd_blocks, fat_stats.max_psd_block, fat_stats.sum_sq, fat_stats.sum_tri_slots))
flushln(@sprintf("  paper (ΔN, K) blocks        : blocks = %d | max = %d | Σ dim² = %d | tri slots = %d",
                 paper_stats.n_psd_blocks, paper_stats.max_psd_block, paper_stats.sum_sq, paper_stats.sum_tri_slots))
flushln("  paper block sizes           : ", paper_stats.psd_block_sizes)

flushln("\n[stage] moment_relax         : route = ", PQG.route_label(CLI.route))
route_moment_relax_time = @timed PQG.build_route_moment_problem(data, CLI.route)
moment_problem = route_moment_relax_time.value
route_mp_stats = PQG.moment_problem_stats(moment_problem)
flushln(@sprintf("[stage] moment_relax         : %.2f s  |  alloc %s",
                 route_moment_relax_time.time, format_bytes(route_moment_relax_time.bytes)))
flushln("Route unique moments   : ", route_mp_stats.n_unique_moments)
flushln("Route total_basis      : ", route_mp_stats.total_basis_len)
flushln("Route PSD blocks       : ", route_mp_stats.psd_block_sizes)

local build_stats
local solve_stats
local term_status = nothing
local primal_status = nothing
local dual_status = nothing
local energy = nothing
local residuals = nothing

if CLI.direct
    flushln("\n[stage] direct build         : assembling sparse COSMO problem ...")
    direct_build_time = @timed build_direct_cosmo_constraints(moment_problem; report_every_fraction = 0.25)
    direct = direct_build_time.value
    build_stats = (
        n_vars = direct.n_vars,
        total_rows = direct.total_rows,
        total_nnz = direct.total_nnz,
        n_constraint_blocks = length(direct.constraints),
        build_time = direct_build_time.time,
        build_bytes = direct_build_time.bytes,
        assemble_time = 0.0,
    )
    flushln(@sprintf("[stage] direct build         : %.2f s  |  alloc %s",
                     direct_build_time.time, format_bytes(direct_build_time.bytes)))
    flushln(@sprintf("  direct variables                : %d", direct.n_vars))
    flushln(@sprintf("  direct rows                     : %d", direct.total_rows))
    flushln(@sprintf("  direct A nnz                    : %d", direct.total_nnz))

    P = spzeros(Float64, direct.n_vars, direct.n_vars)
    model = COSMO.Model{Float64}()
    settings = cosmo_settings(max_iter = CLI.max_iter, eps = CLI.eps, decompose = CLI.decompose)

    flushln("\n[stage] COSMO assemble!      : loading sparse problem into COSMO ...")
    assemble_time = @timed COSMO.assemble!(model, P, direct.q, direct.constraints; settings = settings)
    build_stats = merge(build_stats, (assemble_time = assemble_time.time,))
    flushln(@sprintf("[stage] COSMO assemble!      : %.2f s  |  alloc %s",
                     assemble_time.time, format_bytes(assemble_time.bytes)))
    flushln(@sprintf("  model size m × n                : %d × %d", model.p.model_size[1], model.p.model_size[2]))
    flushln(@sprintf("  assembled KKT A nnz             : %d", nnz(model.p.A)))
    flushln(@sprintf("  assembled cone count            : %d", length(model.p.C.sets)))

    flushln("\n[stage] COSMO optimize!      : running direct solver ...")
    flushln("---------------- begin COSMO output ----------------")
    solve_time = @timed result = COSMO.optimize!(model)
    flushln("---------------- end   COSMO output ----------------")
    res = solve_time.value
    solve_stats = (solve_time = solve_time.time, solve_bytes = solve_time.bytes, result = res)
    flushln(@sprintf("[stage] COSMO optimize!      : %.2f s  |  alloc %s",
                     solve_time.time, format_bytes(solve_time.bytes)))
    flushln("Status                  : ", res.status)
    flushln("Iterations              : ", res.iter)
    if hasproperty(res, :obj_val)
        flushln("Objective               : ", res.obj_val)
    end
else
    optimizer = cosmo_optimizer(max_iter = CLI.max_iter, eps = CLI.eps, decompose = CLI.decompose)
    flushln("\n[stage] JuMP primal build    : assembling JuMP primal model ...")
    primal_build = @timed build_primal_jump_model(moment_problem, optimizer; typed_exprs = false)
    jump = primal_build.value
    summary = jump_model_summary(jump.model)
    build_stats = (
        n_basis = jump.n_basis,
        n_variables = summary.n_variables,
        n_constraints = summary.n_constraints,
        build_time = primal_build.time,
        build_bytes = primal_build.bytes,
        model = jump.model,
    )
    flushln(@sprintf("[stage] JuMP primal build    : %.2f s  |  alloc %s",
                     primal_build.time, format_bytes(primal_build.bytes)))
    flushln(@sprintf("  JuMP variables                : %d", summary.n_variables))
    flushln(@sprintf("  JuMP constraints              : %d", summary.n_constraints))

    flushln("\n[stage] optimize! (COSMO)    : running JuMP primal solver ...")
    flushln("---------------- begin COSMO output ----------------")
    solve_time = @timed optimize!(jump.model)
    flushln("---------------- end   COSMO output ----------------")
    solve_stats = (solve_time = solve_time.time, solve_bytes = solve_time.bytes, model = jump.model)
    flushln(@sprintf("[stage] optimize! (COSMO)    : %.2f s  |  alloc %s",
                     solve_time.time, format_bytes(solve_time.bytes)))

    term_status = termination_status(jump.model)
    primal_status = primal_status(jump.model)
    dual_status = dual_status(jump.model)
    flushln("termination_status      : ", term_status)
    flushln("primal_status           : ", primal_status)
    flushln("dual_status             : ", dual_status)

    if term_status in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.LOCALLY_SOLVED)
        result = solve_moment_problem(moment_problem, optimizer; silent = true)
        residuals = evaluate_constraint_residuals(data, result.monomap)
        energy_summary_vals = energy_summary(data.asset, result.objective)
        energy = (objective = result.objective,
                  hf_shift = energy_summary_vals.hf_shift,
                  total = energy_summary_vals.e_total)
    end
end

if !isempty(CLI.positional)
    out_path = first(CLI.positional)
    report = markdown_report(
        data = data,
        route = CLI.route,
        static_route_stats = paper_stats,
        static_other_stats = fat_stats,
        route_mp_stats = route_mp_stats,
        route_moment_relax_time = route_moment_relax_time,
        direct_mode = CLI.direct,
        build_stats = build_stats,
        solve_stats = solve_stats,
        energy = energy,
        residuals = residuals,
        term_status = term_status,
        primal_status = primal_status,
        dual_status = dual_status,
    )
    mkpath(dirname(out_path))
    open(out_path, "w") do io
        write(io, report)
    end
    flushln("\nReport written to ", out_path)
end
