include(joinpath(@__DIR__, "h4_periodic_probe_common.jl"))

function parse_flags(args)
    max_iter = 10
    eps = 1e-4
    decompose = false
    for arg in args
        if startswith(arg, "--max-iter=")
            max_iter = parse(Int, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--eps=")
            eps = parse(Float64, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--decompose=")
            v = lowercase(split(arg, "="; limit = 2)[2])
            decompose = v in ("1", "true", "yes")
        elseif startswith(arg, "--")
            error("Unknown flag: $(arg)")
        end
    end
    return (; max_iter, eps, decompose)
end

const CLI = parse_flags(ARGS)

prep = build_h4_nk2_moment_problem()
asset = prep.asset
moment_problem = prep.moment_problem

optimizer = optimizer_with_attributes(
    COSMO.Optimizer,
    MOI.Silent() => false,
    "verbose" => true,
    "verbose_timing" => true,
    "max_iter" => CLI.max_iter,
    "eps_abs" => CLI.eps,
    "eps_rel" => CLI.eps,
    "eps_prim_inf" => 1e-5,
    "eps_dual_inf" => 1e-5,
    "check_termination" => 200,
    "decompose" => CLI.decompose,
)

flushln(@sprintf("\nTyped-expression COSMO config: max_iter = %d, eps_abs/rel = %.1e, decompose = %s",
                 CLI.max_iter, CLI.eps, CLI.decompose))

flushln("\n[stage] primal moment build  : assembling JuMP primal moment SDP with typed affine matrices ...")
build_primal = @timed build_primal_jump_model(moment_problem, optimizer; typed_exprs = true)
built = build_primal.value
jump_model = built.model
summary = jump_model_summary(jump_model)
flushln(@sprintf("[stage] primal moment build  : %.2f s  |  alloc %s",
                 build_primal.time, format_bytes(build_primal.bytes)))
flushln(@sprintf("  JuMP num_variables                : %d", summary.n_variables))
flushln(@sprintf("  JuMP num_constraints (all types)  : %d", summary.n_constraints))

flushln("\n[stage] optimize! (COSMO)    : running SDP solver ...")
flushln("---------------- begin COSMO output ----------------")
solve_wall = @timed optimize!(jump_model)
flushln("---------------- end   COSMO output ----------------")
flushln(@sprintf("[stage] optimize! (COSMO)    : %.2f s  |  alloc %s",
                 solve_wall.time, format_bytes(solve_wall.bytes)))

term_status = termination_status(jump_model)
primal_stat = primal_status(jump_model)
dual_stat = dual_status(jump_model)
moi_solve_t = try solve_time(jump_model) catch; NaN end
obj_active = try objective_value(jump_model) catch; NaN end
energy = energy_summary(asset, obj_active)

flushln("\n-- COSMO termination --")
flushln("  termination_status : ", term_status)
flushln("  primal_status      : ", primal_stat)
flushln("  dual_status        : ", dual_stat)
flushln(@sprintf("  MOI solve_time     : %.2f s", moi_solve_t))
flushln(solution_summary(jump_model))

flushln("\n-- Energy (Ha / unit cell) --")
flushln(@sprintf("  active-space electronic (SDP obj.)   : %.8f", obj_active))
flushln(@sprintf("  + HF constant shift                  : %.8f", energy.hf_shift))
flushln(@sprintf("  => total energy                      : %.8f", energy.e_total))
flushln(@sprintf("  paper figure V2RDM[4,8]              : %.3f  (± %.3f)",
                 energy.figure_v2rdm, energy.figure_unc))
flushln(@sprintf("  paper figure HF                      : %.3f", energy.figure_hf))

flushln("\n-- Aggregate --")
flushln(@sprintf("  total wall time (all stages) : %.2f s",
                 prep.build_time.time + prep.sparsity_time.time + prep.mp_time.time +
                 build_primal.time + solve_wall.time))
