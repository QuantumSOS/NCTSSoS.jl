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

flushln(@sprintf("\nDirect COSMO config: max_iter = %d, eps_abs/rel = %.1e, decompose = %s",
                 CLI.max_iter, CLI.eps, CLI.decompose))

flushln("\n[stage] direct data build     : symbolic moment problem -> sparse COSMO matrices ...")
direct_build = @timed build_direct_cosmo_constraints(moment_problem)
direct = direct_build.value
flushln(@sprintf("[stage] direct data build     : %.2f s  |  alloc %s",
                 direct_build.time, format_bytes(direct_build.bytes)))
flushln(@sprintf("  direct objective variables        : %d", direct.n_vars))
flushln(@sprintf("  direct constraint rows            : %d", direct.total_rows))
flushln(@sprintf("  direct A nnz                      : %d", direct.total_nnz))
flushln(@sprintf("  number of COSMO constraints       : %d", length(direct.constraints)))

settings = COSMO.Settings{Float64}(;
    verbose = true,
    verbose_timing = true,
    max_iter = CLI.max_iter,
    eps_abs = CLI.eps,
    eps_rel = CLI.eps,
    eps_prim_inf = 1e-5,
    eps_dual_inf = 1e-5,
    check_termination = 200,
    decompose = CLI.decompose,
)

P = spzeros(Float64, direct.n_vars, direct.n_vars)
model = COSMO.Model{Float64}()

flushln("\n[stage] COSMO assemble!      : loading sparse problem into COSMO ...")
assemble_time = @timed COSMO.assemble!(model, P, direct.q, direct.constraints; settings = settings)
flushln(@sprintf("[stage] COSMO assemble!      : %.2f s  |  alloc %s",
                 assemble_time.time, format_bytes(assemble_time.bytes)))
flushln(@sprintf("  model size m × n                  : %d × %d", model.p.model_size[1], model.p.model_size[2]))
flushln(@sprintf("  assembled KKT A nnz               : %d", nnz(model.p.A)))
flushln(@sprintf("  assembled b nnz                   : %d", count(!iszero, model.p.b)))
flushln(@sprintf("  assembled cone count              : %d", length(model.p.C.sets)))

flushln("\n[stage] COSMO optimize!      : running direct solver ...")
flushln("---------------- begin COSMO output ----------------")
solve_time = @timed result = COSMO.optimize!(model)
flushln("---------------- end   COSMO output ----------------")
res = solve_time.value
flushln(@sprintf("[stage] COSMO optimize!      : %.2f s  |  alloc %s",
                 solve_time.time, format_bytes(solve_time.bytes)))
flushln("\n-- COSMO direct termination --")
flushln("  status              : ", res.status)
flushln("  iterations          : ", res.iter)
flushln(@sprintf("  objective           : %.8f", res.obj_val))
flushln(res.times)

energy = energy_summary(asset, res.obj_val)
flushln("\n-- Energy (Ha / unit cell) --")
flushln(@sprintf("  active-space electronic (SDP obj.)   : %.8f", res.obj_val))
flushln(@sprintf("  + HF constant shift                  : %.8f", energy.hf_shift))
flushln(@sprintf("  => total energy                      : %.8f", energy.e_total))
flushln(@sprintf("  paper figure V2RDM[4,8]              : %.3f  (± %.3f)",
                 energy.figure_v2rdm, energy.figure_unc))
flushln(@sprintf("  paper figure HF                      : %.3f", energy.figure_hf))

flushln("\n-- Aggregate --")
flushln(@sprintf("  total wall time (all stages) : %.2f s",
                 prep.build_time.time + prep.sparsity_time.time + prep.mp_time.time +
                 direct_build.time + assemble_time.time + solve_time.time))
