include(joinpath(@__DIR__, "h4_periodic_probe_common.jl"))
using TOML

function parse_flags(args)
    max_iter = 5_000
    eps = 1e-6
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
const HUBBARD_EXPECTATIONS_PATH = joinpath(pkgdir(NCTSSoS), "test", "data", "expectations", "hubbard.toml")

function expectation_by_id(path::AbstractString, id::AbstractString)
    data = TOML.parsefile(path)
    for case in data["cases"]
        case["id"] == id && return case["expected"]
    end
    error("Expectation not found: $(id)")
end

function hubbard_problem(nsites::Int; t::Real, U::Real, periodic::Bool)
    registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
        ("c_up", 1:nsites),
        ("c_dn", 1:nsites),
    ])

    bonds = periodic ? [(i, mod1(i + 1, nsites)) for i in 1:nsites] : [(i, i + 1) for i in 1:nsites-1]

    hopping = -t * sum(
        c_up_dag[i] * c_up[j] + c_up_dag[j] * c_up[i] +
        c_dn_dag[i] * c_dn[j] + c_dn_dag[j] * c_dn[i]
        for (i, j) in bonds
    )
    interaction = U * sum(
        (c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
        for i in 1:nsites
    )

    return registry, (c_up, c_up_dag), (c_dn, c_dn_dag), hopping + interaction
end

function run_case(label::String, expectation_id::String, nsites::Int; canonical::Bool)
    flushln("\n", "="^78)
    flushln(label)
    flushln("="^78)

    oracle = expectation_by_id(HUBBARD_EXPECTATIONS_PATH, expectation_id)
    registry, (c_up, c_up_dag), (c_dn, c_dn_dag), ham = hubbard_problem(nsites; t = 1.0, U = 4.0, periodic = true)

    build_time = @timed begin
        if canonical
            n_up_total = 1.0 * sum(c_up_dag[i] * c_up[i] for i in 1:nsites)
            n_dn_total = 1.0 * sum(c_dn_dag[i] * c_dn[i] for i in 1:nsites)
            polyopt(ham, registry;
                    moment_eq_constraints = [n_up_total - (nsites / 2) * one(ham),
                                             n_dn_total - (nsites / 2) * one(ham)])
        else
            polyopt(ham, registry)
        end
    end
    pop = build_time.value
    flushln(@sprintf("[stage] build polyopt        : %.2f s  |  alloc %s",
                     build_time.time, format_bytes(build_time.bytes)))

    sparsity_time = @timed compute_sparsity(pop, SolverConfig(optimizer = nothing, order = 2, ts_algo = MMD()))
    sparsity = sparsity_time.value
    flushln(@sprintf("[stage] compute_sparsity     : %.2f s  |  alloc %s",
                     sparsity_time.time, format_bytes(sparsity_time.bytes)))

    blocks = sparsity.cliques_term_sparsities[1][1].block_bases
    sizes = summarize_blocks(blocks)
    stats = block_stats(sizes)
    flushln(@sprintf("  TS PSD blocks                     : %d", stats.n_blocks))
    flushln(@sprintf("  TS largest block                  : %d × %d", stats.max_size, stats.max_size))
    flushln(@sprintf("  sum of upper-tri PSD scalar slots : %d", stats.total_scalar_entries))

    mp_time = @timed NCTSSoS.moment_relax(pop, sparsity.corr_sparsity,
                                          sparsity.cliques_term_sparsities)
    moment_problem = mp_time.value
    flushln(@sprintf("[stage] moment_relax         : %.2f s  |  alloc %s",
                     mp_time.time, format_bytes(mp_time.bytes)))
    flushln(@sprintf("  unique moment variables           : %d",
                     moment_problem.n_unique_moment_matrix_elements))

    optimizer = optimizer_with_attributes(
        COSMO.Optimizer,
        MOI.Silent() => false,
        "verbose" => true,
        "verbose_timing" => true,
        "max_iter" => CLI.max_iter,
        "eps_abs" => CLI.eps,
        "eps_rel" => CLI.eps,
        "eps_prim_inf" => 1e-6,
        "eps_dual_inf" => 1e-6,
        "rho" => 1.0,
        "adaptive_rho" => true,
        "alpha" => 1.0,
        "scaling" => 10,
        "check_termination" => 50,
        "decompose" => CLI.decompose,
    )

    build_primal = @timed build_primal_jump_model(moment_problem, optimizer; typed_exprs = false)
    built = build_primal.value
    jump_model = built.model
    summary = jump_model_summary(jump_model)
    flushln(@sprintf("[stage] primal moment build  : %.2f s  |  alloc %s",
                     build_primal.time, format_bytes(build_primal.bytes)))
    flushln(@sprintf("  JuMP num_variables                : %d", summary.n_variables))
    flushln(@sprintf("  JuMP num_constraints (all types)  : %d", summary.n_constraints))

    flushln("[stage] optimize! (COSMO)    : running SDP solver ...")
    flushln("---------------- begin COSMO output ----------------")
    solve_wall = @timed optimize!(jump_model)
    flushln("---------------- end   COSMO output ----------------")
    flushln(@sprintf("[stage] optimize! (COSMO)    : %.2f s  |  alloc %s",
                     solve_wall.time, format_bytes(solve_wall.bytes)))

    obj = objective_value(jump_model)
    status = termination_status(jump_model)
    flushln("  termination_status                : ", status)
    flushln(@sprintf("  objective                         : %.8f", obj))
    flushln(@sprintf("  oracle objective                  : %.8f", oracle["objective"]))
    flushln(@sprintf("  objective error                   : %+0.3e", obj - oracle["objective"]))

    return (
        label = label,
        expectation_id = expectation_id,
        build_time = build_time.time,
        sparsity_time = sparsity_time.time,
        mp_time = mp_time.time,
        build_primal_time = build_primal.time,
        solve_time = solve_wall.time,
        stats = stats,
        nuniq = moment_problem.n_unique_moment_matrix_elements,
        nvars = summary.n_variables,
        ncons = summary.n_constraints,
        objective = obj,
        oracle = oracle["objective"],
        status = status,
    )
end

flushln(@sprintf("Smaller-system probe config: COSMO max_iter = %d, eps_abs/rel = %.1e, decompose = %s",
                 CLI.max_iter, CLI.eps, CLI.decompose))
flushln("Repo check: `test/H4PeriodicAssets.jl` only vendors the H4 Nk=2 asset, so this probe uses the smallest in-repo periodic fermionic stand-ins instead.")

case_n2 = run_case("Periodic Hubbard N=2 (grand-canonical, order 2)", "periodic_n2_u4_order2", 2; canonical = false)
case_n4 = run_case("Periodic Hubbard N=4 (canonical half-filling, order 2)", "periodic_n4_u4_order2_canonical", 4; canonical = true)

flushln("\n", "="^78)
flushln("Smaller-system summary")
flushln("="^78)
for case in (case_n2, case_n4)
    total = case.build_time + case.sparsity_time + case.mp_time + case.build_primal_time + case.solve_time
    flushln(case.label)
    flushln(@sprintf("  total wall                        : %.2f s", total))
    flushln(@sprintf("  PSD upper-tri slots               : %d", case.stats.total_scalar_entries))
    flushln(@sprintf("  unique moments                    : %d", case.nuniq))
    flushln(@sprintf("  JuMP vars / constraints           : %d / %d", case.nvars, case.ncons))
    flushln(@sprintf("  objective error vs oracle         : %+0.3e", case.objective - case.oracle))
    flushln("  termination_status                : ", case.status)
end
