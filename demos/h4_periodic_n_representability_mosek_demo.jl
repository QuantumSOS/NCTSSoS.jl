# demos/h4_periodic_n_representability_mosek_demo.jl
#
# Periodic H4, Nk = 2, written in the direct N-representability / V2RDM form:
#
#   minimize    Tr(K D)
#   subject to  D ⪰ 0,
#               Q(D) ⪰ 0,
#               G(D) ⪰ 0,
#               N_up = N_dn = 4,
#               Tr(D), Tr(Q), Tr(G) fixed to the fermionic PQG targets.
#
# The free variable is only the 2-RDM D. The Q and G matrices are built as
# affine images of D using the usual contraction identities. This is the sane
# thing to hand to Mosek here; the generic order-2 moment route is bigger and,
# for this benchmark, mostly a good way to burn RAM.
#
# Examples:
#   julia --project demos/h4_periodic_n_representability_mosek_demo.jl
#   julia --project demos/h4_periodic_n_representability_mosek_demo.jl --report-only
#   julia --project demos/h4_periodic_n_representability_mosek_demo.jl --refinement=k_only --time-limit=60

using JuMP
using MosekTools
using NCTSSoS
using Printf

const MOI = JuMP.MOI

include(joinpath(@__DIR__, "H4PeriodicNativeV2RDMHelpers.jl"))
using .H4PeriodicNativeV2RDMHelpers

function parse_refinement(value::AbstractString)
    key = lowercase(value)
    key in ("k", "k_only", "k-only") && return :k_only
    key in ("spin", "spin_resolved", "spin-resolved") && return :spin_resolved
    error("Unknown refinement $(value). Use --refinement=k_only or --refinement=spin_resolved.")
end

function parse_flags(args)
    refinement = :spin_resolved
    time_limit = 30.0
    include_traces = true
    num_threads = 0
    log_level = 10
    report_only = false
    for arg in args
        if startswith(arg, "--refinement=")
            refinement = parse_refinement(split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--time-limit=")
            time_limit = parse(Float64, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--include-traces=")
            include_traces = parse(Bool, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--num-threads=")
            num_threads = parse(Int, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--log-level=")
            log_level = parse(Int, split(arg, "="; limit = 2)[2])
        elseif arg == "--report-only"
            report_only = true
        elseif startswith(arg, "--report-only=")
            report_only = parse(Bool, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--")
            error("Unknown flag: $(arg)")
        end
    end
    return (; refinement, time_limit, include_traces, num_threads, log_level, report_only)
end

function print_size_report(label::AbstractString, size_summary)
    println("  ", rpad(label, 28), "D/Q blocks = ", size_summary.d_only_block_sizes,
            " | G blocks = ", size_summary.g_block_sizes)
    println("  ", rpad("", 28), "native D vars = ", size_summary.jump_variables,
            " | D-only moment real vars = ", size_summary.d_only_moment_vector_real_variables,
            " | PSD rows = ", size_summary.solver_psd_total_rows,
            " | total rows = ", size_summary.solver_total_rows)
end

function print_formulation_summary(refinement::Symbol)
    println("formulation:")
    println("  free variable : Hermitian ²D blocks")
    println("  derived blocks: ¹D(D), ²Q(D), ²G(D)")
    println("  cones         : ²D ⪰ 0, ²Q(D) ⪰ 0, ²G(D) ⪰ 0")
    println("  equalities    : N↑ = 4, N↓ = 4, plus Tr(D), Tr(Q), Tr(G)")
    println("  refinement    : ", refinement_label(refinement))
end

function try_optimize!(model)
    try
        solve = @timed optimize!(model)
        return (; ok = true, time = solve.time, bytes = solve.bytes,
                 error = nothing, bt = nothing)
    catch err
        return (; ok = false, time = nothing, bytes = 0,
                 error = err, bt = catch_backtrace())
    end
end

const CLI = parse_flags(ARGS)
const CURRENT_MOMENT_ROUTE = (
    real_variables = 340_546,
    psd_rows = 1_545_187,
    equality_rows = 16_012,
    total_rows = 1_561_199,
)

optimizer_factory = optimizer_with_attributes(
    Mosek.Optimizer,
    MOI.Silent() => false,
    MOI.TimeLimitSec() => CLI.time_limit,
    "MSK_DPAR_OPTIMIZER_MAX_TIME" => CLI.time_limit,
    "MSK_IPAR_NUM_THREADS" => CLI.num_threads,
    "MSK_IPAR_LOG" => CLI.log_level,
)

println("H4 periodic Nk=2 — direct N-representability / V2RDM demo with Mosek")
println(@sprintf("MosekTools config: time_limit = %.1f s, num_threads = %d, log_level = %d",
                 CLI.time_limit, CLI.num_threads, CLI.log_level))
println("include trace equalities: ", CLI.include_traces)
println("report only            : ", CLI.report_only)
print_formulation_summary(CLI.refinement)

build = @timed build_native_jump_model(;
    optimizer_factory = optimizer_factory,
    include_trace_constraints = CLI.include_traces,
    refinement = CLI.refinement,
)
built = build.value
summary = jump_model_summary(built.model)
size_summary = built.size_summary
size_k = native_size_summary(built.data;
    refinement = :k_only,
    include_trace_constraints = CLI.include_traces)
size_spin = native_size_summary(built.data;
    refinement = :spin_resolved,
    include_trace_constraints = CLI.include_traces)

println(@sprintf("\n[stage] native JuMP build : %.2f s  |  alloc %.2f GiB",
                 build.time, build.bytes / 2.0^30))
println("\nsize comparison before solve:")
println("  current generic order-2 route      : real vars = ", CURRENT_MOMENT_ROUTE.real_variables,
        ", total rows = ", CURRENT_MOMENT_ROUTE.total_rows)
print_size_report("K-only native bridge", size_k)
print_size_report("spin-resolved native bridge", size_spin)

println("\nselected model counts:")
println("  JuMP num_variables                 : ", summary.n_variables)
println("  JuMP num_constraints               : ", summary.n_constraints)
println("  JuMP constraint breakdown          : ", summary.constraint_counts)
println("  Hermitian PSD blocks               : ", size_summary.jump_psd_blocks)
println("  scalar equalities                  : ", size_summary.jump_scalar_equalities)
println("  solver PSD rows after real lift    : ", size_summary.solver_psd_total_rows)
println("  solver total rows (+ equalities)   : ", size_summary.solver_total_rows)
println("  Hermitian block sizes              : ", size_summary.hermitian_block_sizes)
println("  real-lift block sizes              : ", size_summary.solver_real_lift_block_sizes)
println("  build stage times                  : ", built.build_times)
println(@sprintf("  spin-vs-K variable reduction       : %.2fx",
                 size_k.jump_variables / size_spin.jump_variables))
println(@sprintf("  spin-vs-K total-row reduction      : %.2fx",
                 size_k.solver_total_rows / size_spin.solver_total_rows))
println(@sprintf("  spin-vs-current variable reduction : %.2fx",
                 CURRENT_MOMENT_ROUTE.real_variables / size_spin.jump_variables))
println(@sprintf("  spin-vs-current total-row reduction: %.2fx",
                 CURRENT_MOMENT_ROUTE.total_rows / size_spin.solver_total_rows))

if CLI.report_only
    println("\nreport-only requested; skipping optimize!.")
    exit()
end

println("\n[stage] optimize! (Mosek) : running selected N-representability model ...")
println("---------------- begin Mosek output ----------------")
flush(stdout)
solve = try_optimize!(built.model)
println("---------------- end   Mosek output ----------------")

if solve.ok
    println(@sprintf("[stage] optimize! (Mosek) : %.2f s  |  alloc %.2f GiB",
                     solve.time, solve.bytes / 2.0^30))
    term = termination_status(built.model)
    primal = primal_status(built.model)
    dual = dual_status(built.model)
    println("\ntermination_status : ", term)
    println("primal_status      : ", primal)
    println("dual_status        : ", dual)
    if has_values(built.model)
        println(@sprintf("objective_value    : %.8f", objective_value(built.model)))
    end
else
    println("[stage] optimize! (Mosek) : threw before completion")
    println(sprint(showerror, solve.error, solve.bt))
    Base.exit(1)
end
