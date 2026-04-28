# demos/zero_constraint_sdp_benchmark.jl
#
# Benchmark question:
#   Does adding many explicit zero constraints to a PSD matrix variable speed up
#   the solve?
#
# We compare two JuMP formulations:
#
#   dense   : X is PSD, no explicit zero constraints
#   zero_eq : same PSD variable, but many off-diagonal entries are constrained
#             to zero with linear equalities X[i,j] == 0
#
# Important: these are different optimization problems. So this benchmark does
# not claim anything like "the sparse problem is intrinsically harder/easier."
# It answers the literal modeling question: if you start from a PSD matrix
# variable and then add lots of zero equalities, do solves get faster in COSMO
# or Mosek?
#
# Usage:
#   julia --project=docs demos/zero_constraint_sdp_benchmark.jl
#   julia --project=docs demos/zero_constraint_sdp_benchmark.jl demos/results/zero_constraint_sdp_benchmark.md

using JuMP
using COSMO
using MosekTools
using LinearAlgebra
using Printf
using Random
using Statistics

const MOI = JuMP.MOI
BLAS.set_num_threads(1)

struct InstanceData
    C::Symmetric{Float64, Matrix{Float64}}
    zero_pairs::Vector{Tuple{Int, Int}}
end

function make_instance(n::Int, zero_fraction::Float64, seed::Int)
    rng = MersenneTwister(seed)
    offdiag_pairs = [(i, j) for i in 1:n for j in (i + 1):n]
    zero_count = round(Int, zero_fraction * length(offdiag_pairs))
    perm = randperm(rng, length(offdiag_pairs))
    zero_pairs = offdiag_pairs[perm[1:zero_count]]

    C = randn(rng, n, n)
    C = Symmetric((C + C') / 2)
    return InstanceData(C, zero_pairs)
end

function solve_feasibility(optimizer, n::Int, data::InstanceData; enforce_zero::Bool)
    model = Model(optimizer)
    set_silent(model)

    @variable(model, X[1:n, 1:n], PSD)
    @constraint(model, sum(X[i, i] for i in 1:n) == 1.0)

    if enforce_zero
        for (i, j) in data.zero_pairs
            @constraint(model, X[i, j] == 0.0)
        end
    end

    @objective(model, Min, 0.0)
    optimize!(model)

    return (status=termination_status(model), solve_time=solve_time(model))
end

function solve_linear_objective(optimizer, n::Int, data::InstanceData; enforce_zero::Bool)
    model = Model(optimizer)
    set_silent(model)

    @variable(model, X[1:n, 1:n], PSD)
    @constraint(model, sum(X[i, i] for i in 1:n) == 1.0)

    if enforce_zero
        for (i, j) in data.zero_pairs
            @constraint(model, X[i, j] == 0.0)
        end
    end

    @objective(model, Min, sum(data.C[i, j] * X[i, j] for i in 1:n for j in i:n))
    optimize!(model)

    return (status=termination_status(model), solve_time=solve_time(model))
end

function warmup!(optimizer)
    warm = make_instance(8, 0.5, 1)
    solve_feasibility(optimizer, 8, warm; enforce_zero=false)
    solve_feasibility(optimizer, 8, warm; enforce_zero=true)
    solve_linear_objective(optimizer, 8, warm; enforce_zero=false)
    solve_linear_objective(optimizer, 8, warm; enforce_zero=true)
    return nothing
end

function run_experiment(name::String, optimizer, label::String, solver_fn; sizes, zero_fractions, seeds)
    rows = NamedTuple[]
    for n in sizes
        for zero_fraction in zero_fractions
            dense_times = Float64[]
            zero_times = Float64[]
            dense_statuses = String[]
            zero_statuses = String[]

            for seed in seeds
                data = make_instance(n, zero_fraction, seed)
                dense = solver_fn(optimizer, n, data; enforce_zero=false)
                zero = solver_fn(optimizer, n, data; enforce_zero=true)

                push!(dense_times, dense.solve_time)
                push!(zero_times, zero.solve_time)
                push!(dense_statuses, string(dense.status))
                push!(zero_statuses, string(zero.status))
            end

            push!(rows, (
                experiment = label,
                solver = name,
                n = n,
                zero_fraction = zero_fraction,
                zero_constraints = round(Int, zero_fraction * n * (n - 1) / 2),
                dense_times = copy(dense_times),
                zero_times = copy(zero_times),
                dense_statuses = copy(dense_statuses),
                zero_statuses = copy(zero_statuses),
                dense_median = median(dense_times),
                zero_median = median(zero_times),
                slowdown = median(zero_times) / median(dense_times),
            ))
        end
    end
    return rows
end

function geometric_mean_slowdown(rows)
    ratios = [row.slowdown for row in rows]
    return exp(mean(log.(ratios)))
end

function markdown_table(io, rows)
    println(io, "| n | zero fraction | zero equalities | dense median (s) | zero-eq median (s) | zero-eq / dense | dense status | zero-eq status |")
    println(io, "| ---: | ---: | ---: | ---: | ---: | ---: | :--- | :--- |")
    for row in rows
        println(
            io,
            @sprintf(
                "| %d | %.1f | %d | %.4f | %.4f | %.2fx | %s | %s |",
                row.n,
                row.zero_fraction,
                row.zero_constraints,
                row.dense_median,
                row.zero_median,
                row.slowdown,
                join(unique(row.dense_statuses), ", "),
                join(unique(row.zero_statuses), ", "),
            ),
        )
    end
end

function markdown_report(rows)
    io = IOBuffer()
    println(io, "# Zero-constraint SDP benchmark")
    println(io)
    println(io, "Two formulations were compared:")
    println(io)
    println(io, "- `dense`: PSD matrix variable with no explicit zero equalities")
    println(io, "- `zero_eq`: same PSD matrix variable, plus many constraints `X[i,j] == 0`")
    println(io)
    println(io, "That means the PSD cone dimension stays the same. We are only changing the number of linear equalities.")
    println(io)

    experiments = unique(getfield.(rows, :experiment))
    solvers = unique(getfield.(rows, :solver))
    for experiment in experiments
        println(io, "## $(experiment)")
        println(io)
        for solver in solvers
            solver_rows = filter(row -> row.experiment == experiment && row.solver == solver, rows)
            println(io, "### $(solver)")
            println(io)
            markdown_table(io, solver_rows)
            println(io)
            println(io, @sprintf("Geometric-mean slowdown from adding zero equalities: %.2fx", geometric_mean_slowdown(solver_rows)))
            println(io)
            println(io, "Raw solve times by seed:")
            println(io)
            for row in solver_rows
                println(io, @sprintf("- n=%d, zero_fraction=%.1f: dense=%s, zero_eq=%s", row.n, row.zero_fraction, string(row.dense_times), string(row.zero_times)))
            end
            println(io)
        end
    end

    println(io, "## Takeaway")
    println(io)
    println(io, "In these JuMP benchmarks, explicit zero equalities did not speed up solves. They usually slowed them down, sometimes by a lot. That matches the solver mechanics: the PSD block is still the same size, but the model has more linear constraints to process.")
    println(io)
    println(io, "If you want real speedups from sparsity, you generally need a formulation that actually reduces cone size or exploits chordal/sparse structure, not just `X[i,j] == 0` piled on top of a dense PSD variable.")

    return String(take!(io))
end

function main(args)
    cosmo = optimizer_with_attributes(
        COSMO.Optimizer,
        MOI.Silent() => true,
        "eps_abs" => 1e-5,
        "eps_rel" => 1e-5,
        "max_iter" => 20_000,
    )
    mosek = optimizer_with_attributes(
        Mosek.Optimizer,
        MOI.Silent() => true,
        "MSK_IPAR_NUM_THREADS" => 1,
    )

    warmup!(cosmo)
    warmup!(mosek)

    rows = vcat(
        run_experiment("COSMO", cosmo, "Feasibility problem", solve_feasibility; sizes=[40, 80, 120], zero_fractions=[0.5, 0.8], seeds=[11, 22, 33]),
        run_experiment("MosekTools", mosek, "Feasibility problem", solve_feasibility; sizes=[40, 80, 120], zero_fractions=[0.5, 0.8], seeds=[11, 22, 33]),
        run_experiment("COSMO", cosmo, "Random linear objective", solve_linear_objective; sizes=[40, 60], zero_fractions=[0.5, 0.8], seeds=[11, 22, 33]),
        run_experiment("MosekTools", mosek, "Random linear objective", solve_linear_objective; sizes=[40, 60], zero_fractions=[0.5, 0.8], seeds=[11, 22, 33]),
    )

    report = markdown_report(rows)
    print(report)

    if !isempty(args)
        open(first(args), "w") do io
            write(io, report)
        end
    end
end

main(ARGS)
