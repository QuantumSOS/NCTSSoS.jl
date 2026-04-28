# demos/scalar_matrix_psd_benchmark.jl
#
# Compare three modeling styles for an SDP matrix:
#
# 1. `matrix_dense`:
#       @variable(model, X[1:n,1:n], PSD)
#
# 2. `matrix_zeroeq`:
#       same PSD matrix variable, plus many constraints X[i,j] == 0
#
# 3. `scalar_dense` / `scalar_zeroeq`:
#       declare scalar variables for the upper triangle, assemble them into a
#       symmetric matrix M, then constrain M in PSDCone(). Optionally add many
#       z[k] == 0 constraints.
#
# 4. `scalar_structural_zero`:
#       declare scalar variables only for the entries that are allowed to be
#       nonzero, fill the forbidden entries with literal 0.0, then constrain
#       the assembled matrix to be PSD.
#
# This isolates the follow-up question:
#   "What if I declare scalar variables first, then set some to zero, and then
#   require the filled matrix to be PSD?"
#
# Usage:
#   julia --project=docs demos/scalar_matrix_psd_benchmark.jl
#   julia --project=docs demos/scalar_matrix_psd_benchmark.jl demos/results/scalar_matrix_psd_benchmark.md

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
    zero_set::Set{Tuple{Int, Int}}
    pairs::Vector{Tuple{Int, Int}}
end

function make_instance(n::Int, zero_fraction::Float64, seed::Int)
    rng = MersenneTwister(seed)
    pairs = [(i, j) for i in 1:n for j in (i + 1):n]
    zero_count = round(Int, zero_fraction * length(pairs))
    perm = randperm(rng, length(pairs))
    zero_pairs = pairs[perm[1:zero_count]]

    C = randn(rng, n, n)
    C = Symmetric((C + C') / 2)
    return InstanceData(C, zero_pairs, Set(zero_pairs), pairs)
end

function scalar_model_full(model::Model, n::Int, pairs)
    @variable(model, d[1:n])
    @variable(model, z[1:length(pairs)])

    pair_index = Dict(pair => idx for (idx, pair) in enumerate(pairs))
    M = Matrix{AffExpr}(undef, n, n)
    for i in 1:n
        M[i, i] = d[i]
    end
    for (idx, (i, j)) in enumerate(pairs)
        M[i, j] = z[idx]
        M[j, i] = z[idx]
    end
    @constraint(model, M in PSDCone())

    return d, z, pair_index
end

function scalar_model_sparse(model::Model, n::Int, free_pairs)
    @variable(model, d[1:n])
    @variable(model, z[1:length(free_pairs)])

    M = Matrix{AffExpr}(undef, n, n)
    for i in 1:n
        M[i, i] = d[i]
    end
    for i in 1:n, j in (i + 1):n
        M[i, j] = AffExpr(0.0)
        M[j, i] = AffExpr(0.0)
    end
    for (idx, (i, j)) in enumerate(free_pairs)
        M[i, j] = z[idx]
        M[j, i] = z[idx]
    end
    @constraint(model, M in PSDCone())

    return d, z
end

function solve_formulation(optimizer, formulation::Symbol, n::Int, data::InstanceData; objective_kind::Symbol)
    model = Model(optimizer)
    set_silent(model)

    if formulation == :matrix_dense
        @variable(model, X[1:n, 1:n], PSD)
        @constraint(model, sum(X[i, i] for i in 1:n) == 1.0)
        if objective_kind == :linear
            @objective(model, Min, sum(data.C[i, j] * X[i, j] for i in 1:n for j in i:n))
        else
            @objective(model, Min, 0.0)
        end

    elseif formulation == :matrix_zeroeq
        @variable(model, X[1:n, 1:n], PSD)
        @constraint(model, sum(X[i, i] for i in 1:n) == 1.0)
        for (i, j) in data.zero_pairs
            @constraint(model, X[i, j] == 0.0)
        end
        if objective_kind == :linear
            @objective(model, Min, sum(data.C[i, j] * X[i, j] for i in 1:n for j in i:n))
        else
            @objective(model, Min, 0.0)
        end

    elseif formulation == :scalar_dense
        d, z, pair_index = scalar_model_full(model, n, data.pairs)
        @constraint(model, sum(d) == 1.0)
        if objective_kind == :linear
            @objective(
                model,
                Min,
                sum(data.C[i, i] * d[i] for i in 1:n) +
                sum(data.C[i, j] * z[pair_index[(i, j)]] for (i, j) in data.pairs),
            )
        else
            @objective(model, Min, 0.0)
        end

    elseif formulation == :scalar_zeroeq
        d, z, pair_index = scalar_model_full(model, n, data.pairs)
        @constraint(model, sum(d) == 1.0)
        for (i, j) in data.zero_pairs
            @constraint(model, z[pair_index[(i, j)]] == 0.0)
        end
        if objective_kind == :linear
            @objective(
                model,
                Min,
                sum(data.C[i, i] * d[i] for i in 1:n) +
                sum(data.C[i, j] * z[pair_index[(i, j)]] for (i, j) in data.pairs),
            )
        else
            @objective(model, Min, 0.0)
        end

    elseif formulation == :scalar_structural_zero
        free_pairs = [(i, j) for (i, j) in data.pairs if !((i, j) in data.zero_set)]
        d, z = scalar_model_sparse(model, n, free_pairs)
        @constraint(model, sum(d) == 1.0)
        if objective_kind == :linear
            @objective(
                model,
                Min,
                sum(data.C[i, i] * d[i] for i in 1:n) +
                sum(data.C[i, j] * z[idx] for (idx, (i, j)) in enumerate(free_pairs)),
            )
        else
            @objective(model, Min, 0.0)
        end

    else
        error("Unknown formulation: $(formulation)")
    end

    optimize!(model)
    return (status=termination_status(model), solve_time=solve_time(model))
end

function benchmark_case(name::String, optimizer, n::Int, zero_fraction::Float64, objective_kind::Symbol)
    formulations = (
        :matrix_dense,
        :matrix_zeroeq,
        :scalar_dense,
        :scalar_zeroeq,
        :scalar_structural_zero,
    )

    rows = NamedTuple[]
    for formulation in formulations
        times = Float64[]
        statuses = String[]
        for seed in (11, 22, 33)
            data = make_instance(n, zero_fraction, seed)
            result = solve_formulation(optimizer, formulation, n, data; objective_kind=objective_kind)
            push!(times, result.solve_time)
            push!(statuses, string(result.status))
        end
        push!(rows, (
            solver=name,
            n=n,
            zero_fraction=zero_fraction,
            objective_kind=objective_kind,
            formulation=String(formulation),
            median_time=median(times),
            times=copy(times),
            statuses=copy(statuses),
        ))
    end
    return rows
end

function warmup!(optimizer)
    warm = make_instance(8, 0.5, 1)
    for formulation in (:matrix_dense, :matrix_zeroeq, :scalar_dense, :scalar_zeroeq, :scalar_structural_zero)
        solve_formulation(optimizer, formulation, 8, warm; objective_kind=:linear)
    end
    return nothing
end

function markdown_report(rows)
    io = IOBuffer()
    println(io, "# Scalar-vs-matrix PSD benchmark")
    println(io)
    println(io, "This report compares direct PSD matrix variables against scalar-built PSD matrices.")
    println(io)
    println(io, "Modeling styles:")
    println(io, "- `matrix_dense`: `@variable(model, X[1:n,1:n], PSD)`")
    println(io, "- `matrix_zeroeq`: same, plus many `X[i,j] == 0`")
    println(io, "- `scalar_dense`: scalar variables assembled into a symmetric matrix, then `M in PSDCone()`")
    println(io, "- `scalar_zeroeq`: same scalar-built matrix, plus many scalar `== 0` constraints")
    println(io, "- `scalar_structural_zero`: only create scalars for allowed nonzeros, fill forbidden entries with literal zeros, then `M in PSDCone()`")
    println(io)

    cases = unique((row.solver, row.n, row.zero_fraction, row.objective_kind) for row in rows)
    for (solver, n, zero_fraction, objective_kind) in cases
        case_rows = filter(row -> row.solver == solver && row.n == n && row.zero_fraction == zero_fraction && row.objective_kind == objective_kind, rows)
        println(io, "## $(solver), n=$(n), zero_fraction=$(zero_fraction), objective=$(objective_kind)")
        println(io)
        println(io, "| formulation | median solve time (s) | statuses | raw times |")
        println(io, "| :--- | ---: | :--- | :--- |")
        for row in case_rows
            println(
                io,
                "| $(row.formulation) | $(@sprintf("%.4f", row.median_time)) | $(join(unique(row.statuses), ", ")) | `$(row.times)` |",
            )
        end
        println(io)
    end

    println(io, "## Takeaway")
    println(io)
    println(io, "Declaring scalar variables first and then setting many of them to zero does **not** buy speed by itself. In these tests it behaves about the same as the matrix-variable version for COSMO, and it is often worse for Mosek. The main reason is unchanged: the PSD cone dimension is still the same, and explicit zero equalities just add more linear constraints.")
    println(io)
    println(io, "Even the `scalar_structural_zero` version does not magically fix this. It removes some scalar variables, but the PSD cone is still an `n × n` cone, so the win is limited unless you reformulate the problem to expose smaller cones or chordal structure.")
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
        benchmark_case("COSMO", cosmo, 40, 0.8, :linear),
        benchmark_case("COSMO", cosmo, 80, 0.8, :feasibility),
        benchmark_case("MosekTools", mosek, 40, 0.8, :linear),
        benchmark_case("MosekTools", mosek, 80, 0.8, :feasibility),
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
