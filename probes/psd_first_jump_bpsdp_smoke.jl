#!/usr/bin/env julia
# Smoke-test the formulation issue analyzed in MOMENT_SOS_PIPELINE_ANALYSIS.md
# and addressed by LOWERING_REFACTOR_PLAN.md.
#
# This does not implement the MomentProblem lowering. It proves the key BPSDP.jl
# interface fact: a JuMP PSD-first model reaches BPSDP without free-scalar ->
# 1x1-PSD block explosion, while an equivalent variable-first affine-PSD model
# does create scalar cone clutter.

using JuMP
using JSON3
using Printf
using SparseArrays

try
    @eval import BPSDP
catch err
    error("BPSDP.jl is required. Run with Pkg.develop(path=\"/home/ubuntu/BPSDP.jl\"). Original error: $(err)")
end

const MOI = JuMP.MOI
const DEFAULT_OUTPUT = normpath(joinpath(@__DIR__, "..", "output", "phase2", "psd_first_jump_bpsdp_smoke.json"))

function jump_optimizer(; kwargs...)
    return BPSDP.Optimizer(
        ;
        max_iter = 5_000,
        cg_max_iter = 100,
        mu_update_frequency = 25,
        penalty_parameter = 0.1,
        cg_convergence = 1e-12,
        dynamic_cg_convergence = false,
        sdp_objective_convergence = 1e-8,
        sdp_error_convergence = 1e-8,
        guess_type = :zero,
        print_level = 0,
        kwargs...,
    )
end

function inspect_backend(model::JuMP.Model)
    raw = JuMP.unsafe_backend(model)
    block_dims = copy(raw.block_dims)
    block_kinds = string.(raw.block_kinds)
    scalar_blocks = count(==(1), block_dims)
    nonscalar_blocks = count(>(1), block_dims)
    return Dict{String,Any}(
        "termination_status" => string(termination_status(model)),
        "objective" => objective_value(model),
        "jump_num_variables" => JuMP.num_variables(model),
        "jump_num_constraints_total" => JuMP.num_constraints(model; count_variable_in_set_constraints = true),
        "bpsdp_block_dims" => block_dims,
        "bpsdp_block_kinds" => block_kinds,
        "bpsdp_n_blocks" => length(block_dims),
        "bpsdp_scalar_1x1_blocks" => scalar_blocks,
        "bpsdp_nonscalar_blocks" => nonscalar_blocks,
        "bpsdp_reals_count" => raw.reals_count,
        "bpsdp_A_shape" => [size(raw.A, 1), size(raw.A, 2)],
        "bpsdp_A_nnz" => nnz(raw.A),
        "bpsdp_variable_map_length" => length(raw.variable_map),
    )
end

function populate_variable_first!(model)
    @variable(model, y[1:3])
    @constraint(model, [y[1] y[2]; y[2] y[3]] in PSDCone())
    @constraint(model, y[2] == 1.0)
    @objective(model, Min, y[1] + y[3])
    return model
end

function run_variable_first()
    model = Model(jump_optimizer)
    set_silent(model)
    populate_variable_first!(model)
    optimize!(model)
    return inspect_backend(model)
end

function run_late_optimizer_variable_first()
    model = JuMP.GenericModel{Float64}()
    populate_variable_first!(model)
    set_optimizer(model, jump_optimizer)
    set_silent(model)
    optimize!(model)
    return inspect_backend(model)
end

function run_psd_first()
    model = Model(jump_optimizer)
    set_silent(model)

    @variable(model, X[1:2, 1:2], PSD)
    @constraint(model, X[1, 2] == 1.0)
    @objective(model, Min, X[1, 1] + X[2, 2])

    optimize!(model)
    return inspect_backend(model)
end

function main(argv = ARGS)
    output = isempty(argv) ? DEFAULT_OUTPUT : argv[1]
    variable_first = run_variable_first()
    late_variable_first = run_late_optimizer_variable_first()
    psd_first = run_psd_first()

    result = Dict{String,Any}(
        "variable_first" => variable_first,
        "late_optimizer_variable_first" => late_variable_first,
        "psd_first" => psd_first,
        "pass" => Dict{String,Any}(
            "objectives_match" => abs(variable_first["objective"] - psd_first["objective"]) <= 1e-5 && abs(late_variable_first["objective"] - psd_first["objective"]) <= 1e-5,
            "variable_first_has_no_scalar_clutter" => variable_first["bpsdp_scalar_1x1_blocks"] == 0,
            "late_optimizer_variable_first_has_no_scalar_clutter" => late_variable_first["bpsdp_scalar_1x1_blocks"] == 0,
            "psd_first_has_no_scalar_clutter" => psd_first["bpsdp_scalar_1x1_blocks"] == 0,
            "all_one_block" => variable_first["bpsdp_block_dims"] == [2] && late_variable_first["bpsdp_block_dims"] == [2] && psd_first["bpsdp_block_dims"] == [2],
        ),
    )

    mkpath(dirname(output))
    open(output, "w") do io
        JSON3.write(io, result)
        write(io, '\n')
    end

    println("== BPSDP JuMP formulation smoke ==")
    @printf("%-28s %-14s %-10s %-14s %-20s %-12s\n", "formulation", "status", "objective", "#blocks", "block_dims", "1x1 blocks")
    for name in ("variable_first", "late_optimizer_variable_first", "psd_first")
        r = result[name]
        @printf("%-28s %-14s %-10.6f %-14d %-20s %-12d\n",
            name,
            r["termination_status"],
            r["objective"],
            r["bpsdp_n_blocks"],
            string(r["bpsdp_block_dims"]),
            r["bpsdp_scalar_1x1_blocks"],
        )
    end
    println("wrote ", output)

    all(values(result["pass"])) || error("PSD-first smoke test failed: $(result["pass"])")
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
