#!/usr/bin/env julia
# Inspect what the current NCTSSoS variable-first H2/Nk=2 JuMP lowering becomes
# after JuMP/MOI bridges and BPSDP.copy_to. Uses BPSDP max_iter=0: copy/build
# happens, but no real solve is attempted.

using JSON3
using JuMP
using SparseArrays

include(joinpath(@__DIR__, "..", "demos", "h2_periodic_nk2_moment_sos.jl"))

function write_json(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, obj)
        write(io, '\n')
    end
end

function main(argv = ARGS)
    output = normpath(joinpath(@__DIR__, "..", "output", "phase2", "h2_nk2_bpsdp_bridge_inspect.json"))
    formulation_args = String[
        "--blocking=spin",
        "--paper-spin",
        "--output-dir=output/phase2/h2_nk2_bridge_inspect",
        "--bpsdp-max-iter=0",
        "--bpsdp-print-level=0",
    ]
    for arg in argv
        if startswith(arg, "--output-json=")
            output = split(arg, "=", limit = 2)[2]
        else
            push!(formulation_args, arg)
        end
    end

    options = parse_options(formulation_args)
    build_seconds = @elapsed data = build_h2_pqg_moment_problem(options)
    jump_seconds = @elapsed jump_data = build_jump_model(data.moment_problem, data.block_labels)

    set_optimizer(jump_data.model, () -> BPSDP.Optimizer(
        max_iter = 0,
        cg_max_iter = 1,
        mu_update_frequency = 25,
        penalty_parameter = 0.1,
        cg_convergence = 1e-12,
        dynamic_cg_convergence = false,
        sdp_objective_convergence = 1e-8,
        sdp_error_convergence = 1e-8,
        guess_type = :zero,
        print_level = 0,
    ))
    set_silent(jump_data.model)

    optimize_error = nothing
    copy_seconds = @elapsed begin
        try
            optimize!(jump_data.model)
        catch err
            optimize_error = sprint(showerror, err)
        end
    end

    raw = JuMP.unsafe_backend(jump_data.model)
    block_dims = copy(raw.block_dims)
    scalar_blocks = count(==(1), block_dims)
    non_scalar_blocks = count(>(1), block_dims)
    large_dims = sort(filter(>(1), block_dims); rev = true)
    dim_counts = Dict{String,Int}()
    for d in block_dims
        key = string(d)
        dim_counts[key] = get(dim_counts, key, 0) + 1
    end

    result = Dict{String,Any}(
        "build_seconds" => build_seconds,
        "jump_build_seconds" => jump_seconds,
        "copy_optimize_seconds_max_iter_0" => copy_seconds,
        "optimize_error" => optimize_error,
        "termination_status" => try string(termination_status(jump_data.model)) catch err string("unavailable: ", err) end,
        "jump_num_variables" => JuMP.num_variables(jump_data.model),
        "jump_num_constraints_total" => JuMP.num_constraints(jump_data.model; count_variable_in_set_constraints = true),
        "moment_basis_length_complex" => length(jump_data.basis),
        "moment_real_variables_expected" => 2 * length(jump_data.basis),
        "hpsd_block_labels" => data.block_labels,
        "hpsd_block_dims_complex" => [item.dim for item in jump_data.hpsd_refs],
        "bpsdp_n_blocks" => length(block_dims),
        "bpsdp_scalar_1x1_blocks" => scalar_blocks,
        "bpsdp_non_scalar_blocks" => non_scalar_blocks,
        "bpsdp_largest_non_scalar_dims" => large_dims[1:min(end, 25)],
        "bpsdp_dim_counts" => dim_counts,
        "bpsdp_reals_count" => raw.reals_count,
        "bpsdp_A_shape" => [size(raw.A, 1), size(raw.A, 2)],
        "bpsdp_A_nnz" => nnz(raw.A),
        "bpsdp_variable_map_length" => length(raw.variable_map),
        "raw_status" => try raw.raw_status catch err string("unavailable: ", err) end,
    )

    write_json(output, result)
    println("wrote ", output)
    println("BPSDP blocks: ", length(block_dims), " total; ", scalar_blocks, " scalar 1x1; ", non_scalar_blocks, " non-scalar")
    println("A shape: ", result["bpsdp_A_shape"], " nnz=", result["bpsdp_A_nnz"])
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
