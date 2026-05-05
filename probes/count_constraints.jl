#!/usr/bin/env julia
# Reconstruct the H4/Nk=2 paper-spin + 1D MomentProblem and count JuMP
# constraints WITHOUT solving the SDP.
#
# This script does NOT modify anything in demos/. It reuses
# `Options` and `build_h4_pqg_moment_problem` from
# demos/h4_periodic_moment_sos.jl, then mirrors the JuMP-build path inside
# NCTSSoS' `_solve_complex_moment_problem` (src/optimization/moment.jl) up to
# but not including `optimize!`. After the model is built it reports
#
#   * num_variables(model)
#   * num_constraints(model; count_variable_in_set_constraints = true)
#   * one row per (FunctionType, SetType) pair from list_of_constraint_types
#
# Usage from the repository root:
#   julia --project=demos probes/count_constraints.jl \
#       --integrals=output/h4_chain_nk2_figure1_integrals_drop1e-12.txt
#
# Default flags match the paper-faithful run that reproduces Schouten 2025
# (--paper-spin --include-1d).

using Printf
using NCTSSoS
using JuMP

include(joinpath(@__DIR__, "..", "demos", "h4_periodic_moment_sos.jl"))

# -----------------------------------------------------------------------------
# Argument parsing — keep small; default to paper-spin + include-1d.
# -----------------------------------------------------------------------------

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_INTEGRALS_PATH = normpath(joinpath(
    REPO_ROOT, "output", "h4_chain_nk2_figure1_integrals_drop1e-12.txt"
))

function parse_probe_options(argv)
    integrals_path = DEFAULT_INTEGRALS_PATH
    blocking = :momentum
    include_one_d = true
    spin_resolved_trace = true
    singlet_s2 = true

    for arg in argv
        if startswith(arg, "--integrals=")
            integrals_path = String(split(arg, "=", limit = 2)[2])
        elseif startswith(arg, "--blocking=")
            value = split(arg, "=", limit = 2)[2]
            if value in ("momentum", "k")
                blocking = :momentum
            elseif value == "spin"
                blocking = :spin
            elseif value == "none"
                blocking = :none
            else
                throw(ArgumentError("unknown --blocking=$value"))
            end
        elseif arg == "--no-1d"
            include_one_d = false
        elseif arg == "--include-1d"
            include_one_d = true
        elseif arg == "--no-paper-spin"
            spin_resolved_trace = false
            singlet_s2 = false
        elseif arg in ("--paper-spin", "--spin-singlet")
            spin_resolved_trace = true
            singlet_s2 = true
        elseif arg in ("-h", "--help")
            println("Usage: julia --project=demos probes/count_constraints.jl [options]")
            println("Options (defaults match the paper-faithful run):")
            println("  --integrals=PATH         text integral dump (default: figure-1 drop1e-12)")
            println("  --blocking=momentum|spin|none")
            println("  --include-1d / --no-1d   default: --include-1d")
            println("  --paper-spin / --no-paper-spin   default: --paper-spin")
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    return Options(integrals_path, blocking, include_one_d,
                   spin_resolved_trace, singlet_s2)
end

# -----------------------------------------------------------------------------
# JuMP build (mirrors NCTSSoS._solve_complex_moment_problem, no solve).
# -----------------------------------------------------------------------------

function build_jump_no_solve(mp)
    C = real(eltype(coefficients(mp.objective)))
    model = JuMP.GenericModel{C}()

    NM = eltype(mp.total_basis)
    basis = [symmetric_canon(NCTSSoS.expval(m)) for m in mp.total_basis]
    sort!(unique!(basis))                    # NCTSSoS.sorted_unique! inlined

    n_basis = length(basis)
    JuMP.@variable(model, y_re[1:n_basis], set_string_name = false)
    JuMP.@variable(model, y_im[1:n_basis], set_string_name = false)

    one_sym = symmetric_canon(NCTSSoS.expval(one(NM)))
    idx_one = findfirst(==(one_sym), basis)
    idx_one === nothing && error("identity moment missing from basis")
    JuMP.@constraint(model, y_re[idx_one] == 1)
    JuMP.@constraint(model, y_im[idx_one] == 0)

    basis_to_idx = Dict(m => i for (i, m) in enumerate(basis))

    for (cone, mat) in mp.constraints
        dim = size(mat, 1)
        mat_re = Matrix{Any}(undef, dim, dim)
        mat_im = Matrix{Any}(undef, dim, dim)
        for i in 1:dim, j in 1:dim
            re_expr, im_expr = NCTSSoS._substitute_complex_poly(
                mat[i, j], basis_to_idx, y_re, y_im,
            )
            mat_re[i, j] = re_expr
            mat_im[i, j] = im_expr
        end

        if cone == :Zero
            JuMP.@constraint(model, [mat_re[i, j] for i in 1:dim, j in 1:dim] .== 0)
            JuMP.@constraint(model, [mat_im[i, j] for i in 1:dim, j in 1:dim] .== 0)
        elseif cone == :HPSD
            embedded = [
                [mat_re[i, j] for i in 1:dim, j in 1:dim] [-mat_im[i, j] for i in 1:dim, j in 1:dim]
                [mat_im[i, j] for i in 1:dim, j in 1:dim] [mat_re[i, j] for i in 1:dim, j in 1:dim]
            ]
            JuMP.@constraint(model, embedded in PSDCone())
        else
            error("unexpected cone $cone (probe only handles Zero / HPSD)")
        end
    end

    obj_re, _ = NCTSSoS._substitute_complex_poly(
        mp.objective, basis_to_idx, y_re, y_im,
    )
    JuMP.@objective(model, Min, obj_re)

    return model
end

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

function main(argv = ARGS)
    options = parse_probe_options(argv)

    println("== JuMP constraint count probe ==")
    @printf("%-40s %s\n", "integrals", options.integrals_path)
    @printf("%-40s %s\n", "blocking", string(options.blocking))
    @printf("%-40s %s\n", "include 1D", string(options.include_one_d))
    @printf("%-40s %s\n", "spin-resolved D-traces", string(options.spin_resolved_trace))
    @printf("%-40s %s\n", "singlet S^2", string(options.singlet_s2))
    println()

    build_seconds = @elapsed data = build_h4_pqg_moment_problem(options)
    @printf("%-40s %.3f s\n", "MomentProblem build wall", build_seconds)

    mp = data.moment_problem
    n_hpsd = count(cone == :HPSD for (cone, _) in mp.constraints)
    n_zero = count(cone == :Zero for (cone, _) in mp.constraints)
    @printf("%-40s %d\n", "HPSD blocks (mp.constraints)", n_hpsd)
    @printf("%-40s %d\n", "Zero blocks (mp.constraints)", n_zero)
    @printf("%-40s %d\n", "unique moment-matrix monomials", mp.n_unique_moment_matrix_elements)
    println()

    jump_seconds = @elapsed model = build_jump_no_solve(mp)
    @printf("%-40s %.3f s\n", "JuMP model build wall", jump_seconds)
    println()

    println("== JuMP statistics ==")
    @printf("%-40s %d\n", "JuMP.num_variables", JuMP.num_variables(model))
    total = JuMP.num_constraints(model; count_variable_in_set_constraints = true)
    @printf("%-40s %d\n", "JuMP.num_constraints (total)", total)
    println()

    println("Constraint breakdown by (FunctionType, SetType):")
    for (F, S) in JuMP.list_of_constraint_types(model)
        n = JuMP.num_constraints(model, F, S)
        @printf("  %-60s %d\n", string(F, " in ", S), n)
    end

    return model
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
