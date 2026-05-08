#!/usr/bin/env julia
raw"""
H4/Nk=2 periodic PQG V2RDM sector-pruned BPSDP run.

This is a demos-only post-processing fix for the symmetry-leaking orphan
moments created by one-sided moment-equality rows. It builds the frozen H4
MomentProblem, lowers through JuMP with scalar moment variables and
`orphan_policy=:free_variables` (never `:aux_psd_free`), optionally pins
symmetry-forbidden orphan moments to zero, deletes resulting tautological
linear equalities, solves with BPSDP.jl, and writes an audit trail.

Usage from the repository root on HAI:

    easy-ssh run 'cd /home/ubuntu/NCTSSoS.jl-h4-periodic-v2rdm-benchmark
    tmpdir=$(mktemp -d)
    julia --startup-file=no --project="$tmpdir" -e "using Pkg; \
      Pkg.develop(path=\"/home/ubuntu/BPSDP.jl\"); \
      Pkg.develop(path=pwd()); \
      Pkg.add([\"JuMP\",\"JSON3\",\"Printf\",\"Dates\",\"LinearAlgebra\",\"SparseArrays\"]); \
      Pkg.instantiate()"
    julia --startup-file=no --project="$tmpdir" demos/h4_periodic_moment_sos_sector_pruned.jl \
        --pruning=both --output-dir=output/h4_pruning'

The integral text format is inherited from `demos/h4_periodic_moment_sos.jl`:

    h1e k p q Re Im
    eri k1 k2 k3 k4 p r q s Re Im

with ERIs in chemist order (p r | q s). The frozen H4 builder normalizes the
one-electron part by 1/Nk and ERIs by 1/Nk^2 before building the per-cell
Hamiltonian.
"""

using NCTSSoS
using JuMP
const MOI = JuMP.MOI
using LinearAlgebra
using SparseArrays
using Printf
using Dates
using JSON3

try
    import BPSDP
catch err
    error("BPSDP.jl is mandatory for this pruning solve. Run on HAI with Pkg.develop(path=\"/home/ubuntu/BPSDP.jl\") as documented in PHASE2_PLAN.md §12.6. Original load error: $(err)")
end

module H4PeriodicBase
include(joinpath(@__DIR__, "h4_periodic_moment_sos.jl"))
end

const DEFAULT_OUTPUT_ROOT = normpath(joinpath(@__DIR__, "..", "output", "h4_pruning"))
const DEFAULT_TARGET_SECTOR = "N=0,K=0,Sz=0,P=even"

struct TargetSector
    delta_N::Int
    two_Sz::Int
    K::Int
    parity::Int
end

struct SectorPruningOptions
    integrals_path::String
    blocking::Symbol
    include_one_d::Bool
    spin_resolved_trace::Bool
    singlet_s2::Bool
    target_sector::TargetSector
    target_sector_raw::String
    pruning::Symbol
    lowering::Symbol
    output_dir::String
    bpsdp_max_iter::Int
    bpsdp_cg_max_iter::Int
    bpsdp_mu_update_frequency::Int
    bpsdp_penalty::Float64
    bpsdp_cg_tol::Float64
    bpsdp_obj_tol::Float64
    bpsdp_err_tol::Float64
    bpsdp_print_level::Int
end

bytes_to_gib(bytes::Integer) = bytes / 1024.0^3
rss_string() = @sprintf("%.3f GiB", bytes_to_gib(Sys.maxrss()))

function _parse_spin_value(raw::AbstractString)
    value = strip(raw)
    if occursin('/', value)
        parts = split(value, '/'; limit = 2)
        length(parts) == 2 || throw(ArgumentError("bad Sz value $(repr(raw))"))
        return parse(Float64, strip(parts[1])) / parse(Float64, strip(parts[2]))
    end
    return parse(Float64, value)
end

function _parse_parity(raw::AbstractString)
    value = lowercase(strip(raw))
    if value in ("even", "e", "+", "+1", "0")
        return 0
    elseif value in ("odd", "o", "-", "-1", "1")
        return 1
    end
    throw(ArgumentError("bad parity $(repr(raw)); expected even/odd"))
end

function parse_target_sector(raw::AbstractString)
    delta_N = 0
    two_Sz = 0
    K = 0
    parity = 0
    for piece in split(raw, ',')
        isempty(strip(piece)) && continue
        kv = split(piece, '='; limit = 2)
        length(kv) == 2 || throw(ArgumentError("bad --target-sector component $(repr(piece)); expected key=value"))
        key = lowercase(strip(kv[1]))
        value = strip(kv[2])
        if key in ("n", "dn", "deltan", "delta_n", "δn")
            delta_N = parse(Int, value)
        elseif key == "k"
            K = parse(Int, value)
        elseif key in ("2sz", "twosz", "two_sz", "two-sz")
            two_Sz = parse(Int, value)
        elseif key in ("sz", "s_z")
            two_Sz = round(Int, 2 * _parse_spin_value(value))
        elseif key in ("p", "parity")
            parity = _parse_parity(value)
        else
            throw(ArgumentError("unknown target-sector key $(repr(key)); expected N,K,Sz,2Sz,P"))
        end
    end
    return TargetSector(delta_N, two_Sz, K, parity)
end

function target_sector_string(s::TargetSector)
    parity = s.parity == 0 ? "even" : "odd"
    return "N=$(s.delta_N),K=$(s.K),2Sz=$(s.two_Sz),P=$(parity)"
end

function parse_pruning(raw::AbstractString)
    value = lowercase(strip(raw))
    value in ("on", "true", "yes", "1", "pruned") && return :on
    value in ("off", "false", "no", "0", "unpruned") && return :off
    value == "both" && return :both
    throw(ArgumentError("unknown --pruning=$raw; expected on, off, or both"))
end

function _parse_blocking(raw::AbstractString)
    value = lowercase(strip(raw))
    if value in ("momentum", "k")
        return :momentum
    elseif value == "spin"
        return :spin
    elseif value == "none"
        return :none
    end
    throw(ArgumentError("unknown --blocking=$raw; expected momentum, spin, or none"))
end

function print_help()
    println("Usage: julia --project=<tmp-with-BPSDP> demos/h4_periodic_moment_sos_sector_pruned.jl [options]\n")
    println("Options:")
    println("  --integrals=PATH                  H4 text integral dump; default inherited from frozen H4 demo")
    println("  --blocking=momentum|spin|none     1D/PQG block grouping; default momentum")
    println("  --include-1d / --no-1d            include/drop extra ¹D PSD block; default --no-1d")
    println("  --paper-spin, --spin-singlet      spin-resolved ²D traces plus singlet S²")
    println("  --spin-resolved-trace / --no-spin-resolved-trace")
    println("  --singlet-s2 / --no-singlet-s2")
    println("  --target-sector=N=0,K=0,Sz=0,P=even   default closed-shell expectation sector")
    println("  --pruning=on|off|both             default on; both runs unpruned then pruned")
    println("  --lowering=moment_variables|psd_blocks  default moment_variables; psd_blocks is the source-true free-orphan path")
    println("  --output-dir=PATH                 output root; default output/h4_pruning")
    println("  --bpsdp-max-iter=N                default 5000")
    println("  --bpsdp-cg-max-iter=N             default 100")
    println("  --bpsdp-mu-update-frequency=N     default 25")
    println("  --bpsdp-penalty=rho               default 0.1")
    println("  --bpsdp-cg-tol=eps                default 1e-12")
    println("  --bpsdp-obj-tol=eps               default 1e-8")
    println("  --bpsdp-err-tol=eps               default 1e-8")
    println("  --bpsdp-print-level=N             default 1")
    return nothing
end

function parse_sector_pruning_options(argv)
    integrals_path = H4PeriodicBase.DEFAULT_INTEGRALS
    blocking = :momentum
    include_one_d = false
    spin_resolved_trace = false
    singlet_s2 = false
    target_sector_raw = DEFAULT_TARGET_SECTOR
    pruning = :on
    lowering = :moment_variables
    output_dir = DEFAULT_OUTPUT_ROOT
    bpsdp_max_iter = 5_000
    bpsdp_cg_max_iter = 100
    bpsdp_mu_update_frequency = 25
    bpsdp_penalty = 0.1
    bpsdp_cg_tol = 1e-12
    bpsdp_obj_tol = 1e-8
    bpsdp_err_tol = 1e-8
    bpsdp_print_level = 1

    for arg in argv
        if startswith(arg, "--integrals=")
            integrals_path = split(arg, "="; limit = 2)[2]
        elseif startswith(arg, "--blocking=")
            blocking = _parse_blocking(split(arg, "="; limit = 2)[2])
        elseif arg == "--include-1d"
            include_one_d = true
        elseif arg == "--no-1d"
            include_one_d = false
        elseif arg in ("--paper-spin", "--spin-singlet")
            spin_resolved_trace = true
            singlet_s2 = true
        elseif arg == "--spin-resolved-trace"
            spin_resolved_trace = true
        elseif arg == "--no-spin-resolved-trace"
            spin_resolved_trace = false
        elseif arg == "--singlet-s2"
            singlet_s2 = true
        elseif arg == "--no-singlet-s2"
            singlet_s2 = false
        elseif startswith(arg, "--target-sector=")
            target_sector_raw = split(arg, "="; limit = 2)[2]
        elseif startswith(arg, "--pruning=")
            pruning = parse_pruning(split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--lowering=")
            raw = lowercase(strip(split(arg, "="; limit = 2)[2]))
            if raw in ("moment_variables", "moment", "moments")
                lowering = :moment_variables
            elseif raw in ("psd_blocks", "psd", "blocks")
                lowering = :psd_blocks
            else
                throw(ArgumentError("unknown --lowering=$raw; expected moment_variables or psd_blocks"))
            end
        elseif startswith(arg, "--output-dir=")
            output_dir = split(arg, "="; limit = 2)[2]
        elseif startswith(arg, "--bpsdp-max-iter=")
            bpsdp_max_iter = parse(Int, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--bpsdp-cg-max-iter=")
            bpsdp_cg_max_iter = parse(Int, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--bpsdp-mu-update-frequency=")
            bpsdp_mu_update_frequency = parse(Int, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--bpsdp-penalty=")
            bpsdp_penalty = parse(Float64, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--bpsdp-cg-tol=")
            bpsdp_cg_tol = parse(Float64, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--bpsdp-obj-tol=")
            bpsdp_obj_tol = parse(Float64, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--bpsdp-err-tol=")
            bpsdp_err_tol = parse(Float64, split(arg, "="; limit = 2)[2])
        elseif startswith(arg, "--bpsdp-print-level=")
            bpsdp_print_level = parse(Int, split(arg, "="; limit = 2)[2])
        elseif arg in ("-h", "--help")
            print_help()
            exit(0)
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
    end

    target_sector = parse_target_sector(target_sector_raw)
    return SectorPruningOptions(
        integrals_path,
        blocking,
        include_one_d,
        spin_resolved_trace,
        singlet_s2,
        target_sector,
        target_sector_raw,
        pruning,
        lowering,
        output_dir,
        bpsdp_max_iter,
        bpsdp_cg_max_iter,
        bpsdp_mu_update_frequency,
        bpsdp_penalty,
        bpsdp_cg_tol,
        bpsdp_obj_tol,
        bpsdp_err_tol,
        bpsdp_print_level,
    )
end

function h4_base_options(options::SectorPruningOptions)
    return H4PeriodicBase.Options(
        options.integrals_path,
        options.blocking,
        options.include_one_d,
        options.spin_resolved_trace,
        options.singlet_s2,
    )
end

function run_output_dir(root::AbstractString, nk::Integer)
    base = basename(normpath(root))
    startswith(base, "h4_nk") && return normpath(root)
    return normpath(joinpath(root, "h4_nk$(nk)"))
end

function write_json(path::AbstractString, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, obj)
        write(io, '\n')
    end
    return path
end

finite_or_string(x::Real) = isfinite(Float64(x)) ? Float64(x) : string(Float64(x))
finite_or_string(x) = x
finite_or_nothing(x) = x === nothing ? nothing : finite_or_string(x)

has_field(obj, name::Symbol) = name in fieldnames(typeof(obj))
field_or_nothing(obj, name::Symbol) = has_field(obj, name) ? getfield(obj, name) : nothing
function first_field_or_nothing(obj, names::Symbol...)
    for name in names
        has_field(obj, name) && return getfield(obj, name)
    end
    return nothing
end

function raw_optimizer(model::JuMP.Model)
    try
        return JuMP.unsafe_backend(model)
    catch
        return nothing
    end
end

function bpsdp_state_summary(raw)
    state = raw === nothing ? nothing : field_or_nothing(raw, :state)
    state === nothing && return nothing
    return Dict{String,Any}(
        "outer_iterations" => first_field_or_nothing(state, :outer_iterations, :iterations),
        "inner_iterations" => field_or_nothing(state, :inner_iterations),
        "primal_error" => finite_or_nothing(field_or_nothing(state, :primal_error)),
        "dual_error" => finite_or_nothing(field_or_nothing(state, :dual_error)),
        "objective_primal" => finite_or_nothing(first_field_or_nothing(state, :objective_primal, :primal_objective)),
        "objective_dual" => finite_or_nothing(first_field_or_nothing(state, :objective_dual, :dual_objective)),
        "objective_gap" => finite_or_nothing(field_or_nothing(state, :objective_gap)),
        "termination_reason" => string(first_field_or_nothing(state, :termination_reason, :status)),
    )
end

function raw_problem_summary(raw)
    raw === nothing && return nothing
    summary = Dict{String,Any}()
    if has_field(raw, :block_dims)
        dims = getfield(raw, :block_dims)
        summary["n_blocks"] = length(dims)
        summary["largest_block_dim"] = isempty(dims) ? 0 : maximum(dims)
        summary["block_dims_sample"] = collect(dims[1:min(end, 20)])
    end
    has_field(raw, :b) && (summary["dual_rows"] = length(getfield(raw, :b)))
    has_field(raw, :c) && (summary["primal_dim"] = length(getfield(raw, :c)))
    has_field(raw, :A) && (summary["A_nnz"] = nnz(getfield(raw, :A)))
    return summary
end

function model_cone_summary(model::JuMP.Model)
    types = Vector{Dict{String,Any}}()
    for (F, S) in JuMP.list_of_constraint_types(model)
        count = length(JuMP.all_constraints(model, F, S))
        push!(types, Dict(
            "function_type" => string(F),
            "set_type" => string(S),
            "count" => count,
        ))
    end
    return Dict{String,Any}(
        "num_variables" => JuMP.num_variables(model),
        "num_constraints_excluding_variable_sets" => JuMP.num_constraints(model; count_variable_in_set_constraints = false),
        "num_constraints_including_variable_sets" => JuMP.num_constraints(model; count_variable_in_set_constraints = true),
        "constraint_types" => types,
    )
end

function signature_of_monomial(m; nk::Int, norb::Int)
    delta_N = 0
    K = 0
    two_Sz = 0
    for idx in m.word
        raw = Int(idx)
        mode = abs(raw)
        family, _ = divrem(mode - 1, norb)
        k = family ÷ 2
        spin_sign = iseven(family) ? 1 : -1
        particle_sign = raw < 0 ? 1 : -1  # creator -> +1, annihilator -> -1
        delta_N += particle_sign
        K = mod(K + particle_sign * k, nk)
        two_Sz += particle_sign * spin_sign
    end
    return (delta_N = delta_N, two_Sz = two_Sz, K = mod(K, nk), parity = mod(delta_N, 2))
end

function signature_string(sig)
    parity = sig.parity == 0 ? "even" : "odd"
    return "N=$(sig.delta_N),K=$(sig.K),2Sz=$(sig.two_Sz),P=$(parity)"
end

function signature_matches(sig, target::TargetSector; nk::Int)
    return sig.delta_N == target.delta_N &&
        sig.two_Sz == target.two_Sz &&
        mod(sig.K, nk) == mod(target.K, nk) &&
        sig.parity == target.parity
end

function classify_orphans(mp, target_sector::TargetSector; nk::Int, norb::Int)
    L = mp.linear
    expected_zero = eltype(L.free_keys)[]
    physical_unexpected = eltype(L.free_keys)[]
    key_signature = Dict{Any,String}()
    histogram = Dict{String,Int}()
    for key in L.free_keys
        mono = get(L.key_to_monomial, key, nothing)
        mono === nothing && error("free key has no representative monomial in mp.linear.key_to_monomial: $(sprint(show, key))")
        sig = signature_of_monomial(mono; nk, norb)
        sig_s = signature_string(sig)
        key_signature[key] = sig_s
        histogram[sig_s] = get(histogram, sig_s, 0) + 1
        if signature_matches(sig, target_sector; nk)
            push!(physical_unexpected, key)
        else
            push!(expected_zero, key)
        end
    end
    return (; expected_zero, physical_unexpected, key_signature, histogram)
end

function key_sample(keys; limit::Int = 50)
    n = min(length(keys), limit)
    return [sprint(show, keys[i]) for i in 1:n]
end

function classification_json(mp, classification, target_sector::TargetSector; nk::Int, norb::Int)
    return Dict{String,Any}(
        "target_sector" => target_sector_string(target_sector),
        "n_moments" => length(mp.linear.moments),
        "n_orphans" => length(mp.linear.free_keys),
        "n_expected_zero" => length(classification.expected_zero),
        "n_physical_unexpected" => length(classification.physical_unexpected),
        "n_unexpected" => length(classification.physical_unexpected),
        "signature_histogram" => classification.histogram,
        "expected_zero_sample" => key_sample(classification.expected_zero; limit = 50),
        "physical_unexpected_sample" => key_sample(classification.physical_unexpected; limit = 200),
        "note" => "expected_zero are symmetry-forbidden orphan keys pinned to zero; physical_unexpected are target-sector orphan keys left free and are the excess-kernel proxy.",
    )
end

function recover_orphan_handles(model::JuMP.Model, mp)
    L = mp.linear
    all_vars = JuMP.all_variables(model)
    n_moments = length(L.moments)
    n_orphans = length(L.free_keys)
    handles = Dict{Any,Tuple{JuMP.VariableRef,JuMP.VariableRef}}()
    spotcheck = Dict{String,Any}(
        "n_all_variables" => length(all_vars),
        "n_moments" => n_moments,
        "n_orphans" => n_orphans,
    )

    if n_orphans == 0
        spotcheck["layout"] = "none"
        return (; layout = :none, handles, spotcheck)
    end

    if length(all_vars) == 2 * n_moments
        # Current `formulation=:moment_variables, representation=:real` complex
        # lowering declares y_re[1:n_moments], then y_im[1:n_moments]. The
        # `orphan_policy` keyword is accepted by build_jump_model but this path
        # does not call `_declare_free_orphan_moments!`; free keys are scalar
        # moment variables indexed by L.moments. That's the source reality; the
        # plan's trailing-orphan assumption belongs to the PSD-block path.
        moment_index = Dict{Any,Int}(key => idx for (idx, key) in enumerate(L.moments))
        y_re = all_vars[1:n_moments]
        y_im = all_vars[(n_moments + 1):(2 * n_moments)]
        for key in L.free_keys
            idx = get(moment_index, key, 0)
            idx == 0 && error("free key not present in L.moments: $(sprint(show, key))")
            handles[key] = (y_re[idx], y_im[idx])
        end
        spotcheck["layout"] = "moment_variables_y_re_y_im"
        spotcheck["checked_keys"] = key_sample(L.free_keys; limit = min(3, n_orphans))
    elseif length(all_vars) >= 2 * n_orphans
        # PSD-block/free-orphan layout: all PSD variables first, then the two
        # arrays created by `_declare_free_orphan_moments!`. Kept here as a
        # guardrail if the formulation changes back to `:psd_blocks`.
        orphan_re = all_vars[(end - 2 * n_orphans + 1):(end - n_orphans)]
        orphan_im = all_vars[(end - n_orphans + 1):end]
        for (idx, key) in enumerate(L.free_keys)
            handles[key] = (orphan_re[idx], orphan_im[idx])
        end
        spotcheck["layout"] = "trailing_orphan_re_orphan_im"
        spotcheck["checked_keys"] = key_sample(L.free_keys; limit = min(3, n_orphans))
    else
        error("cannot recover orphan handles: JuMP has $(length(all_vars)) variables but $(n_orphans) orphan keys")
    end

    # Cheap pre-solve lock: every recovered handle must belong to this model and
    # the first few keys must map to distinct re/im scalar variables. The stronger
    # value-level check runs after optimize! when extract_monomap has values.
    seen = Set{JuMP.VariableRef}()
    for key in Iterators.take(L.free_keys, min(3, n_orphans))
        re_var, im_var = handles[key]
        JuMP.owner_model(re_var) === model || error("orphan real handle belongs to a different model")
        JuMP.owner_model(im_var) === model || error("orphan imag handle belongs to a different model")
        re_var == im_var && error("orphan real/imag handles alias for key $(sprint(show, key))")
        push!(seen, re_var)
        push!(seen, im_var)
    end
    length(seen) == 2 * min(3, n_orphans) || error("orphan handle order spot-check found duplicate variables")

    return (; layout = Symbol(spotcheck["layout"]), handles, spotcheck)
end

function verify_orphan_value_spotcheck!(extract_monomap, recovered, mp; limit::Int = 3, atol::Float64 = 1e-8)
    isempty(mp.linear.free_keys) && return Dict{String,Any}("checked" => 0)
    monomap = extract_monomap()
    checked = 0
    max_abs_diff = 0.0
    for key in Iterators.take(mp.linear.free_keys, min(limit, length(mp.linear.free_keys)))
        haskey(recovered.handles, key) || error("missing recovered handle for key $(sprint(show, key))")
        re_var, im_var = recovered.handles[key]
        got = ComplexF64(JuMP.value(re_var), JuMP.value(im_var))
        expected = ComplexF64(get(monomap, key, 0.0 + 0.0im))
        diff = abs(got - expected)
        max_abs_diff = max(max_abs_diff, diff)
        diff <= atol || error("orphan handle spot-check failed for $(sprint(show, key)): handle=$(got), monomap=$(expected), diff=$(diff)")
        checked += 1
    end
    return Dict{String,Any}("checked" => checked, "max_abs_diff" => max_abs_diff, "atol" => atol)
end

function _pruned_var_set(recovered, keys)
    vars = Set{JuMP.VariableRef}()
    for key in keys
        re_var, im_var = recovered.handles[key]
        push!(vars, re_var)
        push!(vars, im_var)
    end
    return vars
end

function _objective_scrub_pruned_terms!(model::JuMP.Model, pruned_vars::Set{JuMP.VariableRef}; atol::Float64)
    touched = 0
    removed_l1 = 0.0
    obj = JuMP.objective_function(model)
    if obj isa JuMP.GenericAffExpr
        for (coef, var) in collect(JuMP.linear_terms(obj))
            if var in pruned_vars && abs(coef) > atol
                JuMP.set_objective_coefficient(model, var, zero(coef))
                touched += 1
                removed_l1 += abs(coef)
            end
        end
    end
    return (; touched, removed_l1)
end

function _is_scalar_equal_to(::Type{S}) where {S}
    return S <: MOI.EqualTo
end

function apply_sector_pruning!(model::JuMP.Model, recovered, classification; atol::Float64 = 1e-12)
    expected_zero = classification.expected_zero
    pruned_vars = _pruned_var_set(recovered, expected_zero)

    # Do not materialize `JuMP.fix` here by default. In this BPSDP/MOI stack a
    # fix is copied as one scalar equality per real/imag orphan, i.e. it recreates
    # the 1×1-cone bloat we are trying to remove. We instead substitute the known
    # zero into every scalar equality/objective term and delete any row that turns
    # into 0 == 0. The variables may remain in JuMP as unused scalars, but they no
    # longer enter the SDP handed to the optimizer.
    fixed_scalars = 0

    objective_scrub = _objective_scrub_pruned_terms!(model, pruned_vars; atol)

    scalar_equalities_seen = 0
    scalar_equalities_touched = 0
    scalar_coefficients_zeroed = 0
    tautologies_deleted = 0
    contradictions = String[]

    for (F, S) in JuMP.list_of_constraint_types(model)
        F <: JuMP.GenericAffExpr || continue
        _is_scalar_equal_to(S) || continue
        for cref in collect(JuMP.all_constraints(model, F, S))
            scalar_equalities_seen += 1
            obj = JuMP.constraint_object(cref)
            func = obj.func
            touched = false
            for (coef, var) in collect(JuMP.linear_terms(func))
                if var in pruned_vars && abs(coef) > atol
                    JuMP.set_normalized_coefficient(cref, var, zero(coef))
                    scalar_coefficients_zeroed += 1
                    touched = true
                end
            end
            touched || continue
            scalar_equalities_touched += 1

            # Re-fetch after coefficient edits. If only constants remain, delete
            # 0 == 0, or abort on c == nonzero. This is the load-bearing guard
            # against pinning the wrong sector.
            obj2 = JuMP.constraint_object(cref)
            func2 = obj2.func
            nonzero_terms = 0
            for (coef, _var) in JuMP.linear_terms(func2)
                abs(coef) > atol && (nonzero_terms += 1)
            end
            if nonzero_terms == 0
                residual_constant = func2.constant - MOI.constant(obj2.set)
                if abs(residual_constant) > atol
                    push!(contradictions, "$(cref): constant residual $(residual_constant)")
                else
                    JuMP.delete(model, cref)
                    tautologies_deleted += 1
                end
            end
        end
    end

    if !isempty(contradictions)
        msg = join(contradictions[1:min(end, 10)], "\n")
        error("sector pruning produced nonzero constant equality after pinning; target sector/classifier is wrong. First contradictions:\n$msg")
    end

    deleted_scalars = 0
    batch_delete_error = nothing
    if !isempty(pruned_vars)
        try
            JuMP.delete(model, collect(pruned_vars))
            deleted_scalars = length(pruned_vars)
        catch err
            batch_delete_error = sprint(showerror, err)
        end
    end

    return Dict{String,Any}(
        "materialized_fixed_scalar_variables" => fixed_scalars,
        "logical_zero_orphan_scalar_variables" => length(pruned_vars),
        "deleted_scalar_variables_after_substitution" => deleted_scalars,
        "retained_logical_zero_scalar_variables" => length(pruned_vars) - deleted_scalars,
        "batch_delete_error" => batch_delete_error,
        "fixed_orphan_keys" => length(expected_zero),
        "objective_pruned_terms_zeroed" => objective_scrub.touched,
        "objective_pruned_coeff_l1" => objective_scrub.removed_l1,
        "scalar_equalities_seen" => scalar_equalities_seen,
        "scalar_equalities_touched" => scalar_equalities_touched,
        "scalar_coefficients_zeroed" => scalar_coefficients_zeroed,
        "tautological_equalities_deleted" => tautologies_deleted,
        "atol" => atol,
    )
end

function progress_monitor_factory(cg_iterations::Vector{Int})
    return function (print_level, outer, inner, primal, dual, mu, perr, derr)
        push!(cg_iterations, Int(inner))
        if print_level > 0 && outer % print_level == 0
            @printf("BPSDP %6d %6d primal=% .10e dual=% .10e mu=% .3e perr=% .3e derr=% .3e\n",
                outer, inner, primal, dual, mu, perr, derr)
            flush(stdout)
        end
        return nothing
    end
end

function bpsdp_optimizer_factory(options::SectorPruningOptions, cg_iterations::Vector{Int})
    monitor = progress_monitor_factory(cg_iterations)
    return () -> BPSDP.Optimizer(
        max_iter = options.bpsdp_max_iter,
        cg_max_iter = options.bpsdp_cg_max_iter,
        mu_update_frequency = options.bpsdp_mu_update_frequency,
        penalty_parameter = options.bpsdp_penalty,
        cg_convergence = options.bpsdp_cg_tol,
        dynamic_cg_convergence = true,
        sdp_objective_convergence = options.bpsdp_obj_tol,
        sdp_error_convergence = options.bpsdp_err_tol,
        guess_type = :zero,
        print_level = options.bpsdp_print_level,
        progress_monitor = monitor,
        dependent_rows = :keep,
    )
end

function bpsdp_options_json(options::SectorPruningOptions)
    return Dict{String,Any}(
        "max_iter" => options.bpsdp_max_iter,
        "cg_max_iter" => options.bpsdp_cg_max_iter,
        "mu_update_frequency" => options.bpsdp_mu_update_frequency,
        "penalty_parameter" => options.bpsdp_penalty,
        "cg_convergence" => options.bpsdp_cg_tol,
        "dynamic_cg_convergence" => true,
        "sdp_objective_convergence" => options.bpsdp_obj_tol,
        "sdp_error_convergence" => options.bpsdp_err_tol,
        "guess_type" => "zero",
        "dependent_rows" => "keep",
    )
end

function solve_pass!(data, options::SectorPruningOptions, classification, pass::Symbol)
    pruning_enabled = pass == :pruned
    result = Dict{String,Any}(
        "pass" => string(pass),
        "pruning_enabled" => pruning_enabled,
        "stage" => "building_jump_model",
        "started_at" => string(now()),
        "rss_start" => rss_string(),
    )
    cones = Dict{String,Any}()

    model = nothing
    extract_monomap = nothing
    recovered = nothing
    cg_iterations = Int[]
    optimize_error = nothing

    try
        result["jump_build_seconds"] = @elapsed begin
            if options.lowering == :moment_variables
                model, extract_monomap = build_jump_model(
                    data.moment_problem;
                    formulation = :moment_variables,
                    representation = :real,
                    orphan_policy = :free_variables,
                )
            elseif options.lowering == :psd_blocks
                model, extract_monomap = build_jump_model(
                    data.moment_problem;
                    formulation = :psd_blocks,
                    representation = :complex,
                    orphan_policy = :free_variables,
                )
            else
                error("bad lowering mode $(options.lowering)")
            end
        end
        cones["before_pruning"] = model_cone_summary(model)
        recovered = recover_orphan_handles(model, data.moment_problem)
        result["orphan_recovery"] = recovered.spotcheck

        if pruning_enabled
            result["stage"] = "applying_sector_pruning"
            result["pruning_seconds"] = @elapsed begin
                result["pruning"] = apply_sector_pruning!(model, recovered, classification)
            end
        else
            result["pruning_seconds"] = 0.0
            result["pruning"] = Dict{String,Any}(
                "fixed_scalar_variables" => 0,
                "fixed_orphan_keys" => 0,
                "tautological_equalities_deleted" => 0,
                "note" => "unpruned pass",
            )
        end
        cones["after_pruning"] = model_cone_summary(model)

        result["stage"] = "optimizing_bpsdp"
        set_optimizer(model, bpsdp_optimizer_factory(options, cg_iterations))
        result["solve_wall_seconds_measured"] = @elapsed begin
            try
                optimize!(model)
            catch err
                optimize_error = sprint(showerror, err, catch_backtrace())
            end
        end
        result["optimize_error"] = optimize_error
        result["termination_status"] = try string(termination_status(model)) catch err string("unavailable: ", err) end
        result["raw_status"] = try string(MOI.get(model, MOI.RawStatusString())) catch err string("unavailable: ", err) end
        result["objective_value"] = try finite_or_string(objective_value(model)) catch err string("unavailable: ", err) end
        result["cg_iterations_per_outer"] = cg_iterations
        result["cg_iterations_total_from_monitor"] = sum(cg_iterations)

        raw = raw_optimizer(model)
        result["bpsdp_problem"] = raw_problem_summary(raw)
        state_summary = bpsdp_state_summary(raw)
        state_summary !== nothing && (result["bpsdp_state"] = state_summary)
        cones["bpsdp_problem"] = raw_problem_summary(raw)

        if optimize_error === nothing && result["termination_status"] == "OPTIMAL"
            if pruning_enabled
                result["orphan_value_spotcheck"] = Dict{String,Any}(
                    "checked" => 0,
                    "note" => "pruned orphan variables were deleted after zero substitution; unpruned runs perform the value-level handle check",
                )
            else
                result["orphan_value_spotcheck"] = verify_orphan_value_spotcheck!(extract_monomap, recovered, data.moment_problem)
            end
        end

        result["stage"] = "finished"
    catch err
        result["stage"] = "error"
        result["error"] = sprint(showerror, err, catch_backtrace())
    finally
        result["finished_at"] = string(now())
        result["rss_finish"] = rss_string()
    end

    return (; result, cones, model)
end

function objective_float(value)
    value isa Real && isfinite(Float64(value)) && return Float64(value)
    return nothing
end

function equivalence_report(runs::Dict{String,Any})
    haskey(runs, "unpruned") && haskey(runs, "pruned") ||
        return Dict{String,Any}("checked" => false, "reason" => "need both unpruned and pruned runs")
    un = runs["unpruned"]
    pr = runs["pruned"]
    string(get(un, "termination_status", "")) == "OPTIMAL" || begin
        status = get(un, "termination_status", "missing")
        return Dict{String,Any}("checked" => false, "reason" => "unpruned status $(status)")
    end
    string(get(pr, "termination_status", "")) == "OPTIMAL" || begin
        status = get(pr, "termination_status", "missing")
        return Dict{String,Any}("checked" => false, "reason" => "pruned status $(status)")
    end
    obj_un = objective_float(get(un, "objective_value", nothing))
    obj_pr = objective_float(get(pr, "objective_value", nothing))
    (obj_un === nothing || obj_pr === nothing) &&
        return Dict{String,Any}("checked" => false, "reason" => "objective unavailable")
    diff = abs(obj_pr - obj_un)
    tol = max(1e-9, 1e-9 * abs(obj_un))
    return Dict{String,Any}(
        "checked" => true,
        "objective_unpruned" => obj_un,
        "objective_pruned" => obj_pr,
        "abs_diff" => diff,
        "tolerance" => tol,
        "passed" => diff <= tol,
    )
end

function moment_problem_summary(data, options::SectorPruningOptions)
    mp = data.moment_problem
    hpsd_sizes = [size(mat, 1) for (cone, mat) in mp.constraints if cone == :HPSD]
    zero_sizes = [size(mat) for (cone, mat) in mp.constraints if cone == :Zero]
    return Dict{String,Any}(
        "nk" => data.nk,
        "norb" => data.norb,
        "total_electrons" => data.total_electrons,
        "blocking" => string(options.blocking),
        "include_one_d" => options.include_one_d,
        "spin_resolved_trace" => options.spin_resolved_trace,
        "singlet_s2" => options.singlet_s2,
        "integrals" => options.integrals_path,
        "hf_active" => data.hf_active,
        "n_hpsd_blocks" => length(hpsd_sizes),
        "hpsd_block_sizes" => hpsd_sizes,
        "largest_hpsd_block" => isempty(hpsd_sizes) ? 0 : maximum(hpsd_sizes),
        "zero_constraint_sizes" => [collect(s) for s in zero_sizes],
        "n_linear_moments" => length(mp.linear.moments),
        "n_free_keys" => length(mp.linear.free_keys),
        "n_pivots" => length(mp.linear.pivots),
    )
end

function passes_for_pruning(pruning::Symbol)
    pruning == :on && return [:pruned]
    pruning == :off && return [:unpruned]
    pruning == :both && return [:unpruned, :pruned]
    error("bad pruning mode $pruning")
end

function print_run_header(options::SectorPruningOptions)
    println("== H4/Nk sector-pruned JuMP/BPSDP run ==")
    @printf("%-36s %s\n", "integrals", options.integrals_path)
    @printf("%-36s %s\n", "blocking", string(options.blocking))
    @printf("%-36s %s\n", "include ¹D", string(options.include_one_d))
    @printf("%-36s %s\n", "paper spin", string(options.spin_resolved_trace && options.singlet_s2))
    @printf("%-36s %s\n", "target sector", target_sector_string(options.target_sector))
    @printf("%-36s %s\n", "pruning", string(options.pruning))
    @printf("%-36s %s\n", "lowering", string(options.lowering))
    @printf("%-36s %s\n", "output root", options.output_dir)
    @printf("%-36s %s\n", "max RSS at start", rss_string())
    println()
    flush(stdout)
    return nothing
end

function main(argv = ARGS)
    options = parse_sector_pruning_options(argv)
    print_run_header(options)

    build_seconds = @elapsed data = H4PeriodicBase.build_h4_pqg_moment_problem(h4_base_options(options))
    outdir = run_output_dir(options.output_dir, data.nk)
    mkpath(outdir)

    @printf("%-36s %.3f s\n", "MomentProblem build", build_seconds)
    @printf("%-36s %s\n", "run output", outdir)
    @printf("%-36s %d\n", "linear moments", length(data.moment_problem.linear.moments))
    @printf("%-36s %d\n", "orphan free keys", length(data.moment_problem.linear.free_keys))
    flush(stdout)

    classification = classify_orphans(data.moment_problem, options.target_sector; nk = data.nk, norb = data.norb)
    kernel = classification_json(data.moment_problem, classification, options.target_sector; nk = data.nk, norb = data.norb)
    write_json(joinpath(outdir, "kernel.json"), kernel)
    @printf("%-36s %d\n", "symmetry-forbidden orphans", length(classification.expected_zero))
    @printf("%-36s %d\n", "target-sector unexpected orphans", length(classification.physical_unexpected))
    flush(stdout)

    report = Dict{String,Any}(
        "schema_version" => 1,
        "created_at" => string(now()),
        "argv" => argv,
        "julia_version" => string(VERSION),
        "julia_threads" => Threads.nthreads(),
        "blas_threads" => BLAS.get_num_threads(),
        "build_seconds" => build_seconds,
        "moment_problem" => moment_problem_summary(data, options),
        "target_sector" => target_sector_string(options.target_sector),
        "pruning_mode" => string(options.pruning),
        "jump_lowering" => Dict(
            "requested" => string(options.lowering),
            "formulation" => options.lowering == :moment_variables ? "moment_variables" : "psd_blocks",
            "representation" => options.lowering == :moment_variables ? "real" : "complex",
            "orphan_policy" => "free_variables",
            "note" => "No aux_psd_free is used. In current source, moment_variables maps free_keys to y_re/y_im moment variables; psd_blocks is the path that actually calls _declare_free_orphan_moments! and has trailing orphan_re/orphan_im.",
        ),
        "bpsdp_options" => bpsdp_options_json(options),
        "runs" => Dict{String,Any}(),
    )
    cones = Dict{String,Any}("runs" => Dict{String,Any}())

    for pass in passes_for_pruning(options.pruning)
        println()
        @printf("== pass: %s ==\n", string(pass))
        flush(stdout)
        pass_data = solve_pass!(data, options, classification, pass)
        label = pass == :pruned ? "pruned" : "unpruned"
        report["runs"][label] = pass_data.result
        cones["runs"][label] = pass_data.cones
        write_json(joinpath(outdir, "report.json"), report)
        write_json(joinpath(outdir, "cones.json"), cones)
        @printf("%-36s %s\n", "stage", get(pass_data.result, "stage", "unknown"))
        @printf("%-36s %s\n", "termination", get(pass_data.result, "termination_status", "missing"))
        @printf("%-36s %s\n", "objective", string(get(pass_data.result, "objective_value", "missing")))
        @printf("%-36s %s\n", "walltime", string(get(pass_data.result, "solve_wall_seconds_measured", "missing")))
        flush(stdout)
    end

    eq = equivalence_report(report["runs"])
    report["equivalence"] = eq
    write_json(joinpath(outdir, "report.json"), report)
    write_json(joinpath(outdir, "cones.json"), cones)

    if get(eq, "checked", false)
        @printf("%-36s %.6e <= %.6e (%s)\n", "objective diff", eq["abs_diff"], eq["tolerance"], eq["passed"] ? "PASS" : "FAIL")
        eq["passed"] || error("pruned/unpruned objective mismatch exceeds tolerance")
    else
        @printf("%-36s %s\n", "equivalence", get(eq, "reason", "not checked"))
    end

    println("Wrote ", joinpath(outdir, "report.json"))
    println("Wrote ", joinpath(outdir, "cones.json"))
    println("Wrote ", joinpath(outdir, "kernel.json"))
    return (; data, report, cones, kernel)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
