using LinearAlgebra, Random, TOML
include("cases.jl")

function save_result(path, result)
    open(path * ".tmp", "w") do io
        TOML.print(io, result)
    end
    mv(path * ".tmp", path; force=true)
end

BLAS.set_num_threads(1)
Threads.nthreads(:default) == 1 && Threads.nthreads(:interactive) == 0 &&
    BLAS.get_num_threads() == 1 || error("expected Julia 1,0 and BLAS 1 thread")

if ARGS == ["smoke"]
    @assert isapprox(exact_energy("chain", 4), -0.5; atol=1e-13)
    @assert isapprox(exact_energy("chain", 6), -(2 + sqrt(13))/12; atol=1e-13)
    torus = Set(minmax(i, j) for (i, j) in bonds("square", 4))
    @assert length(torus) == 32 && (1, 4) in torus && (1, 13) in torus && (1, 6) ∉ torus
    @assert all(i -> count(e -> i in e, torus) == 4, 1:16)
    println("Ring ED (L=4,6) and 4×4 torus topology: OK")
    exit()
end

method, path, lattice, extent, extra, pair_range, control, strengthened_arg = ARGS
L = parse(Int, extent)
strengthened = parse(Bool, strengthened_arg)
if method == "ed"
    save_result(path, Dict("exact_energy_per_site" => exact_energy(lattice, L)))
    exit()
end

using JuMP, MosekTools
import Mosek
const MOI = JuMP.MOI
const TOL = 1e-8
const THREADS = 1
const MOSEK_PARAMETERS = (
    "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => TOL,
    "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => TOL,
    "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => TOL,
    "MSK_IPAR_NUM_THREADS" => THREADS)
const OPTIMIZER = optimizer_with_attributes(Mosek.Optimizer, MOSEK_PARAMETERS..., MOI.Silent() => true)

function record_status!(result, model)
    termination = MOI.get(model, MOI.TerminationStatus())
    primal = MOI.get(model, MOI.PrimalStatus())
    dual = MOI.get(model, MOI.DualStatus())
    result["termination_status"] = string(termination)
    result["primal_status"] = string(primal)
    result["dual_status"] = string(dual)
    result["raw_status"] = MOI.get(model, MOI.RawStatusString())
    for (name, expected) in MOSEK_PARAMETERS
        actual = MOI.get(model, MOI.RawOptimizerAttribute(name))
        actual == expected || error("Mosek $name=$actual; expected $expected")
    end
    termination == MOI.OPTIMAL && primal == dual == MOI.FEASIBLE_POINT || return false
    p, d = MOI.get(model, MOI.ObjectiveValue()), MOI.get(model, MOI.DualObjectiveValue())
    if isfinite(p) && isfinite(d)
        result["primal_objective_per_site"] = p
        result["dual_objective_per_site"] = d
        result["solver_objective_relative_gap"] = abs(p - d) / max(1, abs(p), abs(d))
    end
    return true
end

if method == "preflight"
    # An actual conic optimization checks the license and settings, not just import.
    model = Model(OPTIMIZER)
    @variable(model, x)
    @constraint(model, [x 1; 1 x] in PSDCone())
    @objective(model, Min, x)
    optimize!(model)
    result = Dict{String,Any}("julia_version" => string(VERSION),
        "mosek_version" => string(Mosek.getversion()), "mosek_jl_version" => string(pkgversion(Mosek)),
        "mosektools_version" => string(pkgversion(MosekTools)), "jump_version" => string(pkgversion(JuMP)),
        "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "blas_config" => sprint(show, BLAS.get_config()), "mosek_threads" => THREADS, "tolerance" => TOL)
    record_status!(result, unsafe_backend(model)) &&
        abs(value(x) - 1) <= 1e-7 || error("Mosek preflight did not give the expected optimal solution")
    save_result(path, result)
    exit()
end

result = Dict{String,Any}()
if method in ("reference", "reference-build-smoke")
    build_only = method == "reference-build-smoke"
    !build_only || L == (strengthened ? 8 : 4) || error("reference build smoke requires a small control case")
    include("reference.jl")
    install_reference_probe()
    options = reference_options(lattice, L, parse(Int, extra), strengthened)
    settings = QMBCertify.mosek_para(TOL, TOL, TOL, THREADS)
    result["reference_options"] = Dict(string(k) => v for (k,v) in pairs(options))
    order = strengthened ? 3 : 2
    result["reference_order"] = order
    if strengthened
        result["reference_stationarity_seed"] = 20261118
        Random.seed!(result["reference_stationarity_seed"])
    end
    save_result(path, result)
    psd_count(model) = sum(num_constraints(model, F, S) for (F,S) in list_of_constraint_types(model)
        if S in (MOI.PositiveSemidefiniteConeTriangle, MOI.PositiveSemidefiniteConeSquare); init=0)
    started = time_ns()
    function checkpoint(stage, model, elapsed)
        if stage == :stationarity
            # In the SOS model each filtered equality has one free multiplier.
            result["stationarity_equalities"] = elapsed
            return
        elseif stage == :before_rdm
            result["psd_before_rdm"] = psd_count(model)
            return
        elseif stage == :after_rdm
            result["rdm_psd_blocks"] = psd_count(model) - pop!(result, "psd_before_rdm")
            return
        elseif stage == :built
            result["build_seconds"] = (time_ns() - started) / 1e9
            result["reference_variable_count"] = num_variables(model)
            result["reference_constraint_count"] = num_constraints(model; count_variable_in_set_constraints=true)
            save_result(path, result)
            !strengthened || (get(result, "rdm_psd_blocks", 0) > 0 &&
                get(result, "stationarity_equalities", 0) > 0) || error("reference strengthening missing")
            build_only && exit()
        else
            result["solve_seconds"] = elapsed
            record_status!(result, unsafe_backend(model))
            result["optimizer_reported_seconds"] = solve_time(model)
        end
        save_result(path, result)
    end
    energy, _ = reference_run(checkpoint, L, order, settings, options)
else
    method == "prototype" || error("unknown worker method $method")
    include("prototype_adapter.jl")
    energy = prototype_run(result, path, lattice, L, parse(Int, pair_range), control == "true", strengthened)
end
if energy !== nothing
    isfinite(energy) || error("non-finite objective")
    result["energy_per_site"] = energy
end
save_result(path, result)
