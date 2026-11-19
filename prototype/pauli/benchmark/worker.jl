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
    @assert all(n -> abs(bethe_ring(n).energy - exact_energy("chain",n)) < 1e-12, (4,6,8,10))
    torus = Set(minmax(i, j) for (i, j) in bonds("square", 4))
    @assert length(torus) == 32 && (1, 4) in torus && (1, 13) in torus && (1, 6) ∉ torus
    @assert all(i -> count(e -> i in e, torus) == 4, 1:16)
    println("Ring ED (L=4,6) and 4×4 torus topology: OK")
    exit()
end

method, path, lattice, extent, extra, strengthened_arg = ARGS
L = parse(Int, extent)
strengthened = parse(Bool, strengthened_arg)
if method == "ed"
    save_result(path, Dict("exact_energy_per_site" => exact_energy(lattice, L)))
    exit()
elseif method == "bethe"
    lattice == "chain" || error("Bethe oracle is only for the periodic ring")
    ground = bethe_ring(L)
    save_result(path, Dict("exact_energy_per_site" => ground.energy, "bethe_max_residual" => ground.residual))
    exit()
elseif startswith(method,"audit_")
    include("reference.jl")
    check_reference_checkout(lattice)
    include("prototype_adapter.jl")
    save_result(path,prototype_audit(lattice,L,parse(Int,extra),strengthened,Symbol(method[7:end])))
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
    result["result_count"] = MOI.get(model, MOI.ResultCount())
    if result["result_count"] > 0
        for (key, attribute) in (("primal_objective_per_site", MOI.ObjectiveValue()),
                                ("dual_objective_per_site", MOI.DualObjectiveValue()))
            try
                value = MOI.get(model, attribute)
                isfinite(value) && (result[key] = value)
            catch err
                err isa Mosek.MosekError || rethrow()
            end
        end
        if haskey(result,"primal_objective_per_site") && haskey(result,"dual_objective_per_site")
            p, d = result["primal_objective_per_site"], result["dual_objective_per_site"]
            gap = abs(p - d) / max(1, abs(p), abs(d))
            isfinite(gap) && (result["solver_objective_relative_gap"] = gap)
        end
    end
    # These are violations of the returned interior-point solution, not the
    # scaled PFEAS/DFEAS in the stock iteration log. Neither changes acceptance.
    task = MOI.get(model, MOI.RawSolver())
    try
        if Mosek.solutiondef(task, Mosek.MSK_SOL_ITR)
            _, pc, pv, pb, pk, pa, _, _, _, dc, dv, db, dk, da =
                Mosek.getsolutioninfonew(task, Mosek.MSK_SOL_ITR)
            for (key, violations) in (("mosek_primal_max_violation", (pc,pv,pb,pk,pa)),
                                      ("mosek_dual_max_violation", (dc,dv,db,dk,da)))
                all(isfinite, violations) && (result[key] = maximum(violations))
            end
        end
    catch err
        err isa Mosek.MosekError || rethrow()
    end
    return termination == MOI.OPTIMAL && primal == dual == MOI.FEASIBLE_POINT
end

if method == "preflight"
    # An actual conic optimization checks the license and settings, not just import.
    model = Model(OPTIMIZER)
    @variable(model, x)
    @constraint(model, [x 1; 1 x] in PSDCone())
    @objective(model, Min, x)
    optimize!(model)
    result = Dict{String,Any}("julia_version" => string(VERSION),
        "mosek_version" => join(Mosek.getversion(),'.'), "mosek_jl_version" => string(pkgversion(Mosek)),
        "mosektools_version" => string(pkgversion(MosekTools)), "jump_version" => string(pkgversion(JuMP)),
        "julia_threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(),
        "blas_config" => sprint(show, BLAS.get_config()), "mosek_threads" => THREADS, "tolerance" => TOL)
    record_status!(result, unsafe_backend(model)) &&
        abs(value(x) - 1) <= 1e-7 || error("Mosek preflight did not give the expected optimal solution")
    save_result(path, result)
    exit()
end

result = Dict{String,Any}()
if method == "reference"
    include("reference.jl")
    check_reference_checkout(lattice)
    options = reference_options(lattice, L, parse(Int, extra), strengthened)
    settings = QMBCertify.mosek_para(TOL, TOL, TOL, THREADS)
    result["reference_options"] = Dict(string(k) => v for (k,v) in pairs(options))
    order = 2
    result["reference_revision"] = REFERENCE_REVS[lattice]
    result["reference_order"] = order
    if strengthened
        result["reference_stationarity_seed"] = 20261118
        Random.seed!(result["reference_stationarity_seed"])
    end
    save_result(path, result)
    energy, data = QMBCertify.GSB(options.H_supp, options.H_coe, L, order; options..., mosek_setting=settings)
    result["reference_raw_rows"] = length.(data.basis)
    result["reference_unique_rows_by_sector"] = length.(unique.(data.basis))
    result["reference_moment_terms_including_identity"] = length(data.tsupp)
else
    method in ("full", "graph_distance", "local_products") || error("unknown worker method $method")
    include("prototype_adapter.jl")
    energy = prototype_run(result, path, lattice, L, parse(Int, extra), strengthened, Symbol(method))
end
if energy !== nothing && isfinite(energy)
    result["reported_objective_per_site"] = energy
end
save_result(path, result)
