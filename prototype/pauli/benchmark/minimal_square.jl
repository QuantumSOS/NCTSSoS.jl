# Run with a Julia environment containing the prototype, JuMP, Mosek and MosekTools:
# julia --project=<env> minimal_square.jl L max_rounds output_directory [--certify] [--jump]
using JuMP, LinearAlgebra, MosekTools, SparseArrays, TOML
import Mosek, NCTSSoS, PauliWorkflowRedesign as WF
include("cases.jl")
include("square_moment.jl")

function square_native_point(task, trm, form)
    sol = Mosek.MSK_SOL_ITR
    status = get(MosekTools._TERMINATION_STATUS_MAP,trm.value,MOI.OTHER_ERROR)
    Mosek.solutiondef(task,sol) || return (;status,primal=MOI.NO_SOLUTION,dual=MOI.NO_SOLUTION,x=nothing,z=nothing)
    raw = Mosek.getsolsta(task,sol)
    # Use the same continuous solution status mappings as the JuMP reference.
    if trm == Mosek.MSK_RES_OK
        status = raw == Mosek.MSK_SOL_STA_OPTIMAL ? MOI.OPTIMAL :
            MosekTools._PROSTA_STATUS_MAP[Mosek.getprosta(task,sol).value]
    end
    primal = get(MosekTools._PRIMAL_STATUS_MAP,raw.value,MOI.UNKNOWN_RESULT_STATUS)
    dual = get(MosekTools._DUAL_STATUS_MAP,raw.value,MOI.UNKNOWN_RESULT_STATUS)
    x = primal == MOI.NO_SOLUTION ? nothing : Mosek.getxx(task,sol)
    # The native writer appends ACCs in original cone order and uses b-A*x
    # with the minimization sense. doty is already positive-dual svec order.
    z = dual == MOI.NO_SOLUTION ? nothing : Mosek.getaccdotys(task,sol)
    @assert isnothing(x) || length(x) == length(form.c)
    @assert isnothing(z) || length(z) == length(form.b)
    return (;status,primal,dual,x,z)
end

function square_native_dual_blocks(form,z)
    all(c -> c.kind in (:psd,:zero),form.cones) || error("G5 requires real PSD and equality cones")
    sqrt2 = Rational{BigInt}(sqrt(2.0))
    return map(form.cones) do c
        c.kind == :zero && return view(z,c.first:c.first+c.length-1)
        # Undo svec using the exact binary scale from the original form.
        [Rational{BigInt}(z[c.first+NCTSSoS.tri_offset(c.side,max(i,j),min(i,j))]) /
            (i == j ? 1 : sqrt2) for i in 1:c.side, j in 1:c.side]
    end
end

function minimal_square(L, max_rounds, out; certify=false, initial_motifs=(), local_row_limit=256,
        max_rows=200L^2, solver_tol=1e-8, rdm_shape=(2,2), solver_path=:native)
    solver_path in (:native,:jump) || throw(ArgumentError("Expected solver_path=:native or :jump."))
    edges = bonds("square",L)
    N = L^2
    @assert length(edges) == length(Set(minmax(i,j) for (i,j) in edges)) == 2N
    _, (x,y,z) = NCTSSoS.pauli_variables(N)
    H = sum(x[i]*x[j]+y[i]*y[j]+z[i]*z[j] for (i,j) in edges)/(4N)
    # The compiler uses Float64; 1/(4N) is not dyadic at e.g. L=6.
    # Convert a rigorous certificate back to the exact per-site Hamiltonian.
    coefficient = only(unique(real(c) for (c,_) in NCTSSoS.terms(H)))
    exact_scale = (1//big(4N))/Rational{BigInt}(coefficient)
    tol = 1e-8
    options = (;MSK_IPAR_NUM_THREADS=Threads.nthreads(), MSK_DPAR_INTPNT_CO_TOL_PFEAS=solver_tol,
        MSK_DPAR_INTPNT_CO_TOL_DFEAS=solver_tol, MSK_DPAR_INTPNT_CO_TOL_REL_GAP=solver_tol)
    optimizer = optimizer_with_attributes(MosekTools.Optimizer,(String(k)=>v for (k,v) in pairs(options))...)
    backend = NCTSSoS.MosekBackend(;form=:moment,options...)
    # Check licensing on the execution host, not just on the build host.
    m = Model(optimizer); set_silent(m)
    @variable(m,t); @constraint(m,[t 1; 1 t] in PSDCone()); @objective(m,Min,t)
    optimize!(m)
    @assert termination_status(m)==MOI.OPTIMAL && primal_status(m)==dual_status(m)==MOI.FEASIBLE_POINT
    @assert abs(dual_objective_value(m)-1) <= tol
    mkpath(out)
    round = 0
    function solve_motifs(motifs)
        round += 1
        r = Dict{String,Any}("L"=>L,"round"=>round,"motifs"=>repr(motifs),
            "bonds"=>length(edges),"J2"=>0,"unit"=>"spin-1/2 energy per site",
            "preflight_dual"=>dual_objective_value(m),"threads"=>Threads.nthreads(),"solver_tol"=>solver_tol,
            "rdm_shape"=>collect(rdm_shape),"solver_path"=>string(solver_path))
        path = joinpath(out,"round-$round.toml")
        save() = open(io->TOML.print(io,r; sorted=true),path,"w")
        function stage(f,k)
            println("START ",k); flush(stdout)
            timed = @timed f()
            r[k*"_seconds"] = timed.time; r[k*"_bytes_allocated"] = timed.bytes
            r["peak_rss_bytes"] = Sys.maxrss(); save()
            println("END ",k," ",timed.time); flush(stdout)
            timed.value
        end
        plan = stage("build") do
            WF.compile_minimal_square(H,L; motifs, rdm_shape, recover_folded_gram=certify)
        end
        f = plan.standard
        r["rows"] = length(plan.rows); r["variables"] = length(f.c); r["nnz"] = nnz(f.A)
        r["psd_sides"] = [c.side for c in f.cones if c.kind==:psd]
        save()
        model = nothing
        point = if solver_path == :native
            handle = stage(()->NCTSSoS.prepare_standard(f,backend),"transfer")
            try
                trm = stage(()->Mosek.optimize(handle.task),"solve")
                r["native_seconds"] = Mosek.getdouinf(handle.task,Mosek.MSK_DINF_OPTIMIZER_TIME)
                stage("diagnostics") do
                    p = square_native_point(handle.task,trm,f)
                    square_point_diagnostics!(r,f,p.x,p.z)
                    p
                end
            finally
                Mosek.deletetask(handle.task)
            end
        else
            model,vars = stage("transfer") do
                mod,vars = NCTSSoS.to_jump(f)
                set_optimizer(mod,optimizer); MOI.Utilities.attach_optimizer(mod)
                mod,vars
            end
            stage(()->optimize!(model),"solve")
            r["native_seconds"] = solve_time(model)
            stage("diagnostics") do
                x,z = square_jump_points(f,model,vars)
                square_point_diagnostics!(r,f,x,z)
                (;status=termination_status(model),primal=primal_status(model),dual=dual_status(model),x,z)
            end
        end
        status = point.status
        r["termination"] = string(status); r["primal_status"] = string(point.primal)
        r["dual_status"] = string(point.dual)
        r["accepted"] = status==MOI.OPTIMAL && point.primal==point.dual==MOI.FEASIBLE_POINT
        metrics = ("original_equality_inf","original_primal_psd_violation","original_balance_inf",
            "original_dual_psd_violation","original_primal_nonnegative_violation",
            "original_dual_nonnegative_violation","original_objective_relative_gap")
        r["strict_original_1e-8"] = all(k->get(r,k,Inf)<=tol,metrics)
        d = get(r,"original_dual_objective_per_site",NaN)
        if r["accepted"] && isfinite(d)
            r["accepted_numerical_dual_per_spin"] = d
        end
        save()
        # The exact operator identity can certify a weaker value from a
        # non-optimal dual. Never promote that iterate to numerical acceptance.
        if certify && !isnothing(point.z)
            cert = stage("verify") do
                dual = solver_path == :native ? square_native_dual_blocks(f,point.z) : WF.pauli_dual_blocks(f,model)
                WF.verify_pauli_dual(H,f,plan.table,dual,plan.lifts();
                    symmetry_generators=plan.symmetry.symmetry_generators, spin_average=true)
            end
            r["certificate_projection"] = string(cert.spin_projection)
            r["certified_per_spin"] = Float64(exact_scale*cert.bound,RoundDown)
            r["certified_rational"] = string(exact_scale*cert.bound)
            r["residual_nonidentity_l1"] = Float64(exact_scale*cert.residual_nonidentity_l1,RoundUp)
            save()
        end
        r["accepted"] && isfinite(d) || error("No accepted numerical point; see $path")
        (;plan.rows,plan.table,plan.constraint_ids,plan.domain,moments=point.x,bound=d,status)
    end
    if iszero(max_rounds)
        solve_motifs(initial_motifs)
    else
        result = WF.price_square(H,L,solve_motifs; max_rounds,initial_motifs,local_row_limit,max_rows,
            batch_size=1,pricing_mode=:orbit)
        open(joinpath(out,"pricing.toml"),"w") do io
            TOML.print(io,Dict("stop_reason"=>string(result.stop_reason),"rounds"=>[
                Dict("rows"=>h.row_count,"selected"=>repr(h.selected),"priced"=>length(h.scores),
                    "cuts"=>count(s->s.cut!==nothing,h.scores),"seconds"=>h.pricing_seconds,
                    "partial"=>h.pricing.partial,"pool"=>h.pricing.largest_pool) for h in result.history]))
        end
        println("PRICING STOP ", result.stop_reason," motifs=",result.motifs)
    end
end

if abspath(PROGRAM_FILE)==@__FILE__
    3 <= length(ARGS) <= 5 && all(in(("--certify","--jump")),ARGS[4:end]) && allunique(ARGS[4:end]) ||
        error("Expected L rounds out [--certify] [--jump]")
    minimal_square(parse(Int,ARGS[1]),parse(Int,ARGS[2]),ARGS[3]; certify="--certify" in ARGS[4:end],
        solver_path="--jump" in ARGS[4:end] ? :jump : :native)
end
