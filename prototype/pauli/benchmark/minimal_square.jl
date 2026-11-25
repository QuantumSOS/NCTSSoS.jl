# Run with a Julia environment containing the prototype, JuMP and MosekTools:
# julia --project=<env> minimal_square.jl L max_rounds output_directory [--certify]
using JuMP, LinearAlgebra, MosekTools, SparseArrays, TOML
import NCTSSoS, PauliWorkflowRedesign as WF
include("cases.jl")
include("square_moment.jl")

function minimal_square(L, max_rounds, out; certify=false, initial_motifs=(), local_row_limit=256,
        max_rows=200L^2)
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
    optimizer = optimizer_with_attributes(MosekTools.Optimizer,
        "MSK_IPAR_NUM_THREADS"=>Threads.nthreads(),
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS"=>tol, "MSK_DPAR_INTPNT_CO_TOL_DFEAS"=>tol,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP"=>tol)
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
            "preflight_dual"=>dual_objective_value(m),"threads"=>Threads.nthreads())
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
            WF.compile_minimal_square(H,L; motifs, recover_gram_rows=certify)
        end
        f = plan.standard
        r["rows"] = length(plan.rows); r["variables"] = length(f.c); r["nnz"] = nnz(f.A)
        r["psd_sides"] = [c.side for c in f.cones if c.kind==:psd]
        model,vars = stage("transfer") do
            mod,vars = NCTSSoS.to_jump(f)
            set_optimizer(mod,optimizer); MOI.Utilities.attach_optimizer(mod)
            mod,vars
        end
        stage(()->optimize!(model),"solve")
        status = termination_status(model)
        r["termination"] = string(status); r["primal_status"] = string(primal_status(model))
        r["dual_status"] = string(dual_status(model)); r["native_seconds"] = solve_time(model)
        r["accepted"] = status==MOI.OPTIMAL && primal_status(model)==dual_status(model)==MOI.FEASIBLE_POINT
        stage(()->square_original_diagnostics!(r,f,model,vars),"diagnostics")
        metrics = ("original_equality_inf","original_primal_psd_violation","original_balance_inf",
            "original_dual_psd_violation","original_primal_nonnegative_violation",
            "original_dual_nonnegative_violation","original_objective_relative_gap")
        r["strict_original_1e-8"] = all(k->get(r,k,Inf)<=tol,metrics)
        d = get(r,"original_dual_objective_per_site",NaN)
        if r["accepted"] && isfinite(d)
            r["accepted_numerical_dual_per_spin"] = d
        end
        save(); r["accepted"] && isfinite(d) || error("No accepted numerical point; see $path")
        if certify
            cert = stage("verify") do
                WF.verify_pauli_dual(H,f,plan.table,WF.pauli_dual_blocks(f,model),plan.lifts();
                    symmetry_generators=plan.symmetry.symmetry_generators, spin_average=true)
            end
            r["certificate_projection"] = string(cert.spin_projection)
            r["certified_per_spin"] = Float64(exact_scale*cert.bound,RoundDown)
            r["certified_rational"] = string(exact_scale*cert.bound)
            r["residual_nonidentity_l1"] = Float64(exact_scale*cert.residual_nonidentity_l1,RoundUp)
            save()
        end
        (;plan.rows,plan.table,plan.constraint_ids,plan.domain,moments=value.(vars),bound=d,status)
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
    length(ARGS) in (3,4) && (length(ARGS)==3 || ARGS[4]=="--certify") || error("Expected L rounds out [--certify]")
    minimal_square(parse(Int,ARGS[1]),parse(Int,ARGS[2]),ARGS[3]; certify=length(ARGS)==4)
end
