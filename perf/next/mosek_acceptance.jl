# Licensed acceptance cases. Arguments: bell|xxz|tfi|heisenberg, auto|moment|sos|legacy_sos|jump.
# Bell uses the dense I_nn22 n=5, order-3 relaxation on both routes.
# Set MOSEK_THREADS, MOSEK_TOLERANCE and MOSEK_MAX_TIME for a bounded run.
include(joinpath(@__DIR__,"probe.jl"))
import Mosek, MosekTools
BLAS.set_num_threads(1)

function acceptance_input(name)
    if name=="bell"
        _,(a,b)=N.site_variables(N.Projector(),2,5)
        return N.Problem(-bell(a,b)),N.Relaxation(order=3)
    end
    n=name=="xxz" ? 4 : name=="tfi" ? 6 : 8
    _,(x,y,z)=N.pauli_variables(n)
    translation=N.SitePermutation([2:n;1])
    if name=="tfi"
        h=-sum(z[i]*z[mod1(i+1,n)]+0.7x[i] for i in 1:n)
        r=N.Relaxation(order=2,grouping=N.Correlative(),split=N.TermSplit(),strengthen=(N.StateOptimality(),))
        symmetry=(translation,)
    elseif name=="xxz"
        h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+0.5z[i]*z[mod1(i+1,n)] for i in 1:n)
        r=N.Relaxation(order=2,strengthen=(N.DensityMatrix(1:n),))
        symmetry=(translation,)
    elseif name=="heisenberg"
        h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
        r=N.Relaxation(order=2)
        symmetry=(translation,N.SitePermutation(collect(n:-1:1)))
    else
        error("Unknown acceptance case $name")
    end
    return N.Problem(h;symmetry),r
end

function print_solver_metrics(task)
    println("SOLVER_METRICS primal_feasibility=",Mosek.getdouinf(task,Mosek.MSK_DINF_INTPNT_PRIMAL_FEAS),
        " dual_feasibility=",Mosek.getdouinf(task,Mosek.MSK_DINF_INTPNT_DUAL_FEAS),
        " primal_objective=",Mosek.getdouinf(task,Mosek.MSK_DINF_INTPNT_PRIMAL_OBJ),
        " dual_objective=",Mosek.getdouinf(task,Mosek.MSK_DINF_INTPNT_DUAL_OBJ))
end

function run_acceptance(name,route)
    tolerance=parse(Float64,get(ENV,"MOSEK_TOLERANCE","1e-10"))
    options=(MSK_IPAR_NUM_THREADS=parse(Int,get(ENV,"MOSEK_THREADS","1")),
        MSK_DPAR_OPTIMIZER_MAX_TIME=parse(Float64,get(ENV,"MOSEK_MAX_TIME","600")),
        MSK_DPAR_INTPNT_CO_TOL_REL_GAP=tolerance,MSK_DPAR_INTPNT_CO_TOL_PFEAS=tolerance,
        MSK_DPAR_INTPNT_CO_TOL_DFEAS=tolerance)
    factory=JuMP.optimizer_with_attributes(MosekTools.Optimizer,(String(k)=>v for (k,v) in pairs(options))...)
    if route=="legacy_sos"
        name=="bell" || error("The legacy SOS reference here is Bell only")
        t=@timed input_case("bell","legacy",5,3,"dense","")()
        println("BUILD ",t.time," ",t.bytes); flush(stdout)
        w=@timed Old.sos_dualize(t.value).model
        model=w.value
        JuMP.set_optimizer(model,factory)
        println("WRITE ",w.time," ",w.bytes); flush(stdout)
        s=@timed JuMP.optimize!(model)
        println("RESULT ",JuMP.termination_status(model)," ",JuMP.has_values(model) ? JuMP.objective_value(model) : nothing,
            " SECONDS ",s.time," RSS ",Sys.maxrss())
        print_solver_metrics(JuMP.unsafe_backend(model).task)
        return
    end
    p,r=acceptance_input(name)
    backend=route=="jump" ? N.JuMPBackend(factory) : N.MosekBackend(;options...,form=Symbol(route))
    t=@timed N.build(p,r)
    println("BUILD ",t.time," ",t.bytes," SIZES ",t.value.sizes); flush(stdout)
    f=N.moment_standard_form(t.value.program)
    println("FORM ",size(f.A)," NNZ ",nnz(f.A)); flush(stdout)
    if route=="jump"
        w=@timed N.to_jump(f)
        model,_=w.value
        JuMP.set_optimizer(model,factory)
        a=@timed JuMP.MOI.Utilities.attach_optimizer(JuMP.backend(model))
        println("WRITE ",w.time," ",w.bytes," ATTACH ",a.time," ",a.bytes); flush(stdout)
        s=@timed JuMP.optimize!(model)
        println("RESULT ",JuMP.termination_status(model)," ",JuMP.has_values(model) ? JuMP.objective_value(model) : nothing,
            " SECONDS ",s.time," RSS ",Sys.maxrss())
    else
        w=@timed N.prepare_standard(f,backend)
        println("WRITE ",w.time," ",w.bytes," FORM ",w.value.representation); flush(stdout)
        Mosek.deletetask(w.value.task)
        s=@timed N.solve(p,r;solver=backend)
        println("RESULT ",s.value.status," ",s.value.bound," SECONDS ",s.time," RSS ",Sys.maxrss(),
            " DIAGNOSTICS ",s.value.solution.diagnostics)
        print_solver_metrics(s.value.solution.handle.task)
        s.value.bound===nothing || println("MOMENT_OBJECTIVE ",N.expectation(s.value,p.objective))
    end
end

abspath(PROGRAM_FILE)==(@__FILE__) && run_acceptance(ARGS...)
