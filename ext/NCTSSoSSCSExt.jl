module NCTSSoSSCSExt
import NCTSSoS
using NCTSSoS: pauli_variables, site_variables, Unipotent, Problem, Relaxation,
    Correlative, SitePermutation, solve,
    TermSplit, fermionic_variables, Free, Positive, word_type, AntiunitaryPauli,
    pauli_letter, DensityMatrix, InvariantCone, spin_rotations, NumberRotation,
    FermionicConditions, Sector, detect_symmetry
using PrecompileTools: @compile_workload
import SCS
import JuMP
using SparseArrays: spzeros

mutable struct SCSHandle{B,D}
    form::NCTSSoS.MomentStandardForm
    backend::B
    permutation::Vector{Int}
    data::D
    x::Vector{Float64}
    y::Vector{Float64}
    slack::Vector{Float64}
    warm_start::Bool
end

function NCTSSoS.prepare_standard(form::NCTSSoS.MomentStandardForm,backend::NCTSSoS.SCSBackend)
    linear_solver = backend.linear_solver
    isnothing(linear_solver) || linear_solver in (SCS.DirectSolver,SCS.IndirectSolver) ||
        throw(ArgumentError("linear_solver must be SCS.DirectSolver or SCS.IndirectSolver; pass the type, not an instance."))
    order = (:zero,:nonnegative,:soc,:psd,:hpsd)
    all(c->c.kind in order,form.cones) || throw(ArgumentError("SCS writer received an unsupported cone."))
    permutation = Int[]
    for kind in order,cone in form.cones
        cone.kind==kind && append!(permutation,cone.first:cone.first+cone.length-1)
    end
    countrows(kind) = sum(c.length for c in form.cones if c.kind==kind;init=0)
    sides(kind) = [c.side for c in form.cones if c.kind==kind]
    m,n=size(form.A)
    data=(A=form.A[permutation,:],P=spzeros(n,n),b=form.b[permutation],
        zero=countrows(:zero),nonnegative=countrows(:nonnegative),
        soc=sides(:soc),psd=sides(:psd),hpsd=sides(:hpsd))
    return SCSHandle(form,backend,permutation,data,zeros(n),zeros(m),zeros(m),false)
end

NCTSSoS.solve_standard(h::SCSHandle) = _solve_standard(h,false)

function _solve_standard(h::SCSHandle,warm_start::Bool)
    f,d = h.form,h.data
    m,n = size(f.A)
    linear_solver = isnothing(h.backend.linear_solver) ? SCS.DirectSolver : h.backend.linear_solver
    # Only the append path may reuse a solution; repeated solves stay cold.
    options = merge(h.backend.options,(warm_start=warm_start && get(h.backend.options,:warm_start,true),))
    sol = SCS.scs_solve(linear_solver,m,n,d.A,d.P,d.b,f.c,
        d.zero,d.nonnegative,Float64[],Float64[],d.soc,d.psd,d.hpsd,
        0,0,Float64[],Int[],Int[],Int[],Int[],Int[],Int[],h.x,h.y,h.slack;options...)
    status = sol.ret_val==1 ? JuMP.MOI.OPTIMAL : sol.ret_val==2 ? JuMP.MOI.ALMOST_OPTIMAL :
        sol.ret_val==-2 ? JuMP.MOI.INFEASIBLE : sol.ret_val==-1 ? JuMP.MOI.DUAL_INFEASIBLE : JuMP.MOI.OTHER_ERROR
    raw_status = SCS.raw_status(sol.info)
    if sol.ret_val==2
        # Match SCS.jl's MOI wrapper: a resource limit is not an optimum.
        if occursin("reached time_limit_secs",raw_status)
            status = JuMP.MOI.TIME_LIMIT
        elseif occursin("reached max_iters",raw_status)
            status = JuMP.MOI.ITERATION_LIMIT
        end
    end
    bound = sol.ret_val==1 ? Float64(sol.info.pobj+f.offset) : nothing
    h.warm_start = sol.ret_val in (1,2) && all(isfinite,h.x) && all(isfinite,h.y) && all(isfinite,h.slack)
    diagnostics = (primal_residual=sol.info.res_pri,dual_residual=sol.info.res_dual,gap=sol.info.gap,
        iterations=Int(sol.info.iter),setup_seconds=sol.info.setup_time/1000,solve_seconds=sol.info.solve_time/1000,
        warm_started=options.warm_start)
    return NCTSSoS.ConicSolution(bound,status,raw_status,copy(sol.x),h,diagnostics)
end

function NCTSSoS.solve_standard(form::NCTSSoS.MomentStandardForm,h::SCSHandle)
    NCTSSoS._check_append_only(h.form,form)
    # SCS.jl's public one-shot API closes the native workspace after each solve.
    next = NCTSSoS.prepare_standard(form,h.backend)
    if h.warm_start
        next.x .= h.x
        next.slack .= (form.b-form.A*h.x)[next.permutation]
        # Appended cones can move old coordinates in SCS's required cone order.
        position = invperm(next.permutation)
        for (i,row) in enumerate(h.permutation)
            j = position[row]
            next.y[j] = h.y[i]
            next.slack[j] = h.slack[i]
        end
    end
    h.form,h.permutation,h.data,h.x,h.y,h.slack = next.form,next.permutation,next.data,next.x,next.y,next.slack
    return _solve_standard(h,h.warm_start)
end

# This extension is loaded only with SCS. Keep native solves tiny and silent;
# the base package's workload builds programs without loading a solver.
@compile_workload let
    n=4
    _,(x,y,z)=pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+0.7z[i]*z[mod1(i+1,n)] for i in 1:n)
    backend=NCTSSoS.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7)
    solve(Problem(h;symmetry=NCTSSoS.Symmetry(SitePermutation([2:n;1]))),
        Relaxation(order=2,grouping=Correlative());solver=backend)
    solve(Problem(h;symmetry=NCTSSoS.Symmetry(SitePermutation([2:n;1]),SitePermutation([mod1(2-i,n) for i in 1:n]))),
        Relaxation(order=2,grouping=Correlative(),split=TermSplit());solver=backend)
    _,(a,b)=site_variables(Unipotent(),2,2)
    bell=Problem(-(a[1]*b[1]+a[1]*b[2]+a[2]*b[1]-a[2]*b[2]))
    solve(bell,Relaxation(order=1);solver=backend)
end

@compile_workload let
    backend=NCTSSoS.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7)
    _,(a,ad)=fermionic_variables(4)
    h=-(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
    solve(Problem(h),Relaxation(order=2);solver=backend)
    _,ops=site_variables(Free(),1,2)
    x,y=ops[1]
    solve(Problem(x*y+y*x;constraints=(Positive(1-x^2-y^2),)),Relaxation(order=2);solver=backend)
end

@compile_workload let
    backend=NCTSSoS.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7)
    registry,(x,y,z)=pauli_variables(3)
    h=-(x[1]*y[2]+x[1]*y[3]+x[2]*x[3]+y[2]*y[3]+z[2]*z[3])
    W=word_type(registry,4)
    symmetry=NCTSSoS.Symmetry(SitePermutation([1,3,2]);moment_only=(AntiunitaryPauli(pauli_letter(W,1,1)),))
    for kind in (DensityMatrix(1:3),InvariantCone(1:3,"YII"))
        solve(Problem(h;symmetry),Relaxation(order=1,strengthen=(kind,));solver=backend)
    end
end

@compile_workload let
    n=4
    registry,(x,y,z)=pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    symmetry=(SitePermutation([2:n;1]),SitePermutation(collect(n:-1:1)),spin_rotations(registry)...)
    solve(Problem(h;symmetry),Relaxation(order=2);
        solver=NCTSSoS.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7))
end

@compile_workload let
    n=4
    registry,(a,ad)=fermionic_variables(n)
    h=-sum(ad[i]*a[mod1(i+1,n)]+ad[mod1(i+1,n)]*a[i] for i in 1:n)
    number=sum(ad[i]*a[i] for i in 1:n)
    problem=Problem(h;symmetry=(NumberRotation(registry),))
    relaxation=Relaxation(order=1,strengthen=(FermionicConditions(1:n;conditions=(:P,:Q,:G)),Sector(number,2)))
    solve(problem,relaxation;solver=NCTSSoS.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7))
end

@compile_workload let
    n=8
    _,(x,y,z)=pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    symmetry=detect_symmetry(h)
    solve(Problem(h;symmetry),Relaxation(order=1);
        solver=NCTSSoS.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7))
end

@compile_workload let
    backend=NCTSSoS.SCSBackend(verbose=0,eps_abs=1e-7,eps_rel=1e-7)
    _,(a,b)=site_variables(Unipotent(),2,2;functional=NCTSSoS.RealState())
    h=NCTSSoS.state_symbol(a[1]*b[1])-NCTSSoS.state_symbol(a[1])*NCTSSoS.state_symbol(b[1])
    solve(Problem(h),Relaxation(order=2);solver=backend)
    _,(x,)=site_variables(Unipotent(),1,4;functional=NCTSSoS.RealTrace())
    h=NCTSSoS.trace_symbol(x[1]*x[3])-NCTSSoS.trace_symbol(x[1])*NCTSSoS.trace_symbol(x[3])
    solve(Problem(h),Relaxation(order=2);solver=backend)
end
end
