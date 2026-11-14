module NCTSSoSSCSExt
using NCTSSoS: Next
import SCS
import JuMP
using SparseArrays: spzeros

mutable struct SCSHandle{B}
    form::Next.MomentStandardForm
    backend::B
    permutation::Vector{Int}
    x::Vector{Float64}
    y::Vector{Float64}
    slack::Vector{Float64}
    warm_start::Bool
end

function Next.prepare_standard(form::Next.MomentStandardForm,backend::Next.SCSBackend)
    linear_solver = backend.linear_solver
    isnothing(linear_solver) || linear_solver in (SCS.DirectSolver,SCS.IndirectSolver) ||
        throw(ArgumentError("linear_solver must be SCS.DirectSolver or SCS.IndirectSolver; pass the type, not an instance."))
    order = (:zero,:nonnegative,:soc,:psd,:hpsd)
    all(c->c.kind in order,form.cones) || throw(ArgumentError("SCS writer received an unsupported cone."))
    permutation = Int[]
    for kind in order,cone in form.cones
        cone.kind==kind && append!(permutation,cone.first:cone.first+cone.length-1)
    end
    return SCSHandle(form,backend,permutation,zeros(length(form.c)),zeros(length(form.b)),zeros(length(form.b)),false)
end

function Next.solve_standard(h::SCSHandle)
    f,p = h.form,h.permutation
    countrows(kind) = sum(c.length for c in f.cones if c.kind==kind;init=0)
    sides(kind) = [c.side for c in f.cones if c.kind==kind]
    m,n = size(f.A)
    linear_solver = isnothing(h.backend.linear_solver) ? SCS.DirectSolver : h.backend.linear_solver
    options = merge((warm_start=h.warm_start,),h.backend.options)
    sol = SCS.scs_solve(linear_solver,m,n,f.A[p,:],spzeros(n,n),f.b[p],f.c,
        countrows(:zero),countrows(:nonnegative),Float64[],Float64[],sides(:soc),sides(:psd),sides(:hpsd),
        0,0,Float64[],Int[],Int[],Int[],Int[],Int[],Int[],h.x,h.y,h.slack;options...)
    status = sol.ret_val==1 ? JuMP.MOI.OPTIMAL : sol.ret_val==2 ? JuMP.MOI.ALMOST_OPTIMAL :
        sol.ret_val==-2 ? JuMP.MOI.INFEASIBLE : sol.ret_val==-1 ? JuMP.MOI.DUAL_INFEASIBLE : JuMP.MOI.OTHER_ERROR
    bound = sol.ret_val==1 ? Float64(sol.info.pobj+f.offset) : nothing
    h.warm_start = sol.ret_val in (1,2) && all(isfinite,h.x) && all(isfinite,h.y) && all(isfinite,h.slack)
    diagnostics = (primal_residual=sol.info.res_pri,dual_residual=sol.info.res_dual,gap=sol.info.gap,
        iterations=Int(sol.info.iter),setup_seconds=sol.info.setup_time/1000,solve_seconds=sol.info.solve_time/1000,
        warm_started=options.warm_start)
    return Next.ConicSolution(bound,status,SCS.raw_status(sol.info),copy(sol.x),h,diagnostics)
end

function Next.solve_standard(form::Next.MomentStandardForm,h::SCSHandle)
    Next._check_append_only(h.form,form)
    # SCS.jl's public one-shot API closes the native workspace after each solve.
    next = Next.prepare_standard(form,h.backend)
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
    h.form,h.permutation,h.x,h.y,h.slack = next.form,next.permutation,next.x,next.y,next.slack
    return Next.solve_standard(h)
end
end
