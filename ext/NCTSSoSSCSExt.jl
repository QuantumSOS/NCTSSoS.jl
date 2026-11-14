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
end

function Next.prepare_standard(form::Next.MomentStandardForm,backend::Next.SCSBackend)
    order = (:zero,:nonnegative,:soc,:psd,:hpsd)
    all(c->c.kind in order,form.cones) || throw(ArgumentError("SCS writer received an unsupported cone."))
    permutation = Int[]
    for kind in order,cone in form.cones
        cone.kind==kind && append!(permutation,cone.first:cone.first+cone.length-1)
    end
    return SCSHandle(form,backend,permutation,zeros(length(form.c)),zeros(length(form.b)),zeros(length(form.b)))
end

function Next.solve_standard(h::SCSHandle)
    f,p = h.form,h.permutation
    countrows(kind) = sum(c.length for c in f.cones if c.kind==kind;init=0)
    sides(kind) = [c.side for c in f.cones if c.kind==kind]
    m,n = size(f.A)
    sol = SCS.scs_solve(SCS.DirectSolver,m,n,f.A[p,:],spzeros(n,n),f.b[p],f.c,
        countrows(:zero),countrows(:nonnegative),Float64[],Float64[],sides(:soc),sides(:psd),sides(:hpsd),
        0,0,Float64[],Int[],Int[],Int[],Int[],Int[],Int[],h.x,h.y,h.slack;h.backend.options...)
    status = sol.ret_val==1 ? JuMP.MOI.OPTIMAL : sol.ret_val==2 ? JuMP.MOI.ALMOST_OPTIMAL :
        sol.ret_val==-2 ? JuMP.MOI.INFEASIBLE : sol.ret_val==-1 ? JuMP.MOI.DUAL_INFEASIBLE : JuMP.MOI.OTHER_ERROR
    bound = sol.ret_val==1 ? Float64(sol.info.pobj+f.offset) : nothing
    return Next.ConicSolution(bound,status,SCS.raw_status(sol.info),copy(sol.x),h)
end

function Next.solve_standard(form::Next.MomentStandardForm,h::SCSHandle)
    Next._check_append_only(h.form,form)
    # SCS.jl's public one-shot API closes the native workspace after each solve.
    next = Next.prepare_standard(form,h.backend)
    h.form,h.permutation,h.x,h.y,h.slack = next.form,next.permutation,next.x,next.y,next.slack
    return Next.solve_standard(h)
end
end
