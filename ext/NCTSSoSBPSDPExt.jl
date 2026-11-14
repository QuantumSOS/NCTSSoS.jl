module NCTSSoSBPSDPExt
using NCTSSoS: Next
import BPSDP
import JuMP
using SparseArrays: sparse
using LinearAlgebra: adjoint, transpose

mutable struct BPSDPHandle{F,B,S,E}
    form::F
    backend::B
    solver::S
    expansion::E
end

function Next.prepare_standard(form::Union{Next.MomentStandardForm,Next.SOSStandardForm},backend::Next.BPSDPBackend)
    f = form isa Next.MomentStandardForm ? Next.sos_standard_form(form) : form
    T = any(c->c.kind==:hpsd,f.cones) ? ComplexF64 : Float64
    I,J,V = Int[],Int[],T[]
    blocks = Union{BPSDP.PSDBlock,BPSDP.HermitianPSDBlock}[]
    offset = 0
    for cone in f.cones
        if cone.kind in (:psd,:hpsd)
            n,k = cone.side,cone.first
            push!(blocks,cone.kind==:hpsd ? BPSDP.HermitianPSDBlock(n) : BPSDP.PSDBlock(n))
            for j in 1:n,i in j:n
                if i==j
                    push!(I,offset+i+(j-1)*n); push!(J,k); push!(V,1); k+=1
                else
                    append!(I,(offset+i+(j-1)*n,offset+j+(i-1)*n)); append!(J,(k,k)); append!(V,(1/sqrt(2),1/sqrt(2))); k+=1
                    if cone.kind==:hpsd
                        append!(I,(offset+i+(j-1)*n,offset+j+(i-1)*n)); append!(J,(k,k)); append!(V,(im/sqrt(2),-im/sqrt(2))); k+=1
                    end
                end
            end
            offset+=n*n
        elseif cone.kind in (:free,:nonnegative)
            for k in cone.first:cone.first+cone.length-1
                for sign in (cone.kind==:free ? (1,-1) : (1,))
                    offset+=1; push!(blocks,BPSDP.PSDBlock(1))
                    push!(I,offset); push!(J,k); push!(V,sign)
                end
            end
        else
            throw(ArgumentError("BPSDP supports PSD, HPSD, nonnegative and free SOS coordinates; use a different backend for $(cone.kind)."))
        end
    end
    E = sparse(I,J,V,offset,length(f.c))
    problem = BPSDP.Problem(blocks,-E*f.c,f.b,BPSDP.SparseOperator(f.A*transpose(E)))
    return BPSDPHandle(form,backend,BPSDP.Solver(problem;backend.options...),E)
end

function Next.solve_standard(h::BPSDPHandle)
    sol = BPSDP.solve!(h.solver)
    optimal = sol.status==:optimal
    status = optimal ? JuMP.MOI.OPTIMAL : sol.status==:max_iter ? JuMP.MOI.ITERATION_LIMIT : JuMP.MOI.OTHER_ERROR
    bound = optimal ? h.form.offset-sol.primal_objective : nothing
    values = h.form isa Next.MomentStandardForm ? copy(sol.y) : real.(adjoint(h.expansion)*sol.x)
    return Next.ConicSolution(bound,status,string(sol.status),values,h)
end

function Next.solve_standard(form::Next.MomentStandardForm,h::BPSDPHandle{Next.MomentStandardForm})
    Next._check_append_only(h.form,form)
    # BPSDP has no public block-resize operation; rebuild its solver for added cones.
    next = Next.prepare_standard(form,h.backend)
    h.form,h.solver,h.expansion = next.form,next.solver,next.expansion
    return Next.solve_standard(h)
end
end
