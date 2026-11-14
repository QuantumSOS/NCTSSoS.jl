module NCTSSoSMosekExt
using NCTSSoS: Next
import Mosek
import JuMP
using SparseArrays: findnz

mutable struct MosekHandle{B}
    task::Mosek.MSKtask
    form::Next.MomentStandardForm
    real_form::Next.MomentStandardForm
    backend::B
end

function append_constraints!(task,f,firstrow=1,firstcone=1)
    m = length(f.b)
    Mosek.appendafes(task,m-firstrow+1)
    I,J,V = findnz(f.A[firstrow:end,:])
    Mosek.putafefentrylist(task,I .+ (firstrow-1),J,-V)
    Mosek.putafegslice(task,firstrow,m+1,f.b[firstrow:end])
    for cone in f.cones[firstcone:end]
        domain = cone.kind==:zero ? Mosek.appendrzerodomain(task,cone.length) :
            cone.kind==:nonnegative ? Mosek.appendrplusdomain(task,cone.length) :
            cone.kind==:soc ? Mosek.appendquadraticconedomain(task,cone.length) :
            cone.kind==:psd ? Mosek.appendsvecpsdconedomain(task,cone.length) :
            throw(ArgumentError("Mosek writer received unsupported cone $(cone.kind)."))
        Mosek.appendaccseq(task,domain,Int64(cone.first),nothing)
    end
end

function Next.prepare_standard(form::Next.MomentStandardForm,backend::Next.MosekBackend)
    f = Next.real_psd_form(form)
    task = Mosek.maketask()
    n = length(f.c)
    Mosek.appendvars(task,n)
    Mosek.putvarboundsliceconst(task,1,n+1,Mosek.MSK_BK_FR,-Inf,Inf)
    Mosek.putclist(task,collect(1:n),f.c)
    Mosek.putcfix(task,f.offset)
    Mosek.putobjsense(task,Mosek.MSK_OBJECTIVE_SENSE_MINIMIZE)
    for (name,value) in pairs(backend.options)
        key = String(name)
        if startswith(key,"MSK_IPAR_")
            Mosek.putnaintparam(task,key,value)
        elseif startswith(key,"MSK_DPAR_")
            Mosek.putnadouparam(task,key,value)
        elseif startswith(key,"MSK_SPAR_")
            Mosek.putnastrparam(task,key,value)
        else
            throw(ArgumentError("Mosek option $name must use its MSK_IPAR_, MSK_DPAR_ or MSK_SPAR_ name."))
        end
    end
    append_constraints!(task,f)
    return MosekHandle(task,form,f,backend)
end

function Next.solve_standard(h::MosekHandle)
    Mosek.optimize(h.task)
    raw = Mosek.getsolsta(h.task,Mosek.MSK_SOL_ITR)
    optimal = raw==Mosek.MSK_SOL_STA_OPTIMAL
    status = optimal ? JuMP.MOI.OPTIMAL : raw==Mosek.MSK_SOL_STA_PRIM_INFEAS_CER ? JuMP.MOI.INFEASIBLE :
        raw==Mosek.MSK_SOL_STA_DUAL_INFEAS_CER ? JuMP.MOI.DUAL_INFEASIBLE : JuMP.MOI.OTHER_ERROR
    bound = optimal ? Mosek.getprimalobj(h.task,Mosek.MSK_SOL_ITR) : nothing
    values = optimal ? Mosek.getxx(h.task,Mosek.MSK_SOL_ITR) : Float64[]
    return Next.ConicSolution(bound,status,sprint(show,raw),values,h)
end

function Next.solve_standard(form::Next.MomentStandardForm,h::MosekHandle)
    Next._check_append_only(h.form,form)
    f = Next.real_psd_form(form)
    append_constraints!(h.task,f,length(h.real_form.b)+1,length(h.real_form.cones)+1)
    h.form,h.real_form = form,f
    return Next.solve_standard(h)
end
end
