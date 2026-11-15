module NCTSSoSMosekExt
using NCTSSoS: Next
import Mosek
import JuMP
using SparseArrays: sparse, findnz, nzrange, rowvals, nonzeros

struct MosekHandle{B}
    task::Mosek.MSKtask
    form::Next.MomentStandardForm
    backend::B
    representation::Symbol
end

# Realify HPSD as [Re(M) -Im(M); Im(M) Re(M)]. This leaves moment
# variables unchanged, so both primal values and SOS equality duals recover them.
function real_psd_form(form)
    any(c->c.kind==:hpsd,form.cones) || return form
    I,J,V = Int[],Int[],Float64[]
    cones = Next.DataCone[]
    row = 0
    for cone in form.cones
        first = row+1
        if cone.kind!=:hpsd
            for k in cone.first:cone.first+cone.length-1
                row+=1; push!(I,row); push!(J,k); push!(V,1.)
            end
            push!(cones,Next.DataCone(cone.kind,cone.side,first,cone.length,cone.label))
            continue
        end
        n=cone.side
        coordinate(i,j)=cone.first+(j-1)*(2n-j+1)+2(i-j)
        for j in 1:2n,i in j:2n
            row+=1
            a,b=mod1(i,n),mod1(j,n)
            if (i>n)==(j>n)
                push!(I,row); push!(J,coordinate(a,b)-(a==b ? 0 : 1)); push!(V,1.)
            elseif a!=b
                push!(I,row); push!(J,coordinate(max(a,b),min(a,b))); push!(V,a>b ? 1. : -1.)
            end
        end
        push!(cones,Next.DataCone(:psd,2n,first,n*(2n+1),cone.label))
    end
    T=sparse(I,J,V,row,length(form.b))
    return Next.MomentStandardForm(T*form.A,T*form.b,form.c,form.offset,cones)
end

function domain(task,kind,length)
    kind==:zero && return Mosek.appendrzerodomain(task,length)
    kind==:nonnegative && return Mosek.appendrplusdomain(task,length)
    kind==:soc && return Mosek.appendquadraticconedomain(task,length)
    kind==:psd && return Mosek.appendsvecpsdconedomain(task,length)
    throw(ArgumentError("Mosek writer received unsupported cone $kind."))
end

function write_moment!(task,f)
    m,n=size(f.A)
    Mosek.appendvars(task,n)
    Mosek.putvarboundsliceconst(task,1,n+1,Mosek.MSK_BK_FR,-Inf,Inf)
    Mosek.putclist(task,collect(1:n),f.c)
    Mosek.putobjsense(task,Mosek.MSK_OBJECTIVE_SENSE_MINIMIZE)
    Mosek.appendafes(task,m)
    I,J,V=findnz(f.A)
    Mosek.putafefentrylist(task,I,J,-V)
    Mosek.putafegslice(task,1,m+1,f.b)
    for cone in f.cones
        Mosek.appendaccseq(task,domain(task,cone.kind,cone.length),Int64(cone.first),nothing)
    end
end

function write_sos!(task,f)
    # max -b'z, A'z=-c, z in K*. PSD coordinates become native bar variables.
    m,n=size(f.A)
    scalar=zeros(Int,m); bar=zeros(Int,m); row=zeros(Int,m); col=zeros(Int,m)
    scale=ones(m)
    ns,nb=0,0
    for cone in f.cones
        if cone.kind==:psd
            nb+=1
            Mosek.appendbarvars(task,[cone.side])
            k=cone.first
            for j in 1:cone.side,i in j:cone.side
                bar[k]=nb; row[k]=i; col[k]=j
                scale[k]=i==j ? 1. : inv(sqrt(2.))
                k+=1
            end
        else
            for k in cone.first:cone.first+cone.length-1
                ns+=1; scalar[k]=ns
            end
        end
    end
    Mosek.appendvars(task,ns)
    Mosek.putvarboundsliceconst(task,1,ns+1,Mosek.MSK_BK_FR,-Inf,Inf)
    Mosek.appendcons(task,n)
    Mosek.putconboundslice(task,1,n+1,fill(Mosek.MSK_BK_FX,n),-f.c,-f.c)
    Mosek.putobjsense(task,Mosek.MSK_OBJECTIVE_SENSE_MAXIMIZE)
    ai,aj,av=Int[],Int[],Float64[]
    bi,bj,bk,bl,bv=Int[],Int[],Int[],Int[],Float64[]
    cj,ck,cl,cv=Int[],Int[],Int[],Float64[]
    for k in 1:m
        if bar[k]>0
            if !iszero(f.b[k])
                push!(cj,bar[k]); push!(ck,row[k]); push!(cl,col[k]); push!(cv,-f.b[k]*scale[k])
            end
        else
            Mosek.putcj(task,scalar[k],-f.b[k])
        end
    end
    for j in 1:n,t in nzrange(f.A,j)
        k=rowvals(f.A)[t]; v=nonzeros(f.A)[t]
        if bar[k]>0
            push!(bi,j); push!(bj,bar[k]); push!(bk,row[k]); push!(bl,col[k]); push!(bv,v*scale[k])
        else
            push!(ai,j); push!(aj,scalar[k]); push!(av,v)
        end
    end
    Mosek.putaijlist(task,ai,aj,av)
    Mosek.putbarablocktriplet(task,bi,bj,bk,bl,bv)
    Mosek.putbarcblocktriplet(task,cj,ck,cl,cv)
    for cone in f.cones
        cone.kind in (:psd,:zero) && continue
        first=scalar[cone.first]
        if cone.kind==:nonnegative
            Mosek.putvarboundsliceconst(task,first,first+cone.length,Mosek.MSK_BK_LO,0.,Inf)
        elseif cone.kind==:soc
            Mosek.appendcone(task,Mosek.MSK_CT_QUAD,0.,collect(first:first+cone.length-1))
        else
            throw(ArgumentError("Mosek writer received unsupported cone $(cone.kind)."))
        end
    end
end

function Next.prepare_standard(form::Next.MomentStandardForm,backend::Next.MosekBackend)
    f=real_psd_form(form)
    representation=backend.form==:auto ? (length(f.c)<=length(f.b) ? :sos : :moment) : backend.form
    task=Mosek.maketask()
    try
        for (name,value) in pairs(backend.options)
            key=String(name)
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
        representation==:sos ? write_sos!(task,f) : write_moment!(task,f)
        Mosek.putcfix(task,f.offset)
    catch
        Mosek.deletetask(task)
        rethrow()
    end
    return MosekHandle(task,form,backend,representation)
end

function Next.solve_standard(h::MosekHandle)
    termination=Mosek.optimize(h.task)
    raw=Mosek.getsolsta(h.task,Mosek.MSK_SOL_ITR)
    optimal=raw==Mosek.MSK_SOL_STA_OPTIMAL
    sos=h.representation==:sos
    status=optimal ? JuMP.MOI.OPTIMAL :
        raw==Mosek.MSK_SOL_STA_PRIM_INFEAS_CER ? (sos ? JuMP.MOI.DUAL_INFEASIBLE : JuMP.MOI.INFEASIBLE) :
        raw==Mosek.MSK_SOL_STA_DUAL_INFEAS_CER ? (sos ? JuMP.MOI.INFEASIBLE : JuMP.MOI.DUAL_INFEASIBLE) :
        termination==Mosek.MSK_RES_TRM_MAX_TIME ? JuMP.MOI.TIME_LIMIT :
        termination==Mosek.MSK_RES_TRM_MAX_ITERATIONS ? JuMP.MOI.ITERATION_LIMIT :
        termination==Mosek.MSK_RES_TRM_STALL ? JuMP.MOI.SLOW_PROGRESS : JuMP.MOI.OTHER_ERROR
    bound=optimal ? Mosek.getprimalobj(h.task,Mosek.MSK_SOL_ITR) : nothing
    values=optimal ? (sos ? -Mosek.gety(h.task,Mosek.MSK_SOL_ITR) : Mosek.getxx(h.task,Mosek.MSK_SOL_ITR)) : Float64[]
    diagnostics=(representation=h.representation,termination=sprint(show,termination),warm_started=false)
    return Next.ConicSolution(bound,status,sprint(show,raw),values,h,diagnostics)
end

function Next.solve_standard(form::Next.MomentStandardForm,h::MosekHandle)
    Next._check_append_only(h.form,form)
    return Next.solve_standard(form,h.backend)
end
end
