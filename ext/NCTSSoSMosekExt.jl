module NCTSSoSMosekExt
import NCTSSoS
# Internal construction helpers shared by this package's precompile workload.
using NCTSSoS: build, moment_standard_form
import NCTSSoS: _prepare
import Mosek
import JuMP
using SparseArrays: sparse, nzrange, rowvals, nonzeros
using PrecompileTools: @compile_workload

struct MosekHandle{B}
    task::Mosek.MSKtask
    form::NCTSSoS.MomentStandardForm
    backend::B
    representation::Symbol
end

function NCTSSoS.task_sizes(h::MosekHandle)
    task=h.task
    sides=Int[Mosek.getdimbarvarj(task,i) for i in 1:Mosek.getnumbarvar(task)]
    for i in 1:Mosek.getnumacc(task)
        domain=Mosek.getaccdomain(task,i)
        if Mosek.getdomaintype(task,domain)==Mosek.MSK_DOMAIN_SVEC_PSD_CONE
            n=Mosek.getdomainn(task,domain)
            push!(sides,(isqrt(8n+1)-1)÷2)
        end
    end
    return (affine_rows=Int(Mosek.getnumcon(task)+Mosek.getnumafe(task)),
        scalar_variables=Int(Mosek.getnumvar(task)),psd_sides=sides)
end

# Realify HPSD as [Re(M) -Im(M); Im(M) Re(M)]. This leaves moment
# variables unchanged, so both primal values and SOS equality duals recover them.
function real_psd_form(form)
    any(c->c.kind==:hpsd,form.cones) || return form
    I,J,V = Int[],Int[],Float64[]
    cones = NCTSSoS.DataCone[]
    row = 0
    for cone in form.cones
        first = row+1
        if cone.kind!=:hpsd
            for k in cone.first:cone.first+cone.length-1
                row+=1; push!(I,row); push!(J,k); push!(V,1.)
            end
            push!(cones,NCTSSoS.DataCone(cone.kind,cone.side,first,cone.length,cone.label))
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
        push!(cones,NCTSSoS.DataCone(:psd,2n,first,n*(2n+1),cone.label))
    end
    T=sparse(I,J,V,row,length(form.b))
    return NCTSSoS.MomentStandardForm(T*form.A,T*form.b,form.c,form.offset,cones)
end

function domain(task,kind,length)
    kind==:zero && return Mosek.appendrzerodomain(task,length)
    kind==:nonnegative && return Mosek.appendrplusdomain(task,length)
    kind==:soc && return Mosek.appendquadraticconedomain(task,length)
    kind==:psd && return Mosek.appendsvecpsdconedomain(task,length)
    throw(ArgumentError("Mosek writer received unsupported cone $kind."))
end

function write_moment_entries!(task, A)
    rows, values = rowvals(A), nonzeros(A)
    isempty(values) && return nothing
    # Bound the triplet workspace and the wrapper's zero-based index copies,
    # independently of nnz(A). The CSC remains the only full Julia matrix.
    count = min(length(values), 1<<18)
    I, J, V = Vector{Int64}(undef,count), Vector{Int32}(undef,count), Vector{Float64}(undef,count)
    k = 0
    for j in axes(A,2), t in nzrange(A,j)
        k += 1
        I[k], J[k], V[k] = rows[t], j, -values[t]
        if k == count
            Mosek.putafefentrylist(task,I,J,V)
            k = 0
        end
    end
    if k > 0
        resize!(I,k); resize!(J,k); resize!(V,k)
        Mosek.putafefentrylist(task,I,J,V)
    end
    return nothing
end

function write_moment!(task,f)
    m,n=size(f.A)
    Mosek.appendvars(task,n)
    Mosek.putvarboundsliceconst(task,1,n+1,Mosek.MSK_BK_FR,-Inf,Inf)
    Mosek.putclist(task,collect(1:n),f.c)
    Mosek.putobjsense(task,Mosek.MSK_OBJECTIVE_SENSE_MINIMIZE)
    Mosek.appendafes(task,m)
    write_moment_entries!(task,f.A)
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

function NCTSSoS.prepare_standard(form::NCTSSoS.MomentStandardForm,backend::NCTSSoS.MosekBackend)
    f=real_psd_form(form)
    return prepare_native(f,backend,form)
end

function _prepare(built::NCTSSoS.BuiltRelaxation,backend::NCTSSoS.MosekBackend)
    f=get!(built.cache,:mosek_real) do
        form=haskey(built.cache,:moment) ? built.cache[:moment] : moment_standard_form(built.program)
        real_psd_form(form)
    end::NCTSSoS.MomentStandardForm
    return prepare_native(f,backend,f)
end

function prepare_native(f,backend,source)
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
    return MosekHandle(task,source,backend,representation)
end

function NCTSSoS.solve_standard(h::MosekHandle)
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
    diagnostics=(representation=h.representation,termination=sprint(show,termination),warm_started=false,
        solve_seconds=Mosek.getdouinf(h.task,Mosek.MSK_DINF_OPTIMIZER_TIME),
        iterations=Int(Mosek.getintinf(h.task,Mosek.MSK_IINF_INTPNT_ITER)),
        primal_residual=Mosek.getdouinf(h.task,Mosek.MSK_DINF_INTPNT_PRIMAL_FEAS),
        dual_residual=Mosek.getdouinf(h.task,Mosek.MSK_DINF_INTPNT_DUAL_FEAS))
    return NCTSSoS.ConicSolution(bound,status,sprint(show,raw),values,h,diagnostics)
end

function NCTSSoS.solve_standard(form::NCTSSoS.MomentStandardForm,h::MosekHandle)
    NCTSSoS._check_append_only(h.form,form)
    return NCTSSoS.solve_standard(form,h.backend)
end

# Build both native forms without calling the optimizer or requesting a license.
@compile_workload let
    forms=NCTSSoS.MomentStandardForm[]
    n=4
    registry,(x,y,z)=NCTSSoS.pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    symmetry=(NCTSSoS.SitePermutation([2:n;1]),NCTSSoS.SitePermutation(collect(n:-1:1)),NCTSSoS.spin_rotations(registry)...)
    push!(forms,moment_standard_form(build(NCTSSoS.Problem(h;symmetry),NCTSSoS.Relaxation(order=2)).program))
    registry,(a,ad)=NCTSSoS.fermionic_variables(n)
    h=-sum(ad[i]*a[mod1(i+1,n)]+ad[mod1(i+1,n)]*a[i] for i in 1:n)
    number=sum(ad[i]*a[i] for i in 1:n)
    problem=NCTSSoS.Problem(h;symmetry=(NCTSSoS.NumberRotation(registry),))
    relaxation=NCTSSoS.Relaxation(order=1,strengthen=(NCTSSoS.FermionicConditions(1:n;conditions=(:P,:Q,:G)),NCTSSoS.Sector(number,2)))
    push!(forms,moment_standard_form(build(problem,relaxation).program))
    n=8
    _,(x,y,z)=NCTSSoS.pauli_variables(n)
    h=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    push!(forms,moment_standard_form(build(NCTSSoS.Problem(h;symmetry=NCTSSoS.detect_symmetry(h)),NCTSSoS.Relaxation(order=1)).program))
    _,(a,b)=NCTSSoS.site_variables(NCTSSoS.Unipotent(),2,2;functional=NCTSSoS.RealState())
    h=NCTSSoS.state_symbol(a[1]*b[1])-NCTSSoS.state_symbol(a[1])*NCTSSoS.state_symbol(b[1])
    push!(forms,moment_standard_form(build(NCTSSoS.Problem(h),NCTSSoS.Relaxation(order=2)).program))
    _,(x,)=NCTSSoS.site_variables(NCTSSoS.Unipotent(),1,4;functional=NCTSSoS.RealTrace())
    h=NCTSSoS.trace_symbol(x[1]*x[3])-NCTSSoS.trace_symbol(x[1])*NCTSSoS.trace_symbol(x[3])
    push!(forms,moment_standard_form(build(NCTSSoS.Problem(h),NCTSSoS.Relaxation(order=2)).program))
    for standard in forms,form in (:sos,:moment)
        backend=NCTSSoS.MosekBackend(;form)
        handle=NCTSSoS.prepare_standard(standard,backend)
        # Compile the solve/extraction code without executing the optimizer.
        precompile(NCTSSoS.solve_standard,(typeof(standard),typeof(backend)))
        Mosek.deletetask(handle.task)
    end
end
end
