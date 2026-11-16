# Standalone worker. See README.md; run through run.py for resource limits.
import NCTSSoS
const N = NCTSSoS.Next
using LinearAlgebra, SparseArrays, Statistics, Profile, TOML, Serialization
import JuMP
if length(ARGS)>=8 && ARGS[8]=="solve"
    if first(ARGS)=="hydrogen"
        import SCS
    else
        import MosekTools
    end
end

function cg_matrix(n)
    n==3 && return [0 -1 -2 0; 0 1 1 -1; -1 1 1 1; 0 -1 1 0]
    c=zeros(Int,n+1,n+1)
    for i in 1:n-1
        c[1,i+1]=c[i+1,1]=-1
    end
    for i in 1:n,j in 1:n
        c[i+1,j+1]=(i==n || j==n) ? (i==j ? -1 : 1) :
            i+j<n ? 0 : i+j==n ? 1 : -1
    end
    c
end
function bell(a,b)
    n=length(a); c=cg_matrix(n)
    sum(c[i+1,1]*a[i]+c[1,i+1]*b[i] for i in 1:n)+
        sum(c[i+1,j+1]*a[i]*b[j] for i in 1:n,j in 1:n)
end
function pauli_actions(n,permutations)
    W=N.word_type(N.Pauli(),n,8)
    actions=N.AbstractAction[N.SitePermutation(p) for p in permutations]
    append!(actions,[N.CliffordMap([N.pauli_letter(W,i,a) for a in axes for i in 1:n])
        for axes in ((1,0),(1,2))])
    N.Symmetry(actions...;moment_only=(N.AntiunitaryPauli(one(W)),))
end
function geometry(nx,ny)
    site(x,y)=mod(x,nx)+nx*mod(y,ny)+1
    bonds=Set{Tuple{Int,Int}}()
    for x in 0:nx-1,y in 0:ny-1,(dx,dy) in ((1,0),(0,1))
        a,b=site(x,y),site(x+dx,y+dy)
        a==b || push!(bonds,minmax(a,b))
    end
    maps=[[(site(f(x,y)...)) for y in 0:ny-1 for x in 0:nx-1]
        for f in ((x,y)->(x+1,y),(x,y)->(x,y+1),(x,y)->(-x,y),(x,y)->(x,-y))]
    nx==ny && push!(maps,[site(y,x) for y in 0:ny-1 for x in 0:nx-1])
    sort!(collect(bonds)),maps
end
function input_case(target,route,n,order,variant,data)
    route=="next" || error("Historical routes require the pinned checkout; see README.md")
    if target=="bell"
        if route=="next"
            _,(a,b)=N.site_variables(N.Projector(),2,n)
            p=N.Problem(-bell(a,b))
            r=variant=="cs_ts" ? N.Relaxation(;order,grouping=N.Correlative(algorithm=N.MMD()),split=N.TermSplit()) : N.Relaxation(;order)
            return ()->N.build(p,r)
        end
    elseif target in ("chain","torus")
        if target=="chain"
            bonds=[(i,mod1(i+1,n)) for i in 1:n]
            maps=[collect([2:n;1]),collect(n:-1:1)]
        else
            nx,ny=n==6 ? (3,2) : n==9 ? (3,3) : (6,6)
            bonds,maps=geometry(nx,ny)
        end
        if route=="next"
            _,ops=N.pauli_variables(n)
            h=sum(ops[a][i]*ops[a][j] for (i,j) in bonds for a in 1:3)
            p=N.Problem(h;symmetry=pauli_actions(n,maps))
            r=N.Relaxation(;order,row_rule=target=="chain" ? N.Contiguous(periodic=true) : nothing,
                grouping=startswith(variant,"regions") ? N.PhysicsInformed(level=1) : nothing,
                split=endswith(variant,"ts") ? N.TermSplit() : nothing)
            return ()->N.build(p,r)
        end
    elseif target=="hydrogen"
        d=TOML.parsefile(data); m=d["n_spatial_orbitals"]
        h=[Float64(d["h_ij"][i][j]) for i in 1:m,j in 1:m]
        eri=[Float64(d["eri"][i][j][k][l]) for i in 1:m,j in 1:m,k in 1:m,l in 1:m]
        if route=="next"
            registry,(a,ad)=N.fermionic_variables(2m)
            H=d["nuclear_repulsion"]*one(a[1])
            for p in 1:m,q in 1:m,sigma in (0,m)
                abs(h[p,q])>=1e-14 && (H+=h[p,q]*ad[p+sigma]*a[q+sigma])
            end
            for p in 1:m,q in 1:m,r in 1:m,s in 1:m,sigma in (0,m),tau in (0,m)
                c=eri[p,q,r,s]/2
                abs(c)>=1e-14 && (H+=c*ad[p+sigma]*ad[r+tau]*a[s+tau]*a[q+sigma])
            end
            H=(H+H')/2
            number=sum(ad[i]*a[i] for i in 1:2m)
            spin=sum(ad[i]*a[i]-ad[i+m]*a[i+m] for i in 1:m)
            conditions=variant in ("pqg","pqg_number") ? (:P,:Q,:G) : (:P,:Q,:G,:T1,:T2prime)
            symmetry=endswith(variant,"_number") ? N.Symmetry(;continuous=(
                N.NumberRotation(registry;modes=1:m),N.NumberRotation(registry;modes=m+1:2m))) : nothing
            p=N.Problem(H;symmetry)
            r=N.Relaxation(order=1,strengthen=(N.FermionicConditions(1:2m;conditions),
                N.Sector(number,d["n_electrons"]),N.Sector(spin,d["n_alpha"]-d["n_beta"])))
            return ()->N.build(p,r)
        end
    end
    error("unknown target $target")
end

function dimensions(built,route)
    if route=="next"
        s=hasproperty(built.sizes,:program) ? built.sizes.program : built.sizes
        return Dict{String,Any}("rows"=>s.rows,"moments"=>s.moments,"variables"=>s.real_variables,
            "blocks"=>s.blocks,"side"=>s.largest_block,
            "engine"=>(isnothing(built.sizes.closure) ? "plain" : string(built.sizes.closure.route)))
    end
end
function lower_dimensions(lowered,route)
    if route=="next"
        return Dict("lower_rows"=>size(lowered.A,1),"lower_cols"=>size(lowered.A,2),"nonzeros"=>nnz(lowered.A))
    end
end
function checkpoint(result,path)
    result["peak_rss_bytes"]=Sys.maxrss()
    open(path,"w") do io
        TOML.print(io,result;sorted=true)
    end
    println("STAGE ",get(result,"stage","")); flush(stdout)
end
function solve_probe(built,route,target,result,out)
    if target=="hydrogen"
        route=="next" && result["n"]==4 && result["variant"]=="pqg" || error("Only H4 Next PQG solve is authorized")
        form=N.moment_standard_form(built.program)
        result["stage"]="solve"; checkpoint(result,out)
        timed=@timed Base.invokelatest(N.solve_standard,form,N.SCSBackend(verbose=1,eps_abs=1e-5,eps_rel=1e-5,time_limit_secs=600))
        sol=timed.value
        result["status"]=string(sol.status)
        result["bound"]=isnothing(sol.bound) ? "unavailable" : sol.bound
        result["diagnostics"]=string(sol.diagnostics)
    elseif target=="bell"
        factory=JuMP.optimizer_with_attributes(MosekTools.Optimizer,"MSK_IPAR_NUM_THREADS"=>1,"MSK_DPAR_OPTIMIZER_MAX_TIME"=>600.0)
        if route=="next"
            form=N.moment_standard_form(built.program)
            result["stage"]="solve"; checkpoint(result,out)
            timed=@timed Base.invokelatest(N.solve_standard,form,N.JuMPBackend(factory))
            sol=timed.value
            result["status"]=string(sol.status)
            result["bound"]=isnothing(sol.bound) ? "unavailable" : -sol.bound
        end
    else
        error("No solves for spin scale targets")
    end
    result["solve_seconds"]=timed.time
    result["stage"]="complete"; checkpoint(result,out)
    println("| target | route | n | status | bound | solve seconds |")
    println("| $target | $route | $(result["n"]) | $(result["status"]) | $(result["bound"]) | $(timed.time) |")
end
function main()
    target,route,ns,os,variant,data,out,mode=ARGS[1:8]
    result=Dict{String,Any}("target"=>target,"route"=>route,"n"=>parse(Int,ns),"order"=>parse(Int,os),
        "variant"=>variant,"julia"=>string(VERSION),"cpu"=>Sys.cpu_info()[1].model,
        "threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads(),"stage"=>"input")
    checkpoint(result,out)
    t=@timed input_case(target,route,parse(Int,ns),parse(Int,os),variant,data)
    build=t.value; result["input_seconds"]=t.time
    result["stage"]="cold_build"; checkpoint(result,out)
    cold=@timed build()
    result["cold_seconds"]=cold.time; result["cold_bytes"]=cold.bytes
    if length(ARGS)==9 && ARGS[9]=="save_program"
        result["program_path"]=out*".program.jls"
        serialize(result["program_path"],cold.value.program)
    end
    if mode=="solve"
        merge!(result,dimensions(cold.value,route))
        solve_probe(cold.value,route,target,result,out)
        return
    end
    merge!(result,dimensions(cold.value,route)); cold=nothing; GC.gc()
    lower=b->N.moment_standard_form(b.program)
    times=Float64[]; bytes=Int[]; built=nothing
    for i in 1:3
        built=nothing; GC.gc(); result["stage"]="warm_build_$i"; checkpoint(result,out)
        sample=@timed build(); built=sample.value
        push!(times,sample.time); push!(bytes,sample.bytes); sample=nothing
        result["build_seconds"]=times; result["build_bytes"]=bytes; checkpoint(result,out)
    end
    if mode=="profile"
        built=nothing; GC.gc(); result["stage"]="profile"; checkpoint(result,out)
        repetitions=clamp(ceil(Int,3/median(times)),1,100)
        result["profile_builds"]=repetitions
        Profile.init(n=10^7,delay=0.002)
        Profile.clear()
        Profile.@profile for _ in 1:repetitions
            build()
        end
        open(out*".profile","w") do io
            Profile.print(io;format=:flat,sortedby=:count,mincount=5)
        end
        open(out*".tree","w") do io
            Profile.print(io;format=:tree,mincount=5)
        end
    else
        result["stage"]="cold_lower"; checkpoint(result,out)
        lowered=lower(built); dims=lower_dimensions(lowered,route)
        @assert result["variables"]==dims["lower_cols"]
        merge!(result,dims); lowered=nothing; GC.gc()
        lt=Float64[]; lb=Int[]
        for i in 1:3
            result["stage"]="warm_lower_$i"; checkpoint(result,out)
            sample=@timed lower(built)
            push!(lt,sample.time); push!(lb,sample.bytes); sample=nothing; GC.gc()
            result["lower_seconds"]=lt; result["lower_bytes"]=lb; checkpoint(result,out)
        end
    end
    result["stage"]="complete"; checkpoint(result,out)
    println("| target | route | n | order | option | warm s (median) | allocated bytes (median) | peak RSS |")
    println("| $target | $route | $ns | $os | $variant | $(median(times)) | $(median(bytes)) | $(Sys.maxrss()) |")
end
if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
