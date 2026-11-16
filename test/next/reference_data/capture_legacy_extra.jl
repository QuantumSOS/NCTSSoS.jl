# Historical-only capture: run against fa3484e before removal, or its pinned checkout.
using NCTSSoS, JuMP, MosekTools, TOML, SHA
const solver = optimizer_with_attributes(MosekTools.Optimizer,"QUIET"=>true)
const captured = Dict{String,Any}[]

function bell(a,b)
    n=length(a)
    c=if n==3
        [0 -1 -2 0; 0 1 1 -1; -1 1 1 1; 0 -1 1 0]
    else
        matrix=zeros(Int,n+1,n+1)
        for i in 1:n-1
            matrix[1,i+1]=matrix[i+1,1]=-1
        end
        for i in 1:n,j in 1:n
            matrix[i+1,j+1]=(i==n || j==n) ? (i==j ? -1 : 1) :
                i+j<n ? 0 : i+j==n ? 1 : -1
        end
        matrix
    end
    sum(c[i+1,1]*a[i]+c[1,i+1]*b[i] for i in 1:n)+
        sum(c[i+1,j+1]*a[i]*b[j] for i in 1:n,j in 1:n)
end

r,(a,b)=create_projector_variables([("A",1:5),("B",1:5)])
old=NCTSSoS.correlative_sparsity(polyopt(-bell(a,b),r),5,MMD())
row_digest(rows)=bytes2hex(sha256(join(sort!([join(Int.(w.word),",") for w in rows]),";")))
push!(captured,Dict("id"=>"inn22-n5-order5-mmd","expected"=>Dict(
    "letters"=>[Int(only(p.word)) for p in [a;b]],
    "cliques"=>[Int.(c) for c in old.cliques],
    "row_counts"=>length.(old.clq_mom_mtx_bases),
    "row_sha256"=>row_digest.(old.clq_mom_mtx_bases))))

r,(a,b)=create_projector_variables([("A",1:3),("B",1:3)])
old=cs_nctssos(polyopt(-bell(a,b),r),SolverConfig(optimizer=solver,order=3,
    cs_algo=MMD(),ts_algo=MaximalElimination());dualize=false)
termination_status(old.model)==MOI.OPTIMAL || error("I3322 solve did not reach OPTIMAL")
push!(captured,Dict("id"=>"i3322-order3-mmd-term","expected"=>Dict("objective"=>old.objective)))
println("I3322 ",old.objective)

h=ComplexF64[-0.4 -0.7+0.2im;-0.7-0.2im 0.2]
eri=zeros(ComplexF64,2,2,2,2); eri[1,1,1,1]=4; eri[2,2,2,2]=3.5
linear=build_pqg_moment_data(Dict(0=>h),Dict((0,0,0,0)=>eri);
    nk=1,norb=2,nelec_per_cell=2,blocking=:none,spin_resolved_trace=false,
    singlet_s2=false,include_one_d=true,three_index=:T1T2prime)
M=typeof(first(values(linear.key_to_monomial))); T=eltype(linear.identity)
P=Polynomial{FermionicAlgebra,T,ComplexF64}
mp=NCTSSoS._moment_problem_with_linear(FermionicAlgebra,T,M,P,zero(P),Tuple{Symbol,Matrix{P}}[],M[],0,linear)
old=NCTSSoS.sos_dualize(mp)
set_optimizer(old.model,solver); optimize!(old.model)
termination_status(old.model)==MOI.OPTIMAL || error("T2prime solve did not reach OPTIMAL")
objective=[Dict("word"=>Int.(word),"real"=>real(c),"imag"=>imag(c))
    for (word,c) in sort!(collect(linear.objective_lin);by=p->Tuple(first(p)))]
push!(captured,Dict("id"=>"fermionic-t2prime","expected"=>Dict(
    "objective"=>objective_value(old.model),"terms"=>objective)))
println("T2prime ",objective_value(old.model))

for trace in (false,true), (label,signs) in (("docs",[1 1 1;1 1 -1;1 -1 0]),("bench",[1 1 -1;1 1 1;-1 1 0]))
    if trace
        r,(x,)=create_unipotent_variables([("x",1:6)])
        A,B=x[1:3],x[4:6]; symbol=NCTSSoS.tr
    else
        r,(A,B)=create_unipotent_variables([("A",1:3),("B",1:3)])
        symbol=ς
    end
    h=-sum(signs[i,j]*(symbol(A[i]*B[j])-symbol(A[i])*symbol(B[j])) for i in 1:3,j in 1:3 if signs[i,j]!=0)
    old=cs_nctssos(polyopt(h*one(A[1]),r),SolverConfig(optimizer=solver,order=2,
        cs_algo=NoElimination(),ts_algo=NoElimination());dualize=false)
    termination_status(old.model)==MOI.OPTIMAL || error("Covariance solve did not reach OPTIMAL")
    id="covariance-$(trace ? "trace" : "state")-$label"
    push!(captured,Dict("id"=>id,"expected"=>Dict("objective"=>old.objective)))
    println(id," ",old.objective)
end

version=MOI.get(backend(Model(solver)),MOI.SolverVersion())
data=Dict("schema_version"=>1,"old_commit"=>"fa3484eb06d962e78b7b562788e09496fe6497bc",
    "solver_version"=>"Mosek $version; MosekTools $(pkgversion(MosekTools))",
    "command"=>"OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=<solver environment> test/next/reference_data/capture_legacy_extra.jl <output.toml>",
    "notes"=>"Legacy bounds and basis sets. Row digests encode sorted comma-separated legacy integer words joined by semicolons; empty word is an empty string. Frozen T2prime terms preserve the old structured Hamiltonian, including its complex coefficients.",
    "cases"=>captured)
open(io->TOML.print(io,data;sorted=true),only(ARGS),"w")
