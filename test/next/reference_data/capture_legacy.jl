# Run before deleting the legacy route; afterwards use the pinned old checkout.
# julia --startup-file=no --project=<solver environment> test/next/reference_data/capture_legacy.jl <output.toml>
using NCTSSoS, JuMP, MosekTools, TOML

const solver = optimizer_with_attributes(MosekTools.Optimizer,"QUIET"=>true)
const captured = Dict{String,Any}[]
function capture_bound(id,h,r,order;cs=NoElimination(),ts=NoElimination())
    result = cs_nctssos(polyopt(h,r),SolverConfig(optimizer=solver,order=order,cs_algo=cs,ts_algo=ts);dualize=false)
    termination_status(result.model) == MOI.OPTIMAL || error("$id did not reach OPTIMAL")
    push!(captured,Dict("id"=>id,"expected"=>Dict("objective"=>result.objective)))
    println(id," ",result.objective)
end

r,(x,y,z) = create_pauli_variables(1:4)
capture_bound("heisenberg-open4-order1",sum(x[i]*x[i+1]+y[i]*y[i+1]+z[i]*z[i+1] for i in 1:3),r,1)
r,(a,ad) = create_fermionic_variables(1:3)
h = -(ad[1]*a[2]+ad[2]*a[1])-.7*(ad[2]*a[3]+ad[3]*a[2])+1.3*ad[1]*a[1]*ad[2]*a[2]+.4*ad[3]*a[3]
capture_bound("fermion3-order2",h,r,2)
r,(a,ad) = create_fermionic_variables(1:4)
h = -(ad[1]*a[3]+ad[3]*a[1]+ad[2]*a[4]+ad[4]*a[2])+2*(ad[1]*a[1]*ad[2]*a[2]+ad[3]*a[3]*ad[4]*a[4])
capture_bound("hubbard2-order2",h,r,2)
r,(x,y) = create_unipotent_variables([("x",1:2),("y",1:2)])
capture_bound("chsh-order1",-(x[1]*y[1]+x[1]*y[2]+x[2]*y[1]-x[2]*y[2]),r,1)
r,(b,bd) = create_bosonic_variables(1:2)
h = -(bd[1]*b[2]+bd[2]*b[1])+sum(bd[i]^2*b[i]^2-0.5bd[i]*b[i] for i in 1:2)
for order in (2,3)
    capture_bound("boson-bh-dense-order$order",h,r,order)
end
capture_bound("boson-bh-correlative-order2",h,r,2;cs=MF())
capture_bound("boson-bh-split-order2",h,r,2;ts=MaximalElimination())

# Capture the actual clique/row sets, not only their sizes.
r,ops = create_pauli_variables(1:4)
h = heisenberg_chain_hamiltonian(ops;periodic=true,coupling=1)
for algorithm in (MF(),MMD())
    old = NCTSSoS.correlative_sparsity(polyopt(h,r),1,algorithm)
    groups = [sort!(unique([(Int(i)-1)÷3+1 for i in clique])) for clique in old.cliques]
    push!(captured,Dict("id"=>"pauli-ring4-groups-$(nameof(typeof(algorithm)))",
        "expected"=>Dict("groups"=>groups)))
end
basis = get_ncbasis(r,1)
active = NCTSSoS.init_activated_supp(h,typeof(h)[],basis)
split = NCTSSoS.term_sparsities(active,typeof(h)[],basis,typeof(basis)[],MaximalElimination())
push!(captured,Dict("id"=>"pauli-ring4-term-rows",
    "expected"=>Dict("blocks"=>[[Int.(w.word) for w in block] for block in first(split).block_bases])))

model = Model(solver)
version = MOI.get(backend(model),MOI.SolverVersion())
data = Dict("schema_version"=>1,
    "old_commit"=>"fa3484eb06d962e78b7b562788e09496fe6497bc",
    "solver_version"=>"Mosek $version; MosekTools $(pkgversion(MosekTools))",
    "command"=>"OPENBLAS_NUM_THREADS=1 julia --startup-file=no --project=<solver environment> test/next/reference_data/capture_legacy.jl <output.toml>",
    "notes"=>"Legacy-route captures before removal. Runtime src and Project.toml equal fa3484e. Numeric bounds are solver results, not exact certificates. Clique and row sets are implementation reference values.",
    "cases"=>captured)
open(io->TOML.print(io,data;sorted=true),only(ARGS),"w")
