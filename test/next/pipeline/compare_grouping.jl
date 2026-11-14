# Standalone order-2 grouping comparison. Requires a working Mosek license.
using NCTSSoS, JuMP, MosekTools
const Nx=NCTSSoS.Next
optimizer=optimizer_with_attributes(MosekTools.Optimizer,"MSK_IPAR_LOG"=>0)
check=Model(optimizer)
@variable(check,v>=1)
@objective(check,Min,v)
optimize!(check)
@assert termination_status(check)==JuMP.MOI.OPTIMAL && abs(objective_value(check)-1)<1e-8

registry,ops=create_pauli_variables(1:4)
H=heisenberg_chain_hamiltonian(ops;periodic=true,coupling=1)
dense=cs_nctssos(polyopt(H,registry),SolverConfig(optimizer=optimizer,order=2);dualize=false)
println("DENSE current package bound=",dense.objective)
registry,(x,y,z)=Nx.pauli_variables(4)
problem=Nx.Problem(sum(x[i]*x[mod1(i+1,4)]+y[i]*y[mod1(i+1,4)]+z[i]*z[mod1(i+1,4)] for i in 1:4))
for (name,symmetry) in (("SITE MF",Nx.FamilySymmetry()),
    ("SITE MF + TRANSLATION",Nx.FamilySymmetry([Nx.SitePermutation([2,3,4,1])])))
    description=Nx.BlockRelaxation(order=2,grouping=Nx.Correlative(),symmetry=symmetry)
    result=Nx.solve(problem,description;solver=Nx.JuMPBackend(optimizer))
    closure=result.sizes.closure
    println(name," bound=",result.bound," initial=",closure.initial_blocks," closed=",closure.closed_blocks,
        " orbits=",length(closure.block_orbits)," largest=",result.sizes.program.largest_block)
end
