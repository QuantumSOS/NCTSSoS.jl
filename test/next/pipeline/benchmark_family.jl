# Same seven-pair protocol as benchmark_chain.jl, including the model entry
# and row generation. Record cold cost as well as warm cost; no solves.
include("benchmark_chain.jl")

function model_chain_case(n)
    registry,(x,y,z)=Nx.pauli_variables(n)
    H=sum(x[i]*x[mod1(i+1,n)]+y[i]*y[mod1(i+1,n)]+z[i]*z[mod1(i+1,n)] for i in 1:n)
    _,_,translation,point,anti=chain_case(n)
    description=Nx.BlockRelaxation(order=4,row_rule=Nx.Contiguous(periodic=true),
        symmetry=Nx.FamilySymmetry([translation;point];moment_only=(anti,)))
    Nx.Problem(H),description
end

if abspath(PROGRAM_FILE)==@__FILE__
    println("ENV Julia=",VERSION," threads=",Threads.nthreads()," BLAS=",BLAS.get_num_threads()," CPU=",Sys.cpu_info()[1].model)
    for n in (16,32)
        problem,description=model_chain_case(n)
        registry,ops=create_pauli_variables(1:n)
        pop=polyopt(heisenberg_chain_hamiltonian(ops;periodic=true,coupling=1),registry)
        comparison(n,(("family",()->begin
            built=Nx.build(problem,description)
            (built.program,built.table,built.sizes)
        end),("legacy",()->pauli_translation_invariant_moment_relaxation(pop,ops,4;axis_permutation_symmetry=true))))
    end
end
