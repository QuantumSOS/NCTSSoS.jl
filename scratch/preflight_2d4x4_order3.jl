using Dates
using NCTSSoS
import Clarabel

const Lx = 4
const Ly = 4
const N = Lx * Ly
const LI = LinearIndices((1:Lx, 1:Ly))
site(i, j) = LI[CartesianIndex(mod1(i, Lx), mod1(j, Ly))]

function square_lattice_hamiltonian()
    reg, (sx, sy, sz) = create_pauli_variables(1:N)
    bonds = Tuple{Int,Int}[]
    for i in 1:Lx, j in 1:Ly
        u = site(i, j)
        push!(bonds, (u, site(i + 1, j)))
        push!(bonds, (u, site(i, j + 1)))
    end
    H = sum(ComplexF64(1/4) * op[u] * op[v] for (u, v) in bonds for op in (sx, sy, sz))
    return reg, H
end

println("Started: ", Dates.now())
println("4x4 Heisenberg order-3 dense/no-CS/no-TS preflight")
reg, H = square_lattice_hamiltonian()
println("Pauli terms = ", length(terms(H)))
pop = polyopt(H, reg)
config = SolverConfig(
    optimizer=Clarabel.Optimizer,
    order=3,
    cs_algo=NoElimination(),
    ts_algo=NoElimination(),
)
flush(stdout)
elapsed = @elapsed sparsity = compute_sparsity(pop, config)
println("compute_sparsity seconds = ", elapsed)
println("cliques = ", length(sparsity.corr_sparsity.cliques))
println("clique variable sizes = ", length.(sparsity.corr_sparsity.cliques))
println("moment matrix sizes = ", NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities))
println("basis length = ", length(only(sparsity.corr_sparsity.clq_mom_mtx_bases)))
println("Finished: ", Dates.now())
