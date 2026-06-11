using Dates
using SparseArrays
using NCTSSoS
import Clarabel

const Lx = 4
const Ly = 4
const N = Lx * Ly
const LI = LinearIndices((1:Lx, 1:Ly))
site(i, j) = LI[CartesianIndex(mod1(i, Lx), mod1(j, Ly))]
coord(s) = Tuple(CartesianIndices((1:Lx, 1:Ly))[s])

function site_permutation_clifford(sx, sy, sz, perm::Vector{Int})
    M = typeof(sx[1])
    images = Dict{M,Tuple{Int,M}}()
    for i in 1:N
        j = perm[i]
        images[sx[i]] = (1, sx[j])
        images[sy[i]] = (1, sy[j])
        images[sz[i]] = (1, sz[j])
    end
    return CliffordSymmetry(images; nqubits=N)
end

println("Started ", Dates.now())
reg, (sx, sy, sz) = create_pauli_variables(1:N)
bonds = Tuple{Int,Int}[]
for i in 1:Lx, j in 1:Ly
    u = site(i, j)
    push!(bonds, (u, site(i + 1, j)))
    push!(bonds, (u, site(i, j + 1)))
end
H = sum(ComplexF64(1/4) * op[u] * op[v] for (u, v) in bonds for op in (sx, sy, sz))
pop = polyopt(H, reg)
config = SolverConfig(optimizer=Clarabel.Optimizer, order=3, cs_algo=NoElimination(), ts_algo=NoElimination())
sp = compute_sparsity(pop, config)
basis = only(sp.corr_sparsity.clq_mom_mtx_bases)
println("basis len=", length(basis))

tx = [site(coord(s)[1] + 1, coord(s)[2]) for s in 1:N]
ty = [site(coord(s)[1], coord(s)[2] + 1) for s in 1:N]
rot90 = [site(coord(s)[2], Lx + 1 - coord(s)[1]) for s in 1:N]
reflx = [site(Lx + 1 - coord(s)[1], coord(s)[2]) for s in 1:N]
gens = CliffordSymmetry[
    site_permutation_clifford(sx, sy, sz, tx),
    site_permutation_clifford(sx, sy, sz, ty),
    site_permutation_clifford(sx, sy, sz, rot90),
    site_permutation_clifford(sx, sy, sz, reflx),
]
sym = SymmetrySpec(gens; offblock_check=:off)
support_domain = NCTSSoS._symmetry_domain(pop, sp.corr_sparsity, sp.cliques_term_sparsities)
nqubits = max(NCTSSoS._clifford_nqubits_from_domain(support_domain), NCTSSoS._clifford_max_generator_nqubits(sym.clifford_generators))
sw_group = CliffordSymmetryGroup(sym.clifford_generators; nqubits, integer_type=UInt8, domain=support_domain)
println("group order=", length(sw_group))
flush(stdout)

t = @elapsed blocks = NCTSSoS._sw_decompose_half_basis(basis, sw_group)
println("decompose seconds=", t)
println("nblocks=", length(blocks))
for (k,b) in enumerate(blocks)
    m,n = size(b)
    nz = b isa SparseMatrixCSC ? nnz(b) : count(!iszero, b)
    println(k, ": type=", typeof(b), " size=", (m,n), " nnz=", nz, " density=", nz/(m*n))
end
println("Finished ", Dates.now())
