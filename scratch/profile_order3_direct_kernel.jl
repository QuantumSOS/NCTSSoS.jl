using Dates
using Profile
using SparseArrays
using Statistics
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

function build_problem()
    reg, (sx, sy, sz) = create_pauli_variables(1:N)
    bonds = Tuple{Int,Int}[]
    for i in 1:Lx, j in 1:Ly
        u = site(i, j)
        push!(bonds, (u, site(i + 1, j)))
        push!(bonds, (u, site(i, j + 1)))
    end
    H = sum(ComplexF64(1/4) * op[u] * op[v] for (u, v) in bonds for op in (sx, sy, sz))
    pop = polyopt(H, reg)

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
    return pop, SymmetrySpec(gens; offblock_check=:off)
end

function stage(msg)
    println("[", Dates.now(), "] ", msg)
    flush(stdout)
end

function row_nnz(A)
    At = SparseMatrixCSC(transpose(A))
    return [At.colptr[i + 1] - At.colptr[i] for i in 1:size(A, 1)]
end

function main()
    block_idx = parse(Int, get(ENV, "BLOCK_INDEX", "14"))
    profile_path = get(ENV, "PROFILE_OUT", "scratch/order3_direct_kernel_profile.txt")
    profile_flat_path = replace(profile_path, ".txt" => "_flat.txt")

    stage("build problem")
    pop, symmetry = build_problem()
    config = SolverConfig(
        optimizer=Clarabel.Optimizer,
        order=3,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry,
    )

    stage("compute_sparsity start")
    sparsity = compute_sparsity(pop, config)
    stage("compute_sparsity done")
    basis = sparsity.cliques_term_sparsities[1][1].block_bases[1]
    T = typeof(first(basis).word[1])
    MP = Polynomial{PauliAlgebra,T,ComplexF64}
    one_poly = convert(MP, one(pop.objective))

    stage("build Clifford group/reducer")
    support_domain = NCTSSoS._symmetry_domain(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)
    nqubits = max(
        NCTSSoS._clifford_nqubits_from_domain(support_domain),
        NCTSSoS._clifford_max_generator_nqubits(symmetry.clifford_generators),
    )
    lazy_generators = CliffordSymmetry{T}[
        NCTSSoS._convert_clifford_symmetry(T, g; nqubits) for g in symmetry.clifford_generators
    ]
    sw_group = CliffordSymmetryGroup(lazy_generators; nqubits, integer_type=T, domain=support_domain)
    group = CliffordSymmetry{T}[NCTSSoS._clifford_group_value(element) for element in sw_group]
    reducer = NCTSSoS._LazyOrbitReducer(PauliAlgebra, T, group, lazy_generators)

    stage("decompose basis")
    row_bases = NCTSSoS._symmetry_sparse_row_bases(basis, sw_group)
    sizes = [size(B, 1) for B in row_bases]
    nnzs = [nnz(B) for B in row_bases]
    row_stats = [(minimum(row_nnz(B)), median(row_nnz(B)), maximum(row_nnz(B))) for B in row_bases]
    println("basis size = ", length(basis))
    println("group order = ", length(group))
    println("block sizes = ", sizes)
    println("block nnz = ", nnzs)
    println("row nnz min/median/max = ", row_stats)
    flush(stdout)

    B = row_bases[block_idx]
    stage("profile block $block_idx size=$(size(B,1)) nnz=$(nnz(B))")
    Profile.clear()
    Profile.init(n=20_000_000, delay=0.01)
    elapsed = @elapsed begin
        @profile begin
            block = NCTSSoS._transform_constraint_block_from_basis(one_poly, basis, B, reducer, MP; hermitian=true)
            println("constructed block size = ", size(block), ", nonzero entries = ", count(!iszero, block))
        end
    end
    stage("profiled block done elapsed=$elapsed seconds reducer_cache=$(length(reducer.reduction))")

    open(profile_path, "w") do io
        Profile.print(io; sortedby=:count, maxdepth=80, mincount=20)
    end
    open(profile_flat_path, "w") do io
        Profile.print(io; format=:flat, sortedby=:count, maxdepth=80, mincount=20)
    end
    println("profile samples = ", length(Profile.fetch()))
    println("profile_path = ", profile_path)
    println("profile_flat_path = ", profile_flat_path)
    stage("done")
end

main()
