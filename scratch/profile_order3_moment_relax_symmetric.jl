using Dates
using Profile
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

function main()
    println("Started: ", Dates.now())
    pop, symmetry = build_problem()
    config = SolverConfig(
        optimizer=Clarabel.Optimizer,
        order=3,
        cs_algo=NoElimination(),
        ts_algo=NoElimination(),
        symmetry=symmetry,
    )

    println("compute_sparsity start: ", Dates.now())
    sparsity = compute_sparsity(pop, config)
    println("compute_sparsity done: ", Dates.now())
    println("sizes = ", NCTSSoS._compute_moment_matrix_sizes(sparsity.cliques_term_sparsities))
    NCTSSoS._check_symmetry_mvp_support(pop, config, sparsity)

    Profile.clear()
    Profile.init(n=20_000_000, delay=0.01)
    println("profiled moment_relax_symmetric start: ", Dates.now())
    flush(stdout)

    try
        @profile begin
            NCTSSoS.moment_relax_symmetric(
                pop,
                sparsity.corr_sparsity,
                sparsity.cliques_term_sparsities,
                symmetry,
            )
        end
        println("moment_relax_symmetric completed: ", Dates.now())
    catch err
        println("moment_relax_symmetric interrupted/error: ", typeof(err), " ", sprint(showerror, err))
    finally
        println("writing profile: ", Dates.now())
        open("scratch/order3_moment_relax_profile.txt", "w") do io
            Profile.print(io; sortedby=:count, maxdepth=80, mincount=20)
        end
        open("scratch/order3_moment_relax_profile_flat.txt", "w") do io
            Profile.print(io; format=:flat, sortedby=:count, maxdepth=80, mincount=20)
        end
        println("profile samples = ", length(Profile.fetch()))
        println("Finished: ", Dates.now())
    end
end

main()
