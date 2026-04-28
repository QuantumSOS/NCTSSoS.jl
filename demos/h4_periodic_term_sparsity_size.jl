# Demo: periodic H4 (Nk = 2) term-sparsity footprint at order 2
#
# Purpose:
#   Build the exact 32-mode Bloch Hamiltonian for the vendored H4 periodic asset,
#   apply term sparsity at order 2, and print the actual SDP block footprint.
#
# Notes:
#   - This is a symbolic size probe, not an energy solve.
#   - It can take a couple of minutes and allocates heavily (~10 GiB transiently)
#     because the term-sparsity graph itself is still large.
#   - We use `ts_algo = MMD()` because that is the standard chordal TS path in
#     the examples and source.

using NCTSSoS

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: fermionic_order2_basis_size,
                         fermionic_order2_nuniq,
                         load_nk2_asset

function build_h4_nk2_hamiltonian(h1e, eri; nk::Int, norb::Int)
    registry, ((c_up_k0, c_up_k0_dag),
               (c_dn_k0, c_dn_k0_dag),
               (c_up_k1, c_up_k1_dag),
               (c_dn_k1, c_dn_k1_dag)) = create_fermionic_variables([
        ("c_up_k0", 1:norb),
        ("c_dn_k0", 1:norb),
        ("c_up_k1", 1:norb),
        ("c_dn_k1", 1:norb),
    ])

    ann = Dict(
        (0, :up) => c_up_k0,
        (0, :dn) => c_dn_k0,
        (1, :up) => c_up_k1,
        (1, :dn) => c_dn_k1,
    )
    dag = Dict(
        (0, :up) => c_up_k0_dag,
        (0, :dn) => c_dn_k0_dag,
        (1, :up) => c_up_k1_dag,
        (1, :dn) => c_dn_k1_dag,
    )
    spin_channels = ((:up, :up), (:up, :dn), (:dn, :up), (:dn, :dn))

    ham = (0.0 + 0.0im) * (c_up_k0_dag[1] * c_up_k0[1])

    # One-body part: diagonal in k, diagonal in spin, scaled by 1 / nk.
    for k in 0:nk-1
        h_block = h1e[k] / nk
        for p in 1:norb, s in 1:norb, sigma in (:up, :dn)
            coeff = h_block[p, s]
            iszero(coeff) && continue
            ham += coeff * dag[(k, sigma)][p] * ann[(k, sigma)][s]
        end
    end

    # Two-body part: only the 8 momentum-conserving ERI blocks survive.
    for (k1, k2, k3, k4) in sort!(collect(keys(eri)))
        block = eri[(k1, k2, k3, k4)] / nk^2
        for p in 1:norb, r in 1:norb, q in 1:norb, s in 1:norb
            coeff = block[p, r, q, s]
            iszero(coeff) && continue
            for (sigma, tau) in spin_channels
                ham += 0.5 * coeff *
                    dag[(k1, sigma)][p] * dag[(k2, tau)][q] *
                    ann[(k4, tau)][s] * ann[(k3, sigma)][r]
            end
        end
    end

    # Keep roundoff asymmetry out of the symbolic pipeline.
    return registry,
           ((c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
            (c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag)),
           0.5 * (ham + adjoint(ham))
end

@inline mode_k(mode::Int) = mode <= 16 ? 0 : 1

function basis_label(m::NormalMonomial)
    ΔN = 0
    K = 0
    for idx in m.word
        sign = idx < 0 ? 1 : -1
        ΔN += sign
        K = mod(K + sign * mode_k(abs(Int(idx))), 2)
    end
    return (ΔN = ΔN, K = K, deg = length(m.word))
end

function block_labelset(block)
    labels = Set(basis_label(m) for m in block)
    return Tuple(sort!(collect(labels); by = x -> (x.ΔN, x.K, x.deg)))
end

function summarize_blocks(blocks)
    summary = Dict{Tuple, Vector{Int}}()
    for block in blocks
        key = block_labelset(block)
        push!(get!(summary, key, Int[]), length(block))
    end
    return summary
end

function dense_sector_counts(nmodes_total::Int, nk::Int)
    nmodes_per_k = nmodes_total ÷ nk
    @assert nk == 2

    return Dict(
        (ΔN = -2, K = 0, deg = 2) => 2 * binomial(nmodes_per_k, 2),
        (ΔN = -2, K = 1, deg = 2) => nmodes_per_k^2,
        (ΔN = -1, K = 0, deg = 1) => nmodes_per_k,
        (ΔN = -1, K = 1, deg = 1) => nmodes_per_k,
        (ΔN = 0, K = 0, deg = 0) => 1,
        (ΔN = 0, K = 0, deg = 2) => 2 * nmodes_per_k^2,
        (ΔN = 0, K = 1, deg = 2) => 2 * nmodes_per_k^2,
        (ΔN = 1, K = 0, deg = 1) => nmodes_per_k,
        (ΔN = 1, K = 1, deg = 1) => nmodes_per_k,
        (ΔN = 2, K = 0, deg = 2) => 2 * binomial(nmodes_per_k, 2),
        (ΔN = 2, K = 1, deg = 2) => nmodes_per_k^2,
    )
end

function print_sector_line(io, name, dense_count, sizes)
    println(io,
        rpad(name, 34),
        " dense=", lpad(dense_count, 4),
        " | blocks=", lpad(length(sizes), 3),
        " | size range=", lpad(minimum(sizes), 3), "..", lpad(maximum(sizes), 3))
end

function print_overlap_line(io, name, sizes)
    println(io,
        rpad(name, 34),
        " blocks=", lpad(length(sizes), 3),
        " | size range=", lpad(minimum(sizes), 3), "..", lpad(maximum(sizes), 3))
end

asset = load_nk2_asset()
registry, vars, ham = build_h4_nk2_hamiltonian(asset.h1e, asset.eri; nk = asset.nk, norb = asset.n_active_orb)
(c_up_k0, c_up_k0_dag), (c_dn_k0, c_dn_k0_dag),
(c_up_k1, c_up_k1_dag), (c_dn_k1, c_dn_k1_dag) = vars

n_up = 1.0 * sum(c_up_k0_dag[i] * c_up_k0[i] + c_up_k1_dag[i] * c_up_k1[i] for i in 1:asset.n_active_orb)
n_dn = 1.0 * sum(c_dn_k0_dag[i] * c_dn_k0[i] + c_dn_k1_dag[i] * c_dn_k1[i] for i in 1:asset.n_active_orb)

pop = polyopt(ham, registry;
    moment_eq_constraints = [n_up - 4.0 * one(ham),
                             n_dn - 4.0 * one(ham)])

config = SolverConfig(optimizer = nothing, order = 2, ts_algo = MMD())
println("Computing order-2 term sparsity with MMD()...")
@time sparsity = compute_sparsity(pop, config)

blocks = sparsity.cliques_term_sparsities[1][1].block_bases
block_sizes = length.(blocks)
block_summary = summarize_blocks(blocks)
sector_counts = dense_sector_counts(asset.total_spin_orbital_modes, asset.nk)

all_blocks_k_pure = all(length(Set(label.K for label in block_labelset(block))) == 1 for block in blocks)
@assert all_blocks_k_pure

moment_problem = NCTSSoS.moment_relax(pop, sparsity.corr_sparsity, sparsity.cliques_term_sparsities)

println()
println("=== Dense order-2 baseline ===")
println("spin-orbital modes:                 ", asset.total_spin_orbital_modes)
println("dense moment basis size:            ", Int(fermionic_order2_basis_size(asset.total_spin_orbital_modes)))
println("dense PSD block:                    2081 x 2081")
println("dense unique order-≤4 moments:      ", Int(fermionic_order2_nuniq(asset.total_spin_orbital_modes)))
println("paper-style pair slice:             K=0 => 240, K=1 => 256")

println()
println("=== TS first pass (order 2, ts_algo = MMD()) ===")
println("correlative cliques:                ", length(sparsity.corr_sparsity.cliques))
println("term-sparsity PSD blocks:           ", length(blocks))
println("largest PSD block:                  ", maximum(block_sizes), " x ", maximum(block_sizes))
println("sum(block_size^2) proxy:            ", sum(size^2 for size in block_sizes))
println("unique moment variables:            ", moment_problem.n_unique_moment_matrix_elements)
k0_sizes = [size for (size, block) in zip(block_sizes, blocks) if first(block_labelset(block)).K == 0]
k1_sizes = [size for (size, block) in zip(block_sizes, blocks) if first(block_labelset(block)).K == 1]
println("all PSD blocks are K-pure:          ", all_blocks_k_pure)
println("K=0 blocks / largest:               ", length(k0_sizes), " / ", maximum(k0_sizes))
println("K=1 blocks / largest:               ", length(k1_sizes), " / ", maximum(k1_sizes))

println()
println("=== Sector view ===")
print_sector_line(stdout, "ΔN = -2, K = 0", sector_counts[(ΔN = -2, K = 0, deg = 2)], block_summary[((ΔN = -2, K = 0, deg = 2),)])
print_sector_line(stdout, "ΔN = -2, K = 1", sector_counts[(ΔN = -2, K = 1, deg = 2)], block_summary[((ΔN = -2, K = 1, deg = 2),)])
print_sector_line(stdout, "ΔN = -1, K = 0", sector_counts[(ΔN = -1, K = 0, deg = 1)], block_summary[((ΔN = -1, K = 0, deg = 1),)])
print_sector_line(stdout, "ΔN = -1, K = 1", sector_counts[(ΔN = -1, K = 1, deg = 1)], block_summary[((ΔN = -1, K = 1, deg = 1),)])
print_sector_line(stdout, "ΔN =  0, K = 0 (deg 2 only)", sector_counts[(ΔN = 0, K = 0, deg = 2)], block_summary[((ΔN = 0, K = 0, deg = 2),)])
print_overlap_line(stdout, "identity overlaps inside K = 0", block_summary[((ΔN = 0, K = 0, deg = 0), (ΔN = 0, K = 0, deg = 2))])
print_sector_line(stdout, "ΔN =  0, K = 1", sector_counts[(ΔN = 0, K = 1, deg = 2)], block_summary[((ΔN = 0, K = 1, deg = 2),)])
print_sector_line(stdout, "ΔN = +1, K = 0", sector_counts[(ΔN = 1, K = 0, deg = 1)], block_summary[((ΔN = 1, K = 0, deg = 1),)])
print_sector_line(stdout, "ΔN = +1, K = 1", sector_counts[(ΔN = 1, K = 1, deg = 1)], block_summary[((ΔN = 1, K = 1, deg = 1),)])
print_sector_line(stdout, "ΔN = +2, K = 0", sector_counts[(ΔN = 2, K = 0, deg = 2)], block_summary[((ΔN = 2, K = 0, deg = 2),)])
print_sector_line(stdout, "ΔN = +2, K = 1", sector_counts[(ΔN = 2, K = 1, deg = 2)], block_summary[((ΔN = 2, K = 1, deg = 2),)])

println()
println("Interpretation:")
println("- TS does recover the momentum symmetry: no PSD block mixes K = 0 and K = 1.")
println("- The paper's 2-RDM K sectors (240 / 256 pair rows) appear as the ΔN = ±2 slices,")
println("  but the generic TS chordalization shards each sector into many smaller PSD blocks.")
println("- So the right mental model is: TS sees the symmetry in moment-matrix coordinates,")
println("  not as one hand-written 240x240 or 256x256 2-RDM block.")
