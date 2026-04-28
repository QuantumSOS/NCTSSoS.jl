# # [Periodic V2RDM: Why Correlative Sparsity Misses the $k$-Blocks](@id periodic-v2rdm-cs-k-blocks)
#
# This page builds the exact `N_k = 2` H₄ Bloch Hamiltonian that we hand to
# [`polyopt`](@ref), then inspects the support graph that correlative sparsity
# actually sees.
#
# That distinction is the whole point:
#
# - the periodic V2RDM paper blocks the **2-RDM** by total pair momentum `K`,
# - NCTSSoS.jl builds the **second-quantized Hamiltonian** over single fermionic modes,
# - correlative sparsity only sees which single modes co-occur in some monomial,
# - for the ab-initio H₄ asset, that graph is complete.
#
# So yes, momentum-violating coefficients are zero. But omitting those forbidden
# terms is only *correctness*. It is not enough to recover the paper's `K`
# blocks from the Hamiltonian support graph.
#
# !!! note "What this page is and is not"
#     This is a documentation page about the **Hamiltonian that NCTSSoS.jl
#     builds**. It does **not** build the pair-indexed 2-RDM blocks used in a
#     hand-written periodic V2RDM code. That difference is exactly why the
#     correlative-sparsity story changes.

# ## Load the vendored `N_k = 2` H₄ asset
#
# The repo already vendors the reviewed PySCF dump used by the regression tests.
# The tiny loader below lives under `test/`; convenient for docs, but not public API.

using NCTSSoS

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: load_nk2_asset

asset = load_nk2_asset()
nothing #hide
# `asset` is a named tuple with the reviewed `h1e` / `eri` blocks and metadata.

nk = asset.nk
norb = asset.n_active_orb
allowed_tuples = sort!(collect(keys(asset.eri)))
# `allowed_tuples` are the 8 momentum-conserving ERI blocks that survive at `N_k = 2`.

allowed_tuples

@assert allowed_tuples == [
    (0, 0, 0, 0),
    (0, 0, 1, 1),
    (0, 1, 0, 1),
    (0, 1, 1, 0),
    (1, 0, 0, 1),
    (1, 0, 1, 0),
    (1, 1, 0, 0),
    (1, 1, 1, 1),
]

# ## The operator that NCTSSoS.jl actually builds
#
# In NCTSSoS.jl we do **not** hand the solver a block-diagonal 2-RDM variable.
# We hand it a fermionic polynomial
# ```math
# \hat H = \hat H_1 + \hat H_2
# ```
# over 32 single modes:
#
# - 8 active spatial orbitals,
# - 2 spins,
# - 2 `k`-points.
#
# The one-body part is
# ```math
# \hat H_1
# = \frac{1}{N_k}
#   \sum_{k \in \{0,1\}}
#   \sum_{p,s = 1}^{8}
#   \sum_{\sigma \in \{\uparrow,\downarrow\}}
#   h^{(k)}_{ps}\; a^{\dagger}_{p k \sigma} a_{s k \sigma},
# ```
# so one-body terms are diagonal in `k`, but not necessarily in the orbital labels.
#
# The two-body part is
# ```math
# \hat H_2
# = \frac{1}{2 N_k^2}
#   \sum_{(k_1,k_2,k_3,k_4) \in \mathcal A}
#   \sum_{p,q,r,s = 1}^{8}
#   \sum_{\sigma,\tau \in \{\uparrow,\downarrow\}}
#   \langle p k_1, q k_2 \mid r k_3, s k_4 \rangle\;
#   a^{\dagger}_{p k_1 \sigma}
#   a^{\dagger}_{q k_2 \tau}
#   a_{s k_4 \tau}
#   a_{r k_3 \sigma},
# ```
# where the allowed tuple set is exactly
# ```math
# \mathcal A = \{
# (0,0,0,0), (0,0,1,1), (0,1,0,1), (0,1,1,0),
# (1,0,0,1), (1,0,1,0), (1,1,0,0), (1,1,1,1)
# \}.
# ```
#
# For `N_k = 2`, the prefactor on the quartic part is `1 / 8`, so the operator
# is literally the sum of these eight ERI blocks:
# ```math
# \begin{aligned}
# \hat H_2 = \frac{1}{8}\sum_{p q r s}\sum_{\sigma,\tau} (&
# V^{0000}_{p q, r s}\; a^{\dagger}_{p0\sigma} a^{\dagger}_{q0\tau} a_{s0\tau} a_{r0\sigma}
# + V^{0011}_{p q, r s}\; a^{\dagger}_{p0\sigma} a^{\dagger}_{q0\tau} a_{s1\tau} a_{r1\sigma} \\
# &+ V^{1100}_{p q, r s}\; a^{\dagger}_{p1\sigma} a^{\dagger}_{q1\tau} a_{s0\tau} a_{r0\sigma}
# + V^{1111}_{p q, r s}\; a^{\dagger}_{p1\sigma} a^{\dagger}_{q1\tau} a_{s1\tau} a_{r1\sigma} \\
# &+ V^{0101}_{p q, r s}\; a^{\dagger}_{p0\sigma} a^{\dagger}_{q1\tau} a_{s1\tau} a_{r0\sigma}
# + V^{0110}_{p q, r s}\; a^{\dagger}_{p0\sigma} a^{\dagger}_{q1\tau} a_{s0\tau} a_{r1\sigma} \\
# &+ V^{1001}_{p q, r s}\; a^{\dagger}_{p1\sigma} a^{\dagger}_{q0\tau} a_{s1\tau} a_{r0\sigma}
# + V^{1010}_{p q, r s}\; a^{\dagger}_{p1\sigma} a^{\dagger}_{q0\tau} a_{s0\tau} a_{r1\sigma}) .
# \end{aligned}
# ```
#
# That last line is the one people tend to miss. Terms like `0101`, `0110`,
# `1001`, and `1010` are **momentum-conserving**, so they stay in the Hamiltonian.
# They already contain both `k = 0` and `k = 1` modes in the same monomial.
#
# !!! warning "What you must not do"
#     Do **not** zero those cross-sector terms by hand just to make the support
#     graph sparse. They are part of the physical Hamiltonian. Dropping them is
#     not "recovering symmetry"; it is solving a different problem.

# ## Build the Hamiltonian in code
#
# The vendored ERI tensor is stored in chemist order `(p k1, r k3 | q k2, s k4)`.
# To build the physicist-convention Hamiltonian
# ```math
# \hat H = \sum h\, a^{\dagger} a + \tfrac{1}{2} \sum V\, a^{\dagger} a^{\dagger} a a,
# ```
# we therefore read the coefficient as `eri[(k1,k2,k3,k4)][p, r, q, s]` and place
# the annihilators in the order `a_{s k4} a_{r k3}`.

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

    ## One-body part: diagonal in `k`, diagonal in spin, scaled by `1 / nk`.
    for k in 0:nk-1
        h_block = h1e[k] / nk
        for p in 1:norb, s in 1:norb, sigma in (:up, :dn)
            coeff = h_block[p, s]
            iszero(coeff) && continue
            ham += coeff * dag[(k, sigma)][p] * ann[(k, sigma)][s]
        end
    end

    ## Two-body part: only the 8 momentum-conserving ERI blocks survive.
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

    ## The physical Hamiltonian is Hermitian.  The explicit symmetrization below
    ## keeps tiny floating-point asymmetries out of `polyopt`; it does not change
    ## the support graph we care about.
    return registry, 0.5 * (ham + adjoint(ham))
end

registry, ham = build_h4_nk2_hamiltonian(asset.h1e, asset.eri; nk, norb)
nothing #hide
# `registry` stores both creators and annihilators (64 signed symbols total).
# `ham` itself uses 32 physical spin-orbital modes after `abs`-identification.

n_spin_modes = length(variable_indices(ham))

ham_summary = (
    registry_symbols = length(registry),
    physical_modes = n_spin_modes,
    n_terms = length(monomials(ham)),
)

ham_summary

@assert ham_summary.registry_symbols == 64
@assert ham_summary.physical_modes == 32
@assert ham_summary.n_terms == 23_752

# ## What support graph does correlative sparsity see?
#
# The graph is brutally simple:
#
# - one vertex per fermionic mode,
# - an edge whenever two modes appear together in some monomial of `ham`.
#
# That is *not* the same object as the pair-indexed 2-RDM block structure.

function support_edges(poly)
    edges = Set{Tuple{Int, Int}}()
    for mono in monomials(poly)
        idxs = sort!(collect(variable_indices(mono)))
        for i in 1:length(idxs)-1, j in i+1:length(idxs)
            push!(edges, (idxs[i], idxs[j]))
        end
    end
    return edges
end

# The registry order above is
# `(k0,↑), (k0,↓), (k1,↑), (k1,↓)`.
# Collapsing spin therefore maps 32 spin-orbital modes to 16 compound spatial orbitals.
function compound_spatial_id(mode_idx::Int; norb::Int)
    family, orbital0 = divrem(mode_idx - 1, norb)
    k = family ÷ 2
    return k * norb + orbital0 + 1
end

function collapse_spin(edges; norb::Int)
    spatial_edges = Set{Tuple{Int, Int}}()
    for (i, j) in edges
        a = compound_spatial_id(i; norb)
        b = compound_spatial_id(j; norb)
        a == b && continue
        push!(spatial_edges, a < b ? (a, b) : (b, a))
    end
    return spatial_edges
end

spin_edges = support_edges(ham)
spatial_edges = collapse_spin(spin_edges; norb)

support_summary = (
    spin_edges = length(spin_edges),
    complete_spin_graph = n_spin_modes * (n_spin_modes - 1) ÷ 2,
    spatial_edges = length(spatial_edges),
    complete_spatial_graph = (2 * norb) * (2 * norb - 1) ÷ 2,
)

support_summary

@assert support_summary.spin_edges == 496
@assert support_summary.complete_spin_graph == 496
@assert support_summary.spatial_edges == 120
@assert support_summary.complete_spatial_graph == 120

# So the support graph is the complete graph on:
#
# - all 32 spin-orbital modes, and therefore also
# - all 16 compound spatial orbitals after collapsing spin.
#
# Once the graph is complete, a correlative-sparsity routine has nothing left to
# separate. A complete graph has one maximal clique: the whole problem.

# ## Why the paper's 2-RDM blocks and the support graph disagree
#
# The periodic V2RDM paper blocks the **pair-indexed 2-RDM** by total momentum
# `K = (k_1 + k_2) mod N_k`.
#
# Correlative sparsity here works on a different object:
#
# - node = one single mode `a_{p k σ}`,
# - edge = two single modes appear together in some monomial of `H`.
#
# Those are different abstractions. A quartic term can conserve total pair momentum
# and still contain both `k = 0` and `k = 1` single modes in the same monomial.
# The `0101`, `0110`, `1001`, and `1010` blocks do exactly that.
#
# On a model Hamiltonian like pure hopping + on-site Hubbard, the surviving monomial
# support can stay sparse enough that the support graph still decomposes.
# On this ab-initio H₄ asset, the union of the 8 allowed ERI blocks hits every pair,
# so the graph collapses to one clique.

# ## Bottom line
#
# 1. **Build the correct Hamiltonian.** Omit momentum-violating monomials because
#    their coefficients are zero.
# 2. **Keep the momentum-conserving cross-sector terms.** They belong to the real
#    Hamiltonian even though they frustrate graph sparsity.
# 3. **Do not expect the paper's `K`-sector 2-RDM blocks to fall out of the raw
#    Hamiltonian support graph.** That block structure lives on *pairs of modes*,
#    not on single-mode co-occurrence.
# 4. **If you want the paper's blocking, expose the symmetry structurally** — for
#    example with a `k`-blocked variable factory or a user-supplied partition hook.
#
# For the broader benchmark spec, including normalization, chemist/physicist index
# mapping, and the current end-to-end blockers, see
# [Periodic V2RDM Benchmark (H₄ Chain)](@ref periodic-v2rdm-h4).
