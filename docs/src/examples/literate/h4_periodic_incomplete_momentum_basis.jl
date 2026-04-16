# # [H₄ Periodic Chain — Incomplete Momentum-Basis Relaxation](@id h4-periodic-incomplete-momentum-basis)
#
# This example takes the **full Nk=2 periodic H₄ active-space Hamiltonian**
# and solves a **deliberately incomplete** order-2 fermionic relaxation.
#
# The point is not to reproduce the paper's periodic V2RDM number.  It does
# not.  The point is narrower and more useful:
#
# 1. verify that `NCTSSoS.jl` accepts a **hand-built custom basis** for this
#    32-mode periodic problem,
# 2. keep the solve small enough that it actually finishes on ordinary CPU
#    hardware, and
# 3. show what kind of bound this compromise produces.
#
# !!! warning "Scope"
#     This is **not** the faithful periodic V2RDM formulation from the paper.
#     The basis below keeps only the identity, all annihilation singles,
#     and all annihilation pairs.  That is enough to make the relaxation
#     well-defined and solvable, but it drops a lot of order-2 information.
#     So the bound is expected to be **much weaker** than the true periodic
#     V2RDM bound.
#
# Think of this page as the pragmatic follow-up after
# [H₄ Periodic Chain — Correlative-Sparsity Attempt](@ref h4-periodic-correlative-sparsity):
# the full-Nk=2 Hamiltonian stays, but we give up on the dense automatic order-2
# basis and replace it with a smaller user-specified one.
#
# Relative to the other H₄ pages:
#
# - [H₄ Periodic Active-Space Workflow](@ref h4-periodic-active-space) inspects
#   the full Nk=2 asset honestly but does **not** solve it,
# - [H₄ Periodic Chain — Correlative-Sparsity Attempt](@ref h4-periodic-correlative-sparsity)
#   keeps the full order-2 basis and shows where pure correlative sparsity runs
#   out of road,
# - [H₄ Chain — k=0 Moment Relaxation](@ref h4-chain-energy-benchmark) solves a
#   smaller **k=0-only** reduction.
#
# Here we still solve the **full Nk=2 Hamiltonian**, but with a cheaper,
# user-specified basis.
#
# ## Why this custom basis is even legal
#
# `NCTSSoS.jl` lets you replace the automatically generated order-based basis
# with an explicit one through `SolverConfig(moment_basis=...)`.  For fermions,
# the moment matrix entry is `bᵢ' * bⱼ`, so if the basis contains annihilation
# words then their adjoints contribute the corresponding creation words.
#
# The basis used here is:
#
# - the identity,
# - all annihilation singles `aᵢ`,
# - all annihilation pairs `aᵢ aⱼ`.
#
# The pair words naturally split by total crystal momentum
# `K = kᵢ + kⱼ mod Nk`.  For `Nk = 2` that gives two sectors, `K = 0` and
# `K = 1`.
#
# This basis is **not** complete at order 2.  That's the whole trick.
# It is smaller, cheaper, and weaker.
#
# ## Setup
#
# We use the vendored periodic H₄ asset and a silent `COSMO` backend.
# Expect a solve time on the order of **10–15 minutes** on a laptop-scale CPU.

using NCTSSoS
using COSMO
using JuMP
using LinearAlgebra
using Printf

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: load_nk2_asset, hf_energy

const MOI = NCTSSoS.MOI
const SILENT_COSMO = MOI.OptimizerWithAttributes(COSMO.Optimizer, MOI.Silent() => true)

asset = load_nk2_asset();
(; preflight, figure, nk, n_active_orb, n_active_elec, h1e, eri) = asset;

n_spatial = nk * n_active_orb       # 16 spatial orbitals total
n_modes = 2 * n_spatial             # 32 spin-orbital modes total

println("Nk = ", nk)
println("spatial orbitals = ", n_spatial)
println("spin-orbital modes = ", n_modes)

# ## Full periodic mode convention
#
# We flatten the compound orbital label `(k, p, σ)` into one spin-orbital mode:
#
# - `k ∈ {0, 1}` is crystal momentum,
# - `p ∈ {0, …, 7}` is the active orbital inside that k-sector,
# - `σ ∈ {0, 1}` is spin (`↑`, `↓`).
#
# The first 16 modes are spin-up, the next 16 are spin-down.

mode(k::Int, orb::Int, spin::Int) = spin * n_spatial + k * n_active_orb + orb + 1
k_of_mode(m::Int) = ((m - 1) % n_spatial) ÷ n_active_orb
nothing #hide

# ## Building the full Nk=2 Hamiltonian
#
# Unlike the [k=0 companion page](@ref h4-chain-energy-benchmark), this example
# keeps **all eight momentum-conserving ERI blocks**.  We therefore solve the
# full periodic Hamiltonian on all 32 spin-orbital modes.
#
# To stay consistent with the asset bookkeeping, we normalize the one-body and
# two-body integrals per unit cell by `1 / Nk` and `1 / Nk^2` respectively.

function build_full_periodic_hamiltonian(
    a,
    a_dag,
    h1e,
    eri;
    nk::Int,
    n_orb::Int,
)
    h = zero(a_dag[1] * a[1])
    h1e_scale = 1 / nk
    eri_scale = 1 / nk^2

    for k in 0:nk-1, p in 0:n_orb-1, q in 0:n_orb-1
        hpq = h1e_scale * h1e[k][p + 1, q + 1]
        abs(hpq) < 1e-14 && continue
        for σ in 0:1
            h += hpq * a_dag[mode(k, p, σ)] * a[mode(k, q, σ)]
        end
    end

    for ((k1, k2, k3, k4), block) in eri
        for σ in 0:1, τ in 0:1
            for p in 0:n_orb-1, q in 0:n_orb-1, r in 0:n_orb-1, s in 0:n_orb-1
                v = eri_scale * block[p + 1, r + 1, q + 1, s + 1]
                abs(v) < 1e-14 && continue
                h += 0.5v *
                     a_dag[mode(k1, p, σ)] *
                     a_dag[mode(k2, q, τ)] *
                     a[mode(k4, s, τ)] *
                     a[mode(k3, r, σ)]
            end
        end
    end

    return simplify((h + adjoint(h)) / 2)
end

registry, (a, a_dag) = create_fermionic_variables(1:n_modes);
h = build_full_periodic_hamiltonian(a, a_dag, h1e, eri; nk, n_orb=n_active_orb);

println("Hamiltonian terms = ", length(h.terms))

# ## The incomplete momentum-inspired basis
#
# Start from the smallest useful order-2 sector we can get away with:
#
# - `1`,
# - `aᵢ`,
# - `aᵢ aⱼ`.
#
# This keeps the pair sector explicit while avoiding the full dense order-2
# basis.  For 32 modes, the dense order-2 fermionic basis has
# `1 + 64 + binomial(64, 2) = 2081` canonical words; our custom basis has only
# `1 + 32 + binomial(32, 2) = 529` words.

function incomplete_momentum_basis(a; nk::Int)
    NM = typeof(one(a[1]))
    basis = NM[one(a[1])]
    append!(basis, a)

    pair_sector_counts = zeros(Int, nk)
    for i in 1:length(a)-1, j in i+1:length(a)
        push!(basis, only(monomials(a[i] * a[j])))
        total_k = mod(k_of_mode(i) + k_of_mode(j), nk)
        pair_sector_counts[total_k + 1] += 1
    end

    return basis, pair_sector_counts
end

basis, pair_sector_counts = incomplete_momentum_basis(a; nk);
dense_basis = get_ncbasis(registry, 2);

@assert length(dense_basis) == 2081 #src
@assert length(basis) == 529 #src
@assert pair_sector_counts == [240, 256] #src

println("dense order-2 basis size = ", length(dense_basis))
println("custom basis size        = ", length(basis))
println("pair sector K=0 count    = ", pair_sector_counts[1])
println("pair sector K=1 count    = ", pair_sector_counts[2])

# That is the key tradeoff: we keep the full Hamiltonian, but we no longer ask
# for the full order-2 relaxation.

# ## Canonical particle-number sector
#
# As in the [k=0 companion page](@ref h4-chain-energy-benchmark), we impose the
# target particle numbers through `moment_eq_constraints`, not operator
# equalities.  For the full Nk=2 problem the active space contains 8 electrons:
# 4 spin-up and 4 spin-down.

I_poly = one(h);
n_up = 1.0 * sum(a_dag[i] * a[i] for i in 1:n_spatial);
n_dn = 1.0 * sum(a_dag[i] * a[i] for i in (n_spatial + 1):n_modes);

pop = polyopt(h, registry;
    moment_eq_constraints=[n_up - 4.0 * I_poly,
                           n_dn - 4.0 * I_poly]);

# ## What term sparsity does with this basis
#
# Even though the correlative graph still collapses to one clique, the custom
# basis is small enough that `MMD()` can break the moment matrix into many
# smaller term-sparsity blocks.

config = SolverConfig(
    optimizer=SILENT_COSMO,
    moment_basis=basis,
    ts_algo=MMD(),
);

sparsity = compute_sparsity(pop, config);
block_sizes = length.(only(sparsity.cliques_term_sparsities)[1].block_bases);

@assert length(sparsity.corr_sparsity.cliques) == 1 #src
@assert length(block_sizes) == 88 #src
@assert maximum(block_sizes) == 82 #src
@assert minimum(block_sizes) == 1 #src

println("number of cliques = ", length(sparsity.corr_sparsity.cliques))
println("moment blocks     = ", length(block_sizes))
println("largest block     = ", maximum(block_sizes))
println("smallest block    = ", minimum(block_sizes))

# The large dense 2081-word lift is gone.  What remains is still serious work,
# but at least it is work the solver can finish.

# ## Solving the compromise relaxation
#
# Now solve the SDP.  The result is a **lower bound** for this incomplete
# relaxation — not the full periodic V2RDM value.

result = cs_nctssos(pop, config);

hf_active = hf_energy(h1e, eri; nk, n_active_elec);
hf_shift = preflight["hf_constant_shift"];
result_total = result.objective + hf_shift;

@assert result.objective ≈ -6.366838786026133 atol = 5e-3 #src
@assert result_total ≈ -5.089926089471001 atol = 5e-3 #src

println("termination status        = ", JuMP.termination_status(result.model))
@printf("electronic lower bound    = %.12f Ha\n", result.objective)
@printf("total lower bound (+shift)= %.12f Ha\n", result_total)
@printf("HF active-space energy    = %.12f Ha\n", hf_active)
@printf("HF total / cell           = %.12f Ha\n", preflight["hf_total_energy"])
@printf("figure V2RDM / cell       = %.12f Ha\n", figure["V2RDM_4_8"])
println("unique moment variables   = ", result.n_unique_moment_matrix_elements)

# ## Reading the result honestly
#
# The important outcomes are:
#
# 1. **Feasibility**: a hand-built incomplete order-2 basis is fully supported
#    by the current API.
# 2. **Practicality**: the combination of that custom basis with term sparsity
#    does produce a solvable 32-mode periodic relaxation.
# 3. **Weakness**: the resulting bound is far below the physically relevant
#    energies, so this basis is too incomplete to stand in for the true
#    periodic V2RDM relaxation.
#
# That is not failure.  It is information.  The experiment tells us the right
# next step: add back carefully chosen order-2 sectors — most obviously
# particle-hole words — instead of jumping straight to the full dense basis.
#
# ## Summary
#
# This example answers the practical question left open by the full-asset page:
#
# - **Can we specify an incomplete order-2 momentum basis ourselves?** Yes.
# - **Can we run term sparsity on top of it?** Yes.
# - **Do we get a result?** Yes.
# - **Is this already the periodic V2RDM benchmark?** No.
#
# See also:
#
# - [H₄ Periodic Active-Space Workflow](@ref h4-periodic-active-space) for the
#   full Nk=2 asset inspection and the dense-order-2 blocker,
# - [H₄ Periodic Chain — Correlative-Sparsity Attempt](@ref h4-periodic-correlative-sparsity)
#   for the no-custom-basis Mosek attempt that motivates this compromise,
# - [H₄ Chain — k=0 Moment Relaxation](@ref h4-chain-energy-benchmark) for a
#   smaller reduced model that gives a much more meaningful lower bound.
