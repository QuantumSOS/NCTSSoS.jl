# # [H₄ Chain — k=0 Moment Relaxation](@id h4-chain-energy-benchmark)
#
# This example solves a **reduced model** of the periodic H₄ hydrogen chain
# (Nk=2) using an order-2 fermionic moment relaxation with canonical
# particle-number constraints.
#
# We restrict to the **k=0 intra-sector Hamiltonian** (16 spin-orbital modes)
# and impose $N_\uparrow = N_\downarrow = 2$ via `moment_eq_constraints`.
# This produces a valid SDP lower bound on the k=0-sector ground-state energy.
#
# !!! warning "Scope"
#     This is **not** the full periodic V2RDM benchmark from the paper.
#     It omits cross-k couplings, k-blocked structure, and spin adaptation.
#     For the full asset inspection and a discussion of what the k-blocked
#     formulation requires, see
#     [H₄ Periodic Active-Space Workflow](@ref h4-periodic-active-space).
#
# **System** (Schouten, Ewing & Mazziotti, arXiv:2504.02861):
#
# | Property | Value |
# |:---------|:------|
# | Atoms per unit cell | 4 hydrogen atoms |
# | Bond distance | 1.0 Å |
# | Lattice constant | 4.0 Å |
# | Basis / pseudopotentials | gth-tzvp / gth-pade |
# | Active space per cell | [4 electrons, 8 orbitals] |
# | k-points | Nk = 2 |
#
# The full Nk=2 system has 32 spin-orbital modes.  The k=0 sector alone has
# 16 modes (8 spin-up + 8 spin-down), and the intra-sector ERI block
# $(0,0,0,0)$ gives a self-contained quartic fermionic Hamiltonian on those
# 16 modes.

# ## Setup

using NCTSSoS, LinearAlgebra, Printf

include(joinpath(pkgdir(NCTSSoS), "test", "H4PeriodicAssets.jl"))
using .H4PeriodicAssets: load_nk2_asset, hf_energy

asset = load_nk2_asset()
(; preflight, nk, n_active_orb, n_active_elec,
   total_spin_orbital_modes, h1e, eri) = asset

n_modes_k0 = 2 * n_active_orb  ## 16 spin-orbital modes for k=0

println("System: H₄ chain, Nk=$nk, active space [$(n_active_elec),$(n_active_orb)]")
println("Full system: $total_spin_orbital_modes spin-orbital modes")
println("k=0 sector:  $n_modes_k0 spin-orbital modes")

# ## Spin-orbital index convention (k=0 sector)
#
# Inside the k=0 sector, we flatten $(\text{orbital}, \text{spin})$ into a
# single mode index over 16 modes:
#
# | Mode range | Spin |
# |:-----------|:-----|
# |   1 –  8   | ↑    |
# |   9 – 16   | ↓    |

k0_mode(orb::Int, spin::Int) = spin * n_active_orb + orb + 1
nothing #hide

# ## Building the k=0 Hamiltonian
#
# We use only the one-body block $h^{k=0}_{pq}$ and the intra-sector
# two-electron integral block $(0,0,0,0)$.  The second-quantized form is:
#
# ```math
# H_{k=0} = \sum_{pq\sigma} h^0_{pq}\;
#            a^\dagger_{p\sigma}\, a_{q\sigma}
#          + \frac{1}{2}\sum_{pqrs\,\sigma\tau}
#            (p,r \mid q,s)_{0000}\;
#            a^\dagger_{p\sigma}\, a^\dagger_{q\tau}\,
#            a_{s\tau}\, a_{r\sigma}
# ```
#
# This omits all cross-k coupling terms like $(0,1,0,1)$ — it is a reduced
# model, not the full periodic Hamiltonian.
#
# Use the complex integral blocks **as stored**.  In this orbital gauge, the
# $(0,0,0,0)$ ERI block has nontrivial imaginary parts; the right fix is the
# Hermitianization step below, not silently projecting the Hamiltonian onto its
# real part.

function build_k0_hamiltonian(a, a_dag, h1e, eri; n_orb)
    h = zero(a_dag[1] * a[1])

    ## One-body: h^{k=0}
    h0 = h1e[0]
    for p in 0:n_orb-1, q in 0:n_orb-1
        abs(h0[p+1, q+1]) < 1e-14 && continue
        for σ in 0:1
            h += h0[p+1, q+1] * a_dag[k0_mode(p, σ)] * a[k0_mode(q, σ)]
        end
    end

    ## Two-body: ERI block (0,0,0,0)
    if haskey(eri, (0, 0, 0, 0))
        e00 = eri[(0, 0, 0, 0)]
        for σ in 0:1, τ in 0:1
            for p in 0:n_orb-1, q in 0:n_orb-1
                for r in 0:n_orb-1, s in 0:n_orb-1
                    v = e00[p+1, r+1, q+1, s+1]
                    abs(v) < 1e-14 && continue
                    h += 0.5v *
                         a_dag[k0_mode(p, σ)] *
                         a_dag[k0_mode(q, τ)] *
                         a[k0_mode(s, τ)] *
                         a[k0_mode(r, σ)]
                end
            end
        end
    end
    return h
end

registry, (a, a_dag) = create_fermionic_variables(1:n_modes_k0)
h_k0_raw = build_k0_hamiltonian(a, a_dag, h1e, eri; n_orb = n_active_orb)

println("k=0 Hamiltonian: $(length(h_k0_raw.terms)) terms on $n_modes_k0 modes")

# ## Hermitianization
#
# The raw Hamiltonian is only *numerically* Hermitian — floating-point
# rounding in the integral file means $H \neq H^\dagger$ at machine precision.
# We enforce exact Hermiticity before handing the polynomial to the SDP.
#
# !!! tip "Always Hermitianize ab initio Hamiltonians"
#     Use the complex coefficients from the integral file as-is, then enforce
#     exact $H = H^\dagger$ at the polynomial level.  The Hermitianization step
#     removes machine-epsilon anti-Hermitian noise; it is **not** a license to
#     drop legitimate imaginary ERI coefficients with `real.(...)`.

h_k0 = simplify((h_k0_raw + adjoint(h_k0_raw)) / 2)

println("After Hermitianization: $(length(h_k0.terms)) terms")

# ## Canonical particle-number constraints
#
# The k=0 sector of the Nk=2 H₄ system has 4 active electrons per cell,
# split as $N_\uparrow = 2$ and $N_\downarrow = 2$.  We impose this as a
# state constraint: $(\hat{N}_\sigma - 2)|\psi\rangle = 0$.
#
# This uses `moment_eq_constraints` — **not** `eq_constraints` — because
# fixed particle number is a property of the target state, not an operator
# identity.  See the [Hubbard Model](@ref hubbard-model) example for a
# detailed explanation of why the distinction matters.

n_up = 1.0 * sum(a_dag[i] * a[i] for i in 1:n_active_orb)
n_dn = 1.0 * sum(a_dag[i] * a[i] for i in n_active_orb+1:n_modes_k0)

I_poly = one(h_k0)

println("N↑ modes: 1–$n_active_orb (spin-up)")
println("N↓ modes: $(n_active_orb+1)–$n_modes_k0 (spin-down)")
println("Target sector: N↑ = 2, N↓ = 2")

# ## Order-2 SDP relaxation
#
# An order-2 fermionic moment relaxation builds a moment matrix indexed by
# all fermionic monomials up to degree 2.  For $n$ spin-orbital modes this
# gives $2n$ generators (each $a_i$ and $a_i^\dagger$) and a basis of size
# $\binom{2n}{0} + \binom{2n}{1} + \binom{2n}{2}$.  For $n = 16$:
# $1 + 32 + 496 = 529$ elements.
#
# ### Why we don't need explicit PQG constraints
#
# In quantum chemistry, the "PQG conditions" are positivity constraints on
# the particle ($\mathcal{D}$), hole ($\mathcal{Q}$), and particle-hole
# ($\mathcal{G}$) two-body reduced density matrices.  Our order-2 fermionic
# moment matrix already **implicitly** contains PQG-like positivity:
#
# - The $\mathcal{D}$ (particle) block corresponds to products
#   $a_i^\dagger a_j^\dagger$ applied to $|\psi\rangle$,
# - The $\mathcal{Q}$ (hole) block corresponds to products $a_i a_j$
#   applied to $|\psi\rangle$,
# - The $\mathcal{G}$ (particle-hole) block corresponds to products
#   $a_i^\dagger a_j$ applied to $|\psi\rangle$.
#
# These are all degree-2 fermionic monomials, so they appear as principal
# submatrices of the order-2 moment matrix.  The SDP's
# positive-semidefiniteness constraint on the full moment matrix enforces
# positivity of each block simultaneously — no separate PQG syntax is
# needed.
#
# !!! note "SDP Solver"
#     This example uses [Mosek](https://www.mosek.com/) via `MosekTools`.
#     Any SDP-capable solver works: replace `Mosek.Optimizer` with
#     `COSMO.Optimizer` or `Clarabel.Optimizer` for open-source alternatives.

using MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG" => 0, "MSK_IPAR_NUM_THREADS" => 0)

pop = polyopt(h_k0, registry;
    moment_eq_constraints = [n_up - 2.0 * I_poly,
                             n_dn - 2.0 * I_poly])

config = SolverConfig(optimizer = SOLVER, order = 2, ts_algo = MMD())

println("\nSolving order-2 SDP ($n_modes_k0 modes, k=0, canonical sector) …")
result_first = cs_nctssos(pop, config)

@assert abs(result_first.objective - (-3.6611870442)) < 1e-5 #src
@printf("Order-2 electronic energy (first pass): %.10f Ha\n", result_first.objective)

# ### Refinement
#
# One call to [`cs_nctssos_higher`](@ref) tightens the term-sparsity
# decomposition.  For well-conditioned problems this often improves the
# bound by a small amount at low marginal cost.

println("Refining …")
result_refined = cs_nctssos_higher(pop, result_first, config)

@assert abs(result_refined.objective - (-3.6611883935)) < 1e-5 #src
@assert result_refined.objective ≤ result_first.objective + 1e-8 #src
@printf("Order-2 electronic energy (refined):    %.10f Ha\n", result_refined.objective)

# ## Context: how does this compare?
#
# The k=0 intra-sector SDP bound is a lower bound on the **k=0 sector
# ground-state energy**, not on the full Nk=2 system energy.  To relate it
# to published numbers, we need the constant shift that accounts for
# frozen-core and nuclear-repulsion contributions.

hf_elec = hf_energy(h1e, eri; nk, n_active_elec)
hf_total = preflight["hf_total_energy"]
hf_shift = hf_total - hf_elec

@printf("\nHF active-space electronic energy (all k): %.6f Ha\n", hf_elec)
@printf("HF constant shift:                         %.6f Ha\n", hf_shift)
@printf("HF total energy / cell:                    %.6f Ha\n", hf_total)

# The table below puts the numbers in perspective.  The SDP values above
# apply only to the k=0 intra-sector model; they are **not** directly
# comparable to the full-system V2RDM energy from the paper.
#
# | Quantity | Value (Ha) | Scope |
# |:---------|-----------:|:------|
# | k=0 SDP order 2 (refined) | *see output* | reduced k=0 model |
# | HF total / cell | −2.113 | full Nk=2 system |
# | V2RDM total / cell (figure, ±0.002) | −2.188 | full Nk=2 system |
#
# To reproduce the full periodic V2RDM benchmark, the required formulation
# is a **k-blocked, spin-adapted** SDP that:
#
# 1. keeps crystal momentum $k$ as a conserved quantum number,
# 2. includes cross-k coupling terms ($(k_1,k_2,k_3,k_4) \neq (0,0,0,0)$),
# 3. works in a spin-adapted basis (singlet/triplet).
#
# That formulation is future work for `NCTSSoS.jl`.

# ## Summary
#
# ### What this page showed
#
# | Step | What we did | API |
# |:-----|:------------|:----|
# | 1 | Load vendored Nk=2 active-space integrals | `H4PeriodicAssets.load_nk2_asset` |
# | 2 | Build k=0 intra-sector Hamiltonian (16 modes) | `build_k0_hamiltonian` |
# | 3 | Hermitianize the ab initio Hamiltonian | `simplify((h + adjoint(h)) / 2)` |
# | 4 | Impose canonical $N_\uparrow = N_\downarrow = 2$ | `moment_eq_constraints` in [`polyopt`](@ref) |
# | 5 | Solve order-2 SDP + one refinement step | [`cs_nctssos`](@ref) → [`cs_nctssos_higher`](@ref) |
#
# ### What this page did not do
#
# - Full 32-mode all-k Hamiltonian solve (the coupling graph is complete
#   — correlative sparsity collapses to one dense clique)
# - k-blocked or spin-adapted formulation
# - Reproduce the published V2RDM energy per cell
#
# ### Key concepts
#
# | Concept | One-liner |
# |:--------|:----------|
# | **k=0 intra-sector model** | Use only the $(0,0,0,0)$ ERI block — ignores inter-k coupling |
# | **Hermitianization** | Force exact $H = H^\dagger$ on floating-point integrals before the SDP |
# | **Canonical constraints** | Fix particle number per spin species via `moment_eq_constraints` |
# | **Implicit PQG** | Order-2 fermionic moment matrix contains $\mathcal{D}$/$\mathcal{Q}$/$\mathcal{G}$ blocks as principal submatrices |
# | **Term sparsity** | `MMD()` heuristic decomposes the moment matrix into smaller blocks |
#
# ### See also
#
# - [H₄ Periodic Active-Space Workflow](@ref h4-periodic-active-space) —
#   full asset inspection, HF reconstruction, why the naive 32-mode lift
#   is the wrong interface
# - [Hubbard Model](@ref hubbard-model) — canonical constraints with
#   `moment_eq_constraints`, the `eq_constraints` vs. `moment_eq_constraints`
#   distinction
# - [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
#   creation/annihilation operators, CAR, basic SDP workflow
#
# ### References
#
# - T. Schouten, K. Ewing, D. Mazziotti,
#   "Certified ground-state properties of periodic systems
#   from the variational two-electron reduced-density matrix method,"
#   arXiv:2504.02861 (2025).
