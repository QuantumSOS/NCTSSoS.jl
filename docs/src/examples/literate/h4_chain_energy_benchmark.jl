# # [H₄ Chain — k=0 Moment Relaxation](@id h4-chain-energy-benchmark)
#
# This example computes a **certified lower bound** on the ground-state
# energy of four hydrogen atoms arranged in a repeating chain — the simplest
# model of a crystalline solid.  We do this by converting the quantum
# mechanics problem into a **semidefinite program** (SDP).
#
# We work with a **reduced model**: only the k=0 sector of an Nk=2 periodic
# hydrogen chain, solved with an order-2 fermionic moment relaxation and
# canonical particle-number constraints.
#
# !!! warning "Scope"
#     This is **not** the full periodic V2RDM benchmark from the paper.
#     It omits cross-k couplings, k-blocked structure, and spin adaptation.
#     For the full asset inspection and a discussion of what the k-blocked
#     formulation requires, see
#     [H₄ Periodic Active-Space Workflow](@ref h4-periodic-active-space).
#
# ## Background for newcomers
#
# This section defines every concept you need to follow the code below.
# If you already know quantum chemistry, skip ahead to **Setup**.
#
# ### What is H₄?
#
# "H₄" means **four hydrogen atoms**.  Hydrogen is the simplest atom in
# nature — one proton in the nucleus, one electron orbiting it.  Four atoms
# give us four electrons whose behaviour we want to predict.
#
# We arrange the atoms in a **repeating chain**, like beads on a string that
# repeats forever.  That makes this a tiny model of a *crystal* — the kind of
# repeating atomic pattern found in metals, semiconductors, and insulators.
#
# ### The electronic structure problem
#
# Atoms consist of heavy nuclei and light electrons.  The nuclei barely move
# (the *Born–Oppenheimer approximation*), so the interesting question is:
# **what are the electrons doing?**
#
# Quantum mechanics says the answer is encoded in the *Schrödinger equation*:
#
# ```math
# \hat{H}\,|\Psi\rangle \;=\; E\,|\Psi\rangle
# ```
#
# - $|\Psi\rangle$ is the **wave function** — a complete description of all
#   electrons.
# - $\hat{H}$ is the **Hamiltonian** — an operator that encodes the physics
#   (kinetic energy of electrons, attraction to nuclei, repulsion between
#   electrons).
# - $E$ is the **energy**.
#
# The smallest possible $E$ is called the **ground-state energy** — it
# describes the most stable configuration.  Finding this number is the
# central problem of computational quantum chemistry.
#
# **Why it's hard:** For $N$ electrons that can each occupy any of $r$
# available states ("orbitals"), the wave function lives in a space of
# dimension $\binom{r}{N}$.  For 100 orbitals and 50 electrons that's
# $\approx 10^{29}$.  You cannot store or diagonalize a matrix that large.
#
# ### How SDP (semidefinite programming) helps
#
# The Hamiltonian contains only *one-body* terms (each electron on its own)
# and *two-body* terms (pairs of electrons repelling each other) — no
# three-electron interactions.  So the energy depends on electrons only
# through **pairs**.
#
# !!! info "Analogy"
#     Imagine computing the total "social energy" of a party.  If every
#     interaction is a two-person conversation, you don't need to track all
#     possible group configurations — just the pairwise interaction patterns.
#
# The mathematical object that captures all pairwise information is the
# **2-electron reduced density matrix** (2-RDM): a matrix of size roughly
# $r^2 \times r^2$ — polynomial in $r$, not exponential.  The energy turns
# out to be a *linear* function of the 2-RDM:
# $E = \operatorname{Tr}({}^2\!K \;\cdot\; {}^2\!D)$.
#
# The catch: not every matrix is a physically valid 2-RDM.  The constraints
# that enforce physical validity are called **N-representability conditions**
# — and the key ones are *semidefinite* constraints (certain derived matrices
# must be positive semidefinite).  Minimizing energy subject to these
# constraints is an SDP.
#
# Because the N-representability conditions used are *necessary but not
# sufficient*, the SDP minimum is a **lower bound** on the true ground-state
# energy.  In practice, this bound is remarkably tight — recovering 82–99%
# of the correlation energy for real molecules.  This approach is called
# **V2RDM** (variational 2-RDM method).
#
# ### What "periodic" and "k-points" mean
#
# A **crystal** is a material whose atoms repeat in a regular pattern —
# like tiles on a floor.  The smallest repeating unit is the **unit cell**.
#
# Instead of simulating an infinite crystal, physicists exploit the
# repetition using *Bloch's theorem*: every electronic property decomposes
# into contributions labeled by a **crystal momentum** $k$.  Think of $k$ as
# a "frequency channel": $k = 0$ means the electron wave looks identical in
# every cell (no phase shift); other $k$ values correspond to waves that
# shift phase from cell to cell.
#
# The set of discrete $k$-values used in a finite simulation is called
# the set of **k-points**.  Here we use $N_k = 2$ k-points ($k=0$ and
# $k=1$), and **restrict to the $k=0$ sector only** — giving a
# self-contained 16-mode problem instead of the full 32-mode system.
#
# ### What an "active space" is
#
# Real materials have many electrons, most of which sit in low-energy
# "core" orbitals that barely participate in chemistry.  The **active-space
# approximation** freezes those core electrons at a cheap mean-field level
# and focuses the expensive SDP only on the "active" electrons near the
# energy frontier.
#
# Here, **[4, 8]** means: 4 active electrons in 8 active orbitals per unit
# cell.  The energy of the frozen core is added back as a constant after
# optimization.
#
# ### What a "moment relaxation" is
#
# In polynomial optimization, a **moment relaxation** replaces a hard
# problem with a hierarchy of SDPs, each tighter than the last.  At
# order $d$, you build a *moment matrix* whose entries are expected values
# of products of variables up to degree $2d$, and require it to be positive
# semidefinite.
#
# For electrons, the "variables" are **creation and annihilation operators**
# ($a_i^\dagger$ adds an electron to mode $i$; $a_i$ removes one).  An
# **order-2** relaxation considers all products of up to 2 operators —
# which captures exactly the pairwise correlations the 2-RDM encodes.
#
# !!! info "V2RDM = order-2 moment relaxation"
#     The order-2 fermionic moment relaxation and the V2RDM method with
#     PQG conditions are the **same SDP** in different notation.
#     The moment matrix's subblocks contain the particle ($\mathcal{D}$),
#     hole ($\mathcal{Q}$), and particle-hole ($\mathcal{G}$) matrices
#     that chemists enforce separately.  Requiring the full moment matrix
#     to be PSD enforces all three at once.
#
# ---
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
# Each electron has a spatial position (which orbital it occupies) and a
# **spin** — an intrinsic quantum property that comes in two flavours:
# spin-up (↑) and spin-down (↓).  A **spin-orbital** combines both labels
# into one mode index.  Inside the k=0 sector we have 8 spatial orbitals
# × 2 spin species = 16 spin-orbitals, numbered as:
#
# | Mode range | Spin |
# |:-----------|:-----|
# |   1 –  8   | ↑    |
# |   9 – 16   | ↓    |

k0_mode(orb::Int, spin::Int) = spin * n_active_orb + orb + 1
nothing #hide

# ## Building the k=0 Hamiltonian
#
# The **Hamiltonian** is the mathematical operator that encodes all the
# physics of the system.  In **second-quantized** form — the standard
# language for many-electron problems — it is written in terms of creation
# operators $a^\dagger$ (which add an electron to a mode) and annihilation
# operators $a$ (which remove one).
#
# Our Hamiltonian has two pieces:
#
# - **One-body terms** ($h^0_{pq}$): an electron hops from orbital $q$ to
#   orbital $p$ — kinetic energy plus attraction to nuclei.
# - **Two-body terms** (ERI = electron repulsion integrals): two electrons
#   in orbitals $(p, q)$ repel each other into orbitals $(r, s)$.
#
# The second-quantized form is:
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
# We use only the $(0,0,0,0)$ ERI block — the block where all four orbital
# indices belong to the $k=0$ sector.  This omits all cross-k coupling
# terms like $(0,1,0,1)$ — it is a reduced model, not the full periodic
# Hamiltonian.
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
# A Hamiltonian must be **Hermitian** — physically, this means energy is a
# real number, and mathematically, $H = H^\dagger$ (the operator equals its
# conjugate transpose).
#
# The raw Hamiltonian we built is only *numerically* Hermitian — floating-point
# rounding in the integral file means $H \neq H^\dagger$ at machine precision.
# We enforce exact Hermiticity by averaging: $H \leftarrow (H + H^\dagger)/2$.
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
# Electrons are conserved — you can't create or destroy them.  Our system
# has a fixed number of electrons: 2 spin-up and 2 spin-down in the k=0
# sector.  We tell the SDP about this by imposing **particle-number
# constraints**: the state $|\psi\rangle$ we're optimizing over must satisfy
# $\hat{N}_\uparrow |\psi\rangle = 2|\psi\rangle$ and
# $\hat{N}_\downarrow |\psi\rangle = 2|\psi\rangle$.
#
# This uses `moment_eq_constraints` — **not** `eq_constraints` — because
# fixed particle number is a property of the *target state*, not an operator
# identity that holds for all states.  See the [Hubbard Model](@ref hubbard-model)
# example for a detailed explanation of why the distinction matters.

n_up = 1.0 * sum(a_dag[i] * a[i] for i in 1:n_active_orb)
n_dn = 1.0 * sum(a_dag[i] * a[i] for i in n_active_orb+1:n_modes_k0)

I_poly = one(h_k0)

println("N↑ modes: 1–$n_active_orb (spin-up)")
println("N↓ modes: $(n_active_orb+1)–$n_modes_k0 (spin-down)")
println("Target sector: N↑ = 2, N↓ = 2")

# ## Order-2 SDP relaxation
#
# We now solve the SDP.  An order-2 fermionic moment relaxation builds a
# **moment matrix** indexed by all fermionic monomials up to degree 2.
# For $n$ spin-orbital modes this gives $2n$ generators (each $a_i$ and
# $a_i^\dagger$) and a basis of size
# $\binom{2n}{0} + \binom{2n}{1} + \binom{2n}{2}$.  For $n = 16$:
# $1 + 32 + 496 = 529$ elements.
#
# The SDP solver finds the assignment of expected values to these 529
# basis elements that (a) makes the moment matrix positive semidefinite,
# (b) respects the particle-number constraints, and (c) minimizes the
# energy.
#
# ### Why we don't need explicit PQG constraints
#
# In quantum chemistry, the "PQG conditions" are positivity constraints on
# three derived matrices called the particle ($\mathcal{D}$), hole
# ($\mathcal{Q}$), and particle-hole ($\mathcal{G}$) reduced density
# matrices.  Our order-2 fermionic moment matrix already **implicitly**
# contains all three:
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
# One call to [`cs_nctssos_higher`](@ref) tightens the **term-sparsity
# decomposition** — a technique that breaks the large moment matrix into
# smaller blocks, reducing solver cost.  For well-conditioned problems this
# often improves the bound by a small amount at low marginal cost.

println("Refining …")
result_refined = cs_nctssos_higher(pop, result_first, config)

@assert abs(result_refined.objective - (-3.6611883935)) < 1e-5 #src
@assert result_refined.objective ≤ result_first.objective + 1e-8 #src
@printf("Order-2 electronic energy (refined):    %.10f Ha\n", result_refined.objective)

# ## Context: how does this compare?
#
# The result above is a lower bound on the **k=0 sector ground-state
# energy** — the energy of the electrons in the $k=0$ "frequency channel"
# only.  It is *not* the total energy of the full system.
#
# To relate it to published numbers, we need the **constant shift**: the
# energy from frozen-core electrons and nuclear repulsion that we set aside
# when we defined the active-space Hamiltonian.

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
# ### Glossary
#
# | Term | Meaning |
# |:-----|:--------|
# | **Ground-state energy** | Lowest possible energy of the electrons — the most stable configuration |
# | **Hamiltonian** | The mathematical operator that encodes all the physics (kinetic energy, attractions, repulsions) |
# | **Orbital / spin-orbital** | A possible state for one electron.  A spin-orbital specifies both spatial shape and spin (↑ or ↓) |
# | **2-RDM** | 2-electron reduced density matrix — captures all pairwise electron correlations in a polynomial-size matrix |
# | **N-representability** | Constraints ensuring a candidate 2-RDM could actually come from a physical quantum state |
# | **PQG conditions** | Three positive-semidefiniteness constraints (particle, hole, particle-hole) — the standard N-representability conditions for V2RDM |
# | **Active space [n, r]** | Solve the expensive part for $n$ electrons in $r$ orbitals; freeze the rest |
# | **k-point / crystal momentum** | A label for how an electron's wave shifts phase between unit cells.  $k=0$ means no shift |
# | **Moment relaxation** | An SDP that lower-bounds the optimum by constraining a moment matrix to be PSD |
# | **Hermitianization** | Averaging $H$ and $H^\dagger$ to enforce exact self-adjointness on floating-point data |
# | **ERI** | Electron repulsion integrals — the numbers that quantify how much two electrons in given orbitals repel each other |
# | **Term sparsity** | `MMD()` heuristic that decomposes the moment matrix into smaller blocks for computational efficiency |
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
