# # [Kitaev Chain (Topological Superconductor)](@id kitaev-chain)
#
# This example computes **certified lower bounds** on the ground-state energy
# of the **Kitaev chain** — the simplest model that hosts a *topological*
# superconducting phase with **Majorana zero modes** at its ends.
#
# If you haven't seen fermionic operators before, start with the
# [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) example,
# which builds the language of creation/annihilation operators and the CAR
# from scratch.  Here we assume that background and focus on what's new:
# **pairing terms** that break particle-number conservation, **Majorana
# operators**, the **topological vs. trivial phase**, and the difference
# between **open** and **periodic** boundary conditions.
#
# ## The Kitaev chain Hamiltonian
#
# Picture $N$ sites on an open wire (no ring).  Each site can hold zero or
# one spinless fermion.  The Hamiltonian has three terms:
#
# ```math
# H = -\mu \sum_{i=1}^{N} a_i^\dagger a_i
#     \;-\; t \sum_{i=1}^{N-1}\bigl(a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i\bigr)
#     \;+\; \Delta \sum_{i=1}^{N-1}\bigl(a_i\, a_{i+1} + a_{i+1}^\dagger a_i^\dagger\bigr).
# ```
#
# | Term | Physics | Effect on particle number |
# |:-----|:--------|:-------------------------|
# | $-\mu\, a_i^\dagger a_i$ | **Chemical potential** — energy cost of occupying a site | Conserves $N_{\mathrm{tot}}$ |
# | $-t\,(a_i^\dagger a_{i+1} + \text{h.c.})$ | **Hopping** — fermion tunnels between neighbours | Conserves $N_{\mathrm{tot}}$ |
# | $\Delta\,(a_i\, a_{i+1} + \text{h.c.})$ | **Pairing** — creates or destroys fermion *pairs* from a condensate | **Breaks** $N_{\mathrm{tot}}$ conservation |
#
# The pairing term is the key novelty compared to the XY model.  It models
# a *p-wave* superconductor: the pairing amplitude connects neighbouring
# sites (odd under spatial inversion in momentum space).  Because pairs of
# fermions are created and destroyed, particle number $N_{\mathrm{tot}}$
# is no longer conserved.  What *is* conserved is **fermionic parity**:
# $P = (-1)^{N_{\mathrm{tot}}}$.  Every term in $H$ changes the fermion
# count by $0$ or $\pm 2$, so even/odd parity never flips.
#
# !!! note "Open boundary conditions — no parity constraint needed"
#     The [XY model example](@ref fermionic-ground-state) uses **periodic**
#     boundaries (a ring), where the boundary term picks up a parity
#     operator $P$ and we must impose $P = I$ as a constraint.  Here we
#     use **open** boundaries (a wire) — there is no boundary twist, so
#     no parity constraint is required.  The SDP relaxation searches
#     over all parity sectors automatically.

# ---
#
# ## Majorana operators and the topological sweet spot
#
# Every complex fermion can be split into two **Majorana** (real) fermions:
#
# ```math
# c_{2j-1} = a_j + a_j^\dagger, \qquad
# c_{2j} = \frac{a_j - a_j^\dagger}{i}.
# ```
#
# These are self-conjugate ($c_l^\dagger = c_l$) and obey
# $\{c_l, c_m\} = 2\delta_{lm}$.  Think of them as splitting a complex
# number into its real and imaginary parts — except the "parts" are
# operators that anticommute.
#
# ### The sweet spot: $\mu = 0$, $t = \Delta$
#
# Set $\mu = 0$ and $t = \Delta$.  The Hamiltonian simplifies to:
#
# ```math
# H = it \sum_{j=1}^{N-1} c_{2j}\, c_{2j+1}.
# ```
#
# Each term pairs the Majorana $c_{2j}$ from site $j$ with $c_{2j+1}$ from
# site $j{+}1$ — the pairing reaches **across** sites.  But the first
# Majorana $c_1$ and the last $c_{2N}$ appear in *no* term of $H$: they
# are the **unpaired Majorana zero modes**.  Any operator that commutes
# with $H$ and is built from $c_1$ and $c_{2N}$ alone costs zero energy —
# hence "zero mode."  These edge modes are the hallmark of the topological
# phase.
#
# Under the Jordan–Wigner transformation, $ic_{2j}\,c_{2j+1} = X_j X_{j+1}$,
# so the sweet-spot Hamiltonian maps to the open **XX chain**:
#
# ```math
# H_{\text{sweet}} = \sum_{j=1}^{N-1} X_j X_{j+1},
# \qquad E_0 = -(N-1).
# ```
#
# We'll verify this Majorana↔Pauli connection explicitly in the
# [bonus section](#Bonus:-Sweet-spot-via-PauliAlgebra) below.

# ---
#
# ## Exact diagonalisation oracle
#
# To check our SDP bounds we need an exact answer.  We build the full
# $2^N \times 2^N$ Hamiltonian matrix via the **Jordan–Wigner** mapping
# (fermion operators → tensor products of Pauli matrices) and take the
# smallest eigenvalue.  This is brute-force exact diagonalisation — it
# works for small $N$ and serves as our verification oracle.

using LinearAlgebra

## Jordan–Wigner building blocks
const PAULI_Z  = ComplexF64[1 0; 0 -1]
const PAULI_I  = Matrix{ComplexF64}(I, 2, 2)
const JW_RAISE = ComplexF64[0 1; 0 0]   # σ⁺ = |0⟩⟨1|
const JW_LOWER = ComplexF64[0 0; 1 0]   # σ⁻ = |1⟩⟨0|

"""
    jw_fermion_op(mode, dagger, nmodes)

Jordan–Wigner representation of a fermionic operator on `nmodes` sites.
Returns the `2^nmodes × 2^nmodes` matrix for `a†_mode` (if `dagger=true`)
or `a_mode` (if `dagger=false`).
"""
function jw_fermion_op(mode::Int, dagger::Bool, nmodes::Int)
    mats = [site < mode ? PAULI_Z :
            site == mode ? (dagger ? JW_RAISE : JW_LOWER) :
            PAULI_I for site in 1:nmodes]
    return reduce(kron, mats)
end

"""
    kitaev_chain_exact(nsites; μ=0.0, t=1.0, Δ=1.0)

Exact ground-state energy of the `nsites`-site Kitaev chain via full
exact diagonalisation.
"""
function kitaev_chain_exact(nsites::Int; μ::Real = 0.0, t::Real = 1.0, Δ::Real = 1.0)
    a  = [jw_fermion_op(i, false, nsites) for i in 1:nsites]
    a_dag = [jw_fermion_op(i, true,  nsites) for i in 1:nsites]
    dim = 2^nsites
    H = zeros(ComplexF64, dim, dim)
    for i in 1:nsites                                       # chemical potential
        H .-= μ .* (a_dag[i] * a[i])
    end
    for i in 1:(nsites - 1)
        H .-= t .* (a_dag[i] * a[i+1] + a_dag[i+1] * a[i])  # hopping
        H .+= Δ .* (a[i] * a[i+1] + a_dag[i+1] * a_dag[i])  # pairing
    end
    return eigmin(Hermitian(H))
end

# Quick sanity check — the sweet spot should give $E_0 = -(N-1)$:

for N in [3, 4, 5, 6]
    E = kitaev_chain_exact(N; μ = 0.0, t = 1.0, Δ = 1.0)
    @assert abs(E - (-(N - 1))) < 1e-12   #src
    println("Sweet spot  N = $N :  E₀ = $(round(E; digits=6))  (expected $(-(N-1)))")
end

# ---
#
# ## Solving with NCTSSoS.jl
#
# We now reproduce these exact values using **semidefinite programming**
# relaxations.  The solver doesn't know the answer — it searches for the
# lowest expectation value of $H$ consistent with the fermionic algebra
# (CAR).  The result is a **certified lower bound** on the true $E_0$.

# !!! note "SDP Solver"
#     These examples use [Mosek](https://www.mosek.com/) via `MosekTools`.
#     Any SDP-capable solver works: replace `Mosek.Optimizer` with
#     `COSMO.Optimizer` or `Clarabel.Optimizer` for open-source alternatives.

using NCTSSoS, MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG" => 0,
    "MSK_IPAR_NUM_THREADS" => 0)

# ---
# ## Case 1 — Sweet spot ($\mu = 0$, $t = \Delta = 1$)
#
# At the sweet spot, all Majorana pairs across bonds are occupied and the
# ground-state energy is $E_0 = -(N-1)$.  This is the deep topological phase.
# ---

N₁ = 4
registry₁, (a₁, a₁_dag) = create_fermionic_variables(1:N₁)

# Build the Hamiltonian term by term.  With $\mu = 0$ and $t = \Delta = 1$:
# $H = -\sum (a^\dagger_i a_{i+1} + \text{h.c.}) + \sum (a_i a_{i+1} + \text{h.c.})$

hopping₁ = sum(a₁_dag[i] * a₁[i+1] + a₁_dag[i+1] * a₁[i] for i in 1:(N₁ - 1))
pairing₁ = sum(a₁[i] * a₁[i+1] + a₁_dag[i+1] * a₁_dag[i] for i in 1:(N₁ - 1))
ham₁ = pairing₁ - hopping₁

pop₁ = polyopt(ham₁, registry₁)
config₁ = SolverConfig(optimizer = SOLVER, order = 1)
result₁ = cs_nctssos(pop₁, config₁)

exact₁ = kitaev_chain_exact(N₁; μ = 0.0, t = 1.0, Δ = 1.0)
@assert abs(result₁.objective - exact₁) < 1e-5   #src
println("Case 1 (sweet spot):  SDP = $(result₁.objective),  exact = $exact₁,  gap = $(abs(result₁.objective - exact₁))")

# Tight at **order 1** — the SDP relaxation captures the full ground state
# for this small system.

# ---
# ## Case 2 — Trivial phase ($\mu = 3$, $t = 1$, $\Delta = 1$)
#
# When $|\mu| > 2|t|$ the chemical potential dominates and the system enters
# the **trivial** phase — no Majorana zero modes, no topological protection.
# The gap doesn't close; the chain is a boring band insulator.
# ---

N₂ = 4
μ₂, t₂, Δ₂ = 3.0, 1.0, 1.0
registry₂, (a₂, a₂_dag) = create_fermionic_variables(1:N₂)

chemical₂ = -μ₂ * sum(a₂_dag[i] * a₂[i] for i in 1:N₂)
hopping₂  = -t₂ * sum(a₂_dag[i] * a₂[i+1] + a₂_dag[i+1] * a₂[i] for i in 1:(N₂ - 1))
pairing₂  =  Δ₂ * sum(a₂[i] * a₂[i+1] + a₂_dag[i+1] * a₂_dag[i] for i in 1:(N₂ - 1))
ham₂ = chemical₂ + hopping₂ + pairing₂

pop₂ = polyopt(ham₂, registry₂)
config₂ = SolverConfig(optimizer = SOLVER, order = 1)
result₂ = cs_nctssos(pop₂, config₂)

exact₂ = kitaev_chain_exact(N₂; μ = μ₂, t = t₂, Δ = Δ₂)
@assert abs(result₂.objective - exact₂) < 1e-5   #src
println("Case 2 (trivial):    SDP = $(result₂.objective),  exact = $(round(exact₂; digits=6)),  gap = $(abs(result₂.objective - exact₂))")

# ---
# ## Case 3 — Asymmetric pairing ($\mu = 0$, $t = 1$, $\Delta = 0.5$)
#
# With $\Delta \neq t$ but still $|\mu| < 2|t|$, we remain in the
# topological phase.  The ground-state energy interpolates between
# the sweet spot ($\Delta = t$, maximal gap) and the tight-binding
# limit ($\Delta = 0$, no superconductivity).
# ---

N₃ = 4
μ₃, t₃, Δ₃ = 0.0, 1.0, 0.5
registry₃, (a₃, a₃_dag) = create_fermionic_variables(1:N₃)

hopping₃ = -t₃ * sum(a₃_dag[i] * a₃[i+1] + a₃_dag[i+1] * a₃[i] for i in 1:(N₃ - 1))
pairing₃ =  Δ₃ * sum(a₃[i] * a₃[i+1] + a₃_dag[i+1] * a₃_dag[i] for i in 1:(N₃ - 1))
ham₃ = hopping₃ + pairing₃

pop₃ = polyopt(ham₃, registry₃)
config₃ = SolverConfig(optimizer = SOLVER, order = 1)
result₃ = cs_nctssos(pop₃, config₃)

exact₃ = kitaev_chain_exact(N₃; μ = μ₃, t = t₃, Δ = Δ₃)
@assert abs(result₃.objective - exact₃) < 1e-5   #src
println("Case 3 (asymmetric): SDP = $(result₃.objective),  exact = $(round(exact₃; digits=6)),  gap = $(abs(result₃.objective - exact₃))")

# ---
# ## Case 4 — Pure hopping ($\mu = 0$, $t = 1$, $\Delta = 0$)
#
# With no pairing at all, the Hamiltonian is a standard **tight-binding**
# chain — free fermions hopping on a wire.  The ground-state energy for
# $N = 4$ is $E_0 = -\sqrt{5}$.  No superconductivity, no topology;
# this is the trivial metal limit.
# ---

N₄ = 4
registry₄, (a₄, a₄_dag) = create_fermionic_variables(1:N₄)

hopping₄ = -sum(a₄_dag[i] * a₄[i+1] + a₄_dag[i+1] * a₄[i] for i in 1:(N₄ - 1))

pop₄ = polyopt(hopping₄, registry₄)
config₄ = SolverConfig(optimizer = SOLVER, order = 1)
result₄ = cs_nctssos(pop₄, config₄)

exact₄ = kitaev_chain_exact(N₄; μ = 0.0, t = 1.0, Δ = 0.0)
@assert abs(result₄.objective - exact₄) < 1e-5   #src
println("Case 4 (hopping):   SDP = $(result₄.objective),  exact = $(round(exact₄; digits=6)),  gap = $(abs(result₄.objective - exact₄))")

# ---
# ## Results summary
#
# | Case | $\mu$ | $t$ | $\Delta$ | Phase | Exact $E_0$ | SDP bound |
# |:-----|:-----:|:---:|:--------:|:------|:------------|:----------|
# | 1 — Sweet spot | 0 | 1 | 1 | Topological | $-3$ | tight |
# | 2 — Trivial    | 3 | 1 | 1 | Trivial     | $\approx -12.504$ | tight |
# | 3 — Asymmetric | 0 | 1 | 0.5 | Topological | $\approx -2.422$ | tight |
# | 4 — Hopping    | 0 | 1 | 0 | Trivial (metal) | $-\sqrt{5} \approx -2.236$ | tight |
#
# All four cases converge at **order 1** — the Kitaev chain with open
# boundaries and $N = 4$ sites is small enough for the first relaxation to
# capture the full ground state.

# ---
#
# ## Phase diagram
#
# The Kitaev chain has a sharp **topological phase transition** controlled
# by the ratio $|\mu|/(2|t|)$:
#
# | Regime | Phase | Edge modes? |
# |:-------|:------|:------------|
# | $|\mu| < 2|t|$, $\Delta \neq 0$ | **Topological** | Yes — unpaired Majorana zero modes at each end |
# | $|\mu| > 2|t|$ | **Trivial** | No |
# | $|\mu| = 2|t|$ | **Critical** | Gap closes — phase transition |
# | $\Delta = 0$ (any $\mu$, $t$) | **Trivial (metal)** | No — no superconductivity at all |
#
# Our four cases sample each regime.  Case 1 sits at the sweet spot deep
# in the topological phase; Case 2 is well into the trivial phase;
# Case 3 is topological but not at the sweet spot; Case 4 has no pairing,
# so no superconductivity and no topology.
#
# !!! tip "Why topological matters"
#     In the topological phase, the two Majorana zero modes ($c_1$ and
#     $c_{2N}$) encode a single fermionic degree of freedom that is
#     *delocalised* across the entire chain.  This nonlocal storage makes
#     the qubit robust against local perturbations — the foundation of
#     proposals for **topological quantum computing** with Majorana-based
#     qubits.

# ---
#
# ## Bonus: Sweet-spot via PauliAlgebra
#
# We showed above that at $\mu = 0$, $t = \Delta$, the Majorana form of
# the Kitaev chain maps under Jordan–Wigner to the open XX chain:
# $H = \sum_{j=1}^{N-1} X_j X_{j+1}$.  Let's verify this directly using
# the [`create_pauli_variables`](@ref) interface.

N₅ = 4
registry₅, (x₅, _, _) = create_pauli_variables(1:N₅)

ham₅ = sum(1.0 * x₅[i] * x₅[i+1] for i in 1:(N₅ - 1))

pop₅ = polyopt(ham₅, registry₅)
config₅ = SolverConfig(optimizer = SOLVER, order = 1)
result₅ = cs_nctssos(pop₅, config₅)

@assert abs(result₅.objective - exact₁) < 1e-5   #src
println("Pauli XX chain:  SDP = $(result₅.objective),  fermionic sweet spot = $exact₁")

# Same answer — $E_0 = -3$.  The Pauli and fermionic interfaces are two
# windows onto the same physics.  Use whichever is more natural for your
# model.

# ---
#
# ## Summary
#
# ### The physics in one paragraph
#
# The Kitaev chain is a 1D chain of spinless fermions with nearest-neighbour
# hopping and *p*-wave pairing.  The pairing term breaks particle-number
# conservation but preserves fermionic parity.  When $|\mu| < 2|t|$ and
# $\Delta \neq 0$, the chain is in a **topological superconducting phase**
# with unpaired Majorana zero modes at its ends.  At the sweet spot
# ($\mu = 0$, $t = \Delta$), the Majorana representation makes the physics
# transparent: the Hamiltonian pairs Majoranas across bonds, leaving two
# dangling zero modes.  Our SDP relaxation reproduces the exact
# ground-state energy across all four parameter regimes.
#
# ### The code recipe
#
# | Step | What we did | API |
# |:-----|:------------|:----|
# | 1 | Create fermionic operators (CAR built in) | [`create_fermionic_variables`](@ref) |
# | 2 | Write the Kitaev chain Hamiltonian (hopping + pairing + chemical potential) | standard Julia arithmetic |
# | 3 | Solve the SDP relaxation (no parity constraint — open BCs) | [`polyopt`](@ref), [`cs_nctssos`](@ref) |
# | 4 | Verify against exact diagonalisation | `kitaev_chain_exact` (above) |
# | 5 | Cross-check sweet spot via Pauli XX chain | [`create_pauli_variables`](@ref) |
#
# ### Key concepts
#
# | Concept | One-liner |
# |:--------|:----------|
# | **Kitaev chain** | 1D *p*-wave superconductor: hopping + pairing on an open wire |
# | **Pairing term** | $\Delta(a_i a_{i+1} + \text{h.c.})$ — creates/destroys fermion pairs, breaks $N_{\mathrm{tot}}$ conservation |
# | **Fermionic parity** | $P = (-1)^{N_{\mathrm{tot}}}$ — conserved even when particle number is not |
# | **Majorana operator** | $c = a + a^\dagger$ (self-conjugate, real fermion); two per site |
# | **Sweet spot** | $\mu = 0$, $t = \Delta$: Majoranas pair across bonds, zero modes appear |
# | **Topological phase** | $|\mu| < 2|t|$, $\Delta \neq 0$: unpaired Majorana edge modes |
# | **Trivial phase** | $|\mu| > 2|t|$: no edge modes, ordinary insulator/metal |
#
# ### Fermionic vs. Pauli: when to use which
#
# | Use case | Interface |
# |:---------|:----------|
# | Model written in $a$, $a^\dagger$ (Kitaev chain, Hubbard, BCS) | [`create_fermionic_variables`](@ref) |
# | Model written in Pauli $X$, $Y$, $Z$ (Heisenberg, Ising, XX) | [`create_pauli_variables`](@ref) |
# | Sweet-spot check / Jordan–Wigner cross-validation | Both — compare results |
#
# ### See also
#
# - [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
#   periodic BCs, parity sectors, free-fermion spectrum
# - [PBW Algebra Showcase](@ref pbw-algebras-showcase) — CAR details and
#   normal ordering
# - [Ground State Energy](@ref ground-state-energy) — Pauli-level spin
#   chain examples
#
# ### References
#
# - A. Yu. Kitaev, "Unpaired Majorana fermions in quantum wires,"
#   *Physics-Uspekhi* **44**, 131 (2001),
#   [arXiv:cond-mat/0010440](https://arxiv.org/abs/cond-mat/0010440).
# - J. Alicea, "New directions in the pursuit of Majorana fermions in solid
#   state systems," *Reports on Progress in Physics* **75**, 076501 (2012),
#   [arXiv:1202.1293](https://arxiv.org/abs/1202.1293).
