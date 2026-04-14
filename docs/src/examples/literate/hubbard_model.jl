# # [Hubbard Model (Canonical Particle-Number Constraints)](@id hubbard-model)
#
# This example computes **certified lower bounds** on the ground-state energy
# of the **1D Hubbard model** — the hydrogen atom of strongly correlated
# electron systems.  We solve the problem two ways:
#
# 1. **Grand-canonical** — no particle-number constraint (optimise over the
#    entire Fock space).
# 2. **Canonical** — fix the number of spin-up and spin-down electrons
#    (half-filling).
#
# The canonical case introduces a new API keyword, `moment_eq_constraints`,
# and the difference between the two cases teaches a key lesson about
# **state constraints vs. operator identities** in SDP relaxations.
#
# If you haven't seen fermionic operators before, start with the
# [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) example,
# which builds creation/annihilation operators and the CAR from scratch.
# Here we assume that background and focus on what's new: **spin degrees of
# freedom**, **doubled fermionic modes**, and **moment equality constraints**
# for fixing particle number.

# ---
#
# ## The 1D Hubbard Hamiltonian
#
# Picture $N$ sites on a ring.  Each site can hold up to **two** electrons —
# one spin-up ($\uparrow$) and one spin-down ($\downarrow$).  The
# Hamiltonian has two competing terms:
#
# ```math
# H = -t \sum_{\langle i,j \rangle, \sigma}
#     \bigl(c_{i\sigma}^\dagger c_{j\sigma} + \text{h.c.}\bigr)
#     \;+\; U \sum_{i=1}^{N} n_{i\uparrow}\, n_{i\downarrow},
# ```
#
# where $n_{i\sigma} = c_{i\sigma}^\dagger c_{i\sigma}$ counts electrons of
# spin $\sigma$ at site $i$, and $\langle i,j \rangle$ runs over
# nearest-neighbour pairs on the ring.
#
# | Term | Physics | Favours |
# |:-----|:--------|:--------|
# | $-t\,(c_{i\sigma}^\dagger c_{j\sigma} + \text{h.c.})$ | **Hopping** — kinetic energy, electron tunnels between neighbours | Delocalisation (metallic) |
# | $U\, n_{i\uparrow}\, n_{i\downarrow}$ | **On-site repulsion** — Coulomb energy when two electrons share a site | Localisation (insulating) |
#
# The ground-state physics is controlled by the ratio $U/t$:
# - **Small $U/t$**: electrons delocalise freely — a metal.
# - **Large $U/t$**: double occupancy is too costly — electrons freeze
#   in place, forming a **Mott insulator**.
# - At **half-filling** (one electron per site on average) and large
#   $U/t$, virtual hopping generates an effective antiferromagnetic
#   exchange $J \sim 4t^2/U$ — the **Heisenberg limit**.

# ---
#
# ## Spin and doubled fermionic modes
#
# Each physical site holds **two** fermionic modes (spin-up and spin-down).
# In `NCTSSoS.jl` we represent both species in **one** fermionic algebra
# on **disjoint mode ranges**:
#
# | Mode index | Physical meaning |
# |:-----------|:-----------------|
# | $1, 2, \ldots, N$ | Spin-up: $c_{1\uparrow}, c_{2\uparrow}, \ldots, c_{N\uparrow}$ |
# | $N{+}1, N{+}2, \ldots, 2N$ | Spin-down: $c_{1\downarrow}, c_{2\downarrow}, \ldots, c_{N\downarrow}$ |
#
# Because all $2N$ modes share one algebra, they **all anticommute** with
# each other — including cross-species:
# $\{c_{i\uparrow},\, c_{j\downarrow}^\dagger\} = 0$, exactly as
# required by the full fermionic CAR.  The "disjoint" part refers only to
# the mode *indices*, not to the algebra.
#
# Why not use a single mode per site?  Because the on-site repulsion
# $n_{i\uparrow}\, n_{i\downarrow} = (c_{i\uparrow}^\dagger c_{i\uparrow})(c_{i\downarrow}^\dagger c_{i\downarrow})$
# involves **four** operators on **distinct** mode indices — it is
# genuinely **quartic**.  If up and down shared the same mode index, this
# product would collapse via $n_i^2 = n_i$ and we'd lose the two-body
# physics entirely.

# ---
#
# ## Grand-canonical vs. canonical
#
# | Ensemble | Constraint | Optimises over |
# |:---------|:-----------|:---------------|
# | **Grand-canonical** | None on particle number | Entire Fock space ($4^N$ states) |
# | **Canonical** | $N_\uparrow, N_\downarrow$ fixed | Single particle-number sector |
#
# The grand-canonical ground state may live in a **different**
# particle-number sector than the one of physical interest.  For
# our $N = 4$ ring with $U/t = 4$, the global ground state has
# $E_0 \approx -3.419$ in the dilute sector $(N_\uparrow, N_\downarrow) = (1, 1)$,
# while the half-filled sector $(N_\uparrow = N_\downarrow = 2)$ has
# $E_0 \approx -2.103$.  If you care about the half-filled Mott physics,
# you **must** constrain the particle number.
#
# !!! tip "When to use canonical constraints"
#     Real experiments almost always fix the electron count — electrons
#     don't appear or disappear in a solid.  Use canonical constraints
#     whenever you know which particle-number sector the physics lives in.

# ---
#
# ## Exact diagonalisation oracle
#
# To verify our SDP bounds we build the full $2^{2N} \times 2^{2N}$
# Hamiltonian matrix via the **Jordan–Wigner** mapping and diagonalise it.
# This is brute-force exponential scaling — it works for small $N$, not
# for real materials.

using LinearAlgebra

## Jordan–Wigner building blocks
# We use the occupation-basis convention |0⟩ = empty, |1⟩ = occupied,
# so bit counts in the computational basis equal literal particle numbers.
const PAULI_Z = ComplexF64[1 0; 0 -1]
const PAULI_I = Matrix{ComplexF64}(I, 2, 2)
const JW_CREATE = ComplexF64[0 0; 1 0]   # |1⟩⟨0|
const JW_DESTROY = ComplexF64[0 1; 0 0]  # |0⟩⟨1|

function jw_fermion_op(mode::Int, kind::Symbol, nmodes::Int)
    mats = [site < mode ? PAULI_Z :
            site == mode ? (kind === :annihilate ? JW_DESTROY : JW_CREATE) :
            PAULI_I for site in 1:nmodes]
    return reduce(kron, mats)
end

## Build the full Hamiltonian matrix directly from JW operators.
function hubbard_exact_matrix(nsites::Int; t::Real = 1.0, U::Real = 4.0)
    nmodes = 2 * nsites
    a  = [jw_fermion_op(i, :annihilate, nmodes) for i in 1:nmodes]
    ad = [jw_fermion_op(i, :create,     nmodes) for i in 1:nmodes]
    dim = 2^nmodes
    H = zeros(ComplexF64, dim, dim)
    for i in 1:nsites
        j = mod1(i + 1, nsites)
        H .-= t .* (ad[i] * a[j] + ad[j] * a[i])                                     ## spin-up hop
        H .-= t .* (ad[i + nsites] * a[j + nsites] + ad[j + nsites] * a[i + nsites])  ## spin-down hop
    end
    for i in 1:nsites
        H .+= U .* (ad[i] * a[i]) * (ad[i + nsites] * a[i + nsites])  ## on-site n↑n↓
    end
    return Hermitian(H)
end

## Project onto a particle-number sector and return the ground-state energy.
function sector_ground_state_energy(H::AbstractMatrix, nsites::Int; nup::Int, ndn::Int)
    nmodes = 2 * nsites
    upmask = (UInt(1) << nsites) - 1
    keep = [s + 1 for s in 0:(2^nmodes - 1)
            if count_ones(UInt(s) & upmask) == nup &&
               count_ones(UInt(s) >> nsites) == ndn]
    return eigmin(Hermitian(H[keep, keep]))
end;

# Now compute the exact energies for $N = 4$, $t = 1$, $U = 4$:

N = 4
t = 1.0
U = 4.0

H_mat = hubbard_exact_matrix(N; t, U)
E_full = eigmin(H_mat)
E_dilute = sector_ground_state_energy(H_mat, N; nup = 1, ndn = 1)
E_half = sector_ground_state_energy(H_mat, N; nup = 2, ndn = 2)

@assert abs(E_full - (-3.4185507188738447)) < 1e-12   #src
@assert abs(E_dilute - E_full) < 1e-12                #src
@assert abs(E_half - (-2.102748483462074)) < 1e-12    #src
println("Full Fock space GS:    E₀ = $(round(E_full; digits=6))")
println("Dilute (1,1) GS:       E₀ = $(round(E_dilute; digits=6))")
println("Half-filled (2,2) GS:  E₀ = $(round(E_half; digits=6))")

# The global ground state (energy $\approx -3.42$) sits in the dilute
# $(N_\uparrow, N_\downarrow) = (1, 1)$ sector, not at half-filling.  The
# half-filled sector ($E_0 \approx -2.10$) is the physically relevant one
# for Mott physics.

# ---
#
# ## Building the Hamiltonian with NCTSSoS.jl
#
# We now build the *same* Hamiltonian using the polynomial interface.
# [`create_fermionic_variables`](@ref) returns creation and annihilation
# operators that already obey all CAR — including cross-species
# anticommutation.

using NCTSSoS

registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
    ("c_up", 1:N),   ## modes 1..N   (spin-up)
    ("c_dn", 1:N),   ## modes N+1..2N (spin-down)
]);
println("Created $(length(registry)) variables: $(2N) fermionic modes ($(N) up + $(N) down)")

# ### Hopping — kinetic energy
#
# Each spin species hops independently around the periodic ring.

bonds = [(i, mod1(i + 1, N)) for i in 1:N]

hopping = -t * sum(
    c_up_dag[i] * c_up[j] + c_up_dag[j] * c_up[i] +
    c_dn_dag[i] * c_dn[j] + c_dn_dag[j] * c_dn[i]
    for (i, j) in bonds
);

# ### On-site repulsion — Coulomb energy
#
# The product $n_{i\uparrow} n_{i\downarrow}$ uses operators from both
# mode ranges, keeping it genuinely quartic.

interaction = U * sum(
    (c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
    for i in 1:N
)

ham = hopping + interaction   ## quadratic hops + quartic on-site terms

# ---
#
# ## Setting up the SDP solver

# !!! note "SDP Solver"
#     These examples use [Mosek](https://www.mosek.com/) via `MosekTools`.
#     Any SDP-capable solver works: replace `Mosek.Optimizer` with
#     `COSMO.Optimizer` or `Clarabel.Optimizer` for open-source alternatives.

using MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG" => 0,
    "MSK_IPAR_NUM_THREADS" => 0);

# This quartic Hubbard example is a good reminder that sparsity is not magic.
# Plain order-2 dense relaxations leave visible gaps, while dense order 3 is
# essentially exact but much slower.  The practical compromise here is to use
# **term sparsity only** with one refinement step.  The grand-canonical case is
# already tight at order 3, while the canonical half-filled case benefits from
# pushing to order 4.

config_gc = SolverConfig(optimizer = SOLVER, order = 3, ts_algo = MMD())
config_can = SolverConfig(optimizer = SOLVER, order = 4, ts_algo = MMD())

# ---
#
# ## Case 1 — Grand-canonical (no particle-number constraint)
#
# Without constraining particle number, the SDP relaxation searches for the
# lowest energy across **all** particle-number sectors simultaneously.  The
# result is a valid lower bound on the global ground-state energy.

pop_gc = polyopt(ham, registry)
result_gc_first = cs_nctssos(pop_gc, config_gc)
result_gc = cs_nctssos_higher(pop_gc, result_gc_first, config_gc)

@assert result_gc.objective ≤ E_full + 1e-5   #src
@assert abs(result_gc.objective - E_full) < 1e-6   #src
println("Grand-canonical SDP (first sparse pass):  $(round(result_gc_first.objective; digits=6))")
println("Grand-canonical SDP (refined):            $(round(result_gc.objective; digits=6))")
println("Exact full-Fock GS:                       $(round(E_full; digits=6))")
println("Gap after refinement:                     $(round(abs(result_gc.objective - E_full); digits=6))")

# The first sparse pass is still loose, but one call to
# [`cs_nctssos_higher`](@ref) refines the term-sparsity graph and drives the
# bound essentially onto the exact value.

# ---
#
# ## Why canonical constraints need special treatment
#
# To restrict to half-filling, we want to impose
# $(\hat{N}_\uparrow - 2)|\psi\rangle = 0$ and
# $(\hat{N}_\downarrow - 2)|\psi\rangle = 0$.  This looks like an
# equality constraint, so you might reach for `eq_constraints`.
# **That would be wrong.**  Here's why.
#
# !!! warning "`eq_constraints` vs. `moment_eq_constraints`"
#     [`polyopt`](@ref) offers two kinds of equality constraint:
#
#     **`eq_constraints`** treats ``g = 0`` as an **operator identity** —
#     true in *every* quantum state.  Internally it builds a full bilinear
#     localizing matrix, imposing
#     $\langle b_i^\dagger \, g \, b_j \rangle = 0$ for all basis pairs
#     $(b_i, b_j)$.  Use this for identities like the CAR or $P^2 = I$.
#     Do **not** use it for sector-fixing conditions such as
#     $P|\psi\rangle = |\psi\rangle$.
#
#     **`moment_eq_constraints`** treats ``g|\psi\rangle = 0`` as a **state
#     constraint** — true for the *target state*, not universally.  It
#     builds the one-sided localizing form
#     $\langle b_i^\dagger \, g \rangle = 0$ for each basis element $b_i$.
#
#     The difference matters when $g$ changes particle number.  The
#     bilinear form $\langle b_i^\dagger \, g \, b_j \rangle = 0$
#     constrains matrix elements **between different particle-number
#     sectors** — elements that the true state never connects.  This
#     overconstrains the SDP and produces garbage (for these parameters,
#     objective $\approx +2$ instead of $\approx -2$).
#
#     **Rule of thumb:** if ``g = 0`` is a law of the algebra (holds for
#     every state), use `eq_constraints`.  If ``g|\psi\rangle = 0`` is a
#     property of the state you're looking for, use
#     `moment_eq_constraints`.

# ---
#
# ## Case 2 — Canonical half-filling ($N_\uparrow = N_\downarrow = 2$)
#
# We build the total number operators and impose half-filling as moment
# equality constraints.

n_up_total = 1.0 * sum(c_up_dag[i] * c_up[i] for i in 1:N)  ## 1.0 * promotes to float coefficients
n_dn_total = 1.0 * sum(c_dn_dag[i] * c_dn[i] for i in 1:N)

pop_can = polyopt(ham, registry;
    moment_eq_constraints = [n_up_total - 2.0 * one(ham),
                             n_dn_total - 2.0 * one(ham)])
result_can_first = cs_nctssos(pop_can, config_can)
result_can = cs_nctssos_higher(pop_can, result_can_first, config_can)

@assert result_can.objective ≤ E_half + 1e-4   #src
@assert abs(result_can.objective - E_half) < 1e-6   #src
@assert result_can.objective > -3.0             #src
println("Canonical SDP (first sparse pass):        $(round(result_can_first.objective; digits=6))")
println("Canonical SDP (refined):                  $(round(result_can.objective; digits=6))")
println("Exact half-filled GS:                     $(round(E_half; digits=6))")
println("Gap after refinement:                     $(round(abs(result_can.objective - E_half); digits=6))")

# The canonical bound tracks the half-filled sector — it is *not* pulled
# down to the lower full-Fock minimum.  At order 4, one refinement step drives
# the gap down to about ${\sim}10^{-7}$.

# ---
#
# ## Results
#
# | Quantity | Value |
# |:--------|------:|
# | Exact full-Fock $E_0$ | $-3.419$ |
# | Grand-canonical SDP (order 3 TS + higher) | $\approx -3.419$ |
# | Exact half-filled $E_0\;(N_\uparrow{=}N_\downarrow{=}2)$ | $-2.103$ |
# | Canonical SDP (order 4 TS + higher) | $\approx -2.103$ |
#
# Both SDP values are **valid lower bounds** on their respective exact
# energies.  For this example, the practical sweet spot is not plain dense
# order 2, and it is not naive correlative sparsity either.  Using **term
# sparsity plus one refinement step**, the grand-canonical case is already
# essentially exact at order 3, while the canonical case reaches numerical
# precision at order 4.

# ---
#
# ## Summary
#
# ### The physics in one paragraph
#
# The 1D Hubbard model captures the competition between electron hopping
# ($t$) and on-site Coulomb repulsion ($U$).  We solved $N = 4$ sites on a
# periodic ring at $U/t = 4$ — firmly in the correlated regime — both
# without and with particle-number constraints.  The grand-canonical SDP
# finds a lower bound across all Fock sectors; the canonical relaxation,
# using `moment_eq_constraints` to fix $N_\uparrow = N_\downarrow = 2$,
# targets the half-filled Mott sector specifically.  Numerically, this is a
# case where **term sparsity plus one refinement step** works well, but the
# two sectors are not equally hard: the grand-canonical problem is already
# essentially exact at order 3, while the canonical half-filled problem needs
# order 4 to close the gap to numerical precision.
#
# ### The code recipe
#
# | Step | What we did | API |
# |:-----|:------------|:----|
# | 1 | Create spin-up and spin-down operators on disjoint mode ranges | [`create_fermionic_variables`](@ref) with species list |
# | 2 | Build the Hubbard Hamiltonian (hopping + on-site repulsion) | Standard Julia arithmetic |
# | 3 | Grand-canonical solve with order-3 term sparsity refinement | [`polyopt`](@ref) $\to$ [`cs_nctssos`](@ref) $\to$ [`cs_nctssos_higher`](@ref) |
# | 4 | Canonical solve with fixed particle number, order 4, and one refinement | `moment_eq_constraints` in [`polyopt`](@ref), then [`cs_nctssos`](@ref) $\to$ [`cs_nctssos_higher`](@ref) |
# | 5 | Verify against exact diagonalisation | `hubbard_exact_matrix` + `sector_ground_state_energy` |
#
# ### Key concepts
#
# | Concept | One-liner |
# |:--------|:----------|
# | **Hubbard model** | Hopping $t$ vs. on-site repulsion $U$ — simplest model of strongly correlated electrons |
# | **Doubled modes** | Spin-up and spin-down as disjoint mode ranges in one registry ($2N$ modes for $N$ sites) |
# | **Grand-canonical** | No particle-number constraint — optimise over entire Fock space |
# | **Canonical** | Fix $N_\uparrow$ and $N_\downarrow$ — optimise within a single particle-number sector |
# | **`moment_eq_constraints`** | One-sided localizing: $\langle b_i^\dagger g \rangle = 0$ — for state constraints (``g|\psi\rangle = 0``) |
# | **`eq_constraints`** | Full bilinear localizing: $\langle b_i^\dagger g\, b_j \rangle = 0$ — for operator identities (``g = 0`` always) |
# | **Mott insulator** | Large $U/t$ at half-filling: electrons localise, charge transport freezes |
#
# ### See also
#
# - [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
#   creation/annihilation operators, CAR, parity constraints
# - [Kitaev Chain](@ref kitaev-chain) — pairing terms, Majorana
#   operators, topological phases
# - [Hubbard Model (Filling-Sector Scan)](@ref hubbard-filling-scan) —
#   what the unconstrained SDP does, and how it compares with exact
#   diagonalisation across all $(N_\uparrow, N_\downarrow)$ sectors
# - [Ground State Energy](@ref ground-state-energy) — Pauli-level spin
#   chain examples
#
# ### References
#
# - J. Hubbard, "Electron correlations in narrow energy bands,"
#   *Proceedings of the Royal Society A* **276**, 238 (1963).
# - E. H. Lieb and F. Y. Wu, "Absence of Mott transition in an exact
#   solution of the short-range, one-band model in one dimension,"
#   *Physical Review Letters* **20**, 1445 (1968).
# - F. H. L. Essler, H. Frahm, F. Göhmann, A. Klümper, and V. E.
#   Korepin, *The One-Dimensional Hubbard Model*, Cambridge University
#   Press (2005).
