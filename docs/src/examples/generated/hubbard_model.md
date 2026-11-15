```@meta
EditURL = "../literate/hubbard_model.jl"
```

# [Hubbard Model (Canonical Particle-Number Constraints)](@id hubbard-model)

This example computes **certified lower bounds** on the ground-state energy
of the **1D Hubbard model** — the hydrogen atom of strongly correlated
electron systems.  We solve the problem two ways:

1. **Grand-canonical** — no particle-number constraint (optimise over the
   entire Fock space).
2. **Canonical** — fix the number of spin-up and spin-down electrons
   (half-filling).

The canonical case uses `Sector` for each fixed particle-number sector. The
difference between the two cases teaches a key lesson about **state
constraints vs. operator identities** in SDP relaxations.

If you haven't seen fermionic operators before, start with the
[Fermionic Ground State (XY Model)](@ref fermionic-ground-state) example,
which builds creation/annihilation operators and the CAR from scratch.
Here we assume that background and focus on what's new: **spin degrees of
freedom**, **doubled fermionic modes**, and **state constraints**
for fixing particle number.

---

## The 1D Hubbard Hamiltonian

Picture $N$ sites on a ring.  Each site can hold up to **two** electrons —
one spin-up ($\uparrow$) and one spin-down ($\downarrow$).  The
Hamiltonian has two competing terms:

```math
H = -t \sum_{\langle i,j \rangle, \sigma}
    \bigl(c_{i\sigma}^\dagger c_{j\sigma} + \text{h.c.}\bigr)
    \;+\; U \sum_{i=1}^{N} n_{i\uparrow}\, n_{i\downarrow},
```

where $n_{i\sigma} = c_{i\sigma}^\dagger c_{i\sigma}$ counts electrons of
spin $\sigma$ at site $i$, and $\langle i,j \rangle$ runs over
nearest-neighbour pairs on the ring.

| Term | Physics | Favours |
|:-----|:--------|:--------|
| $-t\,(c_{i\sigma}^\dagger c_{j\sigma} + \text{h.c.})$ | **Hopping** — kinetic energy, electron tunnels between neighbours | Delocalisation (metallic) |
| $U\, n_{i\uparrow}\, n_{i\downarrow}$ | **On-site repulsion** — Coulomb energy when two electrons share a site | Localisation (insulating) |

The ground-state physics is controlled by the ratio $U/t$:
- **Small $U/t$**: electrons delocalise freely — a metal.
- **Large $U/t$**: double occupancy is too costly — electrons freeze
  in place, forming a **Mott insulator**.
- At **half-filling** (one electron per site on average) and large
  $U/t$, virtual hopping generates an effective antiferromagnetic
  exchange $J \sim 4t^2/U$ — the **Heisenberg limit**.

---

## Spin and doubled fermionic modes

Each physical site holds **two** fermionic modes (spin-up and spin-down).
In `NCTSSoS.jl` we represent both species in **one** fermionic algebra
on **disjoint mode ranges**:

| Mode index | Physical meaning |
|:-----------|:-----------------|
| $1, 2, \ldots, N$ | Spin-up: $c_{1\uparrow}, c_{2\uparrow}, \ldots, c_{N\uparrow}$ |
| $N{+}1, N{+}2, \ldots, 2N$ | Spin-down: $c_{1\downarrow}, c_{2\downarrow}, \ldots, c_{N\downarrow}$ |

Because all $2N$ modes share one algebra, they **all anticommute** with
each other — including cross-species:
$\{c_{i\uparrow},\, c_{j\downarrow}^\dagger\} = 0$, exactly as
required by the full fermionic CAR.  The "disjoint" part refers only to
the mode *indices*, not to the algebra.

Why not use a single mode per site?  Because the on-site repulsion
$n_{i\uparrow}\, n_{i\downarrow} = (c_{i\uparrow}^\dagger c_{i\uparrow})(c_{i\downarrow}^\dagger c_{i\downarrow})$
involves **four** operators on **distinct** mode indices — it is
genuinely **quartic**.  If up and down shared the same mode index, this
product would collapse via $n_i^2 = n_i$ and we'd lose the two-body
physics entirely.

---

## Grand-canonical vs. canonical

| Ensemble | Constraint | Optimises over |
|:---------|:-----------|:---------------|
| **Grand-canonical** | None on particle number | Entire Fock space ($4^N$ states) |
| **Canonical** | $N_\uparrow, N_\downarrow$ fixed | Single particle-number sector |

The grand-canonical ground state may live in a **different**
particle-number sector than the one of physical interest.  For
our $N = 4$ ring with $U/t = 4$, the global ground state has
$E_0 \approx -3.419$ in the dilute sector $(N_\uparrow, N_\downarrow) = (1, 1)$,
while the half-filled sector $(N_\uparrow = N_\downarrow = 2)$ has
$E_0 \approx -2.103$.  If you care about the half-filled Mott physics,
you **must** constrain the particle number.

!!! tip "When to use canonical constraints"
    Real experiments almost always fix the electron count — electrons
    don't appear or disappear in a solid.  Use canonical constraints
    whenever you know which particle-number sector the physics lives in.

---

## Exact diagonalisation oracle

To verify our SDP bounds we build the full $2^{2N} \times 2^{2N}$
Hamiltonian matrix via the **Jordan–Wigner** mapping and diagonalise it.
This is brute-force exponential scaling — it works for small $N$, not
for real materials.

````julia
using LinearAlgebra

# Jordan–Wigner building blocks
````

We use the occupation-basis convention |0⟩ = empty, |1⟩ = occupied,
so bit counts in the computational basis equal literal particle numbers.

````julia
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

# Build the full Hamiltonian matrix directly from JW operators.
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

# Project onto a particle-number sector and return the ground-state energy.
function sector_ground_state_energy(H::AbstractMatrix, nsites::Int; nup::Int, ndn::Int)
    nmodes = 2 * nsites
    upmask = (UInt(1) << nsites) - 1
    keep = [s + 1 for s in 0:(2^nmodes - 1)
            if count_ones(UInt(s) & upmask) == nup &&
               count_ones(UInt(s) >> nsites) == ndn]
    return eigmin(Hermitian(H[keep, keep]))
end;
````

Now compute the exact energies for $N = 4$, $t = 1$, $U = 4$:

````julia
N = 4
t = 1.0
U = 4.0

H_mat = hubbard_exact_matrix(N; t, U)
E_full = eigmin(H_mat)
E_dilute = sector_ground_state_energy(H_mat, N; nup = 1, ndn = 1)
E_half = sector_ground_state_energy(H_mat, N; nup = 2, ndn = 2)

println("Full Fock space GS:    E₀ = $(round(E_full; digits=6))")
println("Dilute (1,1) GS:       E₀ = $(round(E_dilute; digits=6))")
println("Half-filled (2,2) GS:  E₀ = $(round(E_half; digits=6))")
````

````
Full Fock space GS:    E₀ = -3.418551
Dilute (1,1) GS:       E₀ = -3.418551
Half-filled (2,2) GS:  E₀ = -2.102748

````

The global ground state (energy $\approx -3.42$) sits in the dilute
$(N_\uparrow, N_\downarrow) = (1, 1)$ sector, not at half-filling.  The
half-filled sector ($E_0 \approx -2.10$) is the physically relevant one
for Mott physics.

---

## Building the Hamiltonian with NCTSSoS.jl

We now build the *same* Hamiltonian using the polynomial interface.
`fermionic_variables` returns creation and annihilation
operators that already obey all CAR — including cross-species
anticommutation.

````julia
using NCTSSoS.Next
import NCTSSoS.Next as Next

registry, (c, c_dag) = fermionic_variables(1:2N);
c_up, c_up_dag = c[1:N], c_dag[1:N];
c_dn, c_dn_dag = c[N+1:2N], c_dag[N+1:2N];
println("Created $(length(registry)) variables: $(2N) fermionic modes ($(N) up + $(N) down)")
````

````
Created 16 variables: 8 fermionic modes (4 up + 4 down)

````

### Hopping — kinetic energy

Each spin species hops independently around the periodic ring.

````julia
bonds = [(i, mod1(i + 1, N)) for i in 1:N]

hopping = -t * sum(
    c_up_dag[i] * c_up[j] + c_up_dag[j] * c_up[i] +
    c_dn_dag[i] * c_dn[j] + c_dn_dag[j] * c_dn[i]
    for (i, j) in bonds
);
````

### On-site repulsion — Coulomb energy

The product $n_{i\uparrow} n_{i\downarrow}$ uses operators from both
mode ranges, keeping it genuinely quartic.

````julia
interaction = U * sum(
    (c_up_dag[i] * c_up[i]) * (c_dn_dag[i] * c_dn[i])
    for i in 1:N
)

ham = hopping + interaction   ## quadratic hops + quartic on-site terms
````

````
(-1.0) ad1 a2 + (-1.0) ad1 a4 + (-1.0) ad2 a1 + (-1.0) ad2 a3 + (-1.0) ad3 a2 + (-1.0) ad3 a4 + (-1.0) ad4 a1 + (-1.0) ad4 a3 + … (20 terms)
````

---

## Setting up the SDP solver

!!! note "SDP Solver"
    These examples use [Mosek](https://www.mosek.com/) through `MosekBackend`.
    For other SDP solvers, see [Optimizers](@ref overview-of-optimizers).

````julia
using Mosek

SOLVER = MosekBackend();
````

This quartic Hubbard example is a good reminder that sparsity is not magic.
Plain order-2 dense relaxations leave visible gaps, while dense order 3 is
essentially exact but much slower.  The practical compromise here is to use
**term sparsity only** with one refinement step.  The grand-canonical case is
tight at order 3. With support propagated from the completed graphs, the
canonical half-filled case also reaches the exact energy at order 3.

````julia
relaxation_gc_first = Relaxation(order=3, split=TermSplit(algorithm=Next.MMD(), level=1))
relaxation_gc = Relaxation(order=3, split=TermSplit(algorithm=Next.MMD(), level=2))
````

````
Relaxation(order=3, split=TermSplit(algorithm=MMD(), level=2))
````

---

## Case 1 — Grand-canonical (no particle-number constraint)

Without constraining particle number, the SDP relaxation searches for the
lowest energy across **all** particle-number sectors simultaneously.  The
result is a valid lower bound on the global ground-state energy.

````julia
problem = Problem(ham)
result_gc_first = solve(problem, relaxation_gc_first; solver=SOLVER)
result_gc = solve(problem, relaxation_gc; solver=SOLVER)

println("Grand-canonical SDP (first sparse pass):  $(round(result_gc_first.bound; digits=6))")
println("Grand-canonical SDP (refined):            $(round(result_gc.bound; digits=6))")
println("Exact full-Fock GS:                       $(round(E_full; digits=6))")
println("Gap after refinement:                     $(round(abs(result_gc.bound - E_full); digits=6))")
````

````
Grand-canonical SDP (first sparse pass):  -4.95709
Grand-canonical SDP (refined):            -3.418551
Exact full-Fock GS:                       -3.418551
Gap after refinement:                     0.0

````

The first sparse pass is still loose. `level=2` adds one support-extension
and completion step and drives the bound essentially onto the exact value.
It rebuilds the graph from the problem; it does not read the first solution.
The new and earlier sparse programs need not be identical at the same
level. Here the first bound changed from −4.962560 to −4.957090; the
refined bounds agree within 1e-6. `MMD()` depends on the row order.
The new row order selects a different chordal extension, so the first
pass uses different blocks. This does not change the dense relaxation.
Each sparse block is a principal submatrix of the dense moment matrix.

---

## Why canonical constraints need special treatment

To restrict to half-filling, we want to impose
$(\hat{N}_\uparrow - 2)|\psi\rangle = 0$ and
$(\hat{N}_\downarrow - 2)|\psi\rangle = 0$.  This looks like an
equality constraint, so you might reach for `Zero`.
**That would be wrong.**  Here's why.

!!! warning "`Zero` vs. `Sector`"
    Operator and state constraints have different meanings:

    **`Zero(g)` in `Problem(constraints=...)`** treats ``g = 0`` as an **operator identity** —
    true in *every* quantum state.  Internally it builds a full bilinear
    localizing matrix, imposing
    $\langle b_i^\dagger \, g \, b_j \rangle = 0$ for all basis pairs
    $(b_i, b_j)$.  Use this for identities like the CAR or $P^2 = I$.
    Do **not** use it for sector-fixing conditions such as
    $P|\psi\rangle = |\psi\rangle$.

    **`Sector(Q, q)` in `Relaxation(strengthen=...)`** treats
    ``g|\psi\rangle = 0``, with ``g=Q-q``, as a **state
    constraint** — true for the *target state*, not universally.  It
    imposes zero variance $\langle g^2\rangle=0$ and left-product
    equations $\langle g\,u\rangle=0$ for the represented moment words $u$.

    The difference matters when a basis operator changes particle number. The
    bilinear form $\langle b_i^\dagger \, g \, b_j \rangle = 0$
    constrains matrix elements **between different particle-number
    sectors** — elements that the true state never connects.  This
    overconstrains the SDP. The earlier page reported an incorrect
    objective near $+2$ instead of $-2$ with an operator-equality setup;
    that is not a reference value for a sector solve.

    **Rule of thumb:** if ``g = 0`` is a law of the algebra (holds for
    every state), use `Zero(g)`. If ``(Q-q)|\psi\rangle = 0`` is a
    property of the state you're looking for, use
    `Sector(Q, q)`.

---

## Case 2 — Canonical half-filling ($N_\uparrow = N_\downarrow = 2$)

We impose half-filling with two `Sector` conditions. Each number operator
is a sum of occupations, one for spin up and one for spin down.

````julia
number_up = sum(c_up_dag[i] * c_up[i] for i in 1:N)
number_dn = sum(c_dn_dag[i] * c_dn[i] for i in 1:N)
canonical_constraints = (Sector(number_up, 2), Sector(number_dn, 2))

relaxation_can_first = Relaxation(order=3,
    split=TermSplit(algorithm=Next.MMD(), level=1), strengthen=canonical_constraints)
relaxation_can = Relaxation(order=3,
    split=TermSplit(algorithm=Next.MMD(), level=2), strengthen=canonical_constraints)
result_can_first = solve(problem, relaxation_can_first; solver=SOLVER)
result_can = solve(problem, relaxation_can; solver=SOLVER)

println("Canonical SDP (first sparse pass):        $(round(result_can_first.bound; digits=6))")
println("Canonical SDP (refined):                  $(round(result_can.bound; digits=6))")
println("Exact half-filled GS:                     $(round(E_half; digits=6))")
println("Gap after refinement:                     $(abs(result_can.bound - E_half))")
````

````
Canonical SDP (first sparse pass):        -2.341756
Canonical SDP (refined):                  -2.102748
Exact half-filled GS:                     -2.102748
Gap after refinement:                     7.506466559448199e-10

````

`Sector` supplies the state equations without requiring them in every model. The canonical
bound tracks the half-filled sector — it is *not* pulled down to the lower
full-Fock minimum. At order 3, the measured first bound is about −2.341756;
one refinement gives −2.102748484213, with a gap of about 7.5×10⁻¹⁰ to
the exact half-filled energy. The earlier first bound was −3.631145 and
its refined bound was −2.102748483464. `Sector` adds variance and
left-product equations. The first-pass blocks also change because
`MMD()` depends on the row order, as in the full-Fock case above.
The refined bounds agree within the unchanged 1e-6 tolerance.

---

## Results

| Quantity | Value |
|:--------|------:|
| Exact full-Fock $E_0$ | $-3.419$ |
| Grand-canonical SDP (order 3 TS + higher) | $\approx -3.419$ |
| Exact half-filled $E_0\;(N_\uparrow{=}N_\downarrow{=}2)$ | $-2.103$ |
| Canonical SDP (order 3 TS + higher) | $\approx -2.103$ |

Both SDP values are **valid lower bounds** on their respective exact
energies up to solver tolerance. Using **term sparsity plus one refinement
step**, both sectors reach numerical agreement with exact diagonalisation
at order 3.

---

## Summary

### The physics in one paragraph

The 1D Hubbard model captures the competition between electron hopping
($t$) and on-site Coulomb repulsion ($U$).  We solved $N = 4$ sites on a
periodic ring at $U/t = 4$ — firmly in the correlated regime — both
without and with particle-number constraints.  The grand-canonical SDP
finds a lower bound across all Fock sectors; the canonical relaxation,
using `Sector(number_up, 2)` and `Sector(number_dn, 2)` to fix
$N_\uparrow = N_\downarrow = 2$, targets the half-filled Mott sector
specifically. Here **term sparsity plus one refinement step** at order 3
reaches both exact energies to within numerical tolerance.

### The code recipe

| Step | What we did | API |
|:-----|:------------|:----|
| 1 | Create spin-up and spin-down operators on disjoint mode ranges | `fermionic_variables(1:2N)`, then split the returned vectors |
| 2 | Build the Hubbard Hamiltonian (hopping + on-site repulsion) | Standard Julia arithmetic |
| 3 | Grand-canonical solve with order-3 term sparsity refinement | `Problem`, `Relaxation`, `solve`; `TermSplit(level=1)` then `level=2`, both with `Next.MMD()` |
| 4 | Canonical solve with fixed particle number, order 3, and one refinement | Two `Sector` conditions in `Relaxation(strengthen=...)` |
| 5 | Verify against exact diagonalisation | `hubbard_exact_matrix` + `sector_ground_state_energy` |

### Key concepts

| Concept | One-liner |
|:--------|:----------|
| **Hubbard model** | Hopping $t$ vs. on-site repulsion $U$ — simplest model of strongly correlated electrons |
| **Doubled modes** | Spin-up and spin-down as disjoint mode ranges in one registry ($2N$ modes for $N$ sites) |
| **Grand-canonical** | No particle-number constraint — optimise over entire Fock space |
| **Canonical** | Fix $N_\uparrow$ and $N_\downarrow$ — optimise within a single particle-number sector |
| **`Sector(Q, q)`** | State constraint: ``(Q-q)|\psi\rangle = 0`` |
| **`Zero(g)`** | Full bilinear localizing: $\langle b_i^\dagger g\, b_j \rangle = 0$ — for operator identities (``g = 0`` always) |
| **Mott insulator** | Large $U/t$ at half-filling: electrons localise, charge transport freezes |

### See also

- [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
  creation/annihilation operators, CAR, parity constraints
- [Kitaev Chain](@ref kitaev-chain) — pairing terms, Majorana
  operators, topological phases
- [Hubbard Model (Filling-Sector Scan)](@ref hubbard-filling-scan) —
  what the unconstrained SDP does, and how it compares with exact
  diagonalisation across all $(N_\uparrow, N_\downarrow)$ sectors
- [Ground State Energy](@ref ground-state-energy) — Pauli-level spin
  chain examples

### References

- J. Hubbard, "Electron correlations in narrow energy bands,"
  *Proceedings of the Royal Society A* **276**, 238 (1963).
- E. H. Lieb and F. Y. Wu, "Absence of Mott transition in an exact
  solution of the short-range, one-band model in one dimension,"
  *Physical Review Letters* **20**, 1445 (1968).
- F. H. L. Essler, H. Frahm, F. Göhmann, A. Klümper, and V. E.
  Korepin, *The One-Dimensional Hubbard Model*, Cambridge University
  Press (2005).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

