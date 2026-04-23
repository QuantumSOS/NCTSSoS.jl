```@meta
EditURL = "../literate/h2_fermionic_symmetry.jl"
```

# [H₂/STO-3G with Fermionic Symmetry Adaptation](@id h2-fermionic-symmetry)

This page walks through the **full hydrogen pipeline** — from
molecular-orbital integrals to a certified electronic ground-state
energy — using the **fermionic + symmetry** path of `NCTSSoS.jl`.

H₂ is the smallest non-trivial molecule.  In the minimal STO-3G basis
it has just **two spatial orbitals** ($\sigma_g$, $\sigma_u$) and
**four spin orbitals** ($\sigma_g\!\uparrow$, $\sigma_g\!\downarrow$,
$\sigma_u\!\uparrow$, $\sigma_u\!\downarrow$).  Two electrons live on
this stage.  That is small enough to:

1. write down the *exact* electronic Hamiltonian as a polynomial in
   fermionic creation/annihilation operators;
2. solve it exactly by full diagonalisation, restricted to the physical
   $N = 2$ sector, as a verification oracle;
3. solve the same problem as an SDP relaxation, **once dense** and
   **once symmetry-adapted**, and recover the same energy with much
   smaller PSD blocks.

The point of using H₂ here is **not** chemistry.  The point is that
this is the smallest setting where every piece of the fermionic
symmetry machinery — sector splitting (parity / particle number /
$S_z$ / Abelian orbital irrep) and SU(2) spin adaptation — does
something visible and easy to interpret.

Prerequisites:

- Comfort with creation/annihilation operators and the fermionic CAR;
  see [Fermionic Ground State (XY Model)](@ref fermionic-ground-state)
  for a from-scratch build.
- The [Hubbard Model](@ref hubbard-model) example, which introduces
  `moment_eq_constraints` for fixing particle number.
- The [CHSH with Symmetry Reduction](@ref chsh-symmetry) example and
  the manual page [Symmetry-Adapted Basis](@ref symmetry-adapted-basis),
  which describe the design of NCTSSoS's symmetry path.

Throughout this page, the **reference energy is electronic only** —
we never add the nuclear-repulsion constant.  Tacking that constant on
at the end gives the usual textbook total energy, but it is a
scalar that does not affect any of the operator-algebra discussion.

---

## Why H₂ deserves a fermionic + symmetry treatment

Even at this size, H₂ already has *all* of the structure that makes
molecular electronic structure hard:

| Structure | What it gives us |
|:----------|:------------------|
| **Fermionic CAR** | electrons obey Pauli exclusion; the algebra is non-commutative |
| **Particle-number conservation** | only the $N = 2$ sector is physical |
| **Total $S_z$ conservation** | the Hamiltonian is spin-rotation invariant |
| **Total $S^2$ conservation (SU(2))** | physical eigenstates are singlets ($S = 0$) or triplets ($S = 1$) |
| **Inversion symmetry $i$ (D$_{2h}$ ↓)** | orbitals split into gerade ($\sigma_g$) and ungerade ($\sigma_u$); operators carry an Abelian $\{ \mathrm{Ag}, \mathrm{B1u} \}$ label |

A *dense* moment relaxation ignores all of that: it builds one large PSD
block over every monomial up to some order.  The fermionic **sector
split** discards monomials whose quantum numbers cannot survive in the
physical sector, and the **spin adaptation** further block-diagonalises
the surviving piece into total-$S^2$ blocks.

Concretely, the order-2 dense moment matrix on these four spin orbitals
is a single PSD constraint of side close to $50$.  The symmetry-adapted
relaxation in this page returns the *same* electronic ground-state
energy, but the largest reduced PSD block has side $3$.

---

## Exact diagonalisation oracle (reviewed integrals)

We need an exact answer to compare against.  Build the four-mode
Hamiltonian via the **Jordan–Wigner** mapping (fermion operators →
tensor products of Pauli matrices) and take the smallest eigenvalue
*restricted to* the physically relevant $N = 2$ sector.

````julia
using LinearAlgebra

# Single-mode Jordan–Wigner pieces (occupation-basis convention).
const PAULI_Z   = ComplexF64[1 0; 0 -1]
const PAULI_I   = Matrix{ComplexF64}(I, 2, 2)
const JW_RAISE  = ComplexF64[0 1; 0 0]   ## σ⁺ = |0⟩⟨1|, the annihilator on |1⟩
const JW_LOWER  = ComplexF64[0 0; 1 0]   ## σ⁻ = |1⟩⟨0|, the creator on |0⟩

"""
    jw_op(mode, kind, nmodes)

Jordan–Wigner image of a fermion operator on `nmodes` sites.  `kind` is
`:annihilate` or `:create`.
"""
function jw_op(mode::Int, kind::Symbol, nmodes::Int)
    mats = [site < mode ? PAULI_Z :
            site == mode ? (kind === :annihilate ? JW_RAISE : JW_LOWER) :
            PAULI_I for site in 1:nmodes]
    return reduce(kron, mats)
end
````

````
Main.var"##288".jw_op
````

H₂'s reviewed offline molecular-orbital integrals — one-body matrix
elements `h[p,q]` and antisymmetric two-body tensor `v[p,q,r,s]`
(chemists' notation, electronic-only).  These exact constants come
from `test/problems/fermionic/fermionic_symmetry.jl`.

````julia
function h2_sto3g_integrals()
    h = [
        -1.252477303982   0.0;
         0.0            -0.475934275367;
    ]
    v = zeros(2, 2, 2, 2)
    v[1, 1, 1, 1] = 0.674493166046
    v[2, 2, 2, 2] = 0.69739794957
    v[1, 2, 1, 2] = 0.663472101055
    v[2, 1, 2, 1] = 0.663472101055
    v[1, 2, 2, 1] = 0.181287518306
    v[2, 1, 1, 2] = 0.181287518306
    v[1, 1, 2, 2] = 0.181287518306
    v[2, 2, 1, 1] = 0.181287518306
    return h, v
end
````

````
h2_sto3g_integrals (generic function with 1 method)
````

Mode layout for the 4-mode Jordan–Wigner construction:
modes `1, 2` are $\sigma_g\!\uparrow, \sigma_u\!\uparrow$,
modes `3, 4` are $\sigma_g\!\downarrow, \sigma_u\!\downarrow$.

````julia
function h2_exact_n2_energy()
    h, v = h2_sto3g_integrals()
    nmodes = 4
    a   = [jw_op(i, :annihilate, nmodes) for i in 1:nmodes]
    a_d = [jw_op(i, :create,     nmodes) for i in 1:nmodes]

    up = (1, 2)   ## (σ_g↑, σ_u↑)
    dn = (3, 4)   ## (σ_g↓, σ_u↓)
    dim = 2^nmodes
    H = zeros(ComplexF64, dim, dim)

    for p in 1:2, q in 1:2
        H .+= h[p, q] * (a_d[up[p]] * a[up[q]] + a_d[dn[p]] * a[dn[q]])
    end
    for p in 1:2, q in 1:2, r in 1:2, s in 1:2
        coef = 0.5 * v[p, q, r, s]
        coef == 0 && continue
        for (cd_σ, c_σ) in ((a_d[up[p]], a[up[r]]), (a_d[dn[p]], a[dn[r]])),
            (cd_τ, c_τ) in ((a_d[up[q]], a[up[s]]), (a_d[dn[q]], a[dn[s]]))
            # chemists' notation: ½ Σ v_{pqrs} c†_{pσ} c†_{qτ} c_{sτ} c_{rσ}
            # with the up/down assignment matching (cd_σ, c_σ) / (cd_τ, c_τ).
            H .+= coef * cd_σ * cd_τ * c_τ * c_σ
        end
    end

    # Restrict to the N = 2 occupation sector.
    keep = [idx for idx in 1:dim if count_ones(UInt(idx - 1)) == 2]
    return real(eigmin(Hermitian(H[keep, keep])))
end

E_exact = h2_exact_n2_energy()
println("Exact electronic N = 2 ground state:  E₀ = $E_exact")
````

````
Exact electronic N = 2 ground state:  E₀ = -1.8510462258674931

````

This is the **electronic-only** ground-state energy of H₂/STO-3G with
two electrons.  Add the nuclear-repulsion constant separately if you
want the textbook total energy; the polynomial we feed the SDP does
not contain it, and neither does the value above.

---

## Building H₂ as a fermionic polynomial

We now build the *same* Hamiltonian via the polynomial interface.
[`create_fermionic_variables`](@ref) gives us creation and annihilation
operators that already obey every CAR — including cross-species
(different spin) anticommutation.  Two species, two spatial orbitals:

````julia
using NCTSSoS

registry, ((c_up, c_up_dag), (c_dn, c_dn_dag)) = create_fermionic_variables([
    ("c_up", 1:2),   ## (σ_g↑, σ_u↑)
    ("c_dn", 1:2),   ## (σ_g↓, σ_u↓)
]);
length(registry)
````

````
8
````

The four spin-orbital modes, in the order NCTSSoS assigned them:

````julia
up_modes = [Int(op.word[1]) for op in c_up]
dn_modes = [Int(op.word[1]) for op in c_dn]
(up_modes, dn_modes)
````

````
([1, 2], [3, 4])
````

!!! note "What are `op.word[1]` indices?"
    Every fermionic operator returned by
    [`create_fermionic_variables`](@ref) is a one-letter
    `NormalMonomial`; `op.word[1]` is the underlying registry index of
    that letter.  We use these indices below to tell the symmetry layer
    which physical mode is spin-up vs. spin-down, and which orbital is
    gerade vs. ungerade.

The electronic Hamiltonian is built term by term from the integrals.

````julia
function h2_hamiltonian(c_up, c_up_dag, c_dn, c_dn_dag)
    h, v = h2_sto3g_integrals()
    H = zero(1.0 * one(c_up[1]))

    for p in 1:2, q in 1:2
        H += h[p, q] * (c_up_dag[p] * c_up[q] + c_dn_dag[p] * c_dn[q])
    end

    for p in 1:2, q in 1:2, r in 1:2, s in 1:2
        coef = 0.5 * v[p, q, r, s]
        coef == 0 && continue
        for (cσ, cσ_dag) in ((c_up, c_up_dag), (c_dn, c_dn_dag)),
            (cτ, cτ_dag) in ((c_up, c_up_dag), (c_dn, c_dn_dag))
            H += coef * (cσ_dag[p] * cτ_dag[q] * cτ[s] * cσ[r])
        end
    end
    return H
end

ham = h2_hamiltonian(c_up, c_up_dag, c_dn, c_dn_dag);
````

---

## The fermionic mode layout

Sector splitting and spin adaptation both need to know **which physical
mode is what** — its orbital, its spin label, and which Abelian
irrep its orbital carries.  The [`FermionicModeLayout`](@ref) is the
struct that records this metadata.

For H₂/STO-3G:

| Mode | Spin (2$S_z$) | Orbital | Irrep |
|:-----|:-------------:|:-------:|:------|
| `up_modes[1]` (σ_g↑) | `+1` | 1 | `:Ag`  |
| `up_modes[2]` (σ_u↑) | `+1` | 2 | `:B1u` |
| `dn_modes[1]` (σ_g↓) | `-1` | 1 | `:Ag`  |
| `dn_modes[2]` (σ_u↓) | `-1` | 2 | `:B1u` |

The irrep labels live on **orbitals**, not on modes — both spins of
the same spatial orbital share the same irrep.  This matches the
physics: $\sigma_g$ vs. $\sigma_u$ is a property of the spatial wave
function under inversion ($i$), not of the spin.

The composition rule is the discrete-$\mathbb{Z}_2$ table:
any two equal labels multiply to `:Ag` (totally symmetric, gerade
product); any two unequal labels multiply to `:B1u` (the lone
odd irrep).  We give that table to [`AbelianIrrepTable`](@ref).

````julia
h2_irrep_multiply(a, b) = a === b ? :Ag : :B1u

layout = FermionicModeLayout(
    Dict(vcat(
        [up_modes[i] => i for i in eachindex(up_modes)],
        [dn_modes[i] => i for i in eachindex(dn_modes)],
    ));
    spin2_of = Dict(vcat(
        [up_modes[i] => 1  for i in eachindex(up_modes)],
        [dn_modes[i] => -1 for i in eachindex(dn_modes)],
    )),
    irrep_of    = Dict(1 => :Ag, 2 => :B1u),
    irrep_table = AbelianIrrepTable(:Ag; multiply = h2_irrep_multiply, dual = identity),
)
````

````
NCTSSoS.FermionicModeLayout(Dict(4 => 2, 2 => 2, 3 => 1, 1 => 1), Dict(4 => -1, 2 => 1, 3 => -1, 1 => 1), Dict{Int64, Any}(2 => :B1u, 1 => :Ag), NCTSSoS.AbelianIrrepTable{Symbol, typeof(Main.var"##288".h2_irrep_multiply), typeof(identity)}(:Ag, Main.var"##288".h2_irrep_multiply, identity))
````

Two things to notice in the layout:

- `orbital_of[mode] = orbital_index` (the first positional argument)
  is what links the spin orbital back to its parent spatial orbital.
- `spin2_of[mode] = ±1` records *twice* the spin (in units of $1/2$),
  so up is `+1` and down is `-1`.  The factor of two avoids fractional
  keys.

---

## Sector split + SU(2) spin adaptation

A [`FermionicSectorSpec`](@ref) declares which conserved quantum
numbers should be used to **partition** monomials into independent
sectors.  We turn on all four for H₂:

````julia
sector = FermionicSectorSpec(
    mode_layout    = layout,
    split_parity   = true,    ## fermionic (-1)^N parity, always exact
    split_number   = true,    ## total electron count N
    split_spin     = true,    ## total 2 S_z
    split_irrep    = true,    ## Abelian orbital-irrep label (here Ag/B1u)
)
````

````
NCTSSoS.FermionicSectorSpec(NCTSSoS.FermionicModeLayout(Dict(4 => 2, 2 => 2, 3 => 1, 1 => 1), Dict(4 => -1, 2 => 1, 3 => -1, 1 => 1), Dict{Int64, Any}(2 => :B1u, 1 => :Ag), NCTSSoS.AbelianIrrepTable{Symbol, typeof(Main.var"##288".h2_irrep_multiply), typeof(identity)}(:Ag, Main.var"##288".h2_irrep_multiply, identity)), true, true, true, true)
````

A [`FermionicSpinAdaptationSpec`](@ref) layers on top of that: within
each sector that has well-defined $S_z$, it diagonalises the SU(2)
Casimir $\hat S^2$ on the basis to identify singlets ($S = 0$,
`total_spin2 = 0`), triplets ($S = 1$, `total_spin2 = 2`), etc.

````julia
spin = FermionicSpinAdaptationSpec(mode_layout = layout)
````

````
NCTSSoS.FermionicSpinAdaptationSpec(NCTSSoS.FermionicModeLayout(Dict(4 => 2, 2 => 2, 3 => 1, 1 => 1), Dict(4 => -1, 2 => 1, 3 => -1, 1 => 1), Dict{Int64, Any}(2 => :B1u, 1 => :Ag), NCTSSoS.AbelianIrrepTable{Symbol, typeof(Main.var"##288".h2_irrep_multiply), typeof(identity)}(:Ag, Main.var"##288".h2_irrep_multiply, identity)), 1.0e-8, false)
````

The two are wrapped in a single [`SymmetrySpec`](@ref).  No raw
[`SignedPermutation`](@ref) generators are needed for this example —
all the structure comes from the conserved quantum numbers and from
the Casimir, not from a finite group action on operators.

````julia
symspec = SymmetrySpec(sector = sector, spin_adaptation = spin)
````

````
NCTSSoS.SymmetrySpec(NCTSSoS.SignedPermutation[], NCTSSoS.FermionicModePermutation[], NCTSSoS.FermionicSectorSpec(NCTSSoS.FermionicModeLayout(Dict(4 => 2, 2 => 2, 3 => 1, 1 => 1), Dict(4 => -1, 2 => 1, 3 => -1, 1 => 1), Dict{Int64, Any}(2 => :B1u, 1 => :Ag), NCTSSoS.AbelianIrrepTable{Symbol, typeof(Main.var"##288".h2_irrep_multiply), typeof(identity)}(:Ag, Main.var"##288".h2_irrep_multiply, identity)), true, true, true, true), NCTSSoS.FermionicSpinAdaptationSpec(NCTSSoS.FermionicModeLayout(Dict(4 => 2, 2 => 2, 3 => 1, 1 => 1), Dict(4 => -1, 2 => 1, 3 => -1, 1 => 1), Dict{Int64, Any}(2 => :B1u, 1 => :Ag), NCTSSoS.AbelianIrrepTable{Symbol, typeof(Main.var"##288".h2_irrep_multiply), typeof(identity)}(:Ag, Main.var"##288".h2_irrep_multiply, identity)), 1.0e-8, false), true)
````

---

## Fixing $N = 2$ as a state constraint

H₂ is a two-electron problem.  Any solver should look only at the
$N = 2$ sector of Fock space.  We do that with `moment_eq_constraints`,
which imposes $g\,|\psi\rangle = 0$ on the **target state** rather than
as an operator identity (cf. the warning in [Hubbard
Model](@ref hubbard-model)).

````julia
n_total = sum(c_up_dag[i] * c_up[i] + c_dn_dag[i] * c_dn[i] for i in 1:2)

pop = polyopt(ham, registry; moment_eq_constraints = [n_total - 2.0 * one(ham)]);
````

---

## Solver setup

!!! note "SDP Solver"
    These examples use [Mosek](https://www.mosek.com/) via `MosekTools`.
    Any SDP-capable solver works: replace `Mosek.Optimizer` with
    `COSMO.Optimizer` or `Clarabel.Optimizer` for open-source
    alternatives.

````julia
using MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG"         => 0,
    "MSK_IPAR_NUM_THREADS" => 0)
````

---

## Run 1 — Dense order-2 baseline

No symmetry, no sparsity: just the full order-2 dense moment relaxation.

````julia
dense_config = SolverConfig(
    optimizer = SOLVER,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
)

dense_result = cs_nctssos(pop, dense_config);

println("Dense  : objective = $(dense_result.objective)")
println("Exact  : objective = $E_exact")
println("Gap    : ", abs(dense_result.objective - E_exact))
````

````
Dense  : objective = -1.8510462240012158
Exact  : objective = -1.8510462258674931
Gap    : 1.8662773548783207e-9

````

A single PSD block; its size is the dense order-2 baseline we will
shrink:

````julia
dense_block_sizes = only(dense_result.moment_matrix_sizes)
````

````
1-element Vector{Int64}:
 37
````

(Here `only(...)` extracts the per-clique vector for the single
correlative-sparsity clique that the relaxation has by default.)

````julia
dense_result.symmetry === nothing
````

````
true
````

---

## Run 2 — Symmetry-adapted relaxation

The only change is the `symmetry` keyword.  Same Hamiltonian, same
constraint, same order, same solver.

````julia
sym_config = SolverConfig(
    optimizer = SOLVER,
    order     = 2,
    cs_algo   = NoElimination(),
    ts_algo   = NoElimination(),
    symmetry  = symspec,
)

sym_result = cs_nctssos(pop, sym_config);

println("Sym    : objective = $(sym_result.objective)")
println("Dense  : objective = $(dense_result.objective)")
println("|Δ|    : ", abs(sym_result.objective - dense_result.objective))
````

````
MOSEK warning 705 (MSK_RES_WRN_ZEROS_IN_SPARSE_ROW): 8 (nearly) zero elements are specified in sparse row ''(46) of matrix 'A'.
MOSEK warning 705 (MSK_RES_WRN_ZEROS_IN_SPARSE_ROW): 8 (nearly) zero elements are specified in sparse row ''(75) of matrix 'A'.
Sym    : objective = -1.8510462079263064
Dense  : objective = -1.8510462240012158
|Δ|    : 1.6074909403585025e-8

````

Same answer — to solver tolerance.  The reduction is **exact**, not an
additional relaxation: every isotypic component is a separate PSD
block, and the original PSD constraint is the direct sum of those
blocks.

---

## Inspect the [`SymmetryReport`](@ref)

`result.symmetry` is no longer `nothing`; it is a small struct
summarising what the symmetry path did.

````julia
report = sym_result.symmetry
````

````
SymmetryReport(group_order=1, invariant_moment_count=14, psd_block_sizes=[1, 2, 1, 1, 1, 2, 2, 2, 2, 3, 2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], basis_half_size=37, basis_full_size=163, block_provenance=[:sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin, :sector_spin])
````

The pre-symmetry half-basis (i.e. the basis the dense path would have
used for each side of its single moment block):

````julia
(report.basis_half_size, report.basis_full_size)
````

````
(37, 163)
````

The single dense PSD block of side $\approx 50$ has been replaced by
many small blocks, the largest of side $3$:

````julia
report.psd_block_sizes
````

````
26-element Vector{Int64}:
 1
 2
 1
 1
 1
 2
 2
 2
 2
 3
 2
 1
 2
 2
 2
 1
 1
 1
 1
 1
 1
 1
 1
 1
 1
 1
````

Comparison with the pre-symmetry single-block layout:

````julia
(maximum(report.psd_block_sizes), maximum(dense_block_sizes))
````

````
(3, 37)
````

Each block carries a **label** describing where it came from.  Most
come from the spin-adaptation step and are
[`FermionicSpinBlockLabel`](@ref); the rest are
[`FermionicSectorLabel`](@ref) for sectors with no spin structure
(e.g. the trivial / parity-only blocks).

````julia
unique(typeof.(report.block_labels))
````

````
1-element Vector{DataType}:
 NCTSSoS.FermionicSpinBlockLabel
````

A quick distribution of block provenance — which step in the pipeline
produced each block:

````julia
unique(report.block_provenance)
````

````
1-element Vector{Symbol}:
 :sector_spin
````

---

## What the labels mean for H₂

The two structural facts you should be able to read off the report:

**(a) Orbital irreps `:Ag` vs. `:B1u`.**
Densities like $a^\dagger_{\sigma_g \uparrow}\,a_{\sigma_g \uparrow}$
or $a^\dagger_{\sigma_u \downarrow}\,a_{\sigma_u \downarrow}$ are
*irrep-neutral* — they multiply to `:Ag`.  Transitions like
$a^\dagger_{\sigma_g \uparrow}\,a_{\sigma_u \uparrow}$ flip gerade $↔$
ungerade and carry the lone odd irrep `:B1u`.  The set of irrep labels
that show up among the spin-adapted blocks is exactly $\{ \mathrm{Ag},
\mathrm{B1u} \}$:

````julia
irreps_present = Set(label.sector.irrep
    for label in report.block_labels
    if label isa FermionicSpinBlockLabel)

irreps_present
````

````
Set{Symbol} with 2 elements:
  :Ag
  :B1u
````

**(b) Total spin $S^2$ values.**
Within each sector that has a well-defined $S_z$, the Casimir-based
spin adaptation splits the basis into singlet ($S = 0$,
`total_spin2 = 0`) and triplet ($S = 1$, `total_spin2 = 2`) blocks.
Both show up:

````julia
spins_present = Set(label.total_spin2
    for label in report.block_labels
    if label isa FermionicSpinBlockLabel)

spins_present
````

````
Set{Int64} with 3 elements:
  0
  2
  1
````

!!! note "`total_spin2` is `2 S²`, not `S²`"
    The convention is the same as for `spin2_of` in the layout: we
    store **twice** the conserved spin quantum number to keep
    everything in integers.  So `total_spin2 = 0` is $S = 0$
    (singlet) and `total_spin2 = 2` is $S = 1$ (triplet).

That the H₂ ground state is a singlet is well known.  The relaxation
does not need that knowledge — both singlet and triplet blocks live
in the SDP, and the $S^2$ adaptation is what keeps them from being
tangled together inside one block.

---

## Summary

### Pipeline at a glance

| Step | What we did | API |
|:-----|:------------|:----|
| 1 | Build the exact electronic Hamiltonian as a polynomial in $a$, $a^\dagger$ from reviewed STO-3G integrals | [`create_fermionic_variables`](@ref) + standard Julia arithmetic |
| 2 | Tag each spin-orbital mode with its spin and orbital irrep (Ag/B1u) | [`FermionicModeLayout`](@ref), [`AbelianIrrepTable`](@ref) |
| 3 | Declare which conserved quantum numbers should split the moment basis into sectors | [`FermionicSectorSpec`](@ref) |
| 4 | Layer Casimir-based SU(2) spin adaptation on top | [`FermionicSpinAdaptationSpec`](@ref) |
| 5 | Wrap (3) + (4) in a single symmetry spec | [`SymmetrySpec`](@ref) |
| 6 | Fix the physical sector $N = 2$ as a state constraint | `moment_eq_constraints` of [`polyopt`](@ref) |
| 7 | Solve dense and symmetry-adapted relaxations | [`cs_nctssos`](@ref) with and without `symmetry = symspec` |
| 8 | Read the symmetry report; check Ag/B1u and singlet/triplet labels | `result.symmetry` ([`SymmetryReport`](@ref)) |

### Numerical headline

| Quantity | Dense order 2 | Symmetry-adapted order 2 |
|:---------|:--------------|:--------------------------|
| objective (electronic only)   | matches $-1.85104622\!\dots$ | matches the same value |
| largest PSD block side        | dense (single big block)     | $3$ |
| block labels available?       | no                            | yes — Ag/B1u + singlet/triplet |
| output type                   | [`PolyOptResult`](@ref NCTSSoS.PolyOptResult) | same, plus `result.symmetry` populated |

### Caveats

- The energy reported here is **electronic only**.  If you want the
  physical total energy, add the constant nuclear-repulsion term
  yourself; it does not interact with the operator algebra and so
  it would only shift the SDP objective by a constant.
- This page exercises the *currently supported* fermionic-symmetry
  path: sector splitting and Casimir-based spin adaptation, with no
  raw [`SignedPermutation`](@ref) generators in the same spec.  See
  [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) for the
  precise scope and fail-fast guards.
- The integrals used are reviewed offline molecular-orbital
  constants for STO-3G H₂; they live alongside the regression test
  `test/problems/fermionic/fermionic_symmetry.jl` and the reviewed
  expectation fixture `test/data/expectations/fermionic_symmetry.toml`
  (case `h2_sto3g_n2_spin_adapted`).  Swapping in your own integrals
  for a different molecule means changing only `h2_sto3g_integrals`
  and the orbital → irrep assignment.

### See also

- [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) —
  from-scratch build of the fermionic interface.
- [Hubbard Model](@ref hubbard-model) — particle-number constraints
  via `moment_eq_constraints`.
- [CHSH with Symmetry Reduction](@ref chsh-symmetry) — the original
  walk-through of the symmetry path, on a non-fermionic problem.
- [Symmetry-Adapted Basis](@ref symmetry-adapted-basis) — design and
  supported scope of the symmetry path.
- The companion regression suite
  `test/problems/fermionic/fermionic_symmetry.jl` (case
  `h2_sto3g_n2_spin_adapted`) keeps the numbers in this page honest.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

