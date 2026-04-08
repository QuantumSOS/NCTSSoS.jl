```@meta
EditURL = "../literate/bosonic_ground_state.jl"
```

# [Bosonic Ground State (Bose–Hubbard Model)](@id bosonic-ground-state)

This example computes a **certified lower bound** on the ground-state
energy of the **Bose–Hubbard model** — the fundamental lattice model for
interacting bosons.  We solve the simplest nontrivial case: **2 sites with
2 particles**, where an exact analytic answer exists.

Along the way we meet a key practical issue: unlike fermions (at most one
per state), bosons can pile up without limit.  The bosonic CCR alone do
not choose a particle-number sector for the optimization problem.  Here we
impose $\hat{N} = 2$ with `moment_eq_constraints`, so the SDP solves the
same canonical-sector problem as the exact $3 \times 3$ reference.

If you haven't seen bosonic operators before, the
[PBW Algebra Showcase](@ref pbw-algebras-showcase) introduces
creation/annihilation operators and the CCR.  If you want fermionic
examples first, see the
[Fermionic Ground State (XY Model)](@ref fermionic-ground-state) or the
[Hubbard Model](@ref hubbard-model) for a fermionic system with
canonical particle-number constraints.

---

## The Bose–Hubbard Hamiltonian

Picture $N$ sites on a lattice.  Each site can hold **any number** of
identical bosons (photons in a cavity array, ultracold atoms in an optical
lattice, etc.).  The Hamiltonian has two competing terms:

```math
H = -t \sum_{\langle i,j \rangle}
    \bigl(b_i^\dagger b_j + b_j^\dagger b_i\bigr)
    \;+\; \frac{U}{2} \sum_{i=1}^{N} \hat{n}_i(\hat{n}_i - 1),
```

where $b_i^\dagger, b_i$ are bosonic creation and annihilation operators,
$\hat{n}_i = b_i^\dagger b_i$ counts particles at site $i$, and
$\langle i,j \rangle$ runs over nearest-neighbour pairs.

| Term | Physics | Favours |
|:-----|:--------|:--------|
| $-t\,(b_i^\dagger b_j + \text{h.c.})$ | **Hopping** — a boson tunnels to a neighbour | Delocalisation (superfluid) |
| $\frac{U}{2}\,\hat{n}_i(\hat{n}_i - 1)$ | **On-site repulsion** — energy cost for two or more bosons sharing a site | Localisation (Mott insulator) |

The ratio $U/t$ drives the **superfluid–Mott insulator** quantum phase
transition.  At small $U/t$ bosons spread out coherently; at large $U/t$
they freeze into integer occupation per site.

---

## Why fix particle number?

Fermions obey the Pauli exclusion principle — each site holds 0 or 1
particle, so the Hilbert space is finite-dimensional.  Bosons have **no
such limit**.  In the polynomial optimization formulation, the bosonic
CCR ($[b_i, b_j^\dagger] = \delta_{ij}$) alone do not choose a particle
number sector, even though the Hamiltonian conserves
$\hat{N} = \sum_i \hat{n}_i$.

If we solve the problem with no extra state constraint, the SDP is free to
optimize over **all** particle numbers at once.  That is a different
question from the exact $N = 2$ diagonalization below.  For the
parameters used here ($t = 1$, $U = 2$), the unconstrained problem is
still bounded; it just is not the canonical-sector benchmark we want.

The physically natural fix is to **specify how many particles are in the
system**.  Imposing $\hat{N} = 2$ on the quantum state restricts us to the
3-dimensional sector spanned by $|2,0\rangle$, $|1,1\rangle$, and
$|0,2\rangle$, so the SDP lower bound can be compared directly with the
exact sector energy.

!!! tip "State constraint, not operator identity"
    The constraint $\hat{N}|\psi\rangle = 2|\psi\rangle$ is a property
    of the **target state**, not an algebraic identity.  This is exactly
    the situation for `moment_eq_constraints` (one-sided localizing), not
    `eq_constraints` (bilinear localizing).  See the
    [Hubbard Model](@ref hubbard-model) example for a detailed
    explanation of the difference.

---

## Exact solution for the 2-site, $N = 2$ sector

With 2 particles on 2 sites, the Hilbert space has just **3 basis
states**: $|2,0\rangle$, $|1,1\rangle$, $|0,2\rangle$.  In this basis
the Hamiltonian is a $3 \times 3$ matrix:

```math
H_{3\times3} = \begin{pmatrix}
  U & -\sqrt{2}\,t & 0 \\
  -\sqrt{2}\,t & 0 & -\sqrt{2}\,t \\
  0 & -\sqrt{2}\,t & U
\end{pmatrix}.
```

The off-diagonal factor $\sqrt{2}$ comes from the bosonic enhancement:
$b^\dagger|1\rangle = \sqrt{2}\,|2\rangle$, so the matrix element of the
hopping between $|2,0\rangle$ and $|1,1\rangle$ picks up a $\sqrt{2}$.
The diagonal entries $U$ appear on $|2,0\rangle$ and $|0,2\rangle$
because $\hat{n}(\hat{n}-1)/2 = 1$ when two bosons share one site.

For $t = 1$, $U = 2$, the eigenvalues are $\{1 - \sqrt{5},\; 2,\;
1 + \sqrt{5}\}$.  The ground-state energy is:

```math
E_0 = 1 - \sqrt{5} \approx -1.236.
```

````julia
using LinearAlgebra

t_hop = 1.0
U = 2.0

# Build the 3×3 sector Hamiltonian
hop = -sqrt(2.0) * t_hop
H_exact = Hermitian(Float64[
    U    hop  0.0
    hop  0.0  hop
    0.0  hop  U
])

exact_e0 = eigmin(H_exact)
println("Exact ground-state energy (3×3 diagonalisation): $exact_e0")
println("Analytic value 1 - √5:                           $(1 - sqrt(5))")
````

````
Exact ground-state energy (3×3 diagonalisation): -1.2360679774997898
Analytic value 1 - √5:                           -1.2360679774997898

````

---

## Building the Hamiltonian with NCTSSoS.jl

[`create_bosonic_variables`](@ref) returns creation and annihilation
operators that obey the CCR automatically — every product is
normal-ordered by the library.

!!! note "SDP Solver"
    These examples use [Mosek](https://www.mosek.com/) via `MosekTools`.
    Any SDP-capable solver works: replace `Mosek.Optimizer` with
    `COSMO.Optimizer` or `Clarabel.Optimizer` for open-source alternatives.

````julia
using NCTSSoS, MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG" => 0,
    "MSK_IPAR_NUM_THREADS" => 0)
````

````
MathOptInterface.OptimizerWithAttributes(Mosek.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[MathOptInterface.RawOptimizerAttribute("MSK_IPAR_LOG") => 0, MathOptInterface.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS") => 0])
````

### Step 1 — Create bosonic variables

````julia
N_sites = 2
registry, (b, b_dag) = create_bosonic_variables(1:N_sites)
println("Created $(length(registry)) bosonic variables: b[1..$(N_sites)] and b†[1..$(N_sites)]")
````

````
Created 4 bosonic variables: b[1..2] and b†[1..2]

````

### Step 2 — Number operators and Hamiltonian

The number operator $\hat{n}_i = b_i^\dagger b_i$ counts particles at
site $i$.

````julia
n = [b_dag[k] * b[k] for k in 1:N_sites]

# Hopping: -t (b†₁ b₂ + b†₂ b₁)
ham_hop = -t_hop * (b_dag[1] * b[2] + b_dag[2] * b[1])

# On-site interaction: (U/2) Σᵢ n̂ᵢ(n̂ᵢ - 1)
ham_int = (U / 2) * sum(n[k] * n[k] - n[k] for k in 1:N_sites)

ham = ham_hop + ham_int
````

````
-c⁺₂c₁ + -c⁺₁c₂ + c⁺₂²c₂² + c⁺₁²c₁²
````

### Step 3 — Canonical particle-number constraint

We fix the total particle number to $N_{\mathrm{tot}} = 2$ via a
**moment equality constraint**.  The constraint polynomial
``g = \hat{N} - 2\cdot\mathbf{1}`` satisfies ``g|\psi\rangle = 0`` for any
state in the $N_{\mathrm{tot}} = 2$ sector.

````julia
total_number = 1.0 * sum(n)
canonical_sector = total_number - 2.0 * one(ham)
````

````
-2.0 + c⁺₂c₂ + c⁺₁c₁
````

### Step 4 — Solve the SDP relaxation

An order-2 relaxation is sufficient for this small system.

````julia
pop = polyopt(ham, registry; moment_eq_constraints = [canonical_sector])
config = SolverConfig(optimizer = SOLVER, order = 2)
result = cs_nctssos(pop, config)
````

````
Objective: -1.2360679763077873
Correlative Sparsity (BosonicAlgebra): 

   maximum clique size: 2
   number of cliques: 1
   Clique 1: 
       Variables: [:c₁, :c₂]
       Bases length: 15
       Constraints: 
   Global Constraints: 
Term Sparsity:
Clique 1:
   Moment Matrix Block Sizes: [15]
   Moment Matrix:
Number of Activated supp:   56
Number of Bases Activated in each sub-block[15]

   Localizing Matrix:
Unique Moment Matrix Elements: 70

````

### Step 5 — Compare with the exact answer

````julia
println("SDP objective (order 2):  $(round(result.objective; digits=8))")
println("Exact ground-state energy: $(round(exact_e0; digits=8))")
println("Gap:                       $(round(abs(result.objective - exact_e0); digits=10))")
````

````
SDP objective (order 2):  -1.23606798
Exact ground-state energy: -1.23606798
Gap:                       1.2e-9

````

The SDP lower bound matches the exact sector energy to high precision.

---

## Results

| Quantity | Value |
|:--------|------:|
| Exact $E_0$ (3×3 diagonalisation) | $1 - \sqrt{5} \approx -1.236$ |
| SDP bound (order 2, canonical $N = 2$) | $\approx -1.236$ |
| Moment matrix size | 15 × 15 |

The order-2 relaxation with a single particle-number constraint is already
tight for this 2-site problem.

---

## Summary

### The physics in one paragraph

The Bose–Hubbard model describes interacting bosons hopping on a lattice —
the simplest model for the superfluid-to-Mott-insulator transition seen in
cold-atom experiments.  Unlike fermions, bosonic CCR do not fix the
particle-number sector, so you usually decide which canonical sector you
want to study.  Here we impose total particle number $\hat{N} = 2$ via
`moment_eq_constraints`, making the SDP solve the same $3 \times 3$
canonical-sector problem as the exact diagonalization.  For this minimal
2-site, 2-particle case, the relaxation reproduces $E_0 = 1 - \sqrt{5}$.

### The code recipe

| Step | What we did | API |
|:-----|:------------|:----|
| 1 | Create bosonic operators (CCR built in) | [`create_bosonic_variables`](@ref) |
| 2 | Build the Bose–Hubbard Hamiltonian (hopping + on-site repulsion) | Standard Julia arithmetic |
| 3 | Impose canonical particle number $\hat{N} = 2$ | `moment_eq_constraints` in [`polyopt`](@ref) |
| 4 | Solve the order-2 SDP relaxation | [`cs_nctssos`](@ref) |
| 5 | Verify against exact 3×3 diagonalisation | `eigmin` from `LinearAlgebra` |

### Key concepts

| Concept | One-liner |
|:--------|:----------|
| **Bose–Hubbard model** | Hopping $t$ vs. on-site repulsion $U$ — simplest model of interacting lattice bosons |
| **Bosonic CCR** | $[b_i, b_j^\dagger] = \delta_{ij}$ — commutation (not anticommutation), no exclusion principle |
| **Fixed particle number** | Without a sector constraint the SDP optimizes over all particle numbers; imposing $\hat{N} = 2$ makes it match the exact canonical reference |
| **`moment_eq_constraints`** | One-sided localizing: $\langle b_i^\dagger g \rangle = 0$ — for state constraints (``g|\psi\rangle = 0``) |
| **Canonical sector** | Fix total $\hat{N} = \sum_i \hat{n}_i$ to restrict to a finite-dimensional Hilbert-space sector |
| **Superfluid–Mott transition** | Bosonic quantum phase transition driven by $U/t$ ratio |

### See also

- [PBW Algebra Showcase](@ref pbw-algebras-showcase) — CCR and CAR normal ordering
- [Hubbard Model](@ref hubbard-model) — fermionic Hubbard model with canonical constraints and the `eq_constraints` vs. `moment_eq_constraints` distinction
- [Fermionic Ground State (XY Model)](@ref fermionic-ground-state) — fermionic creation/annihilation operators from scratch

### References

- M. P. A. Fisher, P. B. Weichman, G. Grinstein, and D. S. Fisher,
  "Boson localization and the superfluid-insulator transition,"
  *Physical Review B* **40**, 546 (1989).
- D. Jaksch, C. Bruder, J. I. Cirac, C. W. Gardiner, and P. Zoller,
  "Cold bosonic atoms in optical lattices,"
  *Physical Review Letters* **81**, 3108 (1998).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

