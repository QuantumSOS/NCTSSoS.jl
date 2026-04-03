```@meta
EditURL = "../literate/fermionic_ground_state.jl"
```

# [Fermionic Ground State (XY Model)](@id fermionic-ground-state)

Many quantum many-body problems are most naturally expressed in terms of
**fermionic** creation and annihilation operators.  The **Jordan–Wigner
transformation** maps spin-chain Hamiltonians to fermions, and models that
become quadratic (free) in fermion language are exactly solvable — making
them ideal benchmarks.

In this example we compute certified lower bounds on the ground-state energy
of the **periodic XY spin chain** for $N = 4, 6, 8$ sites using the native
fermionic interface of `NCTSSoS.jl`, and check each bound against the exact
free-fermion answer.

## Background: XY model and Jordan–Wigner

The periodic XY spin chain is
```math
H_{\mathrm{spin}} = \frac{j_c}{4}\sum_{i=1}^{N}
    \bigl(\sigma^x_i \sigma^x_{i+1} + \sigma^y_i \sigma^y_{i+1}\bigr),
    \qquad \sigma_{N+1} \equiv \sigma_1.
```
The Jordan–Wigner transformation replaces spins with spinless fermions
$a_i, a_i^\dagger$ satisfying the canonical anticommutation relations (CAR)
$\{a_i, a_j^\dagger\} = \delta_{ij}$.
Bulk nearest-neighbour terms become simple hoppings, but the **boundary**
term (site $N$ ↔ site $1$) acquires the **fermion parity operator**
```math
P = \prod_{i=1}^{N}(1 - 2\,a_i^\dagger a_i).
```
The full fermionic Hamiltonian is
```math
H = -\frac{j_c}{2}\sum_{i=1}^{N-1}
    \bigl(a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i\bigr)
  + \frac{j_c}{2}\,P\,
    \bigl(a_N^\dagger a_1 + a_1^\dagger a_N\bigr).
```

### Why a parity constraint is needed

$P$ commutes with $H$ and squares to the identity, so the Hilbert space
splits into **even** ($P = +1$) and **odd** ($P = -1$) sectors.
In the even sector the boundary hopping simplifies (the factor $P$ becomes
$+1$), giving the quadratic Hamiltonian
```math
H_{\mathrm{even}} = -\frac{j_c}{2}\sum_{i=1}^{N-1}
    \bigl(a_i^\dagger a_{i+1} + a_{i+1}^\dagger a_i\bigr)
  + \frac{j_c}{2}
    \bigl(a_N^\dagger a_1 + a_1^\dagger a_N\bigr),
```
together with the equality constraint $P = I$.

Without this constraint we would solve a *different* fermionic ring model
and get a looser (more negative) bound.  Open chains need no parity
constraint because there is no boundary twist.

## Exact energy from the free-fermion spectrum

Because $H_{\mathrm{even}}$ is quadratic, its ground-state energy is the
sum of all negative single-particle energies of the hopping matrix.  The
boundary sign flip corresponds to **anti-periodic boundary conditions**.

````julia
using LinearAlgebra

"""
    xy_exact_energy(N; j_c=1.0)

Exact ground-state energy of the periodic N-site XY spin chain via
the free-fermion spectrum in the even-parity sector.
"""
function xy_exact_energy(N::Int; j_c::Real = 1.0)
    # Single-particle hopping matrix (anti-periodic BC from parity twist)
    h = zeros(N, N)
    for i in 1:N-1
        h[i, i+1] = h[i+1, i] = -j_c / 2
    end
    h[N, 1] = h[1, N] = j_c / 2          # anti-periodic boundary
    εs = eigvals(Symmetric(h))
    return sum(ε for ε in εs if ε < -1e-14)
end
````

````
Main.var"##280".xy_exact_energy
````

The exact energies for our three benchmarks are:

````julia
for N in [4, 6, 8]
    println("N = $N :  E₀ = $(xy_exact_energy(N))")
end
````

````
N = 4 :  E₀ = -1.4142135623730951
N = 6 :  E₀ = -1.732050807568877
N = 8 :  E₀ = -2.6131259297527527

````

!!! note "SDP Solver"
    These examples use [Mosek](https://www.mosek.com/) via `MosekTools`.
    Any SDP-capable solver works: replace `Mosek.Optimizer` with
    `COSMO.Optimizer` or `Clarabel.Optimizer` for open-source alternatives.

````julia
using NCTSSoS, MosekTools, JuMP

SOLVER = optimizer_with_attributes(Mosek.Optimizer,
    "MSK_IPAR_LOG" => 0,                               # quiet output
    "MSK_IPAR_NUM_THREADS" => 0)                        # use all threads
````

````
MathOptInterface.OptimizerWithAttributes(Mosek.Optimizer, Pair{MathOptInterface.AbstractOptimizerAttribute, Any}[MathOptInterface.RawOptimizerAttribute("MSK_IPAR_LOG") => 0, MathOptInterface.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS") => 0])
````

-------------------------------------------------------------------
## Case 1 — N = 4 (tight at first iteration)

For $N = 4$ the parity polynomial has degree $2N = 8$.  At **order 4**
the moment matrix includes monomials up to degree 8, which covers the
full fermionic Fock space, and the parity constraint is enforced.
A single call to [`cs_nctssos`](@ref) is sufficient.
-------------------------------------------------------------------

### Step 1 — Create fermionic variables

[`create_fermionic_variables`](@ref) returns a registry and two vectors:
annihilation operators `a[i]` and creation operators `a_dag[i]`.
All CAR relations are encoded automatically.

````julia
N₁ = 4
j_c = 1.0
registry₁, (a₁, a₁_dag) = create_fermionic_variables(1:N₁)
````

````
(VariableRegistry with 8 variables: a⁺₄, a⁺₃, a⁺₂, a⁺₁, a₁, a₂, a₃, a₄, (NCTSSoS.NormalMonomial{NCTSSoS.FermionicAlgebra, Int8}[a₁, a₂, a₃, a₄], NCTSSoS.NormalMonomial{NCTSSoS.FermionicAlgebra, Int8}[a⁺₁, a⁺₂, a⁺₃, a⁺₄]))
````

### Step 2 — Even-parity Hamiltonian

````julia
ham₁ = sum(
    -ComplexF64(j_c / 2) * (a₁_dag[i] * a₁[i+1] + a₁_dag[i+1] * a₁[i])
    for i in 1:N₁-1
) + ComplexF64(j_c / 2) * (a₁_dag[N₁] * a₁[1] + a₁_dag[1] * a₁[N₁])
````

````
0.5 + 0.0im * a⁺₄a₁ + -0.5 - 0.0im * a⁺₄a₃ + -0.5 - 0.0im * a⁺₃a₂ + -0.5 - 0.0im * a⁺₃a₄ + -0.5 - 0.0im * a⁺₂a₁ + -0.5 - 0.0im * a⁺₂a₃ + -0.5 - 0.0im * a⁺₁a₂ + 0.5 + 0.0im * a⁺₁a₄
````

### Step 3 — Parity constraint  $P - I = 0$

````julia
parity₁ = prod(one(ham₁) - 2.0 * a₁_dag[i] * a₁[i] for i in 1:N₁)
````

````
1 + -2.0 + 0.0im * a⁺₄a₄ + -2.0 + 0.0im * a⁺₃a₃ + -2.0 + 0.0im * a⁺₂a₂ + -2.0 + 0.0im * a⁺₁a₁ + 4.0 - 0.0im * a⁺₄a⁺₃a₃a₄ + 4.0 - 0.0im * a⁺₄a⁺₂a₂a₄ + 4.0 - 0.0im * a⁺₄a⁺₁a₁a₄ + 4.0 + 0.0im * a⁺₃a⁺₂a₂a₃ + 4.0 + 0.0im * a⁺₃a⁺₁a₁a₃ + 4.0 + 0.0im * a⁺₂a⁺₁a₁a₂ + -8.0 + 0.0im * a⁺₄a⁺₃a⁺₂a₂a₃a₄ + -8.0 + 0.0im * a⁺₄a⁺₃a⁺₁a₁a₃a₄ + -8.0 + 0.0im * a⁺₄a⁺₂a⁺₁a₁a₂a₄ + -8.0 + 0.0im * a⁺₃a⁺₂a⁺₁a₁a₂a₃ + 16.0 - 0.0im * a⁺₄a⁺₃a⁺₂a⁺₁a₁a₂a₃a₄
````

### Step 4 — Solve

````julia
pop₁ = polyopt(ham₁, registry₁; eq_constraints = [parity₁ - one(ham₁)])
config₁ = SolverConfig(optimizer = SOLVER, order = 4, ts_algo = MMD())
result₁ = cs_nctssos(pop₁, config₁)
````

````
Objective: -1.4142135596021859
Correlative Sparsity (FermionicAlgebra): 

   maximum clique size: 4
   number of cliques: 1
   Clique 1: 
       Variables: [:a₁, :a₂, :a₃, :a₄]
       Bases length: 163
       Constraints: 
           -2.0 + 0.0im * a⁺₄a₄ + -2.0 + 0.0im * a⁺₃a₃ + -2.0 + 0.0im * a⁺₂a₂ + -2.0 + 0.0im * a⁺₁a₁ + 4.0 - 0.0im * a⁺₄a⁺₃a₃a₄ + 4.0 - 0.0im * a⁺₄a⁺₂a₂a₄ + 4.0 - 0.0im * a⁺₄a⁺₁a₁a₄ + 4.0 + 0.0im * a⁺₃a⁺₂a₂a₃ + 4.0 + 0.0im * a⁺₃a⁺₁a₁a₃ + 4.0 + 0.0im * a⁺₂a⁺₁a₁a₂ + -8.0 + 0.0im * a⁺₄a⁺₃a⁺₂a₂a₃a₄ + -8.0 + 0.0im * a⁺₄a⁺₃a⁺₁a₁a₃a₄ + -8.0 + 0.0im * a⁺₄a⁺₂a⁺₁a₁a₂a₄ + -8.0 + 0.0im * a⁺₃a⁺₂a⁺₁a₁a₂a₃ + 16.0 - 0.0im * a⁺₄a⁺₃a⁺₂a⁺₁a₁a₂a₃a₄ :  with 1 basis monomials
   Global Constraints: 
Term Sparsity:
Clique 1:
   Moment Matrix Block Sizes: [4, 5, 6, 5, 5, 4, 5, 6, 7, 5, 5, 6, 4, 5, 4, 5, 6, 7, 8, 9, 5, 7, 5, 7, 7, 7, 5, 7, 7, 10, 10, 5, 4, 4, 4, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 4, 3, 3, 3, 3, 4, 11, 4, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   Moment Matrix:
Number of Activated supp:   68
Number of Bases Activated in each sub-block[4, 5, 6, 5, 5, 4, 5, 6, 7, 5, 5, 6, 4, 5, 4, 5, 6, 7, 8, 9, 5, 7, 5, 7, 7, 7, 5, 7, 7, 10, 10, 5, 4, 4, 4, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 4, 3, 3, 3, 3, 4, 11, 4, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]

   Localizing Matrix:
Number of Activated supp:   15
Number of Bases Activated in each sub-block[1]
Unique Moment Matrix Elements: 68

````

### Step 5 — Verify

````julia
exact₁ = xy_exact_energy(N₁)
println("N = $N₁:  SDP = $(result₁.objective),  exact = $exact₁,  gap = $(abs(result₁.objective - exact₁))")
````

````
N = 4:  SDP = -1.4142135596021859,  exact = -1.4142135623730951,  gap = 2.7709092798033907e-9

````

-------------------------------------------------------------------
## Scaling up: N = 6 and N = 8

For larger systems the first SDP iteration may yield a loose bound
because term sparsity decomposes the moment matrix into blocks that
miss some cross-correlations.  A second call to
[`cs_nctssos_higher`](@ref) refines the sparsity pattern using the
moment values from the first solve, typically tightening the bound to
machine precision.
-------------------------------------------------------------------

### Helper — build and solve for arbitrary N

````julia
"""
    solve_xy(N; j_c=1.0, order=3)

Build the even-parity fermionic XY Hamiltonian for `N` sites,
run one round of `cs_nctssos` and one of `cs_nctssos_higher`,
and return `(first_obj, refined_obj, exact)`.
"""
function solve_xy(N::Int; j_c::Real = 1.0, order::Int = 3)
    registry, (a, a_dag) = create_fermionic_variables(1:N)

    ham = sum(
        -ComplexF64(j_c / 2) * (a_dag[i] * a[i+1] + a_dag[i+1] * a[i])
        for i in 1:N-1
    ) + ComplexF64(j_c / 2) * (a_dag[N] * a[1] + a_dag[1] * a[N])

    parity = prod(one(ham) - 2.0 * a_dag[i] * a[i] for i in 1:N)

    pop = polyopt(ham, registry; eq_constraints = [parity - one(ham)])
    config = SolverConfig(optimizer = SOLVER, order = order, ts_algo = MMD())

    res1 = cs_nctssos(pop, config)
    res2 = cs_nctssos_higher(pop, res1, config)

    return (first = res1.objective, refined = res2.objective,
            exact = xy_exact_energy(N; j_c))
end
````

````
Main.var"##280".solve_xy
````

### N = 6

````julia
r6 = solve_xy(6; order = 3)
println("N = 6:  first = $(r6.first),  refined = $(r6.refined),  exact = $(r6.exact)")
println("  gap after refinement: $(abs(r6.refined - r6.exact))")
````

````
N = 6:  first = -2.037537944197494,  refined = -1.7320508073230187,  exact = -1.732050807568877
  gap after refinement: 2.4585822266942614e-10

````

### N = 8

````julia
r8 = solve_xy(8; order = 3)
println("N = 8:  first = $(r8.first),  refined = $(r8.refined),  exact = $(r8.exact)")
println("  gap after refinement: $(abs(r8.refined - r8.exact))")
````

````
N = 8:  first = -2.7442065306235786,  refined = -2.6131259294966287,  exact = -2.6131259297527527
  gap after refinement: 2.561240108889251e-10

````

-------------------------------------------------------------------
## Results

| ``N`` | Exact ``E_0``             | First iteration | After `cs_nctssos_higher` |
|:-----:|:-------------------------:|:---------------:|:-------------------------:|
| 4     | ``-\sqrt{2} \approx -1.414``  | tight       | —                         |
| 6     | ``-\sqrt{3} \approx -1.732``  | loose       | tight                     |
| 8     | ``\approx -2.613``            | loose       | tight                     |

For $N = 4$ the first SDP iteration already matches the exact energy.
For $N = 6$ and $N = 8$, the sparsity-refined second iteration
([`cs_nctssos_higher`](@ref)) tightens the bound to the exact value.
-------------------------------------------------------------------

## Summary

| Step | What we did |
|------|-------------|
| 1 | Created fermionic variables with [`create_fermionic_variables`](@ref) — CAR built in |
| 2 | Wrote the even-parity-sector quadratic Hamiltonian |
| 3 | Imposed the parity constraint ``P = I`` as an equality |
| 4 | Solved with [`cs_nctssos`](@ref), refined with [`cs_nctssos_higher`](@ref) |
| 5 | Verified against the free-fermion exact energy |

### When to use the fermionic interface

- Models that are naturally quadratic (or low-degree) in fermion operators
- Jordan–Wigner reductions of spin chains
- Problems where parity or particle-number sectors must be fixed

For Pauli-level formulations of ground-state problems, see
[Ground State Energy](@ref ground-state-energy).
For a showcase of fermionic algebra rules (CAR, nilpotency, parity), see
[PBW Algebra Showcase](@ref pbw-algebras-showcase).

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

