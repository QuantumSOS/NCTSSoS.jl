```@meta
EditURL = "../literate/monoid_algebras_showcase.jl"
```

# [Monoid Algebra Showcase (NonCommutative / Projector / Unipotent)](@id monoid-algebras-showcase)

This page showcases the three site rules supported by `NCTSSoS`.
These algebras provide the building blocks for modeling noncommutative polynomial
optimization problems that arise in quantum information and many-body physics.

## Capabilities Overview

Each monoid algebra type provides **automatic simplification** of operator products
according to its defining relations:

| Algebra Type | Simplification Rule | Typical Use Case |
|:-------------|:--------------------|:-----------------|
| `Free()` | None (free words) | General NC polynomial optimization |
| `Projector()` | `P² = P` (idempotent) | Projective measurements, NPA hierarchy |
| `Unipotent()` | `U² = 𝟙` (involution) | Pauli observables, Bell inequalities |

All three algebras share a common feature: **site-based commutation**.
Operators acting on different physical sites automatically commute during
canonicalization, while the order of operators on the same site is preserved.

## Key Public APIs

Use `site_variables(rule, nsites, nletters)` with `Free()`, `Projector()`
or `Unipotent()`. Each site has `nletters` self-adjoint generators.

The result is a `(registry, variable_groups)` tuple. The registry stores
names and sites; each variable group contains one site's polynomials.

## References and Applications

- **`Free()`** (free words): The standard modeling layer for
  noncommutative polynomial optimization and SDP hierarchies.
  See Burgdorf--Klep--Povh [burgdorfOptimizationPolynomialsNonCommuting2016](@cite),
  Wang--Magron [wangExploitingTermSparsity2021](@cite), and the comprehensive
  treatment in [magronSparsePolynomialOptimization2023](@cite).

- **`Projector()`** (`P² = P`): Projection operators and projective measurements.
  In device-independent quantum information, projective measurements appear naturally
  in the Navascués--Pironio--Acín (NPA) SDP hierarchy for bounding quantum correlations
  [navascues2007Bounding](@cite), [navascues2008Convergent](@cite).
  In many-body physics, projectors enforce local constraints and define subspaces;
  see Sachdev [sachdev2011QuantumPhaseTransitions](@cite).

- **`Unipotent()`** (`U² = 𝟙`): Involutory (±1-valued) observables model
  self-adjoint unitaries like Pauli measurements. These are fundamental in quantum
  information [nielsenChuang2010](@cite) and appear in Bell inequality optimization
  [klep2024State](@cite).

## Setup

````julia
using NCTSSoS
````

Helper: extract the word from a *single-term* polynomial.

````julia
mono1(p) = last(only(NCTSSoS.terms(p)))
````

````
mono1 (generic function with 1 method)
````

Pretty-print a detached word using the registry's names.

````julia
repr_with_registry(reg, m) = sprint(show, m; context=:registry => reg)
````

````
repr_with_registry (generic function with 1 method)
````

## 1) Free noncommutative words

The `Free()` rule represents free words with no algebraic relations
other than site-based commutation. This is the most general monoid algebra type.

### Creating Variables

Request two sites with two generators each. Different sites commute.
The local aliases `x` and `y` do not rename the registry's operators.

````julia
reg_nc, (x, y) = site_variables(Free(), 2, 2)
````

````
(Registry(Free, 2 sites, 4 variables), [[s1_1, s1_2], [s2_1, s2_2]])
````

This creates:
- Variables `x[1], x[2]` on **site 1**
- Variables `y[1], y[2]` on **site 2**

A single generator is a coefficient-one polynomial. Inspect its term:

````julia
collect(NCTSSoS.terms(x[1]))
````

````
1-element Vector{Tuple{ComplexF64, NCTSSoS.SiteWord{NCTSSoS.Free, UInt128, 4, 2}}}:
 (1.0 + 0.0im, SiteWord(operator(site=1, label=1)))
````

### Site-Based Commutation

Operators on different sites commute automatically. The canonical form groups
operators by site, preserving intra-site order:

````julia
lhs = x[2] * y[1] * x[1] * y[2]
rhs = x[2] * x[1] * y[1] * y[2]  # same canonical form

@assert lhs == rhs  # true: sites are sorted

repr_with_registry(reg_nc, mono1(lhs))  # "s1_2 s1_1 s2_1 s2_2"
````

````
"s1_2 s1_1 s2_1 s2_2"
````

### Intra-Site Non-Commutativity

Within a single site, operators do NOT commute (no relations are imposed):

````julia
@assert x[1] * x[2] != x[2] * x[1]
````

Show both representations to see they differ:

````julia
(repr_with_registry(reg_nc, mono1(x[1] * x[2])),
 repr_with_registry(reg_nc, mono1(x[2] * x[1])))
````

````
("s1_1 s1_2", "s1_2 s1_1")
````

### Inspecting the Canonical Word

Inspect letters through public calls, without decoding packed storage:

````julia
word = mono1(lhs);
````

Each pair below gives the registry name and physical site.

````julia
letter_sites = [(NCTSSoS.operator_name(reg_nc, letter), only(NCTSSoS.sites(letter)))
                for letter in NCTSSoS.word_letters(word)]
````

````
4-element Vector{Tuple{Symbol, Int64}}:
 (:s1_2, 1)
 (:s1_1, 1)
 (:s2_1, 2)
 (:s2_2, 2)
````

## 2) Projector: Idempotent Operators (P² = P)

The `Projector()` rule models projection operators (idempotents). Any repeated
application of the same projector collapses to a single instance.

**Quantum Information Applications:**
- Projective measurements (PVM elements) in the NPA hierarchy
  [navascues2007Bounding](@cite), [navascues2008Convergent](@cite)
- Local projectors like |0⟩⟨0| and |1⟩⟨1| for qubit measurements

**Many-Body Physics Applications:**
- Local constraint projectors (e.g., Gutzwiller projectors)
- Subspace projections in variational methods [sachdev2011QuantumPhaseTransitions](@cite)

### Creating Projector Variables

````julia
reg_proj, (P, Q) = site_variables(Projector(), 2, 2)
````

````
(Registry(Projector, 2 sites, 4 variables), [[s1_1, s1_2], [s2_1, s2_2]])
````

- `P[1], P[2]` on **site 1** (with P² = P)
- `Q[1], Q[2]` on **site 2** (with Q² = Q)

### Idempotency Simplification

Repeated projectors automatically collapse:

````julia
@assert P[1] * P[1] == P[1]        # P² = P
@assert P[1] * P[1] * P[1] == P[1] # P³ = P
````

### Site Commutation with Idempotency

Cross-site operators commute, and idempotency is applied after canonicalization:

````julia
@assert Q[2] * P[1] == P[1] * Q[2]  # site commutation
````

More complex example: commutation brings same-site operators together,
then idempotency collapses repeated factors:

````julia
result = mono1(Q[1] * P[1] * Q[1] * P[1])
repr_with_registry(reg_proj, result)  # "s1_1 s2_1" (both collapse)
````

````
"s1_1 s2_1"
````

## 3) Unipotent: Involutory Operators (U² = 𝟙)

The `Unipotent()` rule models operators that square to the identity, such as
Pauli matrices σₓ, σᵧ, σᵤ. This is the natural algebra for ±1-valued observables.

**Quantum Information Applications:**
- Pauli observables in quantum computation [nielsenChuang2010](@cite)
- CHSH and other Bell inequalities with dichotomic measurements [klep2024State](@cite)
- Stabilizer formalism and error correction

**Many-Body Physics Applications:**
- Spin-½ systems with Pauli operators
- Ising model interactions

### Creating Unipotent Variables

````julia
reg_unip, (U, V) = site_variables(Unipotent(), 2, 2)
const ID_UNIP = one(U[1]);
````

- `U[1], U[2]` on **site 1** (with U² = 𝟙)
- `V[1], V[2]` on **site 2** (with V² = 𝟙)

### Involution Simplification

Squared operators simplify to the identity:

````julia
@assert U[1] * U[1] == ID_UNIP     # U² = 𝟙
@assert U[1] * U[1] * U[1] == U[1] # U³ = U
````

### Site Commutation with Involution

Cross-site commutation works as expected:

````julia
@assert V[2] * U[1] == U[1] * V[2]
````

### Cancellation After Canonicalization

Site-based canonicalization can bring paired operators together, enabling cancellation:

V₁ U₁ V₁ U₁ → (after site sort) → V₁ V₁ U₁ U₁ = 𝟙

````julia
@assert V[1] * U[1] * V[1] * U[1] == ID_UNIP
````

Simpler case: U₁² U₂ = U₂

````julia
@assert U[1] * U[1] * U[2] == U[2]
````

## 4) Generating Monomial Bases

`NCTSSoS.basis(registry, degree)` iterates coefficient-one polynomials for all
normal words through the specified degree. Collect it for custom moment rows.

````julia
basis_d1 = collect(NCTSSoS.basis(reg_nc, 1))
````

````
[1.0, s1_1, s2_1, s1_2, s2_2]
````

basis_d1: all monomials up to degree 1

Pretty-print the basis using the registry:

````julia
[repr_with_registry(reg_nc, m) for m in basis_d1]
````

````
5-element Vector{String}:
 "1.0"
 "s1_1"
 "s2_1"
 "s1_2"
 "s2_2"
````

````julia
basis_d2 = collect(NCTSSoS.basis(reg_nc, 2))
````

````
[1.0, s1_1, s2_1, s1_2, s2_2, s1_1 s1_1, s1_2 s1_1, s1_1 s2_1, … (17 entries)]
````

basis_d2: all monomials up to degree 2

Show count and a few examples:

````julia
length(basis_d2)  # total count
````

````
17
````

Pretty-print all degree-2 basis elements:

````julia
[repr_with_registry(reg_nc, m) for m in basis_d2]
````

````
17-element Vector{String}:
 "1.0"
 "s1_1"
 "s2_1"
 "s1_2"
 "s2_2"
 "s1_1 s1_1"
 "s1_2 s1_1"
 "s1_1 s2_1"
 "s2_1 s2_1"
 "s1_2 s2_1"
 "s2_2 s2_1"
 "s1_1 s1_2"
 "s1_2 s1_2"
 "s1_1 s2_2"
 "s2_1 s2_2"
 "s1_2 s2_2"
 "s2_2 s2_2"
````

## 5) Practical Example: Building Polynomials

Monomial algebras integrate seamlessly with `NCTSSoS.jl`'s polynomial optimization
pipeline. Here's a simple example of building a polynomial objective:

Coefficients create polynomials from monomials

````julia
poly = 1.0 * x[1] * x[2] + 2.0 * y[1] - 3.0 * x[1] * y[1]
````

````
(2.0) s2_1 + (-3.0) s1_1 s2_1 + s1_1 s1_2
````

Site commutation is automatic in polynomial arithmetic:

````julia
poly2 = 1.0 * y[1] * x[1] - 3.0 * x[1] * y[1]  # y₁ and x₁ commute!
collect(NCTSSoS.terms(poly2))  # terms with the same canonical word combine
````

````
1-element Vector{Tuple{ComplexF64, NCTSSoS.SiteWord{NCTSSoS.Free, UInt128, 4, 2}}}:
 (-2.0 + 0.0im, SiteWord(operator(site=1, label=1) operator(site=2, label=1)))
````

## Summary: Choosing the Right Algebra

| Your operators satisfy... | Use this algebra | Creation function |
|:--------------------------|:-----------------|:------------------|
| No relations (free) | `Free()` | `site_variables` |
| `P² = P` (projectors) | `Projector()` | `site_variables` |
| `U² = 𝟙` (involutions) | `Unipotent()` | `site_variables` |

All three provide:
- **Automatic site-based commutation** for operators on different sites
- **Canonical words** inspected through `NCTSSoS.terms` and `NCTSSoS.word_letters`
- **Registry-based pretty printing** for human-readable output
- **Optimization** through `Problem`, `Relaxation` and `solve`

For more advanced use cases including Pauli algebra (with cyclic products) and
fermionic/bosonic algebras, see the other examples in this documentation.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

