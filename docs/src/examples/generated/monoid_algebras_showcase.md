```@meta
EditURL = "../literate/monoid_algebras_showcase.jl"
```

# [Monoid Algebra Showcase (NonCommutative / Projector / Unipotent)](@id monoid-algebras-showcase)

This page showcases the three `MonoidAlgebra` variants supported by `NCTSSoS.jl`.
These algebras provide the building blocks for modeling noncommutative polynomial
optimization problems that arise in quantum information and many-body physics.

## Capabilities Overview

Each monoid algebra type provides **automatic simplification** of operator products
according to its defining relations:

| Algebra Type | Simplification Rule | Typical Use Case |
|:-------------|:--------------------|:-----------------|
| `NonCommutativeAlgebra` | None (free words) | General NC polynomial optimization |
| `ProjectorAlgebra` | `P² = P` (idempotent) | Projective measurements, NPA hierarchy |
| `UnipotentAlgebra` | `U² = 𝟙` (involution) | Pauli observables, Bell inequalities |

All three algebras share a common feature: **site-based commutation**.
Operators acting on different physical sites automatically commute during
canonicalization, while the order of operators on the same site is preserved.

## Key Public APIs

Variable creation is done through these functions:
- [`create_noncommutative_variables`](@ref): Create free noncommutative generators
- [`create_projector_variables`](@ref): Create idempotent (P² = P) generators
- [`create_unipotent_variables`](@ref): Create involutory (U² = 𝟙) generators

Each function returns a `(registry, variable_groups)` tuple, where the registry
manages symbol ↔ index mappings and the variable groups provide convenient access
to the created monomials.

## References and Applications

- **`NonCommutativeAlgebra`** (free words): The standard modeling layer for
  noncommutative polynomial optimization and SDP hierarchies.
  See Burgdorf--Klep--Povh [burgdorfOptimizationPolynomialsNonCommuting2016](@cite),
  Wang--Magron [wangExploitingTermSparsity2021](@cite), and the comprehensive
  treatment in [magronSparsePolynomialOptimization2023](@cite).

- **`ProjectorAlgebra`** (`P² = P`): Projection operators and projective measurements.
  In device-independent quantum information, projective measurements appear naturally
  in the Navascués--Pironio--Acín (NPA) SDP hierarchy for bounding quantum correlations
  [navascues2007Bounding](@cite), [navascues2008Convergent](@cite).
  In many-body physics, projectors enforce local constraints and define subspaces;
  see Sachdev [sachdev2011QuantumPhaseTransitions](@cite).

- **`UnipotentAlgebra`** (`U² = 𝟙`): Involutory (±1-valued) observables model
  self-adjoint unitaries like Pauli measurements. These are fundamental in quantum
  information [nielsenChuang2010](@cite) and appear in Bell inequality optimization
  [klep2024State](@cite).

## Setup

````julia
using NCTSSoS
using NCTSSoS: decode_site, decode_operator_id
````

Helper: extract the single monomial from a *single-term* polynomial.

````julia
mono1(p) = monomials(p)[1]
````

````
mono1 (generic function with 1 method)
````

Pretty-print a `NormalMonomial` using the registry's symbols.

````julia
repr_with_registry(reg, m) = sprint(show, m; context=:registry => reg)
````

````
repr_with_registry (generic function with 1 method)
````

## 1) NonCommutativeAlgebra: Free Noncommutative Words

The `NonCommutativeAlgebra` represents free words with no algebraic relations
other than site-based commutation. This is the most general monoid algebra type.

### Creating Variables

Use `create_noncommutative_variables` with a vector of `(prefix, subscripts)` tuples.
Each prefix group is assigned to a distinct **site**, enabling automatic commutation
between groups.

````julia
reg_nc, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 1:2)])
````

````
(VariableRegistry with 4 variables: x₁, x₂, y₁, y₂, (NCTSSoS.NormalMonomial{NCTSSoS.NonCommutativeAlgebra, UInt8}[UInt8[0x05], UInt8[0x09]], NCTSSoS.NormalMonomial{NCTSSoS.NonCommutativeAlgebra, UInt8}[UInt8[0x0e], UInt8[0x12]]))
````

This creates:
- Variables `x[1], x[2]` on **site 1**
- Variables `y[1], y[2]` on **site 2**

A single generator is a `NormalMonomial` (immutable canonical form):

````julia
typeof(x[1])
````

````
NCTSSoS.NormalMonomial{NCTSSoS.NonCommutativeAlgebra, UInt8}
````

### Site-Based Commutation

Operators on different sites commute automatically. The canonical form groups
operators by site, preserving intra-site order:

````julia
lhs = x[2] * y[1] * x[1] * y[2]
rhs = x[2] * x[1] * y[1] * y[2]  # same canonical form

@assert lhs == rhs  # true: sites are sorted

repr_with_registry(reg_nc, mono1(lhs))  # "x₂x₁y₁y₂"
````

````
"x₂x₁y₁y₂"
````

### Intra-Site Non-Commutativity

Within a single site, operators do NOT commute (no relations are imposed):

````julia
@assert x[1] * x[2] != x[2] * x[1]
````

Show both representations to see they differ:

````julia
(repr_with_registry(reg_nc, mono1(x[1] * x[2])),  # "x₁x₂"
 repr_with_registry(reg_nc, mono1(x[2] * x[1])))  # "x₂x₁"
````

````
("x₁x₂", "x₂x₁")
````

### Inspecting the Canonical Word

Each index encodes both operator identity and site. You can decode them:

````julia
word = mono1(lhs).word
````

````
4-element Vector{UInt8}:
 0x09
 0x05
 0x0e
 0x12
````

word: the raw indices of the canonical monomial

Decode each index as (operator_id, site):
- operator_id: which variable (1=x₁, 2=x₂, 3=y₁, 4=y₂)
- site: which physical site (1 or 2)

````julia
[(decode_operator_id(i), decode_site(i)) for i in word]
````

````
4-element Vector{Tuple{Int64, Int64}}:
 (2, 1)
 (1, 1)
 (3, 2)
 (4, 2)
````

Output: [(2,1), (1,1), (3,2), (4,2)] means x₂x₁y₁y₂ with x's on site 1, y's on site 2

## 2) ProjectorAlgebra: Idempotent Operators (P² = P)

The `ProjectorAlgebra` models projection operators (idempotents). Any repeated
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
reg_proj, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 1:2)])
````

````
(VariableRegistry with 4 variables: P₁, P₂, Q₁, Q₂, (NCTSSoS.NormalMonomial{NCTSSoS.ProjectorAlgebra, UInt8}[UInt8[0x05], UInt8[0x09]], NCTSSoS.NormalMonomial{NCTSSoS.ProjectorAlgebra, UInt8}[UInt8[0x0e], UInt8[0x12]]))
````

- `P[1], P[2]` on **site 1** (with P² = P)
- `Q[1], Q[2]` on **site 2** (with Q² = Q)

### Idempotency Simplification

Repeated projectors automatically collapse:

````julia
@assert monomials(P[1] * P[1]) == [P[1]]        # P² = P
@assert monomials(P[1] * P[1] * P[1]) == [P[1]] # P³ = P
````

### Site Commutation with Idempotency

Cross-site operators commute, and idempotency is applied after canonicalization:

````julia
@assert monomials(Q[2] * P[1]) == monomials(P[1] * Q[2])  # site commutation
````

More complex example: commutation brings same-site operators together,
then idempotency collapses repeated factors:

````julia
result = mono1(Q[1] * P[1] * Q[1] * P[1])
repr_with_registry(reg_proj, result)  # "P₁Q₁" (both collapse)
````

````
"P₁Q₁"
````

## 3) UnipotentAlgebra: Involutory Operators (U² = 𝟙)

The `UnipotentAlgebra` models operators that square to the identity, such as
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
reg_unip, (U, V) = create_unipotent_variables([("U", 1:2), ("V", 1:2)])
const ID_UNIP = one(NormalMonomial{UnipotentAlgebra, eltype(indices(reg_unip))})
````

````
𝟙
````

- `U[1], U[2]` on **site 1** (with U² = 𝟙)
- `V[1], V[2]` on **site 2** (with V² = 𝟙)

### Involution Simplification

Squared operators simplify to the identity:

````julia
@assert monomials(U[1] * U[1]) == [ID_UNIP]     # U² = 𝟙
@assert monomials(U[1] * U[1] * U[1]) == [U[1]] # U³ = U
````

### Site Commutation with Involution

Cross-site commutation works as expected:

````julia
@assert monomials(V[2] * U[1]) == monomials(U[1] * V[2])
````

### Cancellation After Canonicalization

Site-based canonicalization can bring paired operators together, enabling cancellation:

V₁ U₁ V₁ U₁ → (after site sort) → V₁ V₁ U₁ U₁ = 𝟙

````julia
@assert monomials(V[1] * U[1] * V[1] * U[1]) == [ID_UNIP]
````

Simpler case: U₁² U₂ = U₂

````julia
@assert monomials(U[1] * U[1] * U[2]) == monomials(U[2])
````

## 4) Generating Monomial Bases

The function `get_ncbasis(registry, degree)` returns all distinct monomials up to
the specified degree. This is useful for constructing moment matrices in SDP
relaxations.

````julia
basis_d1 = get_ncbasis(reg_nc, 1)
````

````
5-element Vector{NCTSSoS.NormalMonomial{NCTSSoS.NonCommutativeAlgebra, UInt8}}:
 𝟙
 UInt8[0x05]
 UInt8[0x09]
 UInt8[0x0e]
 UInt8[0x12]
````

basis_d1: all monomials up to degree 1

Pretty-print the basis using the registry:

````julia
[repr_with_registry(reg_nc, m) for m in basis_d1]
````

````
5-element Vector{String}:
 "𝟙"
 "x₁"
 "x₂"
 "y₁"
 "y₂"
````

````julia
basis_d2 = get_ncbasis(reg_nc, 2)
````

````
17-element Vector{NCTSSoS.NormalMonomial{NCTSSoS.NonCommutativeAlgebra, UInt8}}:
 𝟙
 UInt8[0x05]
 UInt8[0x09]
 UInt8[0x0e]
 UInt8[0x12]
 UInt8[0x05, 0x05]
 UInt8[0x05, 0x09]
 UInt8[0x05, 0x0e]
 UInt8[0x05, 0x12]
 UInt8[0x09, 0x05]
 UInt8[0x09, 0x09]
 UInt8[0x09, 0x0e]
 UInt8[0x09, 0x12]
 UInt8[0x0e, 0x0e]
 UInt8[0x0e, 0x12]
 UInt8[0x12, 0x0e]
 UInt8[0x12, 0x12]
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
 "𝟙"
 "x₁"
 "x₂"
 "y₁"
 "y₂"
 "x₁²"
 "x₁x₂"
 "x₁y₁"
 "x₁y₂"
 "x₂x₁"
 "x₂²"
 "x₂y₁"
 "x₂y₂"
 "y₁²"
 "y₁y₂"
 "y₂y₁"
 "y₂²"
````

## 5) Practical Example: Building Polynomials

Monomial algebras integrate seamlessly with `NCTSSoS.jl`'s polynomial optimization
pipeline. Here's a simple example of building a polynomial objective:

Coefficients create polynomials from monomials

````julia
poly = 1.0 * x[1] * x[2] + 2.0 * y[1] - 3.0 * x[1] * y[1]
````

````
2.0 * UInt8[0x0e] + UInt8[0x05, 0x09] + -3.0 * UInt8[0x05, 0x0e]
````

Site commutation is automatic in polynomial arithmetic:

````julia
poly2 = 1.0 * y[1] * x[1] - 3.0 * x[1] * y[1]  # y₁ and x₁ commute!
monomials(poly2)  # simplified: terms combine if they have the same canonical form
````

````
1-element Vector{NCTSSoS.NormalMonomial{NCTSSoS.NonCommutativeAlgebra, UInt8}}:
 UInt8[0x05, 0x0e]
````

## Summary: Choosing the Right Algebra

| Your operators satisfy... | Use this algebra | Creation function |
|:--------------------------|:-----------------|:------------------|
| No relations (free) | `NonCommutativeAlgebra` | `create_noncommutative_variables` |
| `P² = P` (projectors) | `ProjectorAlgebra` | `create_projector_variables` |
| `U² = 𝟙` (involutions) | `UnipotentAlgebra` | `create_unipotent_variables` |

All three provide:
- **Automatic site-based commutation** for operators on different sites
- **Canonical forms** via `NormalMonomial` for efficient comparison and storage
- **Registry-based pretty printing** for human-readable output
- **Seamless integration** with `polyopt` and `cs_nctssos` for optimization

For more advanced use cases including Pauli algebra (with cyclic products) and
fermionic/bosonic algebras, see the other examples in this documentation.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

