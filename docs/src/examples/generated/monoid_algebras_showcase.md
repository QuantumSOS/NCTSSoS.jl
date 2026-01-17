# Monoid Algebra Showcase (NonCommutative / Projector / Unipotent)

This page showcases the three `MonoidAlgebra` variants supported by `NCTSSoS.jl`:

- `NonCommutativeAlgebra`: free words (no relations)
- `ProjectorAlgebra`: idempotents `P^2 = P`
- `UnipotentAlgebra`: involutions `U^2 = 𝟙`

All three use **site-based commutation**: generators acting on different physical sites commute.
Canonicalization is done by sorting sites while preserving intra-site order.

## References / where these algebras show up

- `NonCommutativeAlgebra` (free words, no relations): a standard modeling layer for
  noncommutative polynomial optimization and SDP hierarchies.
  See Burgdorf--Klep--Povh [burgdorfOptimizationPolynomialsNonCommuting2016](@cite),
  Wang--Magron [wangExploitingTermSparsity2021](@cite).

- `ProjectorAlgebra` (`P^2 = P`): projection operators and projective measurements.
  In device-independent quantum information, projective measurements appear naturally in the
  Navascu{\'e}s--Pironio--Ac{\'i}n (NPA) SDP hierarchy
  [navascues2007Bounding](@cite), [navascues2008Convergent](@cite).
  In many-body physics, projectors are ubiquitous (e.g. projecting onto local subspaces / enforcing
  local constraints); see e.g. Sachdev's textbook [sachdev2011QuantumPhaseTransitions](@cite).

- `UnipotentAlgebra` (`U^2 = I`): involutory (±1-valued) observables such as Pauli measurements,
  common throughout quantum information; see Nielsen--Chuang [nielsenChuang2010](@cite).

## Setup

````julia
using NCTSSoS
using NCTSSoS: decode_site, decode_operator_id
````

Helper: extract the single monomial from a *single-term* polynomial.
For words built via multiplication in a `MonoidAlgebra`, `monomials(p)` has length 1.

````julia
mono1(p) = monomials(p)[1]
````

Pretty-print a `NormalMonomial` using the registry's symbols.

````julia
repr_with_registry(reg, m) = sprint(show, m; context=:registry => reg)
````

## 0) Anatomy: sites, indices, and pretty printing

Variable creation APIs assign *sites* by grouping.
Here `x[*]` live on site 1 and `y[*]` live on site 2.

````julia
reg_nc, (x, y) = create_noncommutative_variables([("x", 1:2), ("y", 1:2)])
````

A single generator is already a `NormalMonomial` (not a polynomial container).

````julia
typeof(x[1])
````

Internally, each generator index is a single unsigned integer that encodes:
- `decode_site(idx)`         → physical site (1-indexed)
- `decode_operator_id(idx)`  → the variable's global id in the registry

In most workflows you never need to construct these indices manually because the registry
creation APIs already do it for you.

````julia
idx_x1 = x[1].word[1]
idx_y1 = y[1].word[1]
decode_site(idx_x1), decode_site(idx_y1)
````

By default, a monomial prints as its raw site-encoded indices.

````julia
repr(x[2] * y[1] * x[1])
````

With the registry in the IOContext, we get symbol names (`x₁`, `y₂`, ...).

````julia
m = mono1(x[2] * y[1] * x[1])
repr_with_registry(reg_nc, m)  # e.g. "x₂x₁y₁"
````

## 1) NonCommutativeAlgebra (free words)

**Use in QI/QMB:** Noncommutative words encode operator products (e.g. correlators) in SDP
hierarchies and polynomial optimization formulations; see
[burgdorfOptimizationPolynomialsNonCommuting2016](@cite), [wangExploitingTermSparsity2021](@cite).

Words on different sites commute (site-based canonicalization).

````julia
lhs_nc = x[2] * y[1] * x[1] * y[2]
rhs_nc = x[2] * x[1] * y[1] * y[2]

@assert lhs_nc == rhs_nc

repr_with_registry(reg_nc, mono1(lhs_nc)), repr_with_registry(reg_nc, mono1(rhs_nc))
````

You can inspect the canonical word as (operator_id, site) pairs.

````julia
word = mono1(lhs_nc).word
[(decode_operator_id(i), decode_site(i)) for i in word]
````

Within a site, order is preserved (no additional relations).

````julia
@assert x[1] * x[2] != x[2] * x[1]
````

## 2) ProjectorAlgebra (P² = P)

**Use in QI:** Projectors model projective measurements (PVM elements). In Bell/correlation contexts,
they are central to the NPA hierarchy [navascues2007Bounding](@cite), [navascues2008Convergent](@cite).

**Use in QMB:** Projection operators (idempotents) appear throughout many-body theory; see e.g.
[sachdev2011QuantumPhaseTransitions](@cite).

````julia
reg_proj, (P, Q) = create_projector_variables([("P", 1:2), ("Q", 1:2)])
````

Idempotency: Pᵢ² = Pᵢ (and repeated factors collapse)

````julia
@assert monomials(P[1] * P[1]) == [P[1]]
@assert monomials(P[1] * P[1] * P[1]) == [P[1]]
````

Explicit site commutation demo: Q (site 2) commutes past P (site 1).

````julia
@assert monomials(Q[2] * P[1]) == monomials(P[1] * Q[2])
````

A slightly larger example showing commutation + idempotency interaction:

````julia
repr_with_registry(reg_proj, mono1(Q[1] * P[1] * Q[1] * P[1]))  # → "P₁Q₁"
````

## 3) UnipotentAlgebra (U² = I)

**Use in QI/QMB:** Involutory generators model ±1-valued observables (self-adjoint unitaries),
e.g. Pauli measurements used throughout quantum information; see [nielsenChuang2010](@cite).

````julia
reg_unip, (U, V) = create_unipotent_variables([("U", 1:2), ("V", 1:2)])
const ID_UNIP = one(NormalMonomial{UnipotentAlgebra, eltype(indices(reg_unip))})
````

Involution: Uᵢ² = 𝟙

````julia
@assert monomials(U[1] * U[1]) == [ID_UNIP]
@assert monomials(U[1] * U[1] * U[1]) == [U[1]]
````

Explicit site commutation demo

````julia
@assert monomials(V[2] * U[1]) == monomials(U[1] * V[2])
````

Cancellation can happen after site canonicalization:

````julia
@assert monomials(V[1] * U[1] * V[1] * U[1]) == [ID_UNIP]
````

A simpler cancellation example (same site, adjacent factors):

````julia
@assert monomials(U[1] * U[1] * U[2]) == monomials(U[2])
````

## 4) Tiny basis example

`get_ncbasis(reg, d)` returns all distinct monomials up to (word) degree `d`.
For these MonoidAlgebra variants, the basis elements are `NormalMonomial`s.

````julia
basis_nc = get_ncbasis(reg_nc, 2)
length(basis_nc)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

