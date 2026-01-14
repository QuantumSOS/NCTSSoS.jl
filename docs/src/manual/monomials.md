# [Monomials and Normal Forms](@id monomials)

This page describes the **implementation-level** design of monomials in
`NCTSSoS.jl`. The goal is:

- **Canonical representation on construction**: algebraically equivalent inputs
  normalize to an identical representation immediately.
- **No unnecessary allocations**: simple monomials store as `Vector{T}` words.
- **Type-safe dispatch**: algebra type `A` drives method specialization.

## Core Types

### `NormalMonomial{A,T}`

`NormalMonomial{A,T}` represents an immutable **word** (a product of generators)
stored as a `Vector{T}`:

```julia
struct NormalMonomial{A<:AlgebraType, T<:Integer} <: AbstractMonomial{A,T}
    word::Vector{T}
end
```

**Invariant**: A `NormalMonomial` is always in the algebra-specific **normal form** for `A`.
This is enforced via `_validate_word(A, word)` (extended by the simplification
code in `src/simplification/*.jl`).

The `create_*_variables` functions (e.g., `create_pauli_variables`,
`create_fermionic_variables`) return `NormalMonomial` objects directly.
When you use these in arithmetic operations (multiplication, addition),
they produce `Polynomial` results.

### `Polynomial{A,T,C}`

`Polynomial{A,T,C}` represents a sum of monomials with coefficients:

```julia
struct Polynomial{A<:AlgebraType, T<:Integer, C<:Number} <: AbstractPolynomial{C}
    terms::Vector{Tuple{C, NormalMonomial{A,T}}}
end
```

**Key invariants**:
- `terms` is sorted by the `NormalMonomial` ordering
- Each monomial appears at most once (like terms are combined)
- Zero coefficients are removed

## Algebra Categories

The storage and simplification behavior depends on the algebra category:

### 1) `MonoidAlgebra`: single terms

Algebras where multiplication of monomials yields a single monomial:
- `NonCommutativeAlgebra` (no relations)
- `ProjectorAlgebra` (P² = P)
- `UnipotentAlgebra` (U² = I)

Multiplication returns a `Polynomial` with exactly one term.

### 2) `TwistedGroupAlgebra`: phase × monomial

Algebras where multiplication yields a scalar phase times a monomial:
- `PauliAlgebra` (σᵢσⱼ = iεᵢⱼₖσₖ)

Phases are tracked as `(im)^k` where `k ∈ 0:3`.

### 3) `PBWAlgebra`: sum of monomials

Algebras where multiplication can expand into multiple terms:
- `FermionicAlgebra` (anticommutation: aᵢaⱼ + aⱼaᵢ = δᵢⱼ)
- `BosonicAlgebra` (commutation: [aᵢ, aⱼ†] = δᵢⱼ)

Results require `Vector`-based storage for multiple terms.

## Simplification

To ensure *immediate* canonical forms, the simplification functions:

```julia
simplify(::Type{A}, word::Vector{T}) where {A<:AlgebraType}
```

Build an internal `NormalMonomial` from the raw word, apply algebra-specific
rules, and return a `Polynomial` in canonical form.

This means:
- Two words that are algebraically equal (under the relations of `A`) become
  equal `Polynomial`s right after construction.

## Term Interface

`Polynomial` exposes a lightweight term API:

- `terms(p::Polynomial)` returns an iterable of `(coefficient, NormalMonomial)` pairs
- `monomials(p)` returns a vector of the monomials
- `coefficients(p)` returns a vector of the coefficients

Usage:
```julia
for (coef, mono) in p.terms
    # coef is the coefficient (Number)
    # mono is a NormalMonomial{A,T} in canonical word form
end
```

## Why This Design

- **Julian dispatch**: `A` and `T` drive method specialization
- **Performance**: simple algebra elements stay compact, PBW expands only when
  mathematically necessary
- **Ergonomics**: consistent `zero/one/iszero` semantics across all algebras
