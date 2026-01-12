# [Monomials and Normal Forms](@id monomials)

This page describes the **implementation-level** design of monomials in
`NCTSSoS.jl`. The goal is:

- **Canonical representation on construction**: algebraically equivalent inputs
  normalize to an identical representation immediately.
- **No unnecessary allocations**: `MonoidAlgebra` and `TwistedGroupAlgebra`
  elements store as **0/1 term scalars**, not `Vector`s.
- **One name, many storage backends**: a single `Monomial` type, specialized by
  type parameters rather than distinct wrapper types.

## Layers: `NormalMonomial` vs `Monomial`

### `NormalMonomial{A,T}`

`NormalMonomial{A,T}` represents a **bare word** (a product of generators)
stored as a `Vector{T}`:

```julia
NormalMonomial{A,T} <: AbstractMonomial{A,T}
```

Invariant:

- A `NormalMonomial` is always in the algebra-specific **normal form** for `A`.
  This is enforced via `_validate_word!(A, word)` (extended by the
  simplification code in `src/simplification/*.jl`).

The `create_*_variables` functions (e.g., `create_pauli_variables`,
`create_fermionic_variables`) return `NormalMonomial` objects directly.
When you use these in arithmetic operations (multiplication, addition),
they are automatically wrapped in `Monomial` for simplification.

### `Monomial{A,T,C,W}`

`Monomial` represents a **simplified algebra element** (possibly with a phase or
as a sum of words):

```julia
struct Monomial{A,T,C,W} <: AbstractMonomial{A,T}
    coeffs::C
    words::W
end
```

Key invariant:

- `Monomial` is in **canonical expansion form** for its algebra category.
  The constructor from a raw word calls `simplify`, so equivalent words
  canonicalize immediately.

## Storage specializations (zero allocations where possible)

The storage backend is chosen by the algebra category (encoded in the type `A`):

### 1) `MonoidAlgebra`: scalar (0 or 1 term)

- Representation: `coeffs::UInt8`, `words::NormalMonomial{A,T}`
- Encoding: `0x00 => 0`, `0x01 => 1`

This avoids allocating `Vector`s for the common “just a word” case.

### 2) `TwistedGroupAlgebra` (Pauli): scalar (phase × word)

- Representation: `coeffs::UInt8`, `words::NormalMonomial{A,T}`
- Encoding:
  - `coeffs = k ∈ 0:3` represents the phase `(im)^k`
  - `coeffs = 0xff` is a reserved sentinel representing the **zero element**

So Pauli simplification stays allocation-free for single-term results while
preserving the phase.

### 3) `PBWAlgebra` (fermionic/bosonic): sum (many terms)

- Representation: `coeffs::Vector{Int}`, `words::Vector{NormalMonomial{A,T}}`
- Invariant: parallel vectors, same length, sorted by `words`, with like terms
  combined and zero coefficients removed.

This is the only category that actually needs vectors.

## Canonicalizing constructor

To ensure *immediate* canonical forms, the preferred constructor is:

```julia
Monomial{A}(word::Vector{T}) where {A<:AlgebraType,T<:Integer}
```

It calls:

- `simplify(A, word)` which first builds an internal `NormalMonomial` from the
  raw word (without requiring it to already satisfy the normal-form invariant),
  then runs the algebra-specific simplifier, returning a canonical `Monomial`.

This means:

- Two words that are algebraically equal (under the relations of `A`) become
  equal `Monomial`s right after construction.

## Term interface: `terms` + `iterate`

`Monomial` is **not** an `AbstractVector`. Instead, it exposes a lightweight term
API:

- `terms(m::Monomial)` returns an object you can iterate.
- Iteration yields `(c_internal, mono::NormalMonomial{A,T})` pairs.

Usage:

```julia
for (c_internal, mono) in terms(m)
    # c_internal is UInt8 (Monoid/Pauli) or Int (PBW)
    # mono is a NormalMonomial{A,T} in canonical word form
end
```

Utilities:

- `iszero(m::Monomial)` works for all storage backends.
- `zero(typeof(m))`, `one(typeof(m))`, `isone(m)`, `length(m)` reflect the number
  of terms in the canonical expansion (0/1 for scalar backends).

## Why this design

- **Julian dispatch**: `A` and `T` drive method specialization; `C` and `W` are
  storage-only and inferred by construction.
- **Performance**: monoid and twisted-group elements stay scalar (no heap), PBW
  expands only when mathematically necessary.
- **Ergonomics**: one `Monomial` name, one iteration protocol, consistent
  `zero/one/iszero` semantics.

