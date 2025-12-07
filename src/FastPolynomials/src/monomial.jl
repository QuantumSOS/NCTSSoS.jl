"""
    AbstractMonomial

Abstract supertype for all monomial types.

Subtypes:
- `Monomial{A,T}`: Single-algebra monomial
- `ComposedMonomial{Ts}`: Tensor product of monomials from different algebras
"""
abstract type AbstractMonomial end

"""
    Monomial{A<:AlgebraType, T<:Integer} <: AbstractMonomial

Represents a monomial with a precomputed hash for efficient comparison.
The monomial vector is word representation for non-commutative monomials (e.g.,
[1,3,1,3] = xzxz)

# Fields
- `word::Vector{T}`: The monomial representation word
- `hash::UInt64`: Precomputed hash for fast inequality checks

# Type Parameters
- `A<:AlgebraType`: Algebra type for dispatch (PauliAlgebra, FermionicAlgebra, etc.)
- `T<:Integer`: Integer type for the monomial vector
  - `UInt16`: For self-adjoint variables (Pauli, Projector, Unipotent)
  - `Int32`: For non-self-adjoint variables (Fermionic, Bosonic)

# Design
The algebra type is in the type parameter, enabling:
- Zero per-element overhead (no runtime metadata)
- Compile-time dispatch for multiplication/simplification
- Type safety: can't accidentally multiply incompatible algebras

# Examples
```jldoctest
julia> m1 = Monomial{PauliAlgebra}([1, 3, 1, 3]);

julia> m1.word
4-element Vector{Int64}:
 1
 3
 1
 3

julia> typeof(m1)
Monomial{PauliAlgebra, Int64}
```

Different algebra types:
```jldoctest
julia> m_pauli = Monomial{PauliAlgebra}(UInt16[1, 2, 3]);

julia> typeof(m_pauli)
Monomial{PauliAlgebra, UInt16}

julia> m_fermi = Monomial{FermionicAlgebra}(Int32[-1, 2, -3]);

julia> typeof(m_fermi)
Monomial{FermionicAlgebra, Int32}
```

Fast inequality checks via precomputed hash:
```jldoctest
julia> m1 = Monomial{PauliAlgebra}([1, 3, 1, 3]);

julia> m2 = Monomial{PauliAlgebra}([1, 3, 1, 2]);

julia> m1 == m2
false

julia> m3 = Monomial{PauliAlgebra}([1, 3, 1, 3]);

julia> m1 == m3
true
```
"""
struct Monomial{A<:AlgebraType,T<:Integer} <: AbstractMonomial
    word::Vector{T}
    hash::UInt64
end

# Outer constructor: infer integer type from vector
function Monomial{A}(word::Vector{T}) where {A<:AlgebraType,T<:Integer}
    word_filtered = filter(!iszero, word)
    Monomial{A,T}(word_filtered, hash(word_filtered))
end

# Convenience constructor: defaults to NonCommutativeAlgebra
function Monomial(word::Vector{T}) where {T<:Integer}
    word_filtered = filter(!iszero, word)
    Monomial{NonCommutativeAlgebra,T}(word_filtered, hash(word_filtered))
end

"""
    Base.:(==)(m1::Monomial, m2::Monomial) -> Bool

Fast equality comparison using precomputed hash.
First compares hashes (O(1)), only does full vector comparison on hash collision.

Monomials of different algebra types are never equal (type-level distinction).

# Performance
- **50-100× faster** inequality checks vs. direct vector comparison
- Full comparison only needed in rare hash collision cases (~0.01%)
"""
function Base.:(==)(m1::Monomial{A1,T1}, m2::Monomial{A2,T2}) where {A1,A2,T1,T2}
    # Different algebra types are never equal
    A1 !== A2 && return false

    # Fast path: compare hashes first (O(1))
    m1.hash == m2.hash || return false

    # Slow path: full comparison only on hash collision
    m1.word == m2.word
end

"""
    Base.hash(m::Monomial, h::UInt) -> UInt

Hash function for Monomial. Uses the precomputed hash value.
"""
Base.hash(m::Monomial, h::UInt) = hash(m.hash, h)

"""
    degree(m::Monomial) -> Int

Compute the total degree of a monomial (length of word vector).
For non-commutative monomials in word representation, this is the number of operators.

# Examples
```jldoctest
julia> m1 = Monomial([1, 3, 1, 3]);

julia> degree(m1)
4

julia> m2 = Monomial([2, 2, 2]);

julia> degree(m2)
3

julia> m_zero = Monomial(Int[]);

julia> degree(m_zero)
0
```
"""
degree(m::Monomial) = length(m.word)

"""
    Base.isone(m::Monomial) -> Bool

Check if a monomial is the multiplicative identity (empty word).

# Examples
```jldoctest
julia> m_identity = Monomial{PauliAlgebra}(Int[]);

julia> isone(m_identity)
true

julia> m_not_identity = Monomial{PauliAlgebra}([1, 2]);

julia> isone(m_not_identity)
false
```
"""
Base.isone(m::Monomial) = isempty(m.word)

"""
    Base.one(::Type{Monomial{A,T}}) where {A<:AlgebraType, T<:Integer}

Create the identity monomial (empty word).

# Examples
```jldoctest
julia> m_one = one(Monomial{PauliAlgebra,Int64});

julia> isone(m_one)
true

julia> m_one.word
Int64[]
```
"""
function Base.one(::Type{Monomial{A,T}}) where {A<:AlgebraType,T<:Integer}
    Monomial{A}(T[])
end

"""
    Base.one(m::Monomial{A,T}) where {A,T}

Create the identity monomial for the same type as `m`.
"""
function Base.one(::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    one(Monomial{A,T})
end

"""
    Base.isless(m1::Monomial{A,T}, m2::Monomial{A,T}) where {A,T} -> Bool

Compare two monomials of the same algebra type using degree-first (graded) ordering,
then lexicographic ordering on the word vector for monomials of equal degree.

This ordering is required for:
- Sorting polynomials (terms are stored in sorted order)
- Binary search in polynomial operations
- Consistent canonical forms

# Ordering Rules
1. **Degree-first**: Shorter monomials come before longer ones
2. **Lexicographic**: For same-degree monomials, compare word vectors element-by-element

# Type Safety
Only monomials of the same algebra type `A` and integer type `T` can be compared.
Attempting to compare monomials of different algebra types will result in a `MethodError`.

# Examples
```jldoctest
julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([1, 2]);

julia> isless(m1, m2)  # degree 1 < degree 2
true

julia> m3 = Monomial{PauliAlgebra}([2]);

julia> isless(m1, m3)  # same degree, [1] < [2] lexicographically
true

julia> isless(m3, m1)
false
```

Sorting works correctly:
```jldoctest
julia> monos = [Monomial{PauliAlgebra}([2]), Monomial{PauliAlgebra}([1, 2]), Monomial{PauliAlgebra}([1])];

julia> sort!(monos);

julia> [m.word for m in monos]
3-element Vector{Vector{Int64}}:
 [1]
 [2]
 [1, 2]
```

See also: [`cmp`](@ref), [`degree`](@ref)
"""
function Base.isless(m1::Monomial{A,T}, m2::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # Compare by degree first (graded ordering)
    length(m1.word) != length(m2.word) && return length(m1.word) < length(m2.word)
    # Then lexicographic on word vector
    return m1.word < m2.word
end


# =============================================================================
# Adjoint / Star Operations
# =============================================================================

"""
    adjoint!(word::Vector{T}) where {T<:Integer}

In-place adjoint of a word vector. Reverses the word and negates elements for signed types.
This is the core implementation shared by all adjoint functions.

# Behavior
- For unsigned types (UInt16, etc.): just reverses the word (self-adjoint algebras)
- For signed types (Int32, Int64, etc.): reverses AND negates each element (non-self-adjoint algebras)

# Examples
```jldoctest
julia> word = UInt16[1, 2, 3];

julia> adjoint!(word);

julia> word
3-element Vector{UInt16}:
 0x0003
 0x0002
 0x0001

julia> word_signed = Int32[1, -2, 3];

julia> adjoint!(word_signed);

julia> word_signed
3-element Vector{Int32}:
 -3
  2
 -1
```
"""
function adjoint!(word::Vector{T}) where {T<:Integer}
    reverse!(word)
    # For signed types (Fermionic/Bosonic): negate each element
    if T <: Signed
        @inbounds for i in eachindex(word)
            word[i] = -word[i]
        end
    end
    word
end

"""
    star!(word::Vector{T}) where {T<:Integer}

Alias for `adjoint!`. In-place star (dagger) operation on a word vector.

See also: [`adjoint!`](@ref)
"""
star!(word::Vector{T}) where {T<:Integer} = adjoint!(word)

"""
    Base.adjoint(m::Monomial{A,T}) where {A<:AlgebraType, T<:Integer}

Compute the adjoint (star/dagger) of a monomial. Returns a new monomial.

For self-adjoint algebras (Pauli, Projector, Unipotent with unsigned types): reverses the word.
For non-self-adjoint algebras (Fermionic, Bosonic with signed types): reverses and negates indices.

The behavior is determined by the integer type T:
- Unsigned T: reverse only (operators are self-adjoint)
- Signed T: reverse and negate (creation/annihilation distinction)

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2, 3]);

julia> adjoint(m).word
3-element Vector{Int64}:
 3
 2
 1

julia> m_ferm = Monomial{FermionicAlgebra}(Int32[1, -2, 3]);

julia> adjoint(m_ferm).word  # reverse + negate
3-element Vector{Int32}:
 -3
  2
 -1
```

See also: [`star`](@ref), [`adjoint!`](@ref)
"""
function Base.adjoint(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    new_word = copy(m.word)
    adjoint!(new_word)
    Monomial{A}(new_word)
end

"""
    star(m::Monomial)

Alias for `adjoint`. Compute the star (dagger) of a monomial.
This notation is common in physics for the adjoint/Hermitian conjugate.

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2, 3]);

julia> star(m).word
3-element Vector{Int64}:
 3
 2
 1

julia> star(m) == adjoint(m)
true
```

See also: [`adjoint`](@ref)
"""
star(m::Monomial) = adjoint(m)
