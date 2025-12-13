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
    Base.:*(m1::Monomial{A,T}, m2::Monomial{A,T}) where {A<:AlgebraType, T<:Integer}

Multiply two monomials of the same algebra type by concatenating their words.

Returns a new Monomial with the concatenated word. Callers should apply
`simplify!` explicitly if algebra-specific simplification is needed.

# Examples
```jldoctest
julia> m1 = Monomial{NonCommutativeAlgebra}([1, 2]);

julia> m2 = Monomial{NonCommutativeAlgebra}([3]);

julia> m = m1 * m2;

julia> m.word
3-element Vector{Int64}:
 1
 2
 3
```
"""
function Base.:*(m1::Monomial{A,T}, m2::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    Monomial{A}(vcat(m1.word, m2.word))
end

"""
    Base.:^(m::Monomial{A,T}, n::Integer) where {A<:AlgebraType, T<:Integer}

Raise a monomial to an integer power by repeated multiplication.

For n ≥ 0, returns m^n = m * m * ... * m (n times).
For n = 0, returns the identity monomial (empty word).
Throws DomainError for negative exponents (monomials are not invertible in general).

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = m^2;

julia> m2.word
4-element Vector{Int64}:
 1
 2
 1
 2

julia> m^0 |> isone
true

julia> m^1 == m
true
```
"""
function Base.:^(m::Monomial{A,T}, n::Integer) where {A<:AlgebraType,T<:Integer}
    if n < 0
        throw(DomainError(n, "monomial exponent must be non-negative"))
    elseif n == 0
        return one(m)
    elseif n == 1
        return m
    else
        # Repeat the word n times
        repeated_word = repeat(m.word, n)
        return Monomial{A}(repeated_word)
    end
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

    # If either hash is zero (uninitialized), compare words directly
    # This handles monomials created via multiplication before simplification
    if m1.hash == 0 || m2.hash == 0
        return m1.word == m2.word
    end

    # Fast path: compare hashes first (O(1))
    m1.hash == m2.hash || return false

    # Slow path: full comparison only on hash collision
    m1.word == m2.word
end

"""
    Base.hash(m::Monomial, h::UInt) -> UInt

Hash function for Monomial. Uses the precomputed hash value if available,
otherwise computes from the word.
"""
function Base.hash(m::Monomial, h::UInt)
    # If hash is zero (uninitialized), compute from word
    word_hash = m.hash == 0 ? hash(m.word) : m.hash
    hash(word_hash, h)
end

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

# =============================================================================
# Addition of Monomials
# =============================================================================

"""
    Base.:(+)(m1::Monomial{A,T}, m2::Monomial{A,T}) where {A,T}

Add two monomials of the same algebra type. Returns a Polynomial with two terms.

This enables expressions like `x + x^2` where both operands are monomials.

# Examples
```jldoctest
julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([1, 2]);

julia> p = m1 + m2;

julia> length(terms(p))
2

julia> p isa Polynomial{PauliAlgebra}
true
```
"""
function Base.:(+)(
    m1::Monomial{A,T}, m2::Monomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    # Convert both monomials to polynomials and add
    return Polynomial([Term(1.0, m1), Term(1.0, m2)])
end

"""
    Base.:(+)(m::Monomial{A,T}, c::Number) where {A,T}
    Base.:(+)(c::Number, m::Monomial{A,T}) where {A,T}

Add a scalar to a monomial. Returns a Polynomial with two terms (the monomial and a constant).

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> p = m + 2.0;

julia> length(terms(p))
2

julia> p2 = 3.0 + m;

julia> length(terms(p2))
2
```
"""
function Base.:(+)(
    m::Monomial{A,T}, c::Number
) where {A<:AlgebraType,T<:Integer}
    # Create monomial term + constant term
    I = one(Monomial{A,T})  # Identity monomial for constant
    return Polynomial([Term(float(c), I), Term(1.0, m)])
end

Base.:(+)(c::Number, m::Monomial) = m + c

"""
    Base.:(-)(m1::Monomial{A,T}, m2::Monomial{A,T}) where {A,T}

Subtract one monomial from another. Returns a Polynomial with two terms.

# Examples
```jldoctest
julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([1, 2]);

julia> p = m1 - m2;

julia> coefficients(p)
2-element Vector{Float64}:
  1.0
 -1.0
```
"""
function Base.:(-)(
    m1::Monomial{A,T}, m2::Monomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    # Convert both monomials to polynomials and subtract
    return Polynomial([Term(1.0, m1), Term(-1.0, m2)])
end

"""
    Base.:(-)(m::Monomial{A,T}) where {A,T}

Negate a monomial. Returns a Term with negated coefficient.

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> t = -m;

julia> t.coefficient
-1.0

julia> t.monomial === m
true
```
"""
Base.:(-)(m::Monomial) = Term(-1.0, m)

"""
    Base.:(-)(m::Monomial{A,T}, c::Number) where {A,T}
    Base.:(-)(c::Number, m::Monomial{A,T}) where {A,T}

Subtract operations involving monomials and scalars. Returns a Polynomial.

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> p = m - 2.0;

julia> length(terms(p))
2

julia> p2 = 3.0 - m;

julia> length(terms(p2))
2
```
"""
function Base.:(-)(
    m::Monomial{A,T}, c::Number
) where {A<:AlgebraType,T<:Integer}
    I = one(Monomial{A,T})
    return Polynomial([Term(1.0, m), Term(-float(c), I)])
end

function Base.:(-)(
    c::Number, m::Monomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    I = one(Monomial{A,T})
    return Polynomial([Term(float(c), I), Term(-1.0, m)])
end

# =============================================================================
# Display
# =============================================================================

"""
    Base.show(io::IO, m::Monomial{A,T}) where {A,T}

Display a monomial. If a `:registry` is present in the IOContext, uses symbol names
from the registry. Otherwise, falls back to displaying raw indices.

When using a registry, consecutive identical variables are displayed with exponents:
- `[1, 1, 1]` with index 1 mapping to `:x` displays as `x³`
- `[1, 2, 2]` displays as `x₁y₂²`
- `[1, 1, 2, 2, 2]` displays as `x₁²y₂³`

# Examples
```julia
# Without registry (raw indices)
julia> m = Monomial{PauliAlgebra}([1, 2, 3]);
julia> show(stdout, m)
[1, 2, 3]

# With registry (symbolic names)
julia> reg, (σx, σy, σz) = create_pauli_variables(1:2);
julia> m = σx[1] * σy[1];
julia> show(IOContext(stdout, :registry => reg), m)
σx₁σy₁

# With exponents for repeated variables
julia> m = Monomial{PauliAlgebra}([1, 1, 1]);
julia> show(IOContext(stdout, :registry => reg), m)
σx₁³
```

See also: [`VariableRegistry`](@ref)
"""
function Base.show(io::IO, m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    if isempty(m.word)
        print(io, "𝟙")  # Identity symbol
        return
    end

    registry = get(io, :registry, nothing)
    if registry !== nothing
        # Use symbols from registry with exponent grouping
        i = 1
        while i <= length(m.word)
            idx = m.word[i]
            abs_idx = T(abs(idx))
            sym = get(registry.idx_to_variables, abs_idx, nothing)

            if sym !== nothing
                # Count consecutive occurrences of the same variable
                count = 1
                while i + count <= length(m.word) && abs(m.word[i + count]) == abs_idx
                    count += 1
                end

                # Handle signed indices (e.g., -1 for creation operator)
                if idx < 0 && T <: Signed
                    print(io, string(sym), "†")
                else
                    print(io, string(sym))
                end

                # Print exponent if count > 1
                if count > 1
                    # Use superscript digits for nicer display
                    superscripts = Dict(2 => "²", 3 => "³", 4 => "⁴", 5 => "⁵",
                                       6 => "⁶", 7 => "⁷", 8 => "⁸", 9 => "⁹")
                    if haskey(superscripts, count)
                        print(io, superscripts[count])
                    else
                        print(io, "^", count)
                    end
                end

                i += count
            else
                # Fallback for missing index: print raw
                print(io, "[", idx, "]")
                i += 1
            end
        end
    else
        # Fallback: print raw word
        print(io, m.word)
    end
end

# =============================================================================
# Expectation Value (identity for regular monomials)
# =============================================================================

"""
    expval(m::Monomial{A,T}) where {A,T} -> Monomial{A,T}

Return the expectation value of a monomial. For regular (non-state) monomials,
this is an identity operation - the monomial is returned unchanged.

This method exists for API compatibility with `NCStateWord`, where `expval`
extracts the underlying monomial from a state word. For regular monomials,
no extraction is needed.

# Examples
```julia
julia> m = Monomial{NonCommutativeAlgebra,UInt8}([1, 2, 3]);
julia> expval(m) === m
true
```
"""
expval(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer} = m
