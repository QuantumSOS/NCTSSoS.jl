# =============================================================================
# Monomial Types for Non-Commutative Algebras
# =============================================================================
#
# This file defines the core monomial representations for non-commutative
# polynomial optimization:
#
#   AbstractTensorMonomial{As}      # Abstract base for tensor-product monomials
#   ├── AbstractMonomial{A,T}       # Single-algebra monomials
#   │   ├── NormalMonomial{A,T}     # Canonical word representation (immutable)
#   └── ComposedMonomial{As,Ts}     # Multi-algebra tensor products (defined elsewhere)
#
# Design Principles:
# - Algebra type `A` is a type parameter for zero-overhead dispatch
# - `NormalMonomial` is always in canonical form (enforced by constructors)
# - All simplification goes through `simplify(A, word)` which returns `Monomial`
#
# Coefficient Encodings (internal):
# - MonoidAlgebra: UInt8 (0x00 = zero, 0x01 = one)
# - TwistedGroupAlgebra: UInt8 phase k where coeff = (im)^k
# - PBWAlgebra: Vector{Int} for multi-term expansions
#
# See also: `src/types/algebra.jl` for AlgebraType hierarchy
# =============================================================================

"""
    AbstractTensorMonomial{As<:Tuple}

Abstract supertype for monomials living in an **ordered tensor-product algebra**
`⊗_{i} As[i]`.

The signature `As` is **semantic and ordered**:
`Tuple{PauliAlgebra,FermionicAlgebra}` is different from
`Tuple{FermionicAlgebra,PauliAlgebra}`.

This type exists to enable signature-level multiple dispatch without wrapper
objects.

See also: [`AbstractMonomial`](@ref), [`ComposedMonomial`](@ref)
"""
abstract type AbstractTensorMonomial{As<:Tuple} end

"""
    AbstractMonomial{A<:AlgebraType, T<:Integer}

Abstract supertype for all monomial types, parameterized by algebra type and integer type.

The type hierarchy is:
```
AbstractTensorMonomial{As}
├── AbstractMonomial{A,T}                 # Single-algebra monomials (As == Tuple{A})
│   ├── NormalMonomial{A,T}               # Bare word (always in canonical form)
│   ├── Monomial{A,T,C,W}                 # User-facing wrapper (coeff + word(s))
│   ├── StateSymbol{ST,A,T}               # State expectation symbol (in src/states/)
│   └── StateWord{ST,A,T}                 # Product of state symbols (in src/states/)
└── ComposedMonomial{As,Ts}               # Multi-algebra tensor monomial (in composed.jl)
```

# Interface
All subtypes should implement:
- `degree(m)`: Return total degree (number of operators)
- `variable_indices(m)`: Return set of variable indices present
- `Base.isless(m1, m2)`: Ordering for sorting
- `Base.:(==)(m1, m2)`: Equality comparison
- `Base.hash(m, h)`: Hash function

See also: [`NormalMonomial`](@ref), [`AbstractTensorMonomial`](@ref)
"""
abstract type AbstractMonomial{A<:AlgebraType,T<:Integer} <: AbstractTensorMonomial{Tuple{A}} end

# =============================================================================
# Shared helpers (keep NormalMonomial / Monomial definitions DRY)
# =============================================================================

"""
    _eq_same_algebra(::Type{A1}, ::Type{A2}, x1, x2) -> Bool

Internal helper for type-safe equality comparison across algebra types.

Returns `false` immediately if `A1 !== A2` (different algebra types are never equal),
otherwise delegates to `x1 == x2` for the actual payload comparison.

This ensures monomials from different algebras with identical word vectors
are correctly distinguished (e.g., Pauli `[1,2]` ≠ Fermionic `[1,2]`).
"""
@inline function _eq_same_algebra(::Type{A1}, ::Type{A2}, x1, x2) where {A1<:AlgebraType,A2<:AlgebraType}
    A1 !== A2 && return false
    return x1 == x2
end

"""
    _hash_with_algebra(::Type{A}, payload, h::UInt) -> UInt

Internal helper that incorporates the algebra type `A` into the hash.

This maintains the hash/equality contract: since monomials from different algebras
are never equal (even with identical payloads), their hashes must also differ.
The algebra type is hashed first, then combined with the payload hash.
"""
@inline _hash_with_algebra(::Type{A}, payload, h::UInt) where {A<:AlgebraType} =
    hash(A, hash(payload, h))

"""
    coeff_type(::Type{T}) -> Type{<:Number}
    coeff_type(x) -> Type{<:Number}

Return the coefficient type for a simplify result type.

This enables compile-time determination of coefficient types for type-stable
processing of simplification results. Used by `ComposedMonomial` simplification
to determine appropriate coefficient types for the Cartesian product of terms.

For `NormalMonomial{A,T}`, returns `coeff_type(A)` (the algebra's default).
For `Polynomial{A,T,C}`, returns `C` (the explicit coefficient type).

This method is shared by `NormalMonomial`, simplified `Monomial` expansions, and
state monomials (`StateSymbol`, `StateWord`) that carry the same algebra type parameter.

# Examples
```jldoctest
julia> using NCTSSoS

julia> coeff_type(NormalMonomial{PauliAlgebra,Int64})
ComplexF64

julia> coeff_type(NormalMonomial{FermionicAlgebra,Int32})
Float64

julia> coeff_type(Polynomial{BosonicAlgebra,Int32,Float64})
Float64
```
"""
coeff_type(::Type{<:AbstractMonomial{A}}) where {A<:AlgebraType} = coeff_type(A)
coeff_type(m::AbstractMonomial) = coeff_type(typeof(m))

"""
    algebra_type(::Type{<:AbstractMonomial{A}}) where {A<:AlgebraType} -> Type{A}
    algebra_type(m::AbstractMonomial) -> Type{A}

Return the algebra type parameter for a monomial type or instance.
"""
algebra_type(::Type{<:AbstractMonomial{A}}) where {A<:AlgebraType} = A
algebra_type(m::AbstractMonomial) = algebra_type(typeof(m))

# =============================================================================
# Normal-form validation hook (extended by algebra-specific simplifiers)
# =============================================================================

"""
    _validate_word(::Type{A}, word::Vector{T}) where {A<:AlgebraType,T<:Integer}

Validate that `word` is already in the algebra-specific normal form for `A`.

This is called by the `NormalMonomial{A,T}` constructor to enforce the invariant:
**a `NormalMonomial` is always in normal form**.

Algebra-specific methods are defined in `src/simplification/*.jl`.
The default implementation throws an error for unimplemented algebras.

# Throws
- `ErrorException` if no validation method is defined for algebra `A`
"""
function _validate_word(::Type{A}, ::Vector{T}) where {A<:AlgebraType,T<:Integer}
    error("_validate_word not implemented for algebra type $A")
end

"""
    SUPERSCRIPT_EXPONENTS::Dict{Int,String}

Module-level constant mapping small exponents (2-9) to Unicode superscript characters.

Used by `Base.show` to display repeated variables with exponents (e.g., `x³` instead of `x^3`).
Avoids allocation in the display loop. Exponents outside this range fall back to `^n` notation.
"""
const SUPERSCRIPT_EXPONENTS = Dict(
    2 => "²",
    3 => "³",
    4 => "⁴",
    5 => "⁵",
    6 => "⁶",
    7 => "⁷",
    8 => "⁸",
    9 => "⁹",
)

"""
    NormalMonomial{A<:AlgebraType, T<:Integer} <: AbstractMonomial{A,T}

Represents an immutable monomial in word representation for non-commutative algebras.
The word vector represents a product of operators (e.g., [1,3,1,3] = xzxz).

# Fields
- `word::Vector{T}`: The monomial representation word

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

The struct is immutable for safety - simplification operations return new monomials
rather than mutating existing ones.

# Examples
```jldoctest
julia> m1 = NormalMonomial{PauliAlgebra}([1, 3, 1, 3]);

julia> m1.word
4-element Vector{Int64}:
 1
 3
 1
 3

julia> typeof(m1)
NormalMonomial{PauliAlgebra, Int64}
```

Different algebra types:
```jldoctest
julia> m_pauli = NormalMonomial{PauliAlgebra}(UInt16[1, 2, 3]);

julia> typeof(m_pauli)
NormalMonomial{PauliAlgebra, UInt16}

julia> m_fermi = NormalMonomial{FermionicAlgebra}(Int32[-1, 2, -3]);

julia> typeof(m_fermi)
NormalMonomial{FermionicAlgebra, Int32}
```

Equality comparison:
```jldoctest
julia> m1 = NormalMonomial{PauliAlgebra}([1, 3, 1, 3]);

julia> m2 = NormalMonomial{PauliAlgebra}([1, 3, 1, 2]);

julia> m1 == m2
false

julia> m3 = NormalMonomial{PauliAlgebra}([1, 3, 1, 3]);

julia> m1 == m3
true
```
"""
struct NormalMonomial{A<:AlgebraType,T<:Integer} <: AbstractMonomial{A,T}
    word::Vector{T}

    # Inner constructor: validate no zeros, then validate algebra invariants.
    function NormalMonomial{A,T}(word::Vector{T}) where {A<:AlgebraType,T<:Integer}
        any(iszero, word) && throw(ArgumentError("NormalMonomial word must not contain zeros"))
        _validate_word(A, word)
        new{A,T}(word)
    end

end

"""
    Base.:(==)(m1::NormalMonomial, m2::NormalMonomial) -> Bool

Equality comparison for monomials.
Monomials of different algebra types are never equal (type-level distinction).
"""
function Base.:(==)(m1::NormalMonomial{A1,T1}, m2::NormalMonomial{A2,T2}) where {A1,A2,T1,T2}
    return _eq_same_algebra(A1, A2, m1.word, m2.word)
end

"""
    Base.hash(m::NormalMonomial, h::UInt) -> UInt

Hash function for NormalMonomial. Includes algebra type for consistency with equality.

The hash includes both the algebra type `A` and the word vector, ensuring that
monomials from different algebras with the same word have different hashes.
This maintains the hash/equality contract: if `m1 == m2`, then `hash(m1) == hash(m2)`.
"""
Base.hash(m::NormalMonomial{A}, h::UInt) where {A} = _hash_with_algebra(A, m.word, h)

"""
    degree(m::NormalMonomial) -> Int

Compute the total degree of a monomial (length of word vector).
For non-commutative monomials in word representation, this is the number of operators.

# Examples
```jldoctest
julia> m1 = NormalMonomial([1, 3, 1, 3]);

julia> degree(m1)
4

julia> m2 = NormalMonomial([2, 2, 2]);

julia> degree(m2)
3

julia> m_zero = NormalMonomial(Int[]);

julia> degree(m_zero)
0
```
"""
degree(m::NormalMonomial) = length(m.word)

"""
    Base.isone(m::NormalMonomial) -> Bool

Check if a monomial is the multiplicative identity (empty word).

# Examples
```jldoctest
julia> m_identity = NormalMonomial{PauliAlgebra}(Int[]);

julia> isone(m_identity)
true

julia> m_not_identity = NormalMonomial{PauliAlgebra}([1, 2]);

julia> isone(m_not_identity)
false
```
"""
Base.isone(m::NormalMonomial) = isempty(m.word)

"""
    Base.one(::Type{NormalMonomial{A,T}}) where {A<:AlgebraType, T<:Integer}

Create the identity monomial (empty word).

# Examples
```jldoctest
julia> m_one = one(NormalMonomial{PauliAlgebra,Int64});

julia> isone(m_one)
true

julia> m_one.word
Int64[]
```
"""
function Base.one(::Type{NormalMonomial{A,T}}) where {A<:AlgebraType,T<:Integer}
    return NormalMonomial{A,T}(T[])
end

"""
    Base.one(m::NormalMonomial{A,T}) where {A,T}

Create the identity monomial for the same type as `m`.
"""
function Base.one(::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    return one(NormalMonomial{A,T})
end

"""
    Base.one(::Type{NormalMonomial}) -> NormalMonomial{NonCommutativeAlgebra,UInt}

Create the identity monomial (empty word) for the generic NormalMonomial type.
This fallback is needed for code that uses `one(NormalMonomial)` without type parameters.
"""
Base.one(::Type{NormalMonomial}) = NormalMonomial{NonCommutativeAlgebra,UInt}(UInt[])

"""
    Base.isless(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A,T} -> Bool

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
julia> m1 = NormalMonomial{PauliAlgebra}([1]);

julia> m2 = NormalMonomial{PauliAlgebra}([1, 2]);

julia> isless(m1, m2)  # degree 1 < degree 2
true

julia> m3 = NormalMonomial{PauliAlgebra}([2]);

julia> isless(m1, m3)  # same degree, [1] < [2] lexicographically
true

julia> isless(m3, m1)
false
```

Sorting works correctly:
```jldoctest
julia> monos = [NormalMonomial{PauliAlgebra}([2]), NormalMonomial{PauliAlgebra}([1, 2]), NormalMonomial{PauliAlgebra}([1])];

julia> sort!(monos);

julia> [m.word for m in monos]
3-element Vector{Vector{Int64}}:
 [1]
 [2]
 [1, 2]
```

See also: [`cmp`](@ref), [`degree`](@ref)
"""
function Base.isless(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # Compare by degree first (graded ordering)
    degree(m1) != degree(m2) && return degree(m1) < degree(m2)
    # Then lexicographic on word vector
    return m1.word < m2.word
end

# Fallback for cross-algebra comparison: provide descriptive error
function Base.isless(m1::NormalMonomial{A1}, m2::NormalMonomial{A2}) where {A1<:AlgebraType,A2<:AlgebraType}
    throw(ArgumentError("Cannot compare monomials of different algebras: $A1 vs $A2"))
end

# =============================================================================
# Display
# =============================================================================

"""
    Base.show(io::IO, m::NormalMonomial{A,T}) where {A,T}

Display a monomial. If a `:registry` is present in the IOContext, uses symbol names
from the registry. Otherwise, falls back to displaying raw indices.

When using a registry, consecutive identical variables are displayed with exponents:
- `[1, 1, 1]` with index 1 mapping to `:x` displays as `x³`
- `[1, 2, 2]` displays as `x₁y₂²`
- `[1, 1, 2, 2, 2]` displays as `x₁²y₂³`

# Examples
```julia
# Without registry (raw indices)
julia> m = NormalMonomial{PauliAlgebra}([1, 2, 3]);
julia> show(stdout, m)
[1, 2, 3]

# With registry (symbolic names)
julia> reg, (σx, σy, σz) = create_pauli_variables(1:2);
julia> m = σx[1] * σy[1];
julia> show(IOContext(stdout, :registry => reg), m)
σx₁σy₁

# With exponents for repeated variables
julia> m = NormalMonomial{PauliAlgebra}([1, 1, 1]);
julia> show(IOContext(stdout, :registry => reg), m)
σx₁³
```

See also: [`VariableRegistry`](@ref)
"""
function Base.show(io::IO, m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    if isempty(m.word)
        print(io, "𝟙")  # Identity symbol
        return nothing
    end

    registry = get(io, :registry, nothing)
    if registry !== nothing
        # Use symbols from registry with exponent grouping
        i = 1
        while i <= length(m.word)
            idx = m.word[i]
            sym = get(registry.idx_to_variables, idx, nothing)

            if sym !== nothing
                # Count consecutive occurrences of the same variable
                count = 1
                while i + count <= length(m.word) && m.word[i + count] == idx
                    count += 1
                end

                print(io, string(sym))

                # Print exponent if count > 1
                if count > 1
                    if haskey(SUPERSCRIPT_EXPONENTS, count)
                        print(io, SUPERSCRIPT_EXPONENTS[count])
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

"""
    monomials(m::NormalMonomial{A,T}) where {A,T} -> Vector{NormalMonomial{A,T}}

Return a single-element vector containing the monomial.

This provides compatibility with the `monomials` function for Polynomials,
allowing uniform iteration over monomials regardless of input type.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{PauliAlgebra}([1, 2]);

julia> monomials(m)
1-element Vector{NormalMonomial{PauliAlgebra, Int64}}:
 NormalMonomial{PauliAlgebra, Int64}([1, 2])
```
"""
monomials(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer} = [m]

# =============================================================================
# Iteration Protocol for NormalMonomial
# =============================================================================

"""
    Base.iterate(m::NormalMonomial{A,T}) where {A,T}
    Base.iterate(m::NormalMonomial{A,T}, state) where {A,T}

Iterate a NormalMonomial, yielding a single `(coefficient, NormalMonomial)` pair.

The coefficient is `one(coeff_type(A))` since NormalMonomial represents a bare
word with unit coefficient. This enables uniform iteration over both
NormalMonomial (single term) and Polynomial (multiple terms) in algorithms
that process bases of monomials.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{PauliAlgebra}([1, 2]);

julia> for (coef, mono) in m
           println("Coefficient: \$coef, Word: \$(mono.word)")
       end
Coefficient: 1.0 + 0.0im, Word: [1, 2]
```
"""
function Base.iterate(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    C = coeff_type(A)
    return ((one(C), m), nothing)
end

function Base.iterate(::NormalMonomial{A,T}, ::Nothing) where {A<:AlgebraType,T<:Integer}
    return nothing
end

# Length for iteration (always 1 term)
Base.length(::NormalMonomial) = 1

# =============================================================================
# Simplification (no-op for canonical NormalMonomial)
# =============================================================================

simplify(m::NormalMonomial) = m
simplify!(m::NormalMonomial) = m

# =============================================================================
# Adjoint (†)
# =============================================================================

function LinearAlgebra.adjoint(m::NormalMonomial{A,T}) where {A<:MonoidAlgebra,T<:Integer}
    isempty(m.word) && return m
    return NormalMonomial{A,T}(simplify!(A, reverse(m.word)))
end

LinearAlgebra.adjoint(m::NormalMonomial{PauliAlgebra,T}) where {T<:Unsigned} = copy(m)

LinearAlgebra.adjoint(m::NormalMonomial{A,T}) where {A<:PBWAlgebra,T<:Signed} = throw(ArgumentError("adjoint not implemented for PBWAlgebra"))
