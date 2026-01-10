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
│   ├── StateSymbol{ST,A,T}               # State expectation symbol
│   └── StateWord{ST,A,T}                 # Product of state symbols
└── ComposedMonomial{As,Ts}               # Multi-algebra tensor monomial
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
# Normal-form validation hook (extended by algebra-specific simplifiers)
# =============================================================================

"""
    _validate_word!(::Type{A}, word::Vector{T}) where {A<:AlgebraType,T<:Integer}

Validate that `word` is already in the algebra-specific normal form for `A`.

This is called by the `NormalMonomial{A,T}` constructor to enforce the invariant:
**a `NormalMonomial` is always in normal form**.

Algebra-specific methods are defined in `src/simplification/*.jl`.
The default implementation performs no checks.
"""
function _validate_word!(::Type{A}, ::Vector{T}) where {A<:AlgebraType,T<:Integer}
    return nothing
end

# Internal token for constructing a `NormalMonomial` from a raw (possibly non-normal)
# word without running `_validate_word!`. Kept non-exported on purpose.
struct _UnsafeNormalMonomial end
const _UNSAFE_NORMAL_MONOMIAL = _UnsafeNormalMonomial()

# Internal token for constructing a `NormalMonomial` from an owned, already-filtered
# word without copying. Still validates normal-form invariants.
struct _OwnedNormalMonomial end
const _OWNED_NORMAL_MONOMIAL = _OwnedNormalMonomial()

# Module-level constant for superscript display (avoids allocation in show loop)
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

    # Inner constructor: take ownership of an already-filtered word and validate invariants.
    function NormalMonomial{A,T}(word::Vector{T}, ::_OwnedNormalMonomial) where {A<:AlgebraType,T<:Integer}
        _validate_word!(A, word)
        new{A,T}(word)
    end

    # Inner constructor: filter zeros, validate invariants, then store an owned word.
    function NormalMonomial{A,T}(word::Vector{T}) where {A<:AlgebraType,T<:Integer}
        word_filtered = filter(!iszero, word)
        return NormalMonomial{A,T}(word_filtered, _OWNED_NORMAL_MONOMIAL)
    end

    # Internal escape hatch: construct without normal-form validation.
    # Used by simplification routines to accept raw words and return canonical outputs.
    function NormalMonomial{A,T}(word::Vector{T}, ::_UnsafeNormalMonomial) where {A<:AlgebraType,T<:Integer}
        word_filtered = filter(!iszero, word)
        new{A,T}(word_filtered)
    end
end

# Outer constructor: infer integer type from vector
function NormalMonomial{A}(word::Vector{T}) where {A<:AlgebraType,T<:Integer}
    return NormalMonomial{A,T}(word)
end

# Convenience constructor: defaults to NonCommutativeAlgebra
function NormalMonomial(word::Vector{T}) where {T<:Integer}
    return NormalMonomial{NonCommutativeAlgebra,T}(word)
end

"""
    Base.:*(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A<:MonoidAlgebra, T<:Integer}

Multiply two monomials of the same algebra type by concatenating their words.

For `MonoidAlgebra` types, the product is still a single monomial. This uses the
specialized `NormalMonomial{A}(word)` outer constructor to re-canonicalize the
concatenated word and returns a new `NormalMonomial`.

For algebras whose product introduces a phase (`TwistedGroupAlgebra`, e.g. Pauli)
or expands into multiple terms (`PBWAlgebra`, e.g. fermionic/bosonic), multiplication
is defined in the corresponding `src/simplification/*.jl` file and returns the
canonical expansion as a `Monomial` (iterable as `(c_internal, NormalMonomial)` pairs).

# Examples
```jldoctest
julia> m1 = NormalMonomial{NonCommutativeAlgebra}([1, 2]);

julia> m2 = NormalMonomial{NonCommutativeAlgebra}([3]);

julia> m = m1 * m2;

julia> m.word
3-element Vector{Int64}:
 1
 2
 3
```
"""
function Base.:*(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A<:MonoidAlgebra,T<:Integer}
    return NormalMonomial{A}(vcat(m1.word, m2.word))
end

"""
    Base.:^(m::NormalMonomial{A,T}, n::Integer)

Raise a normal-form monomial to a non-negative integer power.

For `MonoidAlgebra` (closed on monomials), this returns a `NormalMonomial`.
For algebras where multiplication expands to a scalar/phase times a monomial
(`TwistedGroupAlgebra`) or a sum of monomials (`PBWAlgebra`), this returns the
canonical expansion as a `Monomial` (iterable as `(c_internal, NormalMonomial)` pairs).

This never constructs an invalid `NormalMonomial`.
"""
function Base.:^(m::NormalMonomial{A,T}, n::Integer) where {A<:MonoidAlgebra,T<:Integer}
    n < 0 && throw(DomainError(n, "monomial exponent must be non-negative"))
    n == 0 && return one(m)
    n == 1 && return m
    return NormalMonomial{A}(repeat(m.word, n))
end

function Base.:^(m::NormalMonomial{A,T}, n::Integer) where {A<:Union{TwistedGroupAlgebra,PBWAlgebra},T<:Integer}
    n < 0 && throw(DomainError(n, "monomial exponent must be non-negative"))
    n == 0 && return simplify(A, T[])
    n == 1 && return simplify(m)
    return simplify(A, repeat(m.word, n))
end

"""
    Base.:(==)(m1::NormalMonomial, m2::NormalMonomial) -> Bool

Equality comparison for monomials.
Monomials of different algebra types are never equal (type-level distinction).
"""
function Base.:(==)(m1::NormalMonomial{A1,T1}, m2::NormalMonomial{A2,T2}) where {A1,A2,T1,T2}
    # Different algebra types are never equal
    A1 !== A2 && return false
    return m1.word == m2.word
end

"""
    Base.hash(m::NormalMonomial, h::UInt) -> UInt

Hash function for NormalMonomial. Includes algebra type for consistency with equality.

The hash includes both the algebra type `A` and the word vector, ensuring that
monomials from different algebras with the same word have different hashes.
This maintains the hash/equality contract: if `m1 == m2`, then `hash(m1) == hash(m2)`.
"""
Base.hash(m::NormalMonomial{A}, h::UInt) where {A} = hash(A, hash(m.word, h))

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
    return NormalMonomial{A}(T[])
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
Base.one(::Type{NormalMonomial}) = NormalMonomial{NonCommutativeAlgebra}(UInt[])

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
# Adjoint Operation
# =============================================================================

"""
    Base.adjoint(m::NormalMonomial{A,T}) where {A<:AlgebraType, T<:Integer}

Compute the adjoint (Hermitian conjugate) of a monomial.

This returns a new `NormalMonomial` in the algebra's normal form.

!!! note "Physics notation"
    This is the dagger (†) operation. You can use `m'` as shorthand for `adjoint(m)`.

# Behavior
- `PauliAlgebra`: identity on normal-form monomials (σ operators are Hermitian, and normal
  form has ≤1 operator per site).
- `FermionicAlgebra`, `BosonicAlgebra`: reverse word and flip creation/annihilation
  (negate indices).
- Other algebras: reverse word and re-canonicalize via the `NormalMonomial{A}(word)` constructor.
"""
Base.adjoint(m::NormalMonomial{PauliAlgebra}) = m

function Base.adjoint(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    new_word = reverse(m.word)
    if T <: Signed
        new_word .= .-new_word
    end
    return NormalMonomial{A}(new_word)
end

# =============================================================================
# Addition of Monomials
# =============================================================================

"""
    Base.:(+)(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A,T}

Add two monomials of the same algebra type. Returns a Polynomial with two terms.

This enables expressions like `x + x^2` where both operands are monomials.

# Examples
```jldoctest
julia> m1 = NormalMonomial{PauliAlgebra}([1]);

julia> m2 = NormalMonomial{PauliAlgebra}([1, 2]);

julia> p = m1 + m2;

julia> length(terms(p))
2

julia> p isa Polynomial{PauliAlgebra}
true
```
"""
function Base.:(+)(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # Convert both monomials to polynomials and add
    return Polynomial([Term(1.0, m1), Term(1.0, m2)])
end

"""
    Base.:(+)(m::NormalMonomial{A,T}, c::Number) where {A,T}
    Base.:(+)(c::Number, m::NormalMonomial{A,T}) where {A,T}

Add a scalar to a monomial. Returns a Polynomial with two terms (the monomial and a constant).

# Examples
```jldoctest
julia> m = NormalMonomial{PauliAlgebra}([1, 2]);

julia> p = m + 2.0;

julia> length(terms(p))
2

julia> p2 = 3.0 + m;

julia> length(terms(p2))
2
```
"""
function Base.:(+)(m::NormalMonomial{A,T}, c::Number) where {A<:AlgebraType,T<:Integer}
    # Create monomial term + constant term
    I = one(NormalMonomial{A,T})  # Identity monomial for constant
    return Polynomial([Term(float(c), I), Term(1.0, m)])
end

Base.:(+)(c::Number, m::NormalMonomial) = m + c

"""
    Base.:(-)(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A,T}

Subtract one monomial from another. Returns a Polynomial with two terms.

# Examples
```jldoctest
julia> m1 = NormalMonomial{PauliAlgebra}([1]);

julia> m2 = NormalMonomial{PauliAlgebra}([1, 2]);

julia> p = m1 - m2;

julia> coefficients(p)
2-element Vector{Float64}:
  1.0
 -1.0
```
"""
function Base.:(-)(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    # Convert both monomials to polynomials and subtract
    return Polynomial([Term(1.0, m1), Term(-1.0, m2)])
end

"""
    Base.:(-)(m::NormalMonomial{A,T}) where {A,T}

Negate a monomial. Returns a Term with negated coefficient.

# Examples
```jldoctest
julia> m = NormalMonomial{PauliAlgebra}([1, 2]);

julia> t = -m;

julia> t.coefficient
-1.0

julia> t.monomial === m
true
```
"""
Base.:(-)(m::NormalMonomial) = Term(-1.0, m)

"""
    Base.:(-)(m::NormalMonomial{A,T}, c::Number) where {A,T}
    Base.:(-)(c::Number, m::NormalMonomial{A,T}) where {A,T}

Subtract operations involving monomials and scalars. Returns a Polynomial.

# Examples
```jldoctest
julia> m = NormalMonomial{PauliAlgebra}([1, 2]);

julia> p = m - 2.0;

julia> length(terms(p))
2

julia> p2 = 3.0 - m;

julia> length(terms(p2))
2
```
"""
function Base.:(-)(m::NormalMonomial{A,T}, c::Number) where {A<:AlgebraType,T<:Integer}
    I = one(NormalMonomial{A,T})
    return Polynomial([Term(1.0, m), Term(-float(c), I)])
end

function Base.:(-)(c::Number, m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    I = one(NormalMonomial{A,T})
    return Polynomial([Term(float(c), I), Term(-1.0, m)])
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

# =============================================================================
# Expectation Value (identity for regular monomials)
# =============================================================================

"""
    expval(m::NormalMonomial{A,T}) where {A,T} -> NormalMonomial{A,T}

Return the monomial unchanged (identity operation).

For regular (non-state) monomials, `expval` is an identity operation.
This exists for API compatibility with `NCStateWord`, where `expval`
collapses the state word to a `StateWord`.

!!! todo "Refactor to return StateSymbol"
    This should return `StateSymbol{ST}(m)` instead of `m` for consistency
    with `expval(::NCStateWord)` which returns `StateWord`. Key decisions:
    - Which `StateType` to use: `Arbitrary` (involution canon) or `MaxEntangled` (cyclic symmetric canon)
    - Add `symmetric_canon(::StateSymbol)` method
    - Update optimization code: `monomap` dict keys from `NormalMonomial` to `StateSymbol`
    See also: `StateSymbol` in `src/states/word.jl`

# Examples
```julia
julia> m = NormalMonomial{NonCommutativeAlgebra,UInt8}([1, 2, 3]);
julia> expval(m) === m
true
```
"""
expval(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer} = m

"""
    coeff_type(::Type{NormalMonomial{A,T}}) where {A,T} -> Type{<:Number}

Return the default coefficient type for a NormalMonomial of algebra type A.
"""
coeff_type(::Type{NormalMonomial{A,T}}) where {A<:AlgebraType,T<:Integer} = coeff_type(A)

# Instance method
coeff_type(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer} = coeff_type(NormalMonomial{A,T})

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
# User-facing Monomial (normal-form result)
# =============================================================================

"""
    Monomial{A<:AlgebraType,T<:Integer,C,W} <: AbstractMonomial{A,T}

Unified user-facing representation of a simplified monomial *element* in algebra `A`.

This type replaces the need for separate wrappers like:
- `MonoidMonomial` (single normal monomial)
- `TGMonomial` (phase × normal monomial)
- `PBWMonomial` (∑ coeffᵢ × normal_monoᵢ)

`Monomial` stores the canonical expansion in a representation specialized by algebra
category:
- `MonoidAlgebra` / `TwistedGroupAlgebra`: **0 or 1 term**, stored scalarly (no vectors)
- `PBWAlgebra`: **0 or many terms**, stored as parallel vectors `(coeffs, words)`

The key invariant is: **construction and simplification return canonical normal-form
representations**, so algebraically equivalent inputs normalize to identical `Monomial`s.

`c_internal` uses algebra-specific encodings:
- `MonoidAlgebra`: `UInt8` with `0x00 => 0`, `0x01 => 1`
- `TwistedGroupAlgebra`: `UInt8` phase `k` representing `(im)^k` for `k ∈ 0:3`,
  plus `0xff` as a reserved sentinel for the zero element
- `PBWAlgebra`: `Vector{Int}` coefficients from rewriting/normal-ordering
"""
struct Monomial{A<:AlgebraType,T<:Integer,C,W} <: AbstractMonomial{A,T}
    coeffs::C
    words::W
end

# Reserved internal coefficient sentinel for "zero" in TwistedGroupAlgebra monomials.
const _TG_ZERO = typemax(UInt8)  # 0xff

# Concrete storage aliases (internal convenience).
const ScalarMonomial{A,T} = Monomial{A,T,UInt8,NormalMonomial{A,T}}
const PBWMonomial{A,T} = Monomial{A,T,Vector{Int},Vector{NormalMonomial{A,T}}}

# -----------------------------
# Monomial constructors
# -----------------------------

@inline _unit_internal_coeff(::Type{<:MonoidAlgebra}) = UInt8(1)
@inline _unit_internal_coeff(::Type{<:TwistedGroupAlgebra}) = UInt8(0)
@inline _unit_internal_coeff(::Type{<:PBWAlgebra}) = 1

@inline _zero_internal_coeff(::Type{<:MonoidAlgebra}) = UInt8(0)
@inline _zero_internal_coeff(::Type{<:TwistedGroupAlgebra}) = _TG_ZERO

"""
    Monomial{A}(word::Vector{T}) where {A<:AlgebraType,T<:Integer}
    Monomial(::Type{A}, word::Vector{T}) where {A<:AlgebraType,T<:Integer}

Construct a `Monomial` from a raw word by canonicalizing via `simplify(A, word)`.

This is the preferred constructor when you want algebraically equivalent words to
normalize to identical `Monomial` representations immediately.
"""
Monomial(::Type{A}, word::Vector{T}) where {A<:AlgebraType,T<:Integer} = Monomial{A}(word)

function Monomial{A}(word::Vector{T}) where {A<:AlgebraType,T<:Integer}
    return simplify(A, word)
end

function Monomial(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    return Monomial(_unit_internal_coeff(A), m)
end

function Monomial(c::UInt8, m::NormalMonomial{A,T}) where {A<:Union{MonoidAlgebra,TwistedGroupAlgebra},T<:Integer}
    return Monomial{A,T,UInt8,NormalMonomial{A,T}}(c, m)
end

function Monomial(c::Int, m::NormalMonomial{A,T}) where {A<:PBWAlgebra,T<:Integer}
    return PBWMonomial{A,T}(Int[c], [m])
end

function Monomial(coeffs::Vector{Int}, words::Vector{NormalMonomial{A,T}}) where {A<:PBWAlgebra,T<:Integer}
    length(coeffs) == length(words) || throw(ArgumentError(
        "Monomial coeffs/words length mismatch: $(length(coeffs)) vs $(length(words))"
    ))
    return PBWMonomial{A,T}(coeffs, words)
end

function Monomial(pairs::Vector{Tuple{Int,NormalMonomial{A,T}}}) where {A<:PBWAlgebra,T<:Integer}
    coeffs = Vector{Int}(undef, length(pairs))
    words = Vector{NormalMonomial{A,T}}(undef, length(pairs))
    @inbounds for i in eachindex(pairs)
        coeffs[i] = pairs[i][1]
        words[i] = pairs[i][2]
    end
    return Monomial(coeffs, words)
end

# -----------------------------
# Iteration / term interface
# -----------------------------

"""
    terms(m::Monomial)

Iterate the canonical expansion of a monomial element as `(c_internal, NormalMonomial)` pairs.
"""
terms(m::Monomial) = m

function Base.iterate(m::ScalarMonomial{A,T}) where {A<:Union{MonoidAlgebra,TwistedGroupAlgebra},T<:Integer}
    iszero(m) && return nothing
    return ((m.coeffs, m.words), nothing)
end

Base.iterate(::ScalarMonomial{A,T}, ::Nothing) where {A<:Union{MonoidAlgebra,TwistedGroupAlgebra},T<:Integer} = nothing

function Base.iterate(m::PBWMonomial{A,T}, i::Int=1) where {A<:PBWAlgebra,T<:Integer}
    i > length(m.coeffs) && return nothing
    return ((m.coeffs[i], m.words[i]), i + 1)
end

# Let `collect(m::Monomial)` infer a concrete element type (avoid `Vector{Any}`).
Base.IteratorEltype(::Type{<:Monomial}) = Base.HasEltype()
Base.eltype(::Type{ScalarMonomial{A,T}}) where {A<:AlgebraType,T<:Integer} = Tuple{UInt8,NormalMonomial{A,T}}
Base.eltype(::Type{PBWMonomial{A,T}}) where {A<:PBWAlgebra,T<:Integer} = Tuple{Int,NormalMonomial{A,T}}

# Optional indexing support (keeps callers simple; Monomial is not an AbstractVector).
function Base.getindex(m::ScalarMonomial{A,T}, i::Int) where {A<:Union{MonoidAlgebra,TwistedGroupAlgebra},T<:Integer}
    (i == 1 && !iszero(m)) && return (m.coeffs, m.words)
    throw(BoundsError(m, i))
end

function Base.getindex(m::PBWMonomial{A,T}, i::Int) where {A<:PBWAlgebra,T<:Integer}
    return (m.coeffs[i], m.words[i])
end

# -----------------------------
# Zero/one + size
# -----------------------------

Base.iszero(m::PBWMonomial) = isempty(m.coeffs)
Base.iszero(m::ScalarMonomial{A}) where {A<:MonoidAlgebra} = m.coeffs == 0x00
Base.iszero(m::ScalarMonomial{A}) where {A<:TwistedGroupAlgebra} = m.coeffs == _TG_ZERO

Base.isempty(m::Monomial) = iszero(m)

Base.length(m::PBWMonomial) = length(m.coeffs)
Base.length(m::ScalarMonomial) = iszero(m) ? 0 : 1

function Base.zero(::Type{ScalarMonomial{A,T}}) where {A<:MonoidAlgebra,T<:Integer}
    return Monomial{A,T,UInt8,NormalMonomial{A,T}}(_zero_internal_coeff(A), one(NormalMonomial{A,T}))
end

function Base.one(::Type{ScalarMonomial{A,T}}) where {A<:MonoidAlgebra,T<:Integer}
    return Monomial{A,T,UInt8,NormalMonomial{A,T}}(_unit_internal_coeff(A), one(NormalMonomial{A,T}))
end

function Base.zero(::Type{ScalarMonomial{A,T}}) where {A<:TwistedGroupAlgebra,T<:Integer}
    return Monomial{A,T,UInt8,NormalMonomial{A,T}}(_zero_internal_coeff(A), one(NormalMonomial{A,T}))
end

function Base.one(::Type{ScalarMonomial{A,T}}) where {A<:TwistedGroupAlgebra,T<:Integer}
    return Monomial{A,T,UInt8,NormalMonomial{A,T}}(_unit_internal_coeff(A), one(NormalMonomial{A,T}))
end

function Base.zero(::Type{PBWMonomial{A,T}}) where {A<:PBWAlgebra,T<:Integer}
    return PBWMonomial{A,T}(Int[], NormalMonomial{A,T}[])
end

function Base.one(::Type{PBWMonomial{A,T}}) where {A<:PBWAlgebra,T<:Integer}
    return PBWMonomial{A,T}(Int[1], [one(NormalMonomial{A,T})])
end

Base.zero(m::Monomial) = zero(typeof(m))
Base.one(m::Monomial) = one(typeof(m))

# Fallback identity for generic `one(Monomial)` (mirrors `one(NormalMonomial)` fallback).
Base.one(::Type{Monomial}) = Monomial(one(NormalMonomial))

function Base.isone(m::ScalarMonomial{A}) where {A<:MonoidAlgebra}
    return (m.coeffs == 0x01) && isone(m.words)
end

function Base.isone(m::ScalarMonomial{A}) where {A<:TwistedGroupAlgebra}
    return (m.coeffs == 0x00) && isone(m.words)
end

function Base.isone(m::PBWMonomial)
    length(m.coeffs) == 1 || return false
    m.coeffs[1] == 1 || return false
    return isone(m.words[1])
end

# Monomial "monomials" protocol: return all normal-form monomials appearing in the expansion.
monomials(m::PBWMonomial) = m.words
monomials(m::ScalarMonomial{A,T}) where {A<:AlgebraType,T<:Integer} =
    iszero(m) ? NormalMonomial{A,T}[] : [m.words]

# Degree of a Monomial expansion: maximum degree among its terms (0 for the zero monomial).
degree(m::PBWMonomial) = iszero(m) ? 0 : maximum(degree, m.words)
degree(m::ScalarMonomial) = iszero(m) ? 0 : degree(m.words)

# -----------------------------
# Equality / hashing
# -----------------------------

function Base.:(==)(m1::Monomial{A1}, m2::Monomial{A2}) where {A1<:AlgebraType,A2<:AlgebraType}
    A1 !== A2 && return false
    return (m1.coeffs == m2.coeffs) && (m1.words == m2.words)
end

Base.hash(m::Monomial{A}, h::UInt) where {A<:AlgebraType} = hash(A, hash(m.coeffs, hash(m.words, h)))

"""
    coeff_type(::Type{Monomial{A,T,C,W}}) where {A,T,C,W} -> Type{<:Number}

Return the numeric coefficient type for a simplified monomial expansion represented as
`Monomial`.

The stored coefficient `c` may be an internal encoding (e.g. `UInt8` for Pauli phases,
`UInt8` sentinels for monoid/twisted-group monomials, `Int` for PBW expansions), but the numeric
coefficient type is determined only by the algebra type `A`.
"""
coeff_type(::Type{Monomial{A,T,C,W}}) where {C,W,A<:AlgebraType,T<:Integer} = coeff_type(A)

coeff_type(x::Monomial) = coeff_type(typeof(x))

# Legacy: coefficient type for the old `Vector{(c, mono)}` representation.
coeff_type(::Type{Vector{Tuple{C,NormalMonomial{A,T}}}}) where {C,A<:AlgebraType,T<:Integer} = coeff_type(A)
coeff_type(x::Vector{Tuple{C,NormalMonomial{A,T}}}) where {C,A<:AlgebraType,T<:Integer} = coeff_type(typeof(x))

@inline _coeff_to_number(::Type{A}, c::UInt8) where {A<:MonoidAlgebra} =
    iszero(c) ? zero(coeff_type(A)) : one(coeff_type(A))

@inline _coeff_to_number(::Type{A}, ::Val{1}) where {A<:AlgebraType} = one(coeff_type(A))

@inline function _coeff_to_number(::Type{PauliAlgebra}, phase_k::UInt8)
    phase_k == _TG_ZERO && return ComplexF64(0.0, 0.0)
    return _phase_k_to_complex(phase_k)
end

@inline _coeff_to_number(::Type{A}, c::Number) where {A<:AlgebraType} = coeff_type(A)(c)
@inline _coeff_to_number(m::NormalMonomial{A}, c) where {A<:AlgebraType} = _coeff_to_number(A, c)

# Legacy: treat the old `Vector{(c, mono)}` representation as zero iff empty.
Base.iszero(pairs::Vector{Tuple{C,NormalMonomial{A,T}}}) where {C,A<:AlgebraType,T<:Integer} = isempty(pairs)

"""
    simplify(::Type{A}, word::Vector{T}) where {A<:AlgebraType,T<:Integer}

Simplify a raw word in algebra `A` without requiring it to already satisfy the
`NormalMonomial{A,T}` invariants.

This uses the inner `NormalMonomial{A,T}` constructor (no validation) and then
dispatches to `simplify(m::NormalMonomial{A,T})`, which is responsible for returning
the canonical expansion as a `Monomial`.
"""
function simplify(::Type{A}, word::Vector{T}) where {A<:AlgebraType,T<:Integer}
    return simplify(NormalMonomial{A,T}(word, _UNSAFE_NORMAL_MONOMIAL))
end

simplify(pairs::Vector{Tuple{C,NormalMonomial{A,T}}}) where {C,A<:AlgebraType,T<:Integer} = pairs
simplify(m::Monomial) = m

# -----------------------------------------------------------------------------
# Multiplication for simplified expansions: `Monomial` and legacy `Vector{(c, mono)}`
# -----------------------------------------------------------------------------

function Base.:*(
    terms::ScalarMonomial{PauliAlgebra,T},
    m::NormalMonomial{PauliAlgebra,T},
) where {T<:Integer}
    iszero(terms) && return zero(ScalarMonomial{PauliAlgebra,T})
    phase1 = terms.coeffs
    mono1 = terms.words

    prod = simplify(PauliAlgebra, vcat(mono1.word, m.word))
    iszero(prod) && return zero(ScalarMonomial{PauliAlgebra,T})
    phase2 = prod.coeffs
    mono2 = prod.words

    phase = UInt8((Int(phase1) + Int(phase2)) % 4)
    return Monomial(phase, mono2)
end

function Base.:*(
    m::NormalMonomial{PauliAlgebra,T},
    terms::ScalarMonomial{PauliAlgebra,T},
) where {T<:Integer}
    iszero(terms) && return zero(ScalarMonomial{PauliAlgebra,T})
    phase2 = terms.coeffs
    mono2 = terms.words

    prod = simplify(PauliAlgebra, vcat(m.word, mono2.word))
    iszero(prod) && return zero(ScalarMonomial{PauliAlgebra,T})
    phase1 = prod.coeffs
    mono1 = prod.words

    phase = UInt8((Int(phase1) + Int(phase2)) % 4)
    return Monomial(phase, mono1)
end

function Base.:*(
    t1::ScalarMonomial{PauliAlgebra,T},
    t2::ScalarMonomial{PauliAlgebra,T},
) where {T<:Integer}
    (iszero(t1) || iszero(t2)) && return zero(ScalarMonomial{PauliAlgebra,T})
    phase1 = t1.coeffs
    mono1 = t1.words
    phase2 = t2.coeffs
    mono2 = t2.words

    prod = simplify(PauliAlgebra, vcat(mono1.word, mono2.word))
    iszero(prod) && return zero(ScalarMonomial{PauliAlgebra,T})
    phase3 = prod.coeffs
    mono3 = prod.words

    phase = UInt8((Int(phase1) + Int(phase2) + Int(phase3)) % 4)
    return Monomial(phase, mono3)
end

function Base.:*(
    terms::PBWMonomial{A,T},
    m::NormalMonomial{A,T},
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    iszero(terms) && return zero(PBWMonomial{A,T})

    grouped = Dict{NormalMonomial{A,T},Int}()
    for (c1, mono1) in terms
        prod = mono1 * m
        for (cprod, monoprod) in prod
            grouped[monoprod] = get(grouped, monoprod, 0) + c1 * cprod
        end
    end

    out = Tuple{Int,NormalMonomial{A,T}}[]
    for (mono, c) in grouped
        c == 0 && continue
        push!(out, (c, mono))
    end
    sort!(out, by=t -> t[2])
    return Monomial(out)
end

function Base.:*(
    m::NormalMonomial{A,T},
    terms::PBWMonomial{A,T},
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    iszero(terms) && return zero(PBWMonomial{A,T})

    grouped = Dict{NormalMonomial{A,T},Int}()
    for (c2, mono2) in terms
        prod = m * mono2
        for (cprod, monoprod) in prod
            grouped[monoprod] = get(grouped, monoprod, 0) + c2 * cprod
        end
    end

    out = Tuple{Int,NormalMonomial{A,T}}[]
    for (mono, c) in grouped
        c == 0 && continue
        push!(out, (c, mono))
    end
    sort!(out, by=t -> t[2])
    return Monomial(out)
end

function Base.:*(
    t1::PBWMonomial{A,T},
    t2::PBWMonomial{A,T},
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    (iszero(t1) || iszero(t2)) && return zero(PBWMonomial{A,T})

    grouped = Dict{NormalMonomial{A,T},Int}()
    for (c1, mono1) in t1
        for (c2, mono2) in t2
            prod = mono1 * mono2
            for (cprod, monoprod) in prod
                grouped[monoprod] = get(grouped, monoprod, 0) + c1 * c2 * cprod
            end
        end
    end

    out = Tuple{Int,NormalMonomial{A,T}}[]
    for (mono, c) in grouped
        c == 0 && continue
        push!(out, (c, mono))
    end
    sort!(out, by=t -> t[2])
    return Monomial(out)
end

function Base.:*(
    terms::ScalarMonomial{A,T},
    m::NormalMonomial{A,T},
) where {A<:MonoidAlgebra,T<:Integer}
    iszero(terms) && return zero(ScalarMonomial{A,T})
    mono1 = terms.words
    return simplify(A, vcat(mono1.word, m.word))
end

function Base.:*(
    m::NormalMonomial{A,T},
    terms::ScalarMonomial{A,T},
) where {A<:MonoidAlgebra,T<:Integer}
    iszero(terms) && return zero(ScalarMonomial{A,T})
    mono2 = terms.words
    return simplify(A, vcat(m.word, mono2.word))
end

function Base.:*(
    t1::ScalarMonomial{A,T},
    t2::ScalarMonomial{A,T},
) where {A<:MonoidAlgebra,T<:Integer}
    (iszero(t1) || iszero(t2)) && return zero(ScalarMonomial{A,T})
    mono1 = t1.words
    mono2 = t2.words
    return simplify(A, vcat(mono1.word, mono2.word))
end

function Base.:*(
    terms::Vector{Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}},
    m::NormalMonomial{PauliAlgebra,T},
) where {T<:Integer}
    isempty(terms) && return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[]
    phase1, mono1 = terms[1]
    prod = simplify(PauliAlgebra, vcat(mono1.word, m.word))
    isempty(prod) && return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[]
    phase2, mono2 = prod[1]
    phase = UInt8((Int(phase1) + Int(phase2)) % 4)
    return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[(phase, mono2)]
end

function Base.:*(
    m::NormalMonomial{PauliAlgebra,T},
    terms::Vector{Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}},
) where {T<:Integer}
    isempty(terms) && return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[]
    phase2, mono2 = terms[1]
    prod = simplify(PauliAlgebra, vcat(m.word, mono2.word))
    isempty(prod) && return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[]
    phase1, mono1 = prod[1]
    phase = UInt8((Int(phase1) + Int(phase2)) % 4)
    return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[(phase, mono1)]
end

function Base.:*(
    t1::Vector{Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}},
    t2::Vector{Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}},
) where {T<:Integer}
    (isempty(t1) || isempty(t2)) && return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[]
    phase1, mono1 = t1[1]
    phase2, mono2 = t2[1]
    prod = simplify(PauliAlgebra, vcat(mono1.word, mono2.word))
    isempty(prod) && return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[]
    phase3, mono3 = prod[1]
    phase = UInt8((Int(phase1) + Int(phase2) + Int(phase3)) % 4)
    return Tuple{UInt8,NormalMonomial{PauliAlgebra,T}}[(phase, mono3)]
end

function Base.:*(
    terms::Vector{Tuple{Int,NormalMonomial{A,T}}},
    m::NormalMonomial{A,T},
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    isempty(terms) && return Tuple{Int,NormalMonomial{A,T}}[]

    grouped = Dict{NormalMonomial{A,T},Int}()
    for (c1, mono1) in terms
        prod = mono1 * m
        for (cprod, monoprod) in prod
            grouped[monoprod] = get(grouped, monoprod, 0) + c1 * cprod
        end
    end

    out = Tuple{Int,NormalMonomial{A,T}}[]
    for (mono, c) in grouped
        c == 0 && continue
        push!(out, (c, mono))
    end
    sort!(out, by=t -> t[2])
    return out
end

function Base.:*(
    m::NormalMonomial{A,T},
    terms::Vector{Tuple{Int,NormalMonomial{A,T}}},
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    isempty(terms) && return Tuple{Int,NormalMonomial{A,T}}[]

    grouped = Dict{NormalMonomial{A,T},Int}()
    for (c2, mono2) in terms
        prod = m * mono2
        for (cprod, monoprod) in prod
            grouped[monoprod] = get(grouped, monoprod, 0) + c2 * cprod
        end
    end

    out = Tuple{Int,NormalMonomial{A,T}}[]
    for (mono, c) in grouped
        c == 0 && continue
        push!(out, (c, mono))
    end
    sort!(out, by=t -> t[2])
    return out
end

function Base.:*(
    t1::Vector{Tuple{Int,NormalMonomial{A,T}}},
    t2::Vector{Tuple{Int,NormalMonomial{A,T}}},
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    (isempty(t1) || isempty(t2)) && return Tuple{Int,NormalMonomial{A,T}}[]

    grouped = Dict{NormalMonomial{A,T},Int}()
    for (c1, mono1) in t1
        for (c2, mono2) in t2
            prod = mono1 * mono2
            for (cprod, monoprod) in prod
                grouped[monoprod] = get(grouped, monoprod, 0) + c1 * c2 * cprod
            end
        end
    end

    out = Tuple{Int,NormalMonomial{A,T}}[]
    for (mono, c) in grouped
        c == 0 && continue
        push!(out, (c, mono))
    end
    sort!(out, by=t -> t[2])
    return out
end

function Base.:*(
    terms::Vector{Tuple{Val{1},NormalMonomial{A,T}}},
    m::NormalMonomial{A,T},
) where {A<:MonoidAlgebra,T<:Integer}
    isempty(terms) && return Tuple{Val{1},NormalMonomial{A,T}}[]
    _, mono1 = terms[1]
    return simplify(A, vcat(mono1.word, m.word))
end

function Base.:*(
    m::NormalMonomial{A,T},
    terms::Vector{Tuple{Val{1},NormalMonomial{A,T}}},
) where {A<:MonoidAlgebra,T<:Integer}
    isempty(terms) && return Tuple{Val{1},NormalMonomial{A,T}}[]
    _, mono2 = terms[1]
    return simplify(A, vcat(m.word, mono2.word))
end

function Base.:*(
    t1::Vector{Tuple{Val{1},NormalMonomial{A,T}}},
    t2::Vector{Tuple{Val{1},NormalMonomial{A,T}}},
) where {A<:MonoidAlgebra,T<:Integer}
    (isempty(t1) || isempty(t2)) && return Tuple{Val{1},NormalMonomial{A,T}}[]
    _, mono1 = t1[1]
    _, mono2 = t2[1]
    return simplify(A, vcat(mono1.word, mono2.word))
end
