"""
    AlgebraType

Abstract type representing different algebraic structures with their
commutation/anticommutation relations.

Each concrete algebra type is a singleton struct that enables
multiple dispatch for simplification algorithms.

# Subtypes
- `NonCommutativeAlgebra`: Standard non-commutative variables (xy ≠ yx)
- `PauliAlgebra`: Pauli spin matrices satisfying σᵢ² = I and {σᵢ, σⱼ} = 2δᵢⱼ
- `FermionicAlgebra`: Fermionic creation/annihilation operators with {aᵢ, aⱼ†} = δᵢⱼ
- `BosonicAlgebra`: Bosonic creation/annihilation operators with [cᵢ, cⱼ†] = δᵢⱼ
- `ProjectorAlgebra`: Projector operators satisfying Pᵢ² = Pᵢ (idempotent)
- `UnipotentAlgebra`: Unipotent operators satisfying P² = I

# Design
Singleton types enable zero-cost dispatch on algebra operations:
```julia
simplify(::Type{PauliAlgebra}, m::NormalMonomial) = ...
```

# Examples
```jldoctest
julia> PauliAlgebra()
PauliAlgebra()

julia> FermionicAlgebra() isa AlgebraType
true
```
"""
abstract type AlgebraType end

"""
    MonoidAlgebra <: AlgebraType

Algebra category where the **normal form of a monomial is still a single monomial**
(i.e. closed on monomials).

This matches a *monoid algebra / monoid ring* viewpoint: basis elements are indexed
by elements of a monoid, and multiplying basis elements yields another basis element
in the same basis (extended linearly to polynomials). [@monoidRing]

In `NCTSSoS.jl`, this category includes:
- `NonCommutativeAlgebra` (free/noncommutative words)
- `ProjectorAlgebra` (idempotency P² = P)
- `UnipotentAlgebra` (involution U² = I)
"""
abstract type MonoidAlgebra <: AlgebraType end

"""
    TwistedGroupAlgebra <: AlgebraType

Algebra category where the **normal form of a monomial is a scalar/phase times a single monomial**.

This matches a *twisted group algebra* viewpoint: basis elements are indexed by a group,
but multiplication is twisted by a scalar 2-cocycle, so products satisfy
`u_g * u_h = α(g,h) u_{gh}`. [@twistedGroupAlgebra]

In `NCTSSoS.jl`, this category includes:
- `PauliAlgebra` (Pauli products generate phases {±1, ±i})
"""
abstract type TwistedGroupAlgebra <: AlgebraType end

"""
    PBWAlgebra <: AlgebraType

Algebra category where the **normal form of a monomial expands to a sum of monomials**.

This matches the computational noncommutative algebra pattern of *PBW algebras*: there is
an ordered monomial basis and a rewriting/normal-ordering procedure, but reordering introduces
lower terms, so a product may expand into multiple normal-form monomials. [@pbwAlgebraOscar]

In `NCTSSoS.jl`, this category includes:
- `FermionicAlgebra` (CAR / anticommutation)
- `BosonicAlgebra` (CCR / commutation)
"""
abstract type PBWAlgebra <: AlgebraType end

"""
    NonCommutativeAlgebra <: MonoidAlgebra

Generic non-commutative algebra used when no specific algebra is specified.

This algebra is *site-aware*: operators on different sites commute, while operators on the
same site remain non-commutative.

# Algebraic Rules
- **Within a site**: no simplification rules; operator order is preserved exactly.
- **Across sites**: operators commute; the canonical word is **stable-sorted by site**
  (stable sort preserves within-site order).

# Index Encoding
Variable indices must be **unsigned** and include bit-packed site information; see
`encode_index` / `decode_site`.
"""
struct NonCommutativeAlgebra <: MonoidAlgebra end

"""
    PauliAlgebra <: TwistedGroupAlgebra

Pauli spin matrix algebra.

# Algebraic Rules
- σᵢ² = I (involution: squares to identity)
- {σᵢ, σⱼ} = 2δᵢⱼ (anticommutation)
- σₓ σᵧ = i σᵤ, σᵧ σᵤ = i σₓ, σᵤ σₓ = i σᵧ (cyclic products)
- Operators on different sites commute

# Variable Encoding
Variables are ordered by site first: σx₁, σy₁, σz₁, σx₂, σy₂, σz₂, ...
For index `idx`: site = `(idx - 1) ÷ 3 + 1`, pauli_type = `(idx - 1) % 3`
(0=X, 1=Y, 2=Z)

# Integer Type
Typically uses unsigned integer types (self-adjoint operators).
Concrete type determined by VariableRegistry.
"""
struct PauliAlgebra <: TwistedGroupAlgebra end

"""
    FermionicAlgebra <: PBWAlgebra

Fermionic creation/annihilation operator algebra.

# Algebraic Rules
- {aᵢ, aⱼ†} = δᵢⱼ (anticommutation: creation-annihilation)
- {aᵢ, aⱼ} = 0 (anticommutation: annihilation-annihilation)
- {aᵢ†, aⱼ†} = 0 (anticommutation: creation-creation)

Normal ordering places all creation operators (a†) to the LEFT
of all annihilation operators (a).

# Variable Encoding
- Annihilation `aᵢ`: positive index `i`
- Creation `aᵢ†`: negative index `-i`

# Integer Type
Uses signed integer types (sign distinguishes creation/annihilation).
Concrete type determined by VariableRegistry.
"""
struct FermionicAlgebra <: PBWAlgebra end

"""
    BosonicAlgebra <: PBWAlgebra

Bosonic creation/annihilation operator algebra.

# Algebraic Rules
- [cᵢ, cⱼ†] = δᵢⱼ (commutation: creation-annihilation)
- [cᵢ, cⱼ] = 0 (commutation: annihilation-annihilation)
- [cᵢ†, cⱼ†] = 0 (commutation: creation-creation)

Normal ordering places all creation operators (c†) to the LEFT
of all annihilation operators (c).

**Key difference from fermions**: Commutation does not introduce
sign changes, but adds correction terms. Simplification may return
multiple terms.

# Variable Encoding
- Annihilation `cᵢ`: positive index `i`
- Creation `cᵢ†`: negative index `-i`

# Integer Type
Uses signed integer types (sign distinguishes creation/annihilation).
Concrete type determined by VariableRegistry.
"""
struct BosonicAlgebra <: PBWAlgebra end

"""
    ProjectorAlgebra <: MonoidAlgebra

Projector operator algebra.

# Algebraic Rules
- Pᵢ² = Pᵢ (idempotency: projectors square to themselves)
- Commutativity NOT enforced (treat as non-commutative)

# Variable Naming
Variables use symbols P₁, P₂, P₃, ...
Projectors are self-adjoint.

# Integer Type
Typically uses unsigned integer types (self-adjoint operators).
Concrete type determined by VariableRegistry.
"""
struct ProjectorAlgebra <: MonoidAlgebra end

"""
    UnipotentAlgebra <: MonoidAlgebra

Unipotent operator algebra.

# Algebraic Rules
- P² = I (squares to identity)
- No cyclic products or cross-operator interactions

**Note**: This is simpler than Pauli algebra which also has
cyclic product rules. Unipotent only removes consecutive pairs.

# Integer Type
Typically uses unsigned integer types (self-adjoint operators).
Concrete type determined by VariableRegistry.
"""
struct UnipotentAlgebra <: MonoidAlgebra end

# =============================================================================
# Default Coefficient Types
# =============================================================================

"""
    coeff_type(::Type{A}) where {A<:AlgebraType} -> Type{<:Number}

Return the default coefficient type for a given algebra type.

Different algebras naturally work with different coefficient types:
- `PauliAlgebra`: `ComplexF64` (Pauli products generate complex phases)
- All others: `Float64` (real coefficients suffice)

This is used by constructors like `Polynomial(m::NormalMonomial)` to infer
the appropriate coefficient type when not explicitly specified.

# Examples
```jldoctest
julia> using NCTSSoS

julia> coeff_type(PauliAlgebra) == ComplexF64
true

julia> coeff_type(NonCommutativeAlgebra) == Float64
true
```
"""
coeff_type(::Type{<:AlgebraType}) = Float64
coeff_type(::Type{PauliAlgebra}) = ComplexF64

# Show method for singleton algebra types
Base.show(io::IO, a::AlgebraType) = print(io, nameof(typeof(a)), "()")

# =============================================================================
# Site-Based Index Encoding
# =============================================================================
#
# For algebras that support site-based commutation (ProjectorAlgebra,
# UnipotentAlgebra, NonCommutativeAlgebra), we use bit-packed indices:
#
#   index = (operator_id << site_bits) | site  # 1-indexed
#
# The first n/4 bits encode the physical site, where n is the bit width.
# This allows operators on different sites to commute while maintaining
# non-commutative behavior within a site.

"""
    site_bits(::Type{T}) where {T<:Unsigned} -> Int

Number of bits used for site encoding. Fixed at n/4 where n is bit width.
Equivalently: `sizeof(T) * 8 ÷ 4 = sizeof(T) * 2`.

# Examples
```julia
site_bits(UInt8)   # 2
site_bits(UInt16)  # 4
site_bits(UInt32)  # 8
```
"""
@inline site_bits(::Type{T}) where {T<:Unsigned} = sizeof(T) * 8 ÷ 4

"""
    max_sites(::Type{T}) where {T<:Unsigned} -> Int

Maximum number of sites supported by unsigned type T.
Since sites are 1-indexed and stored directly, max is `2^k - 1` where k is site_bits.

# Examples
```julia
max_sites(UInt8)   # 3 (2 bits can store 1,2,3)
max_sites(UInt16)  # 15 (4 bits can store 1-15)
max_sites(UInt32)  # 255 (8 bits can store 1-255)
```
"""
@inline max_sites(::Type{T}) where {T<:Unsigned} = (1 << site_bits(T)) - 1

"""
    max_operators(::Type{T}) where {T<:Unsigned} -> Int

Maximum number of operators per site supported by unsigned type T.
Since operator IDs are 1-indexed, the maximum is `2^(n-k) - 1` where n is bit width and k is site bits.

# Examples
```julia
max_operators(UInt16)  # 4095
max_operators(UInt32)  # 16777215
```
"""
@inline max_operators(::Type{T}) where {T<:Unsigned} =
    (1 << (sizeof(T) * 8 - site_bits(T))) - 1

"""
    encode_index(::Type{T}, operator_id::Int, site::Int) where {T<:Unsigned} -> T

Encode operator_id and site into a single index of type T.
Site is 1-indexed and stored directly in the lower bits.

# Arguments
- `T`: The unsigned integer type to use for encoding
- `operator_id`: Operator identifier (1-indexed, must be ≤ max_operators(T))
- `site`: Physical site (1-indexed, must be ≤ max_sites(T))

# Examples
```julia
idx = encode_index(UInt16, 1, 3)  # 0x0013 (operator 1, site 3)
decode_site(idx)                   # 3
decode_operator_id(idx)            # 1
```
"""
@inline function encode_index(::Type{T}, operator_id::Integer, site::Integer) where {T<:Unsigned}
    k = site_bits(T)
    ms = max_sites(T)
    mo = max_operators(T)
    @boundscheck if site < 1 || site > ms
        throw(ArgumentError("Site $site out of range for $T (valid: 1-$ms)"))
    end
    @boundscheck if operator_id < 1 || operator_id > mo
        throw(ArgumentError("Operator $operator_id out of range for $T (valid: 1-$mo)"))
    end
    return T((T(operator_id) << k) | site)  # 1-indexed site storage
end

"""
    decode_site(idx::T) where {T<:Unsigned} -> Int

Extract site from encoded index (1-indexed).
"""
@inline function decode_site(idx::T) where {T<:Unsigned}
    k = site_bits(T)
    return Int(idx & ((1 << k) - 1))  # 1-indexed site stored directly
end

"""
    decode_operator_id(idx::T) where {T<:Unsigned} -> Int

Extract operator ID from encoded index.
"""
@inline function decode_operator_id(idx::T) where {T<:Unsigned}
    k = site_bits(T)
    return Int(idx >> k)
end


"""
    select_uint_type(n_operators::Integer, n_sites::Integer) -> Type{<:Unsigned}

Select the smallest unsigned integer type that can encode the given number
of operators and sites.

!!! note "Type Stability"
    This function returns a `Type`, not a value. The returned type depends on
    runtime values, so callers should use a function barrier or type annotation
    to maintain type stability in performance-critical code:
    ```julia
    T = select_uint_type(n_ops, n_sites)
    _inner_loop(T, data)  # function barrier
    ```

# Arguments
- `n_operators`: Number of distinct operators per site
- `n_sites`: Number of physical sites

# Returns
The smallest UInt type (UInt8, UInt16, UInt32, or UInt64) that can fit
the encoding.

# Examples
```julia
select_uint_type(10, 4)    # UInt16
select_uint_type(100, 10)  # UInt16
select_uint_type(1000, 100) # UInt32
```
"""
function select_uint_type(n_operators::Integer, n_sites::Integer)
    for T in (UInt8, UInt16, UInt32, UInt64)
        if n_sites <= max_sites(T) && n_operators <= max_operators(T)
            return T
        end
    end
    throw(ArgumentError("Cannot fit $n_operators operators × $n_sites sites in any UInt type"))
end
