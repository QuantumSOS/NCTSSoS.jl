# =============================================================================
# PauliMonomial - Canonical Pauli monomial with phase
# =============================================================================

"""
    PauliMonomial{T<:Integer} <: AbstractMonomial{PauliAlgebra,T}

A Pauli monomial in canonical form with an accumulated algebraic phase.

PauliMonomial stores a canonical `Monomial{PauliAlgebra,T}` (at most one operator per site,
sorted by site) together with a phase encoded as `phase_k ∈ 0:3` representing `(im)^phase_k`.

# Fields
- `mono::Monomial{PauliAlgebra,T}`: Canonical monomial (≤1 op/site, sorted)
- `phase_k::UInt8`: Phase encoding where `(im)^phase_k` gives {1, i, -1, -i} for k ∈ {0,1,2,3}

# Type Parameters
- `T<:Integer`: Integer type for variable indices

# Invariants
- `mono` is always in canonical form (enforced by `_simplify_pauli_word!`)
- `phase_k ∈ 0:3`

# Phase Encoding
| phase_k | Complex Value |
|---------|---------------|
| 0       | 1             |
| 1       | i             |
| 2       | -1            |
| 3       | -i            |

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm = PauliMonomial([1, 2]);  # σx₁ σy₁ = i σz₁

julia> pm.phase_k  # phase = i = (im)^1
0x01

julia> pm.mono.word
1-element Vector{Int64}:
 3

julia> pm2 = PauliMonomial([1, 2, 1, 2]);  # σx₁ σy₁ σx₁ σy₁ = (i σz₁)² = -1

julia> pm2.phase_k  # phase = -1 = (im)^2
0x02

julia> isempty(pm2.mono.word)  # reduced to identity
true
```

See also: [`Monomial`](@ref), [`_simplify_pauli_word!`](@ref)
"""
struct PauliMonomial{T<:Integer} <: AbstractMonomial{PauliAlgebra,T}
    mono::Monomial{PauliAlgebra,T}
    phase_k::UInt8

    # Inner constructor from already-canonical data
    function PauliMonomial{T}(mono::Monomial{PauliAlgebra,T}, phase_k::UInt8) where {T<:Integer}
        @assert phase_k <= 3 "phase_k must be in 0:3, got $phase_k"
        new{T}(mono, phase_k)
    end
end

# =============================================================================
# Constructors
# =============================================================================

"""
    PauliMonomial(word::Vector{T}) where {T<:Integer}

Construct a PauliMonomial from a raw (non-canonical) word.

The word is automatically canonicalized using `_simplify_pauli_word!`, which:
1. Sorts operators by site
2. Applies Pauli product rules (σ²=I, σₓσᵧ=iσz, etc.)
3. Accumulates phase into `phase_k`

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm = PauliMonomial([2, 1]);  # σy₁ σx₁ = -i σz₁

julia> pm.phase_k
0x03

julia> pm.mono.word
1-element Vector{Int64}:
 3
```
"""
function PauliMonomial(word::Vector{T}) where {T<:Integer}
    canonical_word, phase_k = _simplify_pauli_word!(copy(word))
    mono = Monomial{PauliAlgebra,T}(canonical_word)
    PauliMonomial{T}(mono, phase_k)
end

"""
    PauliMonomial(mono::Monomial{PauliAlgebra,T}) where {T}

Construct a PauliMonomial from an existing Monomial (which may not be canonical).

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> pm = PauliMonomial(m);

julia> pm.phase_k
0x01
```
"""
function PauliMonomial(m::Monomial{PauliAlgebra,T}) where {T<:Integer}
    PauliMonomial(copy(m.word))
end

"""
    PauliMonomial{T}(mono::Monomial{PauliAlgebra,T}, phase_k::Integer) where {T}

Construct from a canonical monomial and phase. Used internally.
"""
function PauliMonomial(mono::Monomial{PauliAlgebra,T}, phase_k::Integer) where {T<:Integer}
    PauliMonomial{T}(mono, UInt8(phase_k % 4))
end

# =============================================================================
# Basic Properties
# =============================================================================

"""
    degree(pm::PauliMonomial) -> Int

Return the degree (number of operators) of the canonical monomial.
"""
degree(pm::PauliMonomial) = degree(pm.mono)

"""
    Base.isone(pm::PauliMonomial) -> Bool

Check if the PauliMonomial represents the identity (empty word, phase=1).
"""
Base.isone(pm::PauliMonomial) = isone(pm.mono) && pm.phase_k == 0

"""
    Base.one(::Type{PauliMonomial{T}}) where {T}

Return the identity PauliMonomial (empty word, phase=1).
"""
Base.one(::Type{PauliMonomial{T}}) where {T<:Integer} =
    PauliMonomial{T}(one(Monomial{PauliAlgebra,T}), UInt8(0))

Base.one(::PauliMonomial{T}) where {T} = one(PauliMonomial{T})

# =============================================================================
# Equality, Hash, Ordering
# =============================================================================

"""
    Base.:(==)(pm1::PauliMonomial, pm2::PauliMonomial) -> Bool

Equality includes both the canonical monomial AND the phase.
"""
function Base.:(==)(pm1::PauliMonomial{T1}, pm2::PauliMonomial{T2}) where {T1,T2}
    pm1.mono == pm2.mono && pm1.phase_k == pm2.phase_k
end

"""
    Base.hash(pm::PauliMonomial, h::UInt) -> UInt

Hash includes both the monomial and phase_k.
"""
Base.hash(pm::PauliMonomial, h::UInt) = hash(pm.phase_k, hash(pm.mono, h))

"""
    Base.isless(pm1::PauliMonomial{T}, pm2::PauliMonomial{T}) -> Bool

Ordering: first by monomial (degree-lexicographic), then by phase_k.
"""
function Base.isless(pm1::PauliMonomial{T}, pm2::PauliMonomial{T}) where {T}
    pm1.mono != pm2.mono && return isless(pm1.mono, pm2.mono)
    return pm1.phase_k < pm2.phase_k
end

# =============================================================================
# Multiplication
# =============================================================================

"""
    Base.:*(pm1::PauliMonomial{T}, pm2::PauliMonomial{T}) -> PauliMonomial{T}

Multiply two PauliMonomials. The result is a PauliMonomial (closed under multiplication).

Multiplies the underlying monomials and adds the phases modulo 4.

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm1 = PauliMonomial([1]);  # σx₁

julia> pm2 = PauliMonomial([2]);  # σy₁

julia> pm3 = pm1 * pm2;           # σx₁ σy₁ = i σz₁

julia> pm3.phase_k
0x01

julia> pm3.mono.word
1-element Vector{Int64}:
 3
```
"""
function Base.:*(pm1::PauliMonomial{T}, pm2::PauliMonomial{T}) where {T<:Integer}
    # Concatenate words and re-canonicalize
    combined_word = vcat(pm1.mono.word, pm2.mono.word)
    canonical_word, new_phase_k = _simplify_pauli_word!(combined_word)

    # Total phase: pm1.phase_k + pm2.phase_k + new_phase_k (mod 4)
    total_phase_k = (pm1.phase_k + pm2.phase_k + new_phase_k) % UInt8(4)

    mono = Monomial{PauliAlgebra,T}(canonical_word)
    PauliMonomial{T}(mono, total_phase_k)
end

"""
    Base.:^(pm::PauliMonomial{T}, n::Integer) -> PauliMonomial{T}

Raise a PauliMonomial to an integer power.
"""
function Base.:^(pm::PauliMonomial{T}, n::Integer) where {T<:Integer}
    n < 0 && throw(DomainError(n, "PauliMonomial exponent must be non-negative"))
    n == 0 && return one(PauliMonomial{T})
    n == 1 && return pm

    result = pm
    for _ in 2:n
        result = result * pm
    end
    return result
end

# =============================================================================
# Adjoint
# =============================================================================

"""
    Base.adjoint(pm::PauliMonomial{T}) -> PauliMonomial{T}

Compute the Hermitian adjoint of a PauliMonomial.

For Pauli operators (self-adjoint), this reverses the word and conjugates the phase.
Phase conjugation: (im)^k → (im)^(-k) = (im)^(4-k) for k>0.

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm = PauliMonomial([1, 2]);  # i σz₁

julia> pm_adj = adjoint(pm);        # -i σz₁

julia> pm_adj.phase_k
0x03
```
"""
function Base.adjoint(pm::PauliMonomial{T}) where {T<:Integer}
    # Pauli operators are self-adjoint, so adjoint just reverses the word
    # But we need to re-canonicalize since reversal changes order
    adj_word = reverse(pm.mono.word)
    canonical_word, adj_phase_k = _simplify_pauli_word!(adj_word)

    # Conjugate the original phase: conj((im)^k) = (-im)^k = (im)^(4-k) for k>0
    # More simply: conj(i) = -i, conj(-1) = -1, conj(-i) = i
    # So: phase_k -> (4 - phase_k) % 4
    conj_phase_k = pm.phase_k == 0 ? UInt8(0) : UInt8((4 - pm.phase_k) % 4)

    # Total phase from adjoint operation
    total_phase_k = (conj_phase_k + adj_phase_k) % UInt8(4)

    mono = Monomial{PauliAlgebra,T}(canonical_word)
    PauliMonomial{T}(mono, total_phase_k)
end

# =============================================================================
# Symmetric Canonicalization
# =============================================================================

"""
    symmetric_canon(pm::PauliMonomial{T}) -> PauliMonomial{T}

Return the symmetric canonical form for moment matrix keys.

For PauliMonomials, this returns min(pm, adjoint(pm)) under the isless ordering.
Since Pauli operators are self-adjoint, the monomial part is unchanged, but
the phase may be conjugated.

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm = PauliMonomial([1, 2]);  # i σz₁

julia> pm_sym = symmetric_canon(pm);

julia> pm_sym.phase_k  # min(1, 3) = 1
0x01
```
"""
function symmetric_canon(pm::PauliMonomial{T}) where {T<:Integer}
    pm_adj = adjoint(pm)
    return isless(pm, pm_adj) ? pm : pm_adj
end

# =============================================================================
# Display
# =============================================================================

function Base.show(io::IO, pm::PauliMonomial{T}) where {T}
    # Show phase prefix
    phase_str = if pm.phase_k == 0
        ""
    elseif pm.phase_k == 1
        "i"
    elseif pm.phase_k == 2
        "-"
    else  # pm.phase_k == 3
        "-i"
    end

    if isone(pm.mono)
        # Just the phase (or "1" if identity)
        print(io, phase_str == "" ? "𝟙" : phase_str)
    else
        print(io, phase_str)
        show(io, pm.mono)
    end
end

# =============================================================================
# Iteration Protocol (for unified simplify result processing)
# =============================================================================

"""
    Base.iterate(pm::PauliMonomial)

Iterate a PauliMonomial as a single-term result, yielding one `(coefficient, monomial)` pair.
The coefficient is the complex phase `(im)^phase_k`.
"""
function Base.iterate(pm::PauliMonomial{T}) where {T}
    coef = _phase_k_to_complex(pm.phase_k)
    return ((coef, pm.mono), nothing)
end

Base.iterate(::PauliMonomial, ::Nothing) = nothing
Base.IteratorSize(::Type{<:PauliMonomial}) = Base.HasLength()
Base.length(::PauliMonomial) = 1
Base.eltype(::Type{PauliMonomial{T}}) where {T} = Tuple{ComplexF64,Monomial{PauliAlgebra,T}}

# =============================================================================
# Coefficient Type
# =============================================================================

coeff_type(::Type{PauliMonomial{T}}) where {T} = ComplexF64
coeff_type(::PauliMonomial{T}) where {T} = ComplexF64

# =============================================================================
# Conversion to Complex
# =============================================================================

"""
    phase_to_complex(pm::PauliMonomial) -> ComplexF64

Return the complex phase value `(im)^phase_k`.
"""
phase_to_complex(pm::PauliMonomial) = _phase_k_to_complex(pm.phase_k)

# =============================================================================
# Multiplication Dispatch: Monomial{PauliAlgebra} -> PauliMonomial
# =============================================================================

"""
    Base.:*(m1::Monomial{PauliAlgebra,T}, m2::Monomial{PauliAlgebra,T}) -> PauliMonomial{T}

Multiply two Pauli monomials. Returns a PauliMonomial with the canonicalized
result and accumulated phase.

This overrides the generic Monomial multiplication to return the proper
wrapper type that tracks phase information.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{PauliAlgebra}([1]);  # σx₁ (assuming canonical)

julia> m2 = Monomial{PauliAlgebra}([2]);  # σy₁

julia> result = m1 * m2;

julia> result isa PauliMonomial
true

julia> result.phase_k  # i from σx σy = i σz
0x01
```
"""
function Base.:*(m1::Monomial{PauliAlgebra,T}, m2::Monomial{PauliAlgebra,T}) where {T<:Integer}
    combined_word = vcat(m1.word, m2.word)
    PauliMonomial(combined_word)
end

# =============================================================================
# Polynomial Conversion
# =============================================================================

"""
    Polynomial(pm::PauliMonomial{T}) where {T}

Convert a PauliMonomial to a Polynomial{PauliAlgebra,T,ComplexF64}.

The phase is converted to a complex coefficient.

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm = PauliMonomial([1, 2]);  # σx₁ σy₁ = i σz₁

julia> p = Polynomial(pm);

julia> coefficients(p)[1]
0.0 + 1.0im
```
"""
function Polynomial(pm::PauliMonomial{T}) where {T<:Integer}
    coef = _phase_k_to_complex(pm.phase_k)
    Polynomial(Term(coef, pm.mono))
end
