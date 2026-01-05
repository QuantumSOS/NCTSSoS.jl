# =============================================================================
# PhysicsMonomial - Linear combination of canonical fermionic/bosonic monomials
# =============================================================================

"""
    PhysicsMonomial{A<:Union{FermionicAlgebra,BosonicAlgebra}, T<:Integer} <: AbstractMonomial{A,T}

A linear combination of normal-ordered monomials with integer coefficients.

PhysicsMonomial represents the result of normal-ordering a fermionic or bosonic
monomial. Due to (anti)commutation relations, a single raw operator product can
expand into multiple normal-ordered terms with integer coefficients.

# Fields
- `coeffs::Vector{Int}`: Integer coefficients for each monomial
- `monos::Vector{Monomial{A,T}}`: Normal-ordered monomials (canonical form)

# Type Parameters
- `A`: Either `FermionicAlgebra` or `BosonicAlgebra`
- `T<:Integer`: Integer type for variable indices

# Invariants
- Each monomial in `monos` is in normal order (creators left, annihilators right)
- `coeffs` and `monos` have the same length
- Monomials are sorted and unique (no duplicate terms)
- Zero-coefficient terms are filtered out (unless result is zero)

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm = PhysicsMonomial{FermionicAlgebra}(Int32[1, -1]);  # a₁ a₁†

julia> length(pm.coeffs)  # Two terms: 1 - a₁† a₁
2

julia> pm.coeffs
2-element Vector{Int64}:
  1
 -1

julia> [m.word for m in pm.monos]
2-element Vector{Vector{Int32}}:
 Int32[]
 Int32[-1, 1]
```

See also: [`Monomial`](@ref), [`_simplify_fermionic_word!`](@ref), [`_simplify_bosonic_word!`](@ref)
"""
struct PhysicsMonomial{A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer} <: AbstractMonomial{A,T}
    coeffs::Vector{Int}
    monos::Vector{Monomial{A,T}}

    # Inner constructor from already-canonical data
    function PhysicsMonomial{A,T}(
        coeffs::Vector{Int},
        monos::Vector{Monomial{A,T}}
    ) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
        @assert length(coeffs) == length(monos) "coeffs and monos must have same length"
        new{A,T}(coeffs, monos)
    end
end

# =============================================================================
# Constructors
# =============================================================================

"""
    PhysicsMonomial{A}(word::Vector{T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra}, T<:Integer}

Construct a PhysicsMonomial from a raw (non-normal-ordered) word.

The word is automatically normal-ordered using algebra-specific simplification,
which expands into a sum of integer-coefficient terms.

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm = PhysicsMonomial{FermionicAlgebra}(Int32[1, -1]);  # a₁ a₁†

julia> length(pm.coeffs)
2
```
"""
function PhysicsMonomial{A}(word::Vector{T}) where {A<:FermionicAlgebra,T<:Integer}
    pairs = _simplify_fermionic_word!(copy(word))
    _build_physics_monomial(A, T, pairs)
end

function PhysicsMonomial{A}(word::Vector{T}) where {A<:BosonicAlgebra,T<:Integer}
    pairs = _simplify_bosonic_word!(copy(word))
    _build_physics_monomial(A, T, pairs)
end

# Helper to build from (coef, word) pairs
function _build_physics_monomial(
    ::Type{A}, ::Type{T},
    pairs::Vector{Tuple{Int,Vector{T}}}
) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    # Filter zeros and sort for canonical ordering
    pairs = filter(p -> p[1] != 0, pairs)

    if isempty(pairs)
        # Zero result - return single zero term with identity monomial
        return PhysicsMonomial{A,T}([0], [Monomial{A}(T[])])
    end

    # Sort by monomial for canonical form
    sort!(pairs, by=p -> Monomial{A}(p[2]))

    coeffs = [p[1] for p in pairs]
    monos = [Monomial{A}(p[2]) for p in pairs]

    PhysicsMonomial{A,T}(coeffs, monos)
end

"""
    PhysicsMonomial(m::Monomial{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra}, T}

Construct a PhysicsMonomial from an existing Monomial (which may not be normal-ordered).

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = Monomial{FermionicAlgebra}(Int32[1, -1]);

julia> pm = PhysicsMonomial(m);

julia> length(pm.coeffs)
2
```
"""
function PhysicsMonomial(m::Monomial{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    PhysicsMonomial{A}(copy(m.word))
end

# =============================================================================
# Basic Properties
# =============================================================================

"""
    degree(pm::PhysicsMonomial) -> Int

Return the maximum degree among all constituent monomials.
"""
function degree(pm::PhysicsMonomial)
    isempty(pm.monos) && return 0
    maximum(degree, pm.monos)
end

"""
    Base.isone(pm::PhysicsMonomial) -> Bool

Check if the PhysicsMonomial represents the identity (single term: 1 * identity).
"""
function Base.isone(pm::PhysicsMonomial)
    length(pm.coeffs) == 1 && pm.coeffs[1] == 1 && isone(pm.monos[1])
end

"""
    Base.iszero(pm::PhysicsMonomial) -> Bool

Check if the PhysicsMonomial is zero.
"""
function Base.iszero(pm::PhysicsMonomial)
    all(==(0), pm.coeffs)
end

"""
    Base.one(::Type{PhysicsMonomial{A,T}}) where {A,T}

Return the identity PhysicsMonomial (single term: 1 * identity monomial).
"""
function Base.one(::Type{PhysicsMonomial{A,T}}) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    PhysicsMonomial{A,T}([1], [one(Monomial{A,T})])
end

Base.one(::PhysicsMonomial{A,T}) where {A,T} = one(PhysicsMonomial{A,T})

"""
    Base.zero(::Type{PhysicsMonomial{A,T}}) where {A,T}

Return the zero PhysicsMonomial.
"""
function Base.zero(::Type{PhysicsMonomial{A,T}}) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    PhysicsMonomial{A,T}([0], [one(Monomial{A,T})])
end

Base.zero(::PhysicsMonomial{A,T}) where {A,T} = zero(PhysicsMonomial{A,T})

# =============================================================================
# Equality, Hash, Ordering
# =============================================================================

"""
    Base.:(==)(pm1::PhysicsMonomial{A,T}, pm2::PhysicsMonomial{A,T}) -> Bool

Equality check: same coefficients and monomials.
"""
function Base.:(==)(pm1::PhysicsMonomial{A1,T1}, pm2::PhysicsMonomial{A2,T2}) where {A1,A2,T1,T2}
    A1 !== A2 && return false
    pm1.coeffs == pm2.coeffs && pm1.monos == pm2.monos
end

"""
    Base.hash(pm::PhysicsMonomial, h::UInt) -> UInt

Hash based on coefficients and monomials.
"""
function Base.hash(pm::PhysicsMonomial{A}, h::UInt) where {A}
    h = hash(A, h)
    h = hash(pm.coeffs, h)
    hash(pm.monos, h)
end

"""
    Base.isless(pm1::PhysicsMonomial{A,T}, pm2::PhysicsMonomial{A,T}) -> Bool

Ordering: first by degree, then by leading monomial, then by leading coefficient.
"""
function Base.isless(pm1::PhysicsMonomial{A,T}, pm2::PhysicsMonomial{A,T}) where {A,T}
    d1, d2 = degree(pm1), degree(pm2)
    d1 != d2 && return d1 < d2

    # Compare leading (first) monomial
    !isempty(pm1.monos) && !isempty(pm2.monos) || return isempty(pm1.monos)
    pm1.monos[1] != pm2.monos[1] && return isless(pm1.monos[1], pm2.monos[1])

    # Compare leading coefficient
    pm1.coeffs[1] != pm2.coeffs[1] && return pm1.coeffs[1] < pm2.coeffs[1]

    # Compare remaining terms
    for i in 2:min(length(pm1.monos), length(pm2.monos))
        pm1.monos[i] != pm2.monos[i] && return isless(pm1.monos[i], pm2.monos[i])
        pm1.coeffs[i] != pm2.coeffs[i] && return pm1.coeffs[i] < pm2.coeffs[i]
    end

    return length(pm1.monos) < length(pm2.monos)
end

# =============================================================================
# Multiplication
# =============================================================================

"""
    Base.:*(pm1::PhysicsMonomial{A,T}, pm2::PhysicsMonomial{A,T}) -> PhysicsMonomial{A,T}

Multiply two PhysicsMonomials by distributing and re-normal-ordering.

# Examples
```jldoctest
julia> using NCTSSoS

julia> pm1 = PhysicsMonomial{FermionicAlgebra}(Int32[-1]);  # a₁†

julia> pm2 = PhysicsMonomial{FermionicAlgebra}(Int32[1]);   # a₁

julia> pm3 = pm1 * pm2;  # a₁† a₁ (already normal)

julia> pm3.coeffs
1-element Vector{Int64}:
 1
```
"""
function Base.:*(pm1::PhysicsMonomial{A,T}, pm2::PhysicsMonomial{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    # Distribute: (Σ c1_i m1_i) * (Σ c2_j m2_j) = Σ c1_i*c2_j * (m1_i * m2_j)
    result_pairs = Tuple{Int,Vector{T}}[]

    for (c1, m1) in zip(pm1.coeffs, pm1.monos)
        c1 == 0 && continue
        for (c2, m2) in zip(pm2.coeffs, pm2.monos)
            c2 == 0 && continue
            combined_word = vcat(m1.word, m2.word)

            # Re-simplify the product
            if A <: FermionicAlgebra
                pairs = _simplify_fermionic_word!(combined_word)
            else  # BosonicAlgebra
                pairs = _simplify_bosonic_word!(combined_word)
            end

            # Scale by c1*c2
            coef_prod = c1 * c2
            for (c, w) in pairs
                push!(result_pairs, (coef_prod * c, w))
            end
        end
    end

    # Combine like terms
    grouped = Dict{Vector{T},Int}()
    for (c, w) in result_pairs
        grouped[w] = get(grouped, w, 0) + c
    end

    final_pairs = [(c, w) for (w, c) in grouped if c != 0]

    if isempty(final_pairs)
        return zero(PhysicsMonomial{A,T})
    end

    sort!(final_pairs, by=p -> Monomial{A}(p[2]))
    coeffs = [p[1] for p in final_pairs]
    monos = [Monomial{A}(p[2]) for p in final_pairs]

    PhysicsMonomial{A,T}(coeffs, monos)
end

"""
    Base.:^(pm::PhysicsMonomial{A,T}, n::Integer) -> PhysicsMonomial{A,T}

Raise a PhysicsMonomial to an integer power.
"""
function Base.:^(pm::PhysicsMonomial{A,T}, n::Integer) where {A,T}
    n < 0 && throw(DomainError(n, "PhysicsMonomial exponent must be non-negative"))
    n == 0 && return one(PhysicsMonomial{A,T})
    n == 1 && return pm

    result = pm
    for _ in 2:n
        result = result * pm
    end
    return result
end

# =============================================================================
# Scalar Multiplication
# =============================================================================

"""
    Base.:*(c::Integer, pm::PhysicsMonomial{A,T}) -> PhysicsMonomial{A,T}

Scale a PhysicsMonomial by an integer.
"""
function Base.:*(c::Integer, pm::PhysicsMonomial{A,T}) where {A,T}
    c == 0 && return zero(PhysicsMonomial{A,T})
    c == 1 && return pm
    PhysicsMonomial{A,T}(c .* pm.coeffs, copy(pm.monos))
end

Base.:*(pm::PhysicsMonomial, c::Integer) = c * pm

# =============================================================================
# Adjoint
# =============================================================================

"""
    Base.adjoint(pm::PhysicsMonomial{A,T}) -> PhysicsMonomial{A,T}

Compute the Hermitian adjoint of a PhysicsMonomial.

For each term c_i * m_i, computes c_i* * m_i† where coefficients are real (integers)
and monomial adjoint reverses word and negates indices (a† ↔ a).
"""
function Base.adjoint(pm::PhysicsMonomial{A,T}) where {A<:Union{FermionicAlgebra,BosonicAlgebra},T<:Integer}
    result_pairs = Tuple{Int,Vector{T}}[]

    for (c, m) in zip(pm.coeffs, pm.monos)
        c == 0 && continue
        # Adjoint of monomial: reverse and negate
        adj_word = T[-idx for idx in reverse(m.word)]

        # Re-simplify to normal order
        if A <: FermionicAlgebra
            pairs = _simplify_fermionic_word!(adj_word)
        else
            pairs = _simplify_bosonic_word!(adj_word)
        end

        # Integer coefficients are real, so conjugate is identity
        for (c_new, w) in pairs
            push!(result_pairs, (c * c_new, w))
        end
    end

    # Combine like terms
    grouped = Dict{Vector{T},Int}()
    for (c, w) in result_pairs
        grouped[w] = get(grouped, w, 0) + c
    end

    final_pairs = [(c, w) for (w, c) in grouped if c != 0]

    if isempty(final_pairs)
        return zero(PhysicsMonomial{A,T})
    end

    sort!(final_pairs, by=p -> Monomial{A}(p[2]))
    coeffs = [p[1] for p in final_pairs]
    monos = [Monomial{A}(p[2]) for p in final_pairs]

    PhysicsMonomial{A,T}(coeffs, monos)
end

# =============================================================================
# Display
# =============================================================================

function Base.show(io::IO, pm::PhysicsMonomial{A,T}) where {A,T}
    if iszero(pm)
        print(io, "0")
        return
    end

    first_term = true
    for (c, m) in zip(pm.coeffs, pm.monos)
        c == 0 && continue

        if !first_term
            if c > 0
                print(io, " + ")
            else
                print(io, " - ")
                c = -c
            end
        elseif c < 0
            print(io, "-")
            c = -c
        end
        first_term = false

        if isone(m)
            print(io, c)
        elseif c == 1
            show(io, m)
        else
            print(io, c, "*")
            show(io, m)
        end
    end
end

# =============================================================================
# Iteration Protocol (for unified simplify result processing)
# =============================================================================

"""
    Base.iterate(pm::PhysicsMonomial)

Iterate a PhysicsMonomial, yielding `(coefficient, monomial)` pairs.
Coefficients are converted to Float64 for consistency with Polynomial.
"""
function Base.iterate(pm::PhysicsMonomial{A,T}) where {A,T}
    isempty(pm.coeffs) && return nothing
    return ((Float64(pm.coeffs[1]), pm.monos[1]), 2)
end

function Base.iterate(pm::PhysicsMonomial{A,T}, state::Int) where {A,T}
    state > length(pm.coeffs) && return nothing
    return ((Float64(pm.coeffs[state]), pm.monos[state]), state + 1)
end

Base.IteratorSize(::Type{<:PhysicsMonomial}) = Base.HasLength()
Base.length(pm::PhysicsMonomial) = length(pm.coeffs)
Base.eltype(::Type{PhysicsMonomial{A,T}}) where {A,T} = Tuple{Float64,Monomial{A,T}}

# =============================================================================
# Coefficient Type
# =============================================================================

coeff_type(::Type{PhysicsMonomial{A,T}}) where {A,T} = Float64
coeff_type(::PhysicsMonomial{A,T}) where {A,T} = Float64

# =============================================================================
# Multiplication Dispatch: Monomial{Fermionic/Bosonic} -> PhysicsMonomial
# =============================================================================

"""
    Base.:*(m1::Monomial{FermionicAlgebra,T}, m2::Monomial{FermionicAlgebra,T}) -> PhysicsMonomial{FermionicAlgebra,T}

Multiply two fermionic monomials. Returns a PhysicsMonomial with the
normal-ordered expansion.

This overrides the generic Monomial multiplication to return the proper
wrapper type that handles anticommutation relations.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{FermionicAlgebra}(Int32[-1]);  # a₁†

julia> m2 = Monomial{FermionicAlgebra}(Int32[1]);   # a₁

julia> result = m1 * m2;

julia> result isa PhysicsMonomial
true

julia> result.coeffs  # a₁† a₁ (single term, already normal)
1-element Vector{Int64}:
 1
```
"""
function Base.:*(m1::Monomial{FermionicAlgebra,T}, m2::Monomial{FermionicAlgebra,T}) where {T<:Integer}
    combined_word = vcat(m1.word, m2.word)
    PhysicsMonomial{FermionicAlgebra}(combined_word)
end

"""
    Base.:*(m1::Monomial{BosonicAlgebra,T}, m2::Monomial{BosonicAlgebra,T}) -> PhysicsMonomial{BosonicAlgebra,T}

Multiply two bosonic monomials. Returns a PhysicsMonomial with the
normal-ordered expansion.

This overrides the generic Monomial multiplication to return the proper
wrapper type that handles commutation relations.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = Monomial{BosonicAlgebra}(Int32[1]);   # c₁

julia> m2 = Monomial{BosonicAlgebra}(Int32[-1]);  # c₁†

julia> result = m1 * m2;

julia> result isa PhysicsMonomial
true

julia> length(result.coeffs)  # c₁ c₁† = c₁† c₁ + 1 (two terms)
2
```
"""
function Base.:*(m1::Monomial{BosonicAlgebra,T}, m2::Monomial{BosonicAlgebra,T}) where {T<:Integer}
    combined_word = vcat(m1.word, m2.word)
    PhysicsMonomial{BosonicAlgebra}(combined_word)
end
