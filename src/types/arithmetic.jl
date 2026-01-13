# =============================================================================
# Arithmetic Operations for NormalMonomial and Polynomial
# =============================================================================
#
# This file provides arithmetic operations where:
# - All NormalMonomial operations return Polynomial
# - Polynomial operations return Polynomial
#
# This enables convenient polynomial construction for optimization problems.
# =============================================================================

# =============================================================================
# Helper: Convert internal coefficient encoding to numeric type
# =============================================================================

"""
    _coeff_to_number(::Type{A}, c) where {A<:AlgebraType}
    _coeff_to_number(m::NormalMonomial{A}, c) where {A<:AlgebraType}

Convert internal coefficient encoding to numeric type.

Different algebras use different internal coefficient representations:
- MonoidAlgebra: UInt8 (0x01 = one) → Float64(1.0)
- TwistedGroupAlgebra: UInt8 phase_k where coeff = (im)^k → ComplexF64
- PBWAlgebra: Int → Float64

The monomial-based dispatch extracts the algebra type from the monomial.
"""
function _coeff_to_number(::Type{A}, c::UInt8) where {A<:MonoidAlgebra}
    return one(Float64)
end

function _coeff_to_number(::Type{A}, phase_k::UInt8) where {A<:TwistedGroupAlgebra}
    phase_k == 0x00 && return ComplexF64(1.0, 0.0)
    phase_k == 0x01 && return ComplexF64(0.0, 1.0)
    phase_k == 0x02 && return ComplexF64(-1.0, 0.0)
    phase_k == 0x03 && return ComplexF64(0.0, -1.0)
    # phase_k == 0x04 is zero, shouldn't reach here
    return ComplexF64(0.0, 0.0)
end

function _coeff_to_number(::Type{A}, c::Integer) where {A<:PBWAlgebra}
    return Float64(c)
end

_coeff_to_number(::Type{<:AlgebraType}, c::Number) = c

# Dispatch via monomial (extracts algebra type from type parameter)
_coeff_to_number(::NormalMonomial{A}, c) where {A<:AlgebraType} = _coeff_to_number(A, c)

# =============================================================================
# Helper: Convert simplified result to Polynomial terms
# =============================================================================

"""
    _simplified_to_terms(::Type{A}, simplified, ::Type{T}) where {A,T}

Convert the result of `simplify(A, word)` to Polynomial terms.

Different algebras return different types from simplify:
- MonoidAlgebra: Vector{T} (simplified word)
- TwistedGroupAlgebra: (Vector{T}, UInt8) (word, phase)
- PBWAlgebra: Vector{Tuple{Int, Vector{T}}} (list of (coeff, word) pairs)
"""
function _simplified_to_terms(
    ::Type{A}, simplified::Vector{T}, ::Type{T}
) where {A<:MonoidAlgebra,T<:Integer}
    C = coeff_type(A)
    mono = NormalMonomial{A,T}(simplified)
    return [(one(C), mono)]
end

function _simplified_to_terms(
    ::Type{A}, simplified::Tuple{Vector{T},UInt8}, ::Type{T}
) where {A<:TwistedGroupAlgebra,T<:Integer}
    word, phase_k = simplified
    C = coeff_type(A)
    # phase_k == _TG_ZERO means zero result (shouldn't happen for valid inputs)
    phase_k == 0x04 && return Tuple{C,NormalMonomial{A,T}}[]
    coef = _coeff_to_number(A, phase_k)
    mono = NormalMonomial{A,T}(word)
    return [(coef, mono)]
end

function _simplified_to_terms(
    ::Type{A}, simplified::Vector{Tuple{Int,Vector{T}}}, ::Type{T}
) where {A<:PBWAlgebra,T<:Integer}
    C = coeff_type(A)
    terms = Tuple{C,NormalMonomial{A,T}}[]
    for (c, word) in simplified
        iszero(c) && continue
        mono = NormalMonomial{A,T}(word)
        push!(terms, (C(c), mono))
    end
    return terms
end

# =============================================================================
# NormalMonomial × NormalMonomial → Polynomial
# =============================================================================

"""
    Base.:(*)(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A,T} -> Polynomial

Multiply two normal monomials, returning a Polynomial.

Concatenates the words, applies algebra-specific simplification, and converts
the result to a Polynomial with appropriate coefficients.

# Examples
```julia
julia> reg, (σx, σy, σz) = create_pauli_variables(1);
julia> σx[1] * σx[1]  # σx² = 1
Polynomial with 1 term
```
"""
function Base.:(*)(
    m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    C = coeff_type(A)
    # Handle identity cases
    isempty(m1.word) && return Polynomial([(one(C), m2)])
    isempty(m2.word) && return Polynomial([(one(C), m1)])

    # Concatenate and simplify
    combined = vcat(m1.word, m2.word)
    simplified = simplify(A, combined)
    terms = _simplified_to_terms(A, simplified, T)

    isempty(terms) && return zero(Polynomial{A,T,C})
    return Polynomial(terms)
end

# =============================================================================
# NormalMonomial Power → Polynomial
# =============================================================================

"""
    Base.:(^)(m::NormalMonomial{A,T}, n::Integer) where {A,T} -> Polynomial

Raise a NormalMonomial to a non-negative integer power, returning a Polynomial.

Uses binary exponentiation (power by squaring) for efficiency.
"""
function Base.:(^)(m::NormalMonomial{A,T}, n::Integer) where {A<:AlgebraType,T<:Integer}
    n < 0 && throw(DomainError(n, "monomial exponent must be non-negative"))
    C = coeff_type(A)

    # Base cases
    n == 0 && return one(Polynomial{A,T,C})
    n == 1 && return Polynomial([(one(C), m)])

    # Repeat word n times, then simplify
    repeated_word = repeat(m.word, n)
    simplified = simplify(A, repeated_word)
    terms = _simplified_to_terms(A, simplified, T)

    isempty(terms) && return zero(Polynomial{A,T,C})
    return Polynomial(terms)
end

# =============================================================================
# NormalMonomial Negation → Polynomial
# =============================================================================

"""
    Base.:(-)(m::NormalMonomial{A,T}) where {A,T} -> Polynomial

Negate a NormalMonomial, returning a Polynomial with coefficient -1.
"""
function Base.:(-)(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    C = coeff_type(A)
    return Polynomial([(-one(C), m)])
end

# =============================================================================
# NormalMonomial + NormalMonomial → Polynomial
# =============================================================================

"""
    Base.:(+)(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A,T} -> Polynomial

Add two NormalMonomials, returning a Polynomial.
"""
function Base.:(+)(
    m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    C = coeff_type(A)
    return Polynomial([(one(C), m1), (one(C), m2)])
end

"""
    Base.:(-)(m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}) where {A,T} -> Polynomial

Subtract two NormalMonomials, returning a Polynomial.
"""
function Base.:(-)(
    m1::NormalMonomial{A,T}, m2::NormalMonomial{A,T}
) where {A<:AlgebraType,T<:Integer}
    C = coeff_type(A)
    return Polynomial([(one(C), m1), (-one(C), m2)])
end

# =============================================================================
# NormalMonomial ± Number → Polynomial
# =============================================================================

"""
    Base.:(+)(m::NormalMonomial{A,T}, c::Number) where {A,T} -> Polynomial

Add a scalar to a NormalMonomial, returning a Polynomial.
"""
function Base.:(+)(m::NormalMonomial{A,T}, c::Number) where {A<:AlgebraType,T<:Integer}
    C = promote_type(typeof(c), coeff_type(A))
    iszero(c) && return Polynomial([(one(C), m)])
    identity = one(NormalMonomial{A,T})
    return Polynomial([(one(C), m), (C(c), identity)])
end

Base.:(+)(c::Number, m::NormalMonomial) = m + c

"""
    Base.:(-)(m::NormalMonomial{A,T}, c::Number) where {A,T} -> Polynomial

Subtract a scalar from a NormalMonomial, returning a Polynomial.
"""
function Base.:(-)(m::NormalMonomial{A,T}, c::Number) where {A<:AlgebraType,T<:Integer}
    return m + (-c)
end

"""
    Base.:(-)(c::Number, m::NormalMonomial{A,T}) where {A,T} -> Polynomial

Subtract a NormalMonomial from a scalar, returning a Polynomial.
"""
function Base.:(-)(c::Number, m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    C = promote_type(typeof(c), coeff_type(A))
    iszero(c) && return Polynomial([(-one(C), m)])
    identity = one(NormalMonomial{A,T})
    return Polynomial([(C(c), identity), (-one(C), m)])
end

# =============================================================================
# Number × NormalMonomial → Polynomial
# =============================================================================

"""
    Base.:(*)(c::Number, m::NormalMonomial{A,T}) where {A,T} -> Polynomial

Multiply a scalar by a NormalMonomial, returning a Polynomial.
"""
function Base.:(*)(c::Number, m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    C = promote_type(typeof(c), coeff_type(A))
    iszero(c) && return zero(Polynomial{A,T,C})
    return Polynomial([(C(c), m)])
end

Base.:(*)(m::NormalMonomial, c::Number) = c * m

# =============================================================================
# NormalMonomial × Polynomial → Polynomial
# =============================================================================

"""
    Base.:(*)(m::NormalMonomial{A,T}, p::Polynomial{A,T,C}) where {A,T,C} -> Polynomial

Multiply a NormalMonomial by a Polynomial from the left.
"""
function Base.:(*)(
    m::NormalMonomial{A,T}, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(p.terms) && return zero(Polynomial{A,T,C})
    isempty(m.word) && return p  # m is identity

    # Delegate to NormalMonomial × NormalMonomial → Polynomial, then sum
    DC = coeff_type(A)
    NC = promote_type(C, DC)
    result = zero(Polynomial{A,T,NC})

    for (coef, mono) in p.terms
        iszero(coef) && continue
        # m * mono returns a Polynomial
        prod_poly = m * mono
        result = result + coef * prod_poly
    end

    return result
end

"""
    Base.:(*)(p::Polynomial{A,T,C}, m::NormalMonomial{A,T}) where {A,T,C} -> Polynomial

Multiply a Polynomial by a NormalMonomial from the right.
"""
function Base.:(*)(
    p::Polynomial{A,T,C}, m::NormalMonomial{A,T}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(p.terms) && return zero(Polynomial{A,T,C})
    isempty(m.word) && return p  # m is identity

    # Delegate to NormalMonomial × NormalMonomial → Polynomial, then sum
    DC = coeff_type(A)
    NC = promote_type(C, DC)
    result = zero(Polynomial{A,T,NC})

    for (coef, mono) in p.terms
        iszero(coef) && continue
        # mono * m returns a Polynomial
        prod_poly = mono * m
        result = result + coef * prod_poly
    end

    return result
end

# =============================================================================
# NormalMonomial ± Polynomial → Polynomial
# =============================================================================

"""
    Base.:(+)(m::NormalMonomial{A,T}, p::Polynomial{A,T,C}) where {A,T,C} -> Polynomial

Add a NormalMonomial to a Polynomial.
"""
function Base.:(+)(
    m::NormalMonomial{A,T}, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    return Polynomial(m) + p
end

Base.:(+)(p::Polynomial, m::NormalMonomial) = m + p

"""
    Base.:(-)(m::NormalMonomial{A,T}, p::Polynomial{A,T,C}) where {A,T,C} -> Polynomial

Subtract a Polynomial from a NormalMonomial.
"""
function Base.:(-)(
    m::NormalMonomial{A,T}, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    return Polynomial(m) + (-p)
end

"""
    Base.:(-)(p::Polynomial{A,T,C}, m::NormalMonomial{A,T}) where {A,T,C} -> Polynomial

Subtract a NormalMonomial from a Polynomial.
"""
function Base.:(-)(
    p::Polynomial{A,T,C}, m::NormalMonomial{A,T}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    return p + (-m)
end

# =============================================================================
# Polynomial × Polynomial → Polynomial
# =============================================================================

"""
    Base.:(*)(p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}) where {A,T,C1,C2} -> Polynomial

Multiply two Polynomials using the distributive law.

For each pair of terms (c1, m1) and (c2, m2):
1. Compute product coefficient c1 * c2
2. Concatenate monomials and simplify: m1 * m2
3. Combine like terms in the result
"""
function Base.:(*)(
    p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}
) where {A<:AlgebraType,T<:Integer,C1<:Number,C2<:Number}
    C = promote_type(C1, C2)
    isempty(p1.terms) && return zero(Polynomial{A,T,C})
    isempty(p2.terms) && return zero(Polynomial{A,T,C})

    DC = coeff_type(A)
    NC = promote_type(C, DC)
    result_terms = Tuple{NC,NormalMonomial{A,T}}[]

    for (coef1, mono1) in p1.terms
        iszero(coef1) && continue
        for (coef2, mono2) in p2.terms
            iszero(coef2) && continue
            base_coef = coef1 * coef2

            # Handle identity cases
            if isempty(mono1.word)
                push!(result_terms, (NC(base_coef), mono2))
                continue
            end
            if isempty(mono2.word)
                push!(result_terms, (NC(base_coef), mono1))
                continue
            end

            # Multiply monomials: concatenate and simplify
            combined = vcat(mono1.word, mono2.word)
            simplified = simplify(A, combined)
            prod_terms = _simplified_to_terms(A, simplified, T)

            for (c_prod, m_prod) in prod_terms
                total_coef = NC(base_coef * c_prod)
                iszero(total_coef) && continue
                push!(result_terms, (total_coef, m_prod))
            end
        end
    end

    isempty(result_terms) && return zero(Polynomial{A,T,NC})
    return Polynomial(result_terms)
end

# =============================================================================
# Polynomial ± Polynomial → Polynomial
# =============================================================================

"""
    Base.:(+)(p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}) where {A,T,C1,C2} -> Polynomial

Add two Polynomials.
"""
function Base.:(+)(
    p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}
) where {A<:AlgebraType,T<:Integer,C1<:Number,C2<:Number}
    C = promote_type(C1, C2)
    terms1 = Tuple{C,NormalMonomial{A,T}}[(C(c), m) for (c, m) in p1.terms]
    terms2 = Tuple{C,NormalMonomial{A,T}}[(C(c), m) for (c, m) in p2.terms]
    return Polynomial(vcat(terms1, terms2))
end

"""
    Base.:(-)(p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}) where {A,T,C1,C2} -> Polynomial

Subtract two Polynomials.
"""
function Base.:(-)(
    p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}
) where {A<:AlgebraType,T<:Integer,C1<:Number,C2<:Number}
    return p1 + (-p2)
end

"""
    Base.:(-)(p::Polynomial{A,T,C}) where {A,T,C} -> Polynomial

Negate a Polynomial (negate all coefficients).
"""
function Base.:(-)(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    negated_terms = Tuple{C,NormalMonomial{A,T}}[(-c, m) for (c, m) in p.terms]
    return Polynomial{A,T,C}(negated_terms)
end

# =============================================================================
# Polynomial ± Number → Polynomial
# =============================================================================

"""
    Base.:(+)(p::Polynomial{A,T,C}, c::Number) where {A,T,C} -> Polynomial

Add a scalar to a Polynomial (adds to constant term).
"""
function Base.:(+)(
    p::Polynomial{A,T,C}, c::Number
) where {A<:AlgebraType,T<:Integer,C<:Number}
    NC = promote_type(typeof(c), C)
    iszero(c) && return Polynomial(Tuple{NC,NormalMonomial{A,T}}[(NC(coef), m) for (coef, m) in p.terms])
    const_term = (NC(c), one(NormalMonomial{A,T}))
    const_poly = Polynomial([const_term])
    p_converted = Polynomial(Tuple{NC,NormalMonomial{A,T}}[(NC(coef), m) for (coef, m) in p.terms])
    return p_converted + const_poly
end

Base.:(+)(c::Number, p::Polynomial) = p + c

"""
    Base.:(-)(p::Polynomial{A,T,C}, c::Number) where {A,T,C} -> Polynomial

Subtract a scalar from a Polynomial.
"""
Base.:(-)(p::Polynomial, c::Number) = p + (-c)

"""
    Base.:(-)(c::Number, p::Polynomial{A,T,C}) where {A,T,C} -> Polynomial

Subtract a Polynomial from a scalar.
"""
function Base.:(-)(
    c::Number, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    NC = promote_type(typeof(c), C)
    const_term = (NC(c), one(NormalMonomial{A,T}))
    const_poly = Polynomial([const_term])
    return const_poly + (-p)
end

# =============================================================================
# Number × Polynomial → Polynomial
# =============================================================================

"""
    Base.:(*)(c::Number, p::Polynomial{A,T,C}) where {A,T,C} -> Polynomial

Multiply a Polynomial by a scalar.
"""
function Base.:(*)(
    c::Number, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    NC = promote_type(typeof(c), C)
    iszero(c) && return zero(Polynomial{A,T,NC})
    scaled_terms = Tuple{NC,NormalMonomial{A,T}}[(NC(c * coef), m) for (coef, m) in p.terms]
    return Polynomial{A,T,NC}(scaled_terms)
end

Base.:(*)(p::Polynomial, c::Number) = c * p

"""
    Base.:(/)(p::Polynomial{A,T,C}, c::Number) where {A,T,C} -> Polynomial

Divide a Polynomial by a scalar.
"""
function Base.:(/)(
    p::Polynomial{A,T,C}, c::Number
) where {A<:AlgebraType,T<:Integer,C<:Number}
    iszero(c) && throw(DivideError())
    return inv(c) * p
end

# =============================================================================
# Polynomial Power → Polynomial
# =============================================================================

"""
    Base.:(^)(p::Polynomial{A,T,C}, n::Int) where {A,T,C} -> Polynomial

Raise a Polynomial to a non-negative integer power.

Uses exponentiation by squaring for efficiency.
"""
function Base.:(^)(p::Polynomial{A,T,C}, n::Int) where {A<:AlgebraType,T<:Integer,C<:Number}
    n < 0 && throw(DomainError(n, "Polynomial exponent must be non-negative"))
    n == 0 && return one(p)
    n == 1 && return copy(p)
    # Exponentiation by squaring
    result = one(p)
    base = copy(p)
    while n > 0
        if isodd(n)
            result = result * base
        end
        base = base * base
        n >>= 1
    end
    return result
end
