"""
    Polynomial{A<:AlgebraType, T<:Integer, C<:Number}

A polynomial represented as a sum of terms with coefficients.
Maintains sorted, unique monomials with non-zero coefficients.

# Fields
- `terms::Vector{Term{Monomial{A,T}, C}}`: Sorted terms with unique monomials and non-zero coefficients

# Type Parameters
- `A<:AlgebraType`: Algebra type for dispatch (PauliAlgebra, FermionicAlgebra, etc.)
- `T<:Integer`: Integer type for monomial word representation
- `C<:Number`: Coefficient type (Float64, ComplexF64, etc.)

# Invariants
- Terms are sorted by monomial (using isless)
- No duplicate monomials (combined during construction)
- All coefficients are non-zero (zeros removed during construction)

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([3]);

julia> t1 = Term(1.0 + 0.0im, m1);

julia> t2 = Term(2.0 + 0.0im, m2);

julia> p = Polynomial([t1, t2]);

julia> length(terms(p))
2

julia> degree(p)
2
```

Construction with automatic deduplication:
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1]);

julia> p = Polynomial([Term(1.0+0im, m), Term(2.0+0im, m)]);  # Same monomial twice

julia> length(terms(p))  # Combined into one term
1

julia> coefficients(p)[1]  # Coefficients added
3.0 + 0.0im
```

See also: [`Term`](@ref), [`Monomial`](@ref), [`coefficients`](@ref), [`monomials`](@ref)
"""
struct Polynomial{A<:AlgebraType,T<:Integer,C<:Number}
    terms::Vector{Term{Monomial{A,T},C}}

    # Inner constructor enforces invariants: sorted, deduplicated, non-zero coefficients
    function Polynomial{A,T,C}(
        input_terms::Vector{Term{Monomial{A,T},C}}
    ) where {A<:AlgebraType,T<:Integer,C<:Number}
        processed = _process_terms(input_terms, C)
        return new{A,T,C}(processed)
    end
end

"""
    _process_terms(input_terms, C) -> Vector{Term}

Process a vector of terms: sort by monomial, combine duplicates, remove zeros.
This is the core algorithm for maintaining polynomial invariants.

# Algorithm
1. If empty, return empty vector
2. Sort terms by monomial (degree-first, then lexicographic)
3. Iterate through sorted terms, combining coefficients for duplicate monomials
4. Filter out terms with zero coefficients
"""
function _process_terms(
    input_terms::Vector{Term{Monomial{A,T},C}}, ::Type{C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(input_terms) && return Term{Monomial{A,T},C}[]

    # Sort by monomial
    sorted = sort(input_terms; by=t -> t.monomial)

    # Combine duplicates and filter zeros
    result = Term{Monomial{A,T},C}[]
    current_mono = sorted[1].monomial
    current_coef = sorted[1].coefficient

    for i in 2:length(sorted)
        if sorted[i].monomial == current_mono
            current_coef += sorted[i].coefficient
        else
            if !iszero(current_coef)
                push!(result, Term(current_coef, current_mono))
            end
            current_mono = sorted[i].monomial
            current_coef = sorted[i].coefficient
        end
    end

    # Don't forget the last accumulated term
    if !iszero(current_coef)
        push!(result, Term(current_coef, current_mono))
    end

    return result
end

# =============================================================================
# Constructors
# =============================================================================

"""
    Polynomial(terms::Vector{Term{Monomial{A,T}, C}}) where {A,T,C}

Construct a polynomial from a vector of terms.
Automatically sorts, combines duplicates, and removes zero coefficients.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1]);

julia> t1 = Term(2.0 + 0.0im, m1);

julia> p = Polynomial([t1]);

julia> coefficients(p)
1-element Vector{ComplexF64}:
 2.0 + 0.0im
```
"""
function Polynomial(
    input_terms::Vector{Term{Monomial{A,T},C}}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    return Polynomial{A,T,C}(input_terms)
end

"""
    Polynomial(t::Term{Monomial{A,T}, C}) where {A,T,C}

Construct a polynomial from a single term.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> t = Term(3.0 + 0.0im, m);

julia> p = Polynomial(t);

julia> length(terms(p))
1
```
"""
function Polynomial(t::Term{Monomial{A,T},C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    if iszero(t.coefficient)
        return Polynomial{A,T,C}(Term{Monomial{A,T},C}[])
    end
    return Polynomial{A,T,C}([t])
end

"""
    Polynomial(m::Monomial{A,T}) where {A,T}

Construct a polynomial from a monomial with coefficient 1.0.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> p = Polynomial(m);

julia> coefficients(p)
1-element Vector{Float64}:
 1.0
```
"""
function Polynomial(m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    return Polynomial(Term(1.0, m))
end

"""
    Polynomial{A,T,C}(c::Number) where {A,T,C}

Construct a constant polynomial (coefficient times identity monomial).

# Examples
```jldoctest
julia> using FastPolynomials

julia> p = Polynomial{PauliAlgebra,Int64,Float64}(5.0);

julia> coefficients(p)
1-element Vector{Float64}:
 5.0

julia> isone(monomials(p)[1])
true
```
"""
function Polynomial{A,T,C}(c::Number) where {A<:AlgebraType,T<:Integer,C<:Number}
    if iszero(c)
        return Polynomial{A,T,C}(Term{Monomial{A,T},C}[])
    end
    return Polynomial{A,T,C}([Term(C(c), Monomial{A}(T[]))])
end

# =============================================================================
# zero and one
# =============================================================================

"""
    Base.zero(::Type{Polynomial{A,T,C}}) where {A,T,C}

Create the zero polynomial (no terms).

# Examples
```jldoctest
julia> using FastPolynomials

julia> p = zero(Polynomial{PauliAlgebra,Int64,Float64});

julia> iszero(p)
true

julia> length(terms(p))
0
```
"""
function Base.zero(::Type{Polynomial{A,T,C}}) where {A<:AlgebraType,T<:Integer,C<:Number}
    return Polynomial{A,T,C}(Term{Monomial{A,T},C}[])
end

"""
    Base.zero(p::Polynomial{A,T,C}) where {A,T,C}

Create the zero polynomial for the same type as `p`.
"""
function Base.zero(::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    return zero(Polynomial{A,T,C})
end

"""
    Base.one(::Type{Polynomial{A,T,C}}) where {A,T,C}

Create the multiplicative identity polynomial (1 times identity monomial).

# Examples
```jldoctest
julia> using FastPolynomials

julia> p = one(Polynomial{PauliAlgebra,Int64,Float64});

julia> isone(p)
true

julia> coefficients(p)
1-element Vector{Float64}:
 1.0
```
"""
function Base.one(::Type{Polynomial{A,T,C}}) where {A<:AlgebraType,T<:Integer,C<:Number}
    return Polynomial{A,T,C}([Term(one(C), Monomial{A}(T[]))])
end

"""
    Base.one(p::Polynomial{A,T,C}) where {A,T,C}

Create the multiplicative identity polynomial for the same type as `p`.
"""
function Base.one(::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    return one(Polynomial{A,T,C})
end

"""
    Base.copy(p::Polynomial{A,T,C}) where {A,T,C}

Create a shallow copy of the polynomial.
"""
function Base.copy(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Create a new Polynomial with copied terms
    # Since terms are mutable, we need to create new Term objects
    copied_terms = [Term(t.coefficient, t.monomial) for t in p.terms]
    return Polynomial{A,T,C}(copied_terms)
end

"""
    Base.iszero(p::Polynomial) -> Bool

Check if a polynomial is the zero polynomial (has no terms).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1]);

julia> p = Polynomial([Term(1.0+0im, m), Term(-1.0+0im, m)]);  # Cancels out

julia> iszero(p)
true
```
"""
Base.iszero(p::Polynomial) = isempty(p.terms)

"""
    Base.isone(p::Polynomial) -> Bool

Check if a polynomial is the multiplicative identity (single term with coefficient 1 and identity monomial).

# Examples
```jldoctest
julia> using FastPolynomials

julia> p = one(Polynomial{PauliAlgebra,Int64,Float64});

julia> isone(p)
true

julia> p2 = Polynomial{PauliAlgebra,Int64,Float64}(2.0);

julia> isone(p2)
false
```
"""
function Base.isone(p::Polynomial)
    length(p.terms) == 1 || return false
    return isone(p.terms[1])
end

# =============================================================================
# Accessors
# =============================================================================

"""
    terms(p::Polynomial) -> Vector{Term}

Get the terms of the polynomial. Terms are sorted by monomial and have unique monomials.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([2, 3]);

julia> p = Polynomial([Term(1.0+0im, m1), Term(2.0+0im, m2)]);

julia> length(terms(p))
2
```
"""
terms(p::Polynomial) = p.terms

"""
    coefficients(p::Polynomial) -> Vector{C}

Extract the coefficients of all terms in the polynomial.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([2]);

julia> p = Polynomial([Term(1.0+0im, m1), Term(2.0+0im, m2)]);

julia> coefficients(p)
2-element Vector{ComplexF64}:
 1.0 + 0.0im
 2.0 + 0.0im
```
"""
coefficients(p::Polynomial) = [t.coefficient for t in p.terms]

"""
    monomials(p::Polynomial) -> Vector{Monomial}

Extract the monomials of all terms in the polynomial.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([2]);

julia> p = Polynomial([Term(1.0+0im, m1), Term(2.0+0im, m2)]);

julia> ms = monomials(p);

julia> length(ms)
2
```
"""
monomials(p::Polynomial) = [t.monomial for t in p.terms]

"""
    degree(p::Polynomial) -> Int

Compute the maximum degree of all monomials in the polynomial.
Returns -1 for the zero polynomial (convention).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([2, 3, 4]);

julia> p = Polynomial([Term(1.0+0im, m1), Term(2.0+0im, m2)]);

julia> degree(p)
3
```
"""
function degree(p::Polynomial)
    isempty(p.terms) && return -1  # Convention for zero polynomial
    return maximum(degree(t.monomial) for t in p.terms)
end

"""
    maxdegree(p::Polynomial) -> Int

Alias for `degree(p)`. Compute the maximum degree of all monomials.
"""
maxdegree(p::Polynomial) = degree(p)

"""
    support(p::Polynomial) -> Vector{Monomial}

Get the support of the polynomial (all unique monomials with non-zero coefficients).
For a valid polynomial, this is the same as `monomials(p)`.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> p = Polynomial(Term(1.0+0im, m));

julia> support(p)
1-element Vector{Monomial{PauliAlgebra, Int64}}:
 Monomial{PauliAlgebra, Int64}([1, 2], 0xf28b9e193dfe87c0)
```
"""
support(p::Polynomial) = monomials(p)

"""
    variable_indices(p::Polynomial) -> Set

Get the set of all variable indices used in the polynomial's monomials.
Returns a Set of integer indices.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1, 2]);

julia> m2 = Monomial{PauliAlgebra}([2, 3]);

julia> p = Polynomial([Term(1.0+0im, m1), Term(1.0+0im, m2)]);

julia> variable_indices(p)
Set{Int64} with 3 elements:
  2
  3
  1
```
"""
function variable_indices(p::Polynomial{A,T,C}) where {A,T,C}
    result = Set{T}()
    for t in p.terms
        for idx in t.monomial.word
            push!(result, abs(idx))  # abs for fermionic/bosonic (negative = annihilation)
        end
    end
    return result
end

"""
    variable_indices(m::Monomial{A,T}) -> Set{T}

Get the set of all variable indices used in a monomial.
Returns a Set of integer indices.

For signed index types (Fermionic/Bosonic), uses `abs(idx)` to normalize
since creation (+) and annihilation (-) refer to the same mode.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2, 1]);

julia> variable_indices(m)
Set{Int64} with 2 elements:
  2
  1
```
"""
function variable_indices(m::Monomial{A,T}) where {A<:AlgebraType, T<:Integer}
    result = Set{T}()
    for idx in m.word
        push!(result, T(abs(idx)))  # abs for fermionic/bosonic (negative = annihilation)
    end
    return result
end

# =============================================================================
# Equality and Hashing
# =============================================================================

"""
    Base.:(==)(p1::Polynomial, p2::Polynomial) -> Bool

Check if two polynomials are equal. Polynomials are equal if they have
the same terms (same monomials with same coefficients).

Since polynomials are maintained in canonical form (sorted, deduplicated),
equality is a straightforward comparison of the terms vectors.
"""
function Base.:(==)(
    p1::Polynomial{A1,T1,C1}, p2::Polynomial{A2,T2,C2}
) where {A1,A2,T1,T2,C1,C2}
    # Different algebra types are never equal
    A1 !== A2 && return false
    T1 !== T2 && return false

    # Compare terms (already sorted and deduplicated)
    length(p1.terms) != length(p2.terms) && return false

    for (t1, t2) in zip(p1.terms, p2.terms)
        t1.monomial != t2.monomial && return false
        t1.coefficient != t2.coefficient && return false
    end
    return true
end

"""
    Base.hash(p::Polynomial, h::UInt) -> UInt

Hash function for polynomial. Combines hashes of all terms.
"""
function Base.hash(p::Polynomial, h::UInt)
    h = hash(:Polynomial, h)
    for t in p.terms
        h = hash(t.coefficient, hash(t.monomial, h))
    end
    return h
end

# =============================================================================
# Display
# =============================================================================

"""
    Base.show(io::IO, p::Polynomial)

Display a polynomial as a sum of terms.
"""
function Base.show(io::IO, p::Polynomial{A,T,C}) where {A,T,C}
    if isempty(p.terms)
        print(io, "0")
        return nothing
    end

    first_term = true
    for t in p.terms
        if !first_term
            print(io, " + ")
        end
        print(io, t)
        first_term = false
    end
end

# =============================================================================
# Arithmetic Operations
# =============================================================================

"""
    Base.:(+)(p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}) where {A,T,C1,C2}

Add two polynomials. Concatenates terms and processes to combine like terms.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{PauliAlgebra}([1]);

julia> m2 = Monomial{PauliAlgebra}([2]);

julia> p1 = Polynomial(Term(1.0+0im, m1));

julia> p2 = Polynomial(Term(2.0+0im, m2));

julia> p_sum = p1 + p2;

julia> length(terms(p_sum))
2
```
"""
function Base.:(+)(
    p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}
) where {A<:AlgebraType,T<:Integer,C1<:Number,C2<:Number}
    C = promote_type(C1, C2)

    # Convert terms to promoted coefficient type and concatenate
    terms1 = [Term(C(t.coefficient), t.monomial) for t in p1.terms]
    terms2 = [Term(C(t.coefficient), t.monomial) for t in p2.terms]

    return Polynomial(vcat(terms1, terms2))
end

"""
    Base.:(-)(p1::Polynomial, p2::Polynomial)

Subtract two polynomials.
"""
function Base.:(-)(
    p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}
) where {A<:AlgebraType,T<:Integer,C1<:Number,C2<:Number}
    return p1 + (-p2)
end

"""
    Base.:(-)(p::Polynomial)

Negate a polynomial (negate all coefficients).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1]);

julia> p = Polynomial(Term(2.0+0im, m));

julia> p_neg = -p;

julia> coefficients(p_neg)[1]
-2.0 - 0.0im
```
"""
function Base.:(-)(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    negated_terms = [Term(-t.coefficient, t.monomial) for t in p.terms]
    # No need to re-process since just negating coefficients preserves invariants
    return Polynomial{A,T,C}(negated_terms)
end

"""
    Base.:(*)(p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}) where {A,T,C1,C2}

Multiply two polynomials. Computes the Cartesian product of terms and
processes the result (handling monomial simplification which may return
single Term or Vector{Term}).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m1 = Monomial{NonCommutativeAlgebra}(UInt16[1]);

julia> m2 = Monomial{NonCommutativeAlgebra}(UInt16[2]);

julia> p1 = Polynomial(Term(2.0, m1));

julia> p2 = Polynomial(Term(3.0, m2));

julia> p_prod = p1 * p2;

julia> coefficients(p_prod)[1]
6.0
```
"""
function Base.:(*)(
    p1::Polynomial{A,T,C1}, p2::Polynomial{A,T,C2}
) where {A<:AlgebraType,T<:Integer,C1<:Number,C2<:Number}
    # Handle zero cases
    isempty(p1.terms) && return zero(Polynomial{A,T,promote_type(C1, C2)})
    isempty(p2.terms) && return zero(Polynomial{A,T,promote_type(C1, C2)})

    C = promote_type(C1, C2)
    result_terms = Term{Monomial{A,T},C}[]

    for t1 in p1.terms
        for t2 in p2.terms
            coef = C(t1.coefficient * t2.coefficient)
            iszero(coef) && continue

            # Monomial multiplication - may return Term or Vector{Term}
            simplified = t1.monomial * t2.monomial
            _add_simplified_terms!(result_terms, coef, simplified)
        end
    end

    return Polynomial(result_terms)
end

"""
    _add_simplified_terms!(result, coef, simplified)

Helper function to add simplified monomial multiplication results.
Handles Monomial (raw concatenation), Term, or Vector{Term} returns from monomial operations.
"""
function _add_simplified_terms!(
    result::Vector{Term{Monomial{A,T},C}},
    coef::C,
    simplified::Monomial{A,T},
) where {A,T,C}
    # Monomial multiplication now returns raw Monomial (no simplification)
    # Simplify it here to get Term or Vector{Term}
    term_result = simplify!(simplified)
    return _add_simplified_terms!(result, coef, term_result)
end

function _add_simplified_terms!(
    result::Vector{Term{Monomial{A,T},C}},
    coef::C,
    simplified::Term{Monomial{A,T},SC},
) where {A,T,C,SC}
    return _add_simplified_terms!(result, coef, [simplified])
end

function _add_simplified_terms!(
    result::Vector{Term{Monomial{A,T},C}},
    coef::C,
    simplified::Vector{Term{Monomial{A,T},SC}},
) where {A,T,C,SC}
    for term in simplified
        combined_coef = C(coef * term.coefficient)
        if !iszero(combined_coef)
            push!(result, Term(combined_coef, term.monomial))
        end
    end
end

"""
    Base.:(*)(c::Number, p::Polynomial)

Scalar multiplication (scalar on left).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1]);

julia> p = Polynomial(Term(2.0+0im, m));

julia> p2 = 3.0 * p;

julia> coefficients(p2)[1]
6.0 + 0.0im
```
"""
function Base.:(*)(
    c::Number, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    iszero(c) && return zero(Polynomial{A,T,promote_type(typeof(c), C)})

    NC = promote_type(typeof(c), C)
    scaled_terms = [Term(NC(c * t.coefficient), t.monomial) for t in p.terms]
    # Scaling preserves invariants (sorted, unique, non-zero if original was)
    return Polynomial{A,T,NC}(scaled_terms)
end

"""
    Base.:(*)(p::Polynomial, c::Number)

Scalar multiplication (scalar on right).
"""
Base.:(*)(p::Polynomial, c::Number) = c * p

"""
    Base.:(*)(c::Number, m::Monomial{A,T}) where {A,T}

Scalar multiplication of a monomial. Returns a Polynomial with a single term.
This is useful for constructing polynomials from scalar-monomial expressions.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> p = 3.0 * m;

julia> coefficients(p)[1]
3.0
```
"""
function Base.:(*)(c::Number, m::Monomial{A,T}) where {A<:AlgebraType,T<:Integer}
    iszero(c) && return Polynomial{A,T,typeof(c)}(Term{Monomial{A,T},typeof(c)}[])
    return Polynomial([Term(c, m)])
end

"""
    Base.:(*)(m::Monomial, c::Number)

Scalar multiplication of a monomial (monomial on left).
"""
Base.:(*)(m::Monomial, c::Number) = c * m

"""
    Base.:(/)(p::Polynomial, c::Number)

Divide polynomial by scalar.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1]);

julia> p = Polynomial(Term(6.0+0im, m));

julia> p2 = p / 2.0;

julia> coefficients(p2)[1]
3.0 + 0.0im
```
"""
function Base.:(/)(
    p::Polynomial{A,T,C}, c::Number
) where {A<:AlgebraType,T<:Integer,C<:Number}
    iszero(c) && throw(DivideError())
    inv_c = inv(c)
    return inv_c * p
end

"""
    Base.:(+)(p::Polynomial, c::Number)

Add a scalar to a polynomial (adds constant term).
"""
function Base.:(+)(
    p::Polynomial{A,T,C}, c::Number
) where {A<:AlgebraType,T<:Integer,C<:Number}
    NC = promote_type(typeof(c), C)
    # Create constant term with identity monomial
    const_term = Term(NC(c), one(Monomial{A,T}))
    const_poly = Polynomial([const_term])
    # Convert p to new coefficient type and add
    p_converted = Polynomial([Term(NC(t.coefficient), t.monomial) for t in p.terms])
    return p_converted + const_poly
end
Base.:(+)(c::Number, p::Polynomial) = p + c

"""
    Base.:(+)(p::Polynomial{A,T,C}, m::Monomial{A,T}) where {A,T,C}
    Base.:(+)(m::Monomial{A,T}, p::Polynomial{A,T,C}) where {A,T,C}

Add a monomial to a polynomial. Converts the monomial to a single-term polynomial and adds.

# Examples
```jldoctest
julia> using FastPolynomials

julia> p = Polynomial([Term(1.0, Monomial{NonCommutativeAlgebra}([1]))]);

julia> m = Monomial{NonCommutativeAlgebra}([2, 3]);

julia> result = p + m;

julia> length(terms(result))
2

julia> result2 = m + p;

julia> length(terms(result2))
2
```
"""
function Base.:(+)(
    p::Polynomial{A,T,C}, m::Monomial{A,T}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Convert monomial to polynomial and add
    m_poly = Polynomial([Term(one(C), m)])
    return p + m_poly
end

Base.:(+)(m::Monomial{A,T}, p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number} = p + m

"""
    Base.:(-)(p::Polynomial{A,T,C}, m::Monomial{A,T}) where {A,T,C}
    Base.:(-)(m::Monomial{A,T}, p::Polynomial{A,T,C}) where {A,T,C}

Subtract a monomial from a polynomial or vice versa.

# Examples
```jldoctest
julia> using FastPolynomials

julia> p = Polynomial([Term(1.0, Monomial{NonCommutativeAlgebra}([1]))]);

julia> m = Monomial{NonCommutativeAlgebra}([2, 3]);

julia> result = p - m;

julia> length(terms(result))
2

julia> result2 = m - p;

julia> length(terms(result2))
2
```
"""
function Base.:(-)(
    p::Polynomial{A,T,C}, m::Monomial{A,T}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Convert monomial to polynomial with negative coefficient and add
    m_poly = Polynomial([Term(-one(C), m)])
    return p + m_poly
end

function Base.:(-)(
    m::Monomial{A,T}, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Convert monomial to polynomial and subtract p
    m_poly = Polynomial([Term(one(C), m)])
    return m_poly - p
end

"""
    Base.:(-)(c::Number, p::Polynomial)

Subtract a polynomial from a scalar.
"""
function Base.:(-)(
    c::Number, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    NC = promote_type(typeof(c), C)
    const_term = Term(NC(c), one(Monomial{A,T}))
    const_poly = Polynomial([const_term])
    return const_poly + (-p)
end

"""
    Base.:(-)(p::Polynomial, c::Number)

Subtract a scalar from a polynomial.
"""
Base.:(-)(p::Polynomial, c::Number) = p + (-c)

"""
    Base.:*(m::Monomial{A,T}, p::Polynomial{A,T,C}) where {A,T,C}

Multiply a monomial by a polynomial. Returns a Polynomial.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{NonCommutativeAlgebra}([1]);

julia> p = Polynomial([Term(1.0, Monomial{NonCommutativeAlgebra}([2])), Term(2.0, Monomial{NonCommutativeAlgebra}([3]))]);

julia> result = m * p;

julia> length(terms(result))
2
```
"""
function Base.:*(
    m::Monomial{A,T}, p::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Convert monomial to polynomial and multiply
    return Polynomial(m) * p
end

"""
    Base.:*(p::Polynomial{A,T,C}, m::Monomial{A,T}) where {A,T,C}

Multiply a polynomial by a monomial. Returns a Polynomial.

# Examples
```jldoctest
julia> using FastPolynomials

julia> p = Polynomial([Term(1.0, Monomial{NonCommutativeAlgebra}([1])), Term(2.0, Monomial{NonCommutativeAlgebra}([2]))]);

julia> m = Monomial{NonCommutativeAlgebra}([3]);

julia> result = p * m;

julia> length(terms(result))
2
```
"""
function Base.:*(
    p::Polynomial{A,T,C}, m::Monomial{A,T}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Convert monomial to polynomial and multiply
    return p * Polynomial(m)
end

"""
    Base.:(^)(p::Polynomial, n::Int)

Raise polynomial to integer power using binary exponentiation (power by squaring).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{NonCommutativeAlgebra}(UInt16[1]);

julia> p = Polynomial(Term(2.0, m));

julia> p0 = p^0;

julia> isone(p0)
true

julia> p1 = p^1;

julia> p1 == p
true
```
"""
function Base.:(^)(p::Polynomial{A,T,C}, n::Int) where {A<:AlgebraType,T<:Integer,C<:Number}
    n < 0 && throw(DomainError(n, "Polynomial exponent must be non-negative"))
    return Base.power_by_squaring(p, n)
end

# =============================================================================
# Adjoint / Star Operations
# =============================================================================

"""
    Base.adjoint(p::Polynomial{A,T,C}) where {A,T,C}

Compute the adjoint (star/dagger) of a polynomial.

The adjoint of a polynomial is computed by:
1. Taking the adjoint of each monomial (reverse word, negate for signed types)
2. Taking the complex conjugate of each coefficient

This ensures that for polynomial operators, `adjoint(AB) = adjoint(B) * adjoint(A)`.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2, 3]);

julia> p = Polynomial(Term(1.0 + 2.0im, m));

julia> p_adj = adjoint(p);

julia> coefficients(p_adj)[1]
1.0 - 2.0im

julia> monomials(p_adj)[1].word
3-element Vector{Int64}:
 3
 2
 1
```

See also: [`star`](@ref), [`adjoint(::Monomial)`](@ref)
"""
function Base.adjoint(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(p.terms) && return zero(Polynomial{A,T,C})

    # Adjoint of polynomial: adjoint each monomial, conjugate each coefficient
    new_terms = [Term(conj(t.coefficient), adjoint(t.monomial)) for t in p.terms]
    return Polynomial(new_terms)
end

"""
    star(p::Polynomial)

Alias for `adjoint`. Compute the star (dagger) of a polynomial.
This notation is common in physics for the adjoint/Hermitian conjugate.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> p = Polynomial(Term(1.0 + 1.0im, m));

julia> star(p) == adjoint(p)
true
```

See also: [`adjoint`](@ref)
"""
star(p::Polynomial) = adjoint(p)

"""
    is_symmetric(p::Polynomial) -> Bool

Check if a polynomial is symmetric (Hermitian), i.e., p == adjoint(p).

A polynomial is symmetric if it equals its adjoint. This is equivalent to
checking if the polynomial represents a Hermitian operator.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> p_real = Polynomial(Term(1.0 + 0.0im, m));

julia> is_symmetric(p_real)  # Real coefficient, symmetric word
true

julia> p_complex = Polynomial(Term(1.0 + 1.0im, m));

julia> is_symmetric(p_complex)  # Complex coefficient breaks symmetry
false
```

Self-adjoint monomials with real coefficients are symmetric:
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1]);  # Single Pauli is self-adjoint

julia> p = Polynomial(Term(2.0 + 0.0im, m));

julia> is_symmetric(p)
true
```

See also: [`adjoint`](@ref), [`star`](@ref)
"""
is_symmetric(p::Polynomial) = p == adjoint(p)
