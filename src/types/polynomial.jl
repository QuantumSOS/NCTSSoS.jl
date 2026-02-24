"""
    AbstractPolynomial{C<:Number}

Abstract supertype for all polynomial types in NCTSSoS.

# Type Parameters
- `C<:Number`: Coefficient type (Float64, ComplexF64, etc.)

# Subtypes
- `Polynomial{A,T,C}`: Standard polynomial over a non-commutative algebra
- `StatePolynomial{C,ST,A,T}`: Polynomial of state expectations (commutative)
- `NCStatePolynomial{C,ST,A,T}`: Polynomial of state-operator products

# Interface
All subtypes should implement:
- `coefficients(p)`: Return vector of coefficients
- `monomials(p)`: Return vector of monomials/terms
- `degree(p)`: Return maximum degree
- `Base.zero(::Type{P})`, `Base.one(::Type{P})`: Identity elements
- `Base.:(+)`, `Base.:(-)`, `Base.:(*)`: Arithmetic operations
"""
abstract type AbstractPolynomial{C<:Number} end

"""
    Polynomial{A<:AlgebraType, T<:Integer, C<:Number}

A polynomial represented as a sum of `(coefficient, monomial)` pairs.
Maintains sorted, unique monomials with non-zero coefficients.

# Fields
- `terms::Vector{Tuple{C,NormalMonomial{A,T}}}`: Sorted terms with unique monomials and non-zero coefficients

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
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}([1, 4]);

julia> m2 = NormalMonomial{PauliAlgebra}([3]);

julia> p = Polynomial([(1.0 + 0.0im, m1), (2.0 + 0.0im, m2)]);

julia> length(terms(p))
2

julia> degree(p)
2
```

Construction with automatic deduplication:
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{PauliAlgebra}([1]);

julia> p = Polynomial([(1.0+0im, m), (2.0+0im, m)]);  # Same monomial twice

julia> length(terms(p))  # Combined into one term
1

julia> coefficients(p)[1]  # Coefficients added
3.0 + 0.0im
```

See also: [`NormalMonomial`](@ref), [`coefficients`](@ref), [`monomials`](@ref)
"""
struct Polynomial{A<:AlgebraType,T<:Integer,C<:Number} <: AbstractPolynomial{C}
    terms::Vector{Tuple{C,NormalMonomial{A,T}}}

    # Inner constructor enforces invariants: sorted, deduplicated, non-zero coefficients
    function Polynomial{A,T,C}(
        input_terms::Vector{Tuple{C,NormalMonomial{A,T}}}
    ) where {A<:AlgebraType,T<:Integer,C<:Number}
        return new{A,T,C}(input_terms)
    end
end

# Delegate explicit type calls to the inner constructor
function Polynomial(
    input_terms::Vector{Tuple{C,NormalMonomial{A,T}}}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    processed = _process_terms(input_terms, C)
    return Polynomial{A,T,C}(processed)
end

"""
    _process_terms(input_terms, C) -> Vector{Tuple{C,NormalMonomial}}

Process a vector of terms: sort by monomial, combine duplicates, remove zeros.
This is the core algorithm for maintaining polynomial invariants.

# Algorithm
1. If empty, return empty vector
2. Sort terms by monomial (degree-first, then lexicographic)
3. Iterate through sorted terms, combining coefficients for duplicate monomials
4. Filter out terms with zero coefficients
"""
function _process_terms(
    input_terms::Vector{Tuple{C,NormalMonomial{A,T}}}, ::Type{C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(input_terms) && return input_terms

    # Sort by monomial
    sorted = sort(input_terms; by=t -> t[2])

    # Combine duplicates and filter zeros
    result = Tuple{C,NormalMonomial{A,T}}[]
    current_coef, current_mono = sorted[1]

    for i in 2:length(sorted)
        coef_i, mono_i = sorted[i]
        if mono_i == current_mono
            current_coef += coef_i
        else
            if !iszero(current_coef)
                push!(result, (current_coef, current_mono))
            end
            current_mono = mono_i
            current_coef = coef_i
        end
    end

    # Don't forget the last accumulated term
    if !iszero(current_coef)
        push!(result, (current_coef, current_mono))
    end

    return result
end

# =============================================================================
# Constructors
# =============================================================================

"""
    Polynomial((c, m)::Tuple{C,NormalMonomial{A,T}}) where {A,T,C}

Construct a polynomial from a single term.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{PauliAlgebra}([1, 4]);

julia> p = Polynomial((3.0 + 0.0im, m));

julia> length(terms(p))
1
```
"""
function Polynomial(t::Tuple{C,NormalMonomial{A,T}}) where {A<:AlgebraType,T<:Integer,C<:Number}
    coef, _ = t
    if iszero(coef)
        return Polynomial{A,T,C}(Tuple{C,NormalMonomial{A,T}}[])
    end
    return Polynomial{A,T,C}([t])
end


"""
    Polynomial(m::NormalMonomial{A,T}) where {A,T}

Construct a polynomial from a monomial with coefficient 1.

The coefficient type is determined by `coeff_type(A)`:
- `PauliAlgebra`: uses `ComplexF64` (Pauli products generate complex phases)
- All others: uses `Float64`

!!! note "No automatic simplification"
    This constructor does NOT call `simplify` on the monomial. If the monomial
    is not in canonical form (e.g., a Pauli product like `σx₁ * σx₁` that should
    simplify to identity), the resulting polynomial will contain the unsimplified
    monomial. Use `simplify(Polynomial(m))` if canonical form is required.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{PauliAlgebra}([1, 4]);

julia> p = Polynomial(m);

julia> coefficients(p)[1]
1.0 + 0.0im
```
"""
function Polynomial(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    C = coeff_type(A)
    return Polynomial((one(C), m))
end

function Polynomial(
    pairs::Vector{Tuple{CIn,NormalMonomial{A,T}}},
) where {CIn,A<:AlgebraType,T<:Integer}
    return Polynomial{A,T,CIn}(pairs)
end

"""
    Polynomial(t::Tuple{Vector{T},UInt8}) where {T<:Integer}

Construct a Polynomial from a (word, phase) tuple returned by Pauli simplification.

The phase encoding: 0=1, 1=i, 2=-1, 3=-i (representing (im)^phase).

# Examples
```jldoctest
julia> using NCTSSoS

julia> result = (Int64[3], UInt8(1));  # word with phase i

julia> p = Polynomial(result);

julia> coefficients(p)[1]
0.0 + 1.0im
```
"""
function Polynomial(t::Tuple{Vector{T},UInt8}) where {T<:Integer}
    word, phase = t
    # Convert phase encoding to complex coefficient
    coef = _phase_to_complex(phase)
    mono = NormalMonomial{PauliAlgebra,T}(word)
    if iszero(coef)
        return Polynomial{PauliAlgebra,T,ComplexF64}(Tuple{ComplexF64,NormalMonomial{PauliAlgebra,T}}[])
    end
    return Polynomial{PauliAlgebra,T,ComplexF64}([(coef, mono)])
end

# _phase_to_complex helper is in util/helpers.jl

"""
    Polynomial(p::Polynomial{A,T,C}) where {A,T,C}

Identity constructor: returns the polynomial unchanged.
Useful for generic code that may receive NormalMonomial or Polynomial.

# Examples
```jldoctest
julia> using NCTSSoS

julia> p1 = Polynomial([(2.0 + 0.0im, NormalMonomial{PauliAlgebra}([1]))]);

julia> p2 = Polynomial(p1);

julia> p1 === p2
true
```
"""
function Polynomial(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    return p
end

"""
    Polynomial{A,T,C}(c::Number) where {A,T,C}

Construct a constant polynomial (coefficient times identity monomial).

# Examples
```jldoctest
julia> using NCTSSoS

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
        return Polynomial{A,T,C}(Tuple{C,NormalMonomial{A,T}}[])
    end
    return Polynomial{A,T,C}([(C(c), one(NormalMonomial{A,T}))])
end

# =============================================================================
# Type Conversion
# =============================================================================

"""
    Base.convert(::Type{Polynomial{A,T,C2}}, p::Polynomial{A,T,C1}) where {A,T,C1,C2}

Convert a polynomial to a different coefficient type.

This is needed for operations like `power_by_squaring` which may require
converting between different coefficient types during computation.

# Examples
```jldoctest
julia> using NCTSSoS

julia> using NCTSSoS: encode_index

julia> m = NormalMonomial{NonCommutativeAlgebra}(UInt8[encode_index(UInt8, 1, 1)]);

julia> p_int = Polynomial([(2, m)]);  # Int coefficients

julia> p_float = convert(Polynomial{NonCommutativeAlgebra,UInt8,Float64}, p_int);

julia> coefficients(p_float)[1]
2.0
```
"""
function Base.convert(::Type{Polynomial{A,T,C2}}, p::Polynomial{A,T,C1}) where {A<:AlgebraType,T<:Integer,C1<:Number,C2<:Number}
    C1 === C2 && return p
    new_terms = [(C2(c), m) for (c, m) in p.terms]
    return Polynomial{A,T,C2}(new_terms)
end

function Base.convert(
    ::Type{Polynomial{A,T2,C2}}, p::Polynomial{A,T1,C1}
) where {A<:AlgebraType,T1<:Integer,T2<:Integer,C1<:Number,C2<:Number}
    (T1 === T2 && C1 === C2) && return p
    new_terms = [(C2(c), NormalMonomial{A,T2}(T2.(m.word))) for (c, m) in p.terms]
    return Polynomial{A,T2,C2}(new_terms)
end

# Identity conversion (no-op)
Base.convert(::Type{Polynomial{A,T,C}}, p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number} = p

# =============================================================================
# zero and one
# =============================================================================

"""
    Base.zero(::Type{Polynomial{A,T,C}}) where {A,T,C}

Create the zero polynomial (no terms).

# Examples
```jldoctest
julia> using NCTSSoS

julia> p = zero(Polynomial{PauliAlgebra,Int64,Float64});

julia> iszero(p)
true

julia> length(terms(p))
0
```
"""
function Base.zero(::Type{Polynomial{A,T,C}}) where {A<:AlgebraType,T<:Integer,C<:Number}
    return Polynomial{A,T,C}(Tuple{C,NormalMonomial{A,T}}[])
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
julia> using NCTSSoS

julia> p = one(Polynomial{PauliAlgebra,Int64,Float64});

julia> isone(p)
true

julia> coefficients(p)
1-element Vector{Float64}:
 1.0
```
"""
function Base.one(::Type{Polynomial{A,T,C}}) where {A<:AlgebraType,T<:Integer,C<:Number}
    return Polynomial{A,T,C}([(one(C), one(NormalMonomial{A,T}))])
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
    # Shallow copy: coefficients are immutable; monomials share underlying word storage.
    return Polynomial{A,T,C}(copy(p.terms))
end

"""
    Base.iszero(p::Polynomial) -> Bool

Check if a polynomial is the zero polynomial (has no terms).

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{PauliAlgebra}([1]);

julia> p = Polynomial([(1.0+0im, m), (-1.0+0im, m)]);  # Cancels out

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
julia> using NCTSSoS

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
    coef, mono = p.terms[1]
    return isone(coef) && isone(mono)
end

# =============================================================================
# Accessors
# =============================================================================

"""
    terms(p::Polynomial)

Iterate the polynomial as `(coefficient, monomial)` pairs.

The returned monomials are `NormalMonomial`s (canonical form).

# Examples
```jldoctest
julia> using NCTSSoS

julia> reg, (σx, σy, σz) = create_pauli_variables(1:1);

julia> m1 = σx[1];

julia> m2 = σy[1];

julia> p = Polynomial([(1.0+0im, m1), (2.0+0im, m2)]);

julia> length(terms(p))
2
```
"""
terms(p::Polynomial) = p

"""
    coefficients(p::Polynomial) -> Vector{C}

Extract the coefficients of all terms in the polynomial.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}([1]);

julia> m2 = NormalMonomial{PauliAlgebra}([2]);

julia> p = Polynomial([(1.0+0im, m1), (2.0+0im, m2)]);

julia> coefficients(p) == [1.0 + 0.0im, 2.0 + 0.0im]
true
```
"""
coefficients(p::Polynomial) = [c for (c, _) in p.terms]

"""
    monomials(p::Polynomial) -> Vector{NormalMonomial}

Extract the monomials of all terms in the polynomial.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}([1]);

julia> m2 = NormalMonomial{PauliAlgebra}([2]);

julia> p = Polynomial([(1.0+0im, m1), (2.0+0im, m2)]);

julia> ms = monomials(p);

julia> length(ms)
2
```
"""
monomials(p::Polynomial{A,T}) where {A<:AlgebraType,T<:Integer} = [m for (_, m) in p.terms]

"""
    degree(p::Polynomial) -> Union{Int, Float64}

Compute the maximum degree of all monomials in the polynomial.
Returns `-Inf` for the zero polynomial to preserve the algebraic identity
`deg(p * q) = deg(p) + deg(q)`.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}([1]);

julia> m2 = NormalMonomial{PauliAlgebra}([2, 5, 8]);

julia> p = Polynomial([(1.0+0im, m1), (2.0+0im, m2)]);

julia> degree(p)
3

julia> degree(zero(Polynomial{PauliAlgebra,Int64,ComplexF64}))
-Inf
```
"""
function degree(p::Polynomial)
    isempty(p.terms) && return -Inf
    return maximum(degree(m) for (_, m) in p.terms)
end



"""
    variable_indices(p::Polynomial) -> Set

Get the set of all variable indices used in the polynomial's monomials.
Returns a Set of integer indices.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}([1]);

julia> m2 = NormalMonomial{PauliAlgebra}([2]);

julia> m3 = NormalMonomial{PauliAlgebra}([3]);

julia> p = Polynomial([(1.0+0im, m1), (1.0+0im, m2), (1.0+0im, m3)]);

julia> variable_indices(p) == Set([1, 2, 3])
true
```
"""
function variable_indices(p::Polynomial{A,T,C}) where {A,T,C}
    result = Set{T}()
    for (_, mono) in p.terms
        for idx in mono.word
            # For signed types (Fermionic/Bosonic), use abs() to treat creation (-) and
            # annihilation (+) of the same mode as the same variable for sparsity analysis
            push!(result, T(abs(idx)))
        end
    end
    return result
end

"""
    variable_indices(m::NormalMonomial{A,T}) -> Set{T}

Get the set of all variable indices used in a monomial.
Returns a Set of integer indices.

For signed index types (Fermionic/Bosonic), uses `abs(idx)` to normalize indices.
This treats creation operators (negative indices like -1 for a₁†) and annihilation
operators (positive indices like 1 for a₁) as referring to the same physical mode,
which is the correct behavior for correlative sparsity analysis.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = NormalMonomial{FermionicAlgebra}(Int8[-2, -1, 1]);

julia> variable_indices(m) == Set(Int8[1, 2])
true
```
"""
function variable_indices(m::NormalMonomial{A,T}) where {A<:AlgebraType,T<:Integer}
    result = Set{T}()
    for idx in m.word
        # For signed types (Fermionic/Bosonic), use abs() to treat creation (-) and
        # annihilation (+) of the same mode as the same variable for sparsity analysis
        push!(result, T(abs(idx)))
    end
    return result
end

function variable_indices(
    pairs::Vector{Tuple{C,NormalMonomial{A,T}}},
) where {C,A<:AlgebraType,T<:Integer}
    result = Set{T}()
    for (_, mono) in pairs
        for idx in mono.word
            push!(result, T(abs(idx)))
        end
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
        t1 != t2 && return false
    end
    return true
end

"""
    Base.hash(p::Polynomial, h::UInt) -> UInt

Hash function for polynomial. Combines hashes of all terms.
"""
function Base.hash(p::Polynomial, h::UInt)
    h = hash(:Polynomial, h)
    for (c, m) in p.terms
        h = hash(c, hash(m, h))
    end
    return h
end

"""
    Base.isless(p1::Polynomial{A,T,C}, p2::Polynomial{A,T,C}) where {A,T,C} -> Bool

Compare two polynomials for sorting. Uses degree-first (graded) ordering:
1. Compare by highest-degree monomial first
2. For equal monomials, compare by coefficient magnitude
3. Continue to next-highest degree monomial if still equal
4. Polynomial with fewer terms is "less" if all compared terms are equal
"""
function Base.isless(
    p1::Polynomial{A,T,C}, p2::Polynomial{A,T,C}
) where {A<:AlgebraType,T<:Integer,C<:Number}
    # Empty polynomials
    isempty(p1.terms) && isempty(p2.terms) && return false
    isempty(p1.terms) && return true   # zero < nonzero
    isempty(p2.terms) && return false

    # Compare from highest-degree terms (terms are sorted low→high, so iterate from end)
    i1, i2 = length(p1.terms), length(p2.terms)
    while i1 > 0 && i2 > 0
        c1, m1 = p1.terms[i1]
        c2, m2 = p2.terms[i2]
        # Compare monomial first
        m1 != m2 && return isless(m1, m2)
        # Same monomial: compare coefficient magnitude
        abs(c1) != abs(c2) && return abs(c1) < abs(c2)
        i1 -= 1
        i2 -= 1
    end
    # Exhausted one polynomial - fewer terms is "less"
    return i1 < i2
end

function _show_term(io::IO, coef::Number, mono)
    if iszero(coef)
        print(io, "0")
    elseif isone(coef) && isone(mono)
        print(io, "1")
    elseif isempty(mono.word)
        print(io, coef)
    elseif coef == one(typeof(coef))
        show(io, mono)
    elseif coef == -one(typeof(coef))
        print(io, "-")
        show(io, mono)
    else
        if !isreal(coef)
            print(io, "(" * string(coef) * ") * ")
        else
            print(io, coef, " * ")
        end
        show(io, mono)
    end
    return nothing
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
    for (coef, mono) in p.terms
        if !first_term
            print(io, " + ")
        end
        _show_term(io, coef, mono)
        first_term = false
    end
end

# =============================================================================
# Iteration Protocol
# =============================================================================

"""
    Base.iterate(p::Polynomial{A,T,C})
    Base.iterate(p::Polynomial{A,T,C}, state::Int)

Iterate a Polynomial, yielding `(coefficient, NormalMonomial)` pairs for each term.

# Examples
```jldoctest
julia> using NCTSSoS

julia> m1 = NormalMonomial{PauliAlgebra}([1]);

julia> m2 = NormalMonomial{PauliAlgebra}([2]);

julia> p = Polynomial([(1.0+0im, m1), (2.0+0im, m2)]);

julia> length(collect(p))
2

julia> for (coef, mono) in p
           println("Coefficient: \$coef, Degree: \$(degree(mono))")
       end
Coefficient: 1.0 + 0.0im, Degree: 1
Coefficient: 2.0 + 0.0im, Degree: 1
```
"""
Base.eltype(::Type{Polynomial{A,T,C}}) where {A<:AlgebraType,T<:Integer,C<:Number} =
    Tuple{C,NormalMonomial{A,T}}

function Base.iterate(p::Polynomial{A,T,C}) where {A<:AlgebraType,T<:Integer,C<:Number}
    isempty(p.terms) && return nothing
    coef, mono = p.terms[1]
    return ((coef, mono), 2)
end

function Base.iterate(
    p::Polynomial{A,T,C}, state::Int
) where {A<:AlgebraType,T<:Integer,C<:Number}
    state > length(p.terms) && return nothing
    coef, mono = p.terms[state]
    return ((coef, mono), state + 1)
end

# Length for iteration (number of terms)
Base.length(p::Polynomial) = length(p.terms)

"""
    coeff_type(::Type{<:AbstractPolynomial{C}}) where {C} -> Type{<:Number}
    coeff_type(p::AbstractPolynomial) -> Type{<:Number}

Return the coefficient type `C` for any `AbstractPolynomial` subtype.

Works for all polynomial types including `Polynomial{A,T,C}` and
`NCStatePolynomial{C,ST,A,T}` since both subtype `AbstractPolynomial{C}`.

# Examples
```jldoctest
julia> coeff_type(Polynomial{PauliAlgebra,Int64,ComplexF64}) == ComplexF64
true

julia> coeff_type(Polynomial{NonCommutativeAlgebra,Int64,Float64}) == Float64
true
```
"""
coeff_type(::Type{<:AbstractPolynomial{C}}) where {C<:Number} = C

# Instance method
coeff_type(p::AbstractPolynomial{C}) where {C<:Number} = C

simplify(p::Polynomial) = p
simplify!(p::Polynomial) = p
