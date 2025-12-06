"""
    Term{M<:AbstractMonomial, C<:Number}

Represents a coefficient paired with a monomial, the result of algebraic simplification.

# Fields
- `coefficient::C`: Numeric coefficient (Float64 or ComplexF64)
- `monomial::M`: The simplified monomial (Monomial or ComposedMonomial)

# Type Parameters
- `M<:AbstractMonomial`: Monomial type (Monomial{A,T} or ComposedMonomial{Ts})
- `C<:Number`: Coefficient type (Float64, ComplexF64, etc.)

# Design
The Term struct is the output of simplification algorithms. Each algebra type
returns a Term with an appropriate coefficient type:
- Pauli: `ComplexF64` (due to i phases from cyclic products)
- Fermionic/Bosonic/Projector/Unipotent: `Float64` (real coefficients only)

Term is mutable to allow coefficient modification during polynomial accumulation.

# Examples
```jldoctest
julia> m = Monomial{PauliAlgebra}([1, 2, 3]);

julia> t = Term(1.0 + 0.0im, m);

julia> t.coefficient
1.0 + 0.0im

julia> t.monomial.word
3-element Vector{Int64}:
 1
 2
 3
```

Identity term (empty monomial with coefficient 1):
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{UnipotentAlgebra}(Int[]);

julia> t = Term(1.0, m);

julia> isone(t)
true

julia> t.coefficient
1.0
```
"""
mutable struct Term{M<:AbstractMonomial,C<:Number}
    coefficient::C
    monomial::M
end

"""
    Base.isone(t::Term) -> Bool

Check if a term is the multiplicative identity (coefficient 1, empty monomial).

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{UnipotentAlgebra}(Int[]);

julia> t = Term(1.0, m);

julia> isone(t)
true

julia> t2 = Term(2.0, Monomial{UnipotentAlgebra}([1]));

julia> isone(t2)
false
```
"""
function Base.isone(t::Term{M,C}) where {M<:AbstractMonomial,C}
    isone(t.coefficient) && isone(t.monomial)
end

"""
    Base.iszero(t::Term) -> Bool

Check if a term is zero (coefficient is zero).

# Examples
```jldoctest
julia> m = Monomial{FermionicAlgebra}(Int32[]);

julia> t = Term(0.0, m);

julia> iszero(t)
true

julia> t2 = Term(1.0, Monomial{FermionicAlgebra}(Int32[]));

julia> iszero(t2)
false
```
"""
Base.iszero(t::Term) = iszero(t.coefficient)

"""
    Base.:(==)(t1::Term, t2::Term) -> Bool

Equality comparison for terms. Two terms are equal if they have the same
coefficient and monomial.
"""
function Base.:(==)(t1::Term{M1,C1}, t2::Term{M2,C2}) where {M1,M2,C1,C2}
    M1 !== M2 && return false
    t1.coefficient == t2.coefficient && t1.monomial == t2.monomial
end

# Show method for clean output
function Base.show(io::IO, t::Term{M,C}) where {M<:Monomial,C}
    if iszero(t)
        print(io, "0")
    elseif isone(t)
        print(io, "1")
    elseif isempty(t.monomial.word)
        print(io, t.coefficient)
    elseif t.coefficient == one(C)
        print(io, t.monomial.word)
    elseif t.coefficient == -one(C)
        print(io, "-", t.monomial.word)
    else
        print(io, t.coefficient, " * ", t.monomial.word)
    end
end

"""
    Base.one(::Type{Term{Monomial{A,T},C}}) where {A<:AlgebraType, T<:Integer, C<:Number}

Create an identity term (coefficient 1 with empty monomial).

# Examples
```jldoctest
julia> using FastPolynomials

julia> one(Term{Monomial{UnipotentAlgebra,Int},Float64})
1

julia> one(Term{Monomial{PauliAlgebra,UInt16},ComplexF64})
1
```
"""
# TODO: shouldn't need to implement this for every specific `AbstractTerm` subtype
function Base.one(::Type{Term{Monomial{A,T},C}}) where {A<:AlgebraType,T<:Integer,C<:Number}
    Term{Monomial{A,T},C}(one(C), Monomial{A}(T[]))
end

"""
    Base.zero(::Type{Term{Monomial{A,T},C}}) where {A<:AlgebraType, T<:Integer, C<:Number}

Create a zero term (coefficient 0 with empty monomial).

# Examples
```jldoctest
julia> using FastPolynomials

julia> zero(Term{Monomial{FermionicAlgebra,Int32},Float64})
0

julia> iszero(zero(Term{Monomial{FermionicAlgebra,Int32},Float64}))
true
```
"""
# TODO: shouldn't need to implement this for every specific `AbstractTerm` subtype
function Base.zero(::Type{Term{Monomial{A,T},C}}) where {A<:AlgebraType,T<:Integer,C<:Number}
    Term{Monomial{A,T},C}(zero(C), Monomial{A}(T[]))
end

# Scalar multiplication
"""
    Base.:*(c::Number, t::Term) -> Term
    Base.:*(t::Term, c::Number) -> Term

Scalar multiplication of a term. Creates a new Term.
"""
function Base.:*(c::Number, t::Term{M,C}) where {M,C}
    NC = promote_type(typeof(c), C)
    Term{M,NC}(NC(c * t.coefficient), t.monomial)
end

Base.:*(t::Term, c::Number) = c * t

"""
    Base.:-(t::Term) -> Term

Negation of a term.
"""
Base.:-(t::Term{M,C}) where {M,C} = Term{M,C}(-t.coefficient, t.monomial)

"""
    Base.iterate(t::Term) -> Tuple
    Base.iterate(t::Term, state) -> Tuple or nothing

Enable tuple-like iteration over a Term, yielding (coefficient, monomial).
This allows destructuring syntax like `(coef, mono) = term`.

# Examples
```jldoctest
julia> using FastPolynomials

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> t = Term(2.5 + 0.0im, m);

julia> coef, mono = t;

julia> coef
2.5 + 0.0im

julia> mono.word
2-element Vector{Int64}:
 1
 2
```
"""
function Base.iterate(t::Term)
    (t.coefficient, 1)
end

function Base.iterate(t::Term, state::Int)
    if state == 1
        (t.monomial, 2)
    else
        nothing
    end
end

# Also define length for completeness
Base.length(::Term) = 2
