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

Term is immutable for performance (enables stack allocation and inlining).
The monomial field stores a reference, so no deep copying occurs.

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
julia> using NCTSSoS

julia> m = Monomial{UnipotentAlgebra}(Int[]);

julia> t = Term(1.0, m);

julia> isone(t)
true

julia> t.coefficient
1.0
```
"""
struct Term{M<:AbstractMonomial,C<:Number}
    coefficient::C
    monomial::M
end

"""
    Base.isone(t::Term) -> Bool

Check if a term is the multiplicative identity (coefficient 1, empty monomial).

# Examples
```jldoctest
julia> using NCTSSoS

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
    return isone(t.coefficient) && isone(t.monomial)
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
    return t1.coefficient == t2.coefficient && t1.monomial == t2.monomial
end

"""
    Base.hash(t::Term, h::UInt) -> UInt

Hash function for Term. Combines hashes of coefficient and monomial.

This ensures the hash/equality contract: if `t1 == t2`, then `hash(t1) == hash(t2)`.
"""
function Base.hash(t::Term, h::UInt)
    h = hash(t.coefficient, h)
    h = hash(t.monomial, h)
    return h
end

# Show method for clean output
# Passes IOContext to monomial show for registry-aware display
function Base.show(io::IO, t::Term{M,C}) where {M<:Monomial,C}
    if iszero(t)
        print(io, "0")
    elseif isone(t)
        print(io, "1")
    elseif isempty(t.monomial.word)
        print(io, t.coefficient)
    elseif t.coefficient == one(C)
        show(io, t.monomial)
    elseif t.coefficient == -one(C)
        print(io, "-")
        show(io, t.monomial)
    else
        if !isreal(t.coefficient)
            print(io, "(" * string(t.coefficient) * ") * ")
        else
            print(io, t.coefficient, " * ")
        end
        show(io, t.monomial)
    end
end

"""
    Base.one(::Type{Term{M,C}}) where {M<:AbstractMonomial, C<:Number}

Create an identity term (coefficient 1 with identity monomial).

Works for any `AbstractMonomial` subtype (`Monomial` or `ComposedMonomial`)
as long as the monomial type implements `one(::Type{M})`.

# Examples
```jldoctest
julia> using NCTSSoS

julia> one(Term{Monomial{UnipotentAlgebra,Int},Float64})
1

julia> one(Term{Monomial{PauliAlgebra,UInt16},ComplexF64})
1
```
"""
function Base.one(::Type{Term{M,C}}) where {M<:AbstractMonomial,C<:Number}
    return Term{M,C}(one(C), one(M))
end

"""
    Base.zero(::Type{Term{M,C}}) where {M<:AbstractMonomial, C<:Number}

Create a zero term (coefficient 0 with identity monomial).

Works for any `AbstractMonomial` subtype (`Monomial` or `ComposedMonomial`)
as long as the monomial type implements `one(::Type{M})`.

# Examples
```jldoctest
julia> using NCTSSoS

julia> zero(Term{Monomial{FermionicAlgebra,Int32},Float64})
0

julia> iszero(zero(Term{Monomial{FermionicAlgebra,Int32},Float64}))
true
```
"""
function Base.zero(::Type{Term{M,C}}) where {M<:AbstractMonomial,C<:Number}
    return Term{M,C}(zero(C), one(M))
end

# =============================================================================
# Scalar Multiplication
# =============================================================================

"""
    Base.:*(c::Number, t::Term) -> Term
    Base.:*(t::Term, c::Number) -> Term

Scalar multiplication of a term. Returns a new Term with promoted coefficient type.
"""
function Base.:*(c::Number, t::Term{M,C}) where {M,C}
    NC = promote_type(typeof(c), C)
    return Term{M,NC}(NC(c * t.coefficient), t.monomial)
end

Base.:*(t::Term, c::Number) = c * t

# =============================================================================
# Negation
# =============================================================================

"""
    Base.:-(t::Term) -> Term

Negation of a term. Returns a new Term.
"""
Base.:-(t::Term{M,C}) where {M,C} = Term{M,C}(-t.coefficient, t.monomial)

# =============================================================================
# Iteration Protocol
# =============================================================================

"""
    Base.iterate(t::Term{M,C}) -> Tuple{Tuple{C, M}, Nothing}
    Base.iterate(::Term, ::Nothing) -> Nothing

Iterate a Term as a single `(coefficient, monomial)` pair.

This enables uniform processing of simplify results across all algebra types.
A Term yields exactly one `(coefficient, monomial)` tuple when iterated.

!!! warning "Breaking Change"
    Previous versions yielded coefficient then monomial separately for destructuring
    (`coef, mono = term`). Now yields a single (coef, mono) tuple. Update code like:
    - Old: `coef, mono = term`
    - New: `(coef, mono), = term` or `coef, mono = term.coefficient, term.monomial`

# Examples
```jldoctest
julia> using NCTSSoS

julia> m = Monomial{PauliAlgebra}([1, 2]);

julia> t = Term(2.5 + 0.0im, m);

julia> collect(t)
1-element Vector{Tuple{ComplexF64, Monomial{PauliAlgebra, Int64}}}:
 (2.5 + 0.0im, Monomial{PauliAlgebra, Int64}([1, 2]))

julia> (coef, mono), = t;

julia> coef
2.5 + 0.0im

julia> mono.word
2-element Vector{Int64}:
 1
 2
```
"""
function Base.iterate(t::Term{M,C}) where {M<:AbstractMonomial,C<:Number}
    return ((t.coefficient, t.monomial), nothing)
end

Base.iterate(::Term, ::Nothing) = nothing

Base.length(::Term) = 1

Base.eltype(::Type{Term{M,C}}) where {M<:AbstractMonomial,C<:Number} = Tuple{C,M}

"""
    coeff_type(::Type{Term{M,C}}) where {M,C} -> Type{<:Number}

Return the coefficient type C for a Term type.
"""
coeff_type(::Type{Term{M,C}}) where {M<:AbstractMonomial,C<:Number} = C

# Instance method
coeff_type(t::Term{M,C}) where {M<:AbstractMonomial,C<:Number} = C
